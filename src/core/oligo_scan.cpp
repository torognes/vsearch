/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.


  GNU General Public License version 3

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The BSD 2-Clause License

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "core/oligo_scan.hpp"
#include "utils/maps/four_bit.hpp"  // vsearch::maps::four_bit::map
#include <algorithm>  // std::min, std::sort, std::fill, std::transform, std::reverse
#include <cassert>  // assert
#include <cmath>  // std::ceil
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <limits>  // std::numeric_limits
#include <string>
#include <tuple>  // std::tie
#include <vector>


namespace vsearch {
  namespace oligo {

    namespace {

      constexpr std::size_t code_count = 16;  /* 4-bit nucleotide codes */
      constexpr unsigned char code_of_N = 15;
      constexpr int64_t no_candidate = -1;
      constexpr auto word_bits = static_cast<int>(max_oligo_length);

      /* The four bits of a code are A, C, G and T: complementing a base set
         swaps A with T and C with G. */
      constexpr unsigned char bit_a = 1U;
      constexpr unsigned char bit_c = 2U;
      constexpr unsigned char bit_g = 4U;
      constexpr unsigned char bit_t = 8U;

      constexpr auto complement(unsigned char const code) noexcept -> unsigned char
      {
        return static_cast<unsigned char>(
          (((code & bit_a) != 0U) ? bit_t : 0U) |
          (((code & bit_c) != 0U) ? bit_g : 0U) |
          (((code & bit_g) != 0U) ? bit_c : 0U) |
          (((code & bit_t) != 0U) ? bit_a : 0U));
      }

      /* the lowest `count` bits set */
      auto low_bits(int const count) noexcept -> uint64_t
      {
        assert(count >= 0);
        assert(count <= word_bits);
        if (count == word_bits) { return ~uint64_t{0}; }
        return (uint64_t{1} << static_cast<unsigned int>(count)) - 1U;
      }

      auto bit_is_set(uint64_t const word, int const position) noexcept -> bool
      {
        assert(position >= 0);
        assert(position < word_bits);
        return ((word >> static_cast<unsigned int>(position)) & 1U) != 0U;
      }

      /* The number of pattern bases that may hang off a read end: an
         occurrence must still align at least --target_cov of the oligo, and
         at least one base. */
      auto free_overhang(int const length, double const target_cov) noexcept -> int
      {
        assert(length > 0);
        auto const needed = static_cast<int>(std::ceil(target_cov * length));
        auto const aligned = std::min(length, std::max(1, needed));
        return length - aligned;
      }

      auto is_match(unsigned char const read_code,
                    unsigned char const pattern_code,
                    NPositions const n_positions) noexcept -> bool
      {
        if ((n_positions == NPositions::mismatch) and
            ((read_code == code_of_N) or (pattern_code == code_of_N)))
          {
            return false;
          }
        return (read_code & pattern_code) != 0U;
      }


      /* The alignment of one occurrence, recomputed over a window of the read:
         a dynamic programming whose cost is the pair (differences, gap
         openings) compared lexicographically, so that among the alignments
         of fewest differences the one with the fewest gap openings wins.
         Both are packed into one integer. */
      constexpr int64_t opening_weight = 1;
      constexpr int64_t difference_weight = 1024;  /* > any count of openings */
      constexpr int64_t unreachable = std::numeric_limits<int64_t>::max() / 2;

      /* the three ways a cell is entered, and the start of an alignment */
      enum struct Move : unsigned char { diagonal, read_gap, pattern_gap, start };
      constexpr std::size_t move_count = 3;  /* the moves stored per cell */

      auto move_index(Move const move) noexcept -> std::size_t
      {
        return static_cast<std::size_t>(move);
      }

    }  // anonymous namespace


    LaneSet::LaneSet(std::vector<View<char>> const & oligos,
                     Strands const strands,
                     NPositions const n_positions,
                     Limits const & limits)
      : limits_(limits)
    {
      for (std::size_t oligo = 0; oligo < oligos.size(); ++oligo)
        {
          auto const & sequence = oligos[oligo];
          assert(not sequence.empty());
          assert(sequence.size() <= max_oligo_length);
          Lane lane;
          lane.oligo = oligo;
          lane.strand = Strand::plus;
          lane.pattern.resize(sequence.size());
          std::transform(sequence.cbegin(), sequence.cend(), lane.pattern.begin(),
                         [](char const nucleotide) -> char {
                           return static_cast<char>(maps::four_bit::map(nucleotide));
                         });
          assert(std::none_of(lane.pattern.cbegin(), lane.pattern.cend(),
                              [](char const code) -> bool { return code == 0; }));
          lane.free_overhang = free_overhang(static_cast<int>(sequence.size()),
                                             limits.target_cov);
          lanes_.push_back(lane);
          if (strands == Strands::plus_only) { continue; }
          lane.strand = Strand::minus;
          std::reverse(lane.pattern.begin(), lane.pattern.end());
          std::transform(lane.pattern.cbegin(), lane.pattern.cend(), lane.pattern.begin(),
                         [](char const code) -> char {
                           return static_cast<char>(complement(static_cast<unsigned char>(code)));
                         });
          lanes_.push_back(lane);
        }

      masks_.assign(code_count * lanes_.size(), 0U);
      for (std::size_t code = 0; code < code_count; ++code)
        {
          for (std::size_t lane = 0; lane < lanes_.size(); ++lane)
            {
              auto const & pattern = lanes_[lane].pattern;
              uint64_t matching = 0;
              for (std::size_t position = 0; position < pattern.size(); ++position)
                {
                  if (is_match(static_cast<unsigned char>(code),
                               static_cast<unsigned char>(pattern[position]),
                               n_positions))
                    {
                      matching |= uint64_t{1} << position;
                    }
                }
              masks_[(code * lanes_.size()) + lane] = matching;
            }
        }
    }


    Scanner::Scanner(LaneSet const & lane_set)
      : lane_set_(lane_set)
    {
      auto const & limits = lane_set.limits();
      threshold_ = std::min(limits.max_diffs, static_cast<int64_t>(max_oligo_length));
      auto const lane_count = lane_set.lanes().size();
      positive_.resize(lane_count);
      negative_.resize(lane_count);
      score_.resize(lane_count);
      high_bit_.resize(lane_count);
      open_candidate_.resize(lane_count);
      last_end_.resize(lane_count);
      for (std::size_t lane = 0; lane < lane_count; ++lane)
        {
          auto const length = lane_set.lanes()[lane].pattern.size();
          high_bit_[lane] = uint64_t{1} << (length - 1);
        }
      if (limits.model == Model::substitutions_only)
        {
          states_.resize(static_cast<std::size_t>(threshold_ + 1) * lane_count);
        }
    }


    /* An end position within the threshold. Neighbouring end positions of
       one lane belong to the same occurrence: an end closer than half the
       pattern length to the previous one joins it, and the occurrence keeps
       its lowest cost, then its leftmost end. Lanes are independent, so
       overlapping occurrences of different oligos are all kept. */
    auto Scanner::note(std::size_t const lane,
                       int64_t const end,
                       int64_t const cost,
                       int const rows) -> void
    {
      auto const length = static_cast<int64_t>(lane_set_.lanes()[lane].pattern.size());
      auto const same_site = std::max(int64_t{1}, length / 2);
      auto const open = open_candidate_[lane];
      if ((open != no_candidate) and (end - last_end_[lane] <= same_site))
        {
          auto & candidate = candidates_[static_cast<std::size_t>(open)];
          if (cost < candidate.cost)
            {
              candidate.end = end;
              candidate.cost = cost;
              candidate.rows = rows;
            }
        }
      else
        {
          open_candidate_[lane] = static_cast<int64_t>(candidates_.size());
          candidates_.push_back(Candidate{lane, end, cost, rows});
        }
      last_end_[lane] = end;
    }


    /* Myers's bit-vector scan, one word per lane. positive_ and negative_
       hold the vertical differences of the current column (bit i: from row i
       to row i + 1), score_ the cost of the whole pattern ending at the
       current read position. The first row is zero in every column (the read
       start is free anywhere). The first column lets the pattern skip up to
       free_overhang bases for nothing, then charges one per base: a pattern
       truncated by the read start. */
    auto Scanner::scan_edit_distance() -> void
    {
      auto const & lanes = lane_set_.lanes();
      auto const & masks = lane_set_.masks();
      auto const lane_count = lanes.size();
      auto const read_length = codes_.size();

      for (std::size_t lane = 0; lane < lane_count; ++lane)
        {
          auto const overhang = lanes[lane].free_overhang;
          positive_[lane] = ~uint64_t{0} << static_cast<unsigned int>(overhang);
          negative_[lane] = 0;
          score_[lane] = static_cast<int64_t>(lanes[lane].pattern.size()) - overhang;
        }

      for (std::size_t position = 0; position < read_length; ++position)
        {
          auto const row = static_cast<std::size_t>(codes_[position]) * lane_count;
          auto best = std::numeric_limits<int64_t>::max();
          /* no branch and no early exit in this loop: it has to vectorise */
          for (std::size_t lane = 0; lane < lane_count; ++lane)
            {
              auto const equal = masks[row + lane];
              auto const vertical_positive = positive_[lane];
              auto const vertical_negative = negative_[lane];
              auto const high = high_bit_[lane];
              auto const vertical_any = equal | vertical_negative;
              auto const horizontal_any =
                (((equal & vertical_positive) + vertical_positive) ^ vertical_positive) | equal;
              auto horizontal_positive = vertical_negative | ~(horizontal_any | vertical_positive);
              auto horizontal_negative = vertical_positive & horizontal_any;
              score_[lane] += static_cast<int64_t>((horizontal_positive & high) != 0U)
                - static_cast<int64_t>((horizontal_negative & high) != 0U);
              horizontal_positive <<= 1U;
              horizontal_negative <<= 1U;
              positive_[lane] = horizontal_negative | ~(vertical_any | horizontal_positive);
              negative_[lane] = horizontal_positive & vertical_any;
              best = std::min(best, score_[lane]);
            }
          if (best > threshold_) { continue; }  /* the common case, far from any hit */
          for (std::size_t lane = 0; lane < lane_count; ++lane)
            {
              if (score_[lane] > threshold_) { continue; }
              note(lane, static_cast<int64_t>(position), score_[lane],
                   static_cast<int>(lanes[lane].pattern.size()));
            }
        }

      /* A pattern truncated by the read end: in the last column, the cost of
         a pattern prefix of i bases is the sum of the first i vertical
         differences. The longest prefix wins a tie. */
      for (std::size_t lane = 0; lane < lane_count; ++lane)
        {
          auto const length = static_cast<int>(lanes[lane].pattern.size());
          auto const shortest = length - lanes[lane].free_overhang;
          auto cost = int64_t{0};
          auto best_cost = std::numeric_limits<int64_t>::max();
          auto best_rows = 0;
          for (auto rows = 1; rows < length; ++rows)
            {
              cost += static_cast<int64_t>(bit_is_set(positive_[lane], rows - 1))
                - static_cast<int64_t>(bit_is_set(negative_[lane], rows - 1));
              if ((rows >= shortest) and (cost <= best_cost))
                {
                  best_cost = cost;
                  best_rows = rows;
                }
            }
          if (best_cost > threshold_) { continue; }
          note(lane, static_cast<int64_t>(read_length) - 1, best_cost, best_rows);
        }
    }


    /* Substitutions only (--maxgaps 0): the shift-and automaton of Wu and
       Manber (1992) with one word per number of mismatches, states_[d][lane]
       bit i: the first i + 1 pattern bases match the read up to the current
       position with at most d mismatches. Starting with the lowest
       free_overhang bits set lets the pattern start that far into itself at
       the read start. */
    auto Scanner::scan_substitutions() -> void
    {
      auto const & lanes = lane_set_.lanes();
      auto const & masks = lane_set_.masks();
      auto const lane_count = lanes.size();
      auto const read_length = codes_.size();
      auto const levels = static_cast<std::size_t>(threshold_) + 1;
      auto const top = static_cast<std::size_t>(threshold_) * lane_count;

      for (std::size_t level = 0; level < levels; ++level)
        {
          for (std::size_t lane = 0; lane < lane_count; ++lane)
            {
              states_[(level * lane_count) + lane] = low_bits(lanes[lane].free_overhang);
            }
        }

      for (std::size_t position = 0; position < read_length; ++position)
        {
          auto const row = static_cast<std::size_t>(codes_[position]) * lane_count;
          /* from the highest level down, so that level - 1 still holds the
             previous position's state when level reads it */
          for (auto level = levels - 1; level > 0; --level)
            {
              auto const here = level * lane_count;
              auto const below = here - lane_count;
              for (std::size_t lane = 0; lane < lane_count; ++lane)
                {
                  states_[here + lane] =
                    (((states_[here + lane] << 1U) | 1U) & masks[row + lane])
                    | ((states_[below + lane] << 1U) | 1U);
                }
            }
          for (std::size_t lane = 0; lane < lane_count; ++lane)
            {
              states_[lane] = ((states_[lane] << 1U) | 1U) & masks[row + lane];
            }

          auto any = uint64_t{0};
          for (std::size_t lane = 0; lane < lane_count; ++lane)
            {
              any |= states_[top + lane] & high_bit_[lane];
            }
          if (any == 0U) { continue; }
          for (std::size_t lane = 0; lane < lane_count; ++lane)
            {
              if ((states_[top + lane] & high_bit_[lane]) == 0U) { continue; }
              auto level = std::size_t{0};
              while ((states_[(level * lane_count) + lane] & high_bit_[lane]) == 0U) { ++level; }
              note(lane, static_cast<int64_t>(position), static_cast<int64_t>(level),
                   static_cast<int>(lanes[lane].pattern.size()));
            }
        }

      /* a pattern truncated by the read end: prefixes, longest wins a tie */
      for (std::size_t lane = 0; lane < lane_count; ++lane)
        {
          auto const length = static_cast<int>(lanes[lane].pattern.size());
          auto const shortest = length - lanes[lane].free_overhang;
          auto best_cost = std::numeric_limits<int64_t>::max();
          auto best_rows = 0;
          for (auto rows = shortest; rows < length; ++rows)
            {
              for (std::size_t level = 0; level < levels; ++level)
                {
                  if (not bit_is_set(states_[(level * lane_count) + lane], rows - 1)) { continue; }
                  if (static_cast<int64_t>(level) <= best_cost)
                    {
                      best_cost = static_cast<int64_t>(level);
                      best_rows = rows;
                    }
                  break;
                }
            }
          if (best_cost > threshold_) { continue; }
          note(lane, static_cast<int64_t>(read_length) - 1, best_cost, best_rows);
        }
    }


    auto Scanner::align_substitutions(Candidate const & candidate,
                                      Occurrence & occurrence) -> bool
    {
      auto const & lane = lane_set_.lanes()[candidate.lane];
      auto const rows = static_cast<int64_t>(candidate.rows);
      auto const skipped = std::max(int64_t{0}, rows - (candidate.end + 1));
      if (skipped > lane.free_overhang) { return false; }
      occurrence.lane = candidate.lane;
      occurrence.read_end = candidate.end;
      occurrence.read_start = candidate.end - (rows - skipped) + 1;
      occurrence.pattern_start = static_cast<int>(skipped);
      occurrence.pattern_end = candidate.rows;
      occurrence.columns.assign(static_cast<std::size_t>(rows - skipped), 'M');
      occurrence.matches = 0;
      occurrence.mismatches = 0;
      occurrence.gap_columns = 0;
      occurrence.gap_openings = 0;
      for (auto position = skipped; position < rows; ++position)
        {
          auto const read_position = occurrence.read_start + (position - skipped);
          auto const code = codes_[static_cast<std::size_t>(read_position)];
          if (bit_is_set(lane_set_.mask(code, candidate.lane), static_cast<int>(position)))
            {
              ++occurrence.matches;
            }
          else
            {
              ++occurrence.mismatches;
            }
        }
      assert(occurrence.mismatches == candidate.cost);
      return true;
    }


    /* The window ends at the candidate's end position and is as wide as the
       widest alignment within the threshold. Rows are pattern bases, columns
       read bases. An alignment starts on the first row anywhere, or, when the
       window reaches the read start, on the first column after skipping up
       to free_overhang pattern bases. It ends on the last column: on the
       last row, or, at the read end, on any row that leaves at most
       free_overhang pattern bases unaligned. An alignment neither starts nor
       ends with a read base facing a gap (the read is free at both ends), and
       it does not start or end with a pattern base facing a gap where that
       base could hang off the read instead: those columns would be terminal
       gaps, not differences. */
    auto Scanner::align_edit_distance(Candidate const & candidate,
                                      Occurrence & occurrence) -> bool
    {
      auto const & lane = lane_set_.lanes()[candidate.lane];
      auto const length = static_cast<int64_t>(lane.pattern.size());
      auto const read_length = static_cast<int64_t>(codes_.size());
      auto const width = std::min(candidate.end + 1, length + threshold_);
      auto const window_start = candidate.end + 1 - width;
      auto const at_read_start = (window_start == 0);
      auto const at_read_end = (candidate.end == read_length - 1);
      auto const columns = static_cast<std::size_t>(width) + 1;
      auto const cells = (static_cast<std::size_t>(length) + 1) * columns;

      cell_costs_.assign(cells * move_count, unreachable);
      cell_moves_.assign(cells * move_count, static_cast<unsigned char>(Move::start));

      auto const cell = [columns](int64_t const row, int64_t const column, Move const move) -> std::size_t {
        return (((static_cast<std::size_t>(row) * columns) + static_cast<std::size_t>(column)) * move_count)
          + move_index(move);
      };
      /* an alignment may start here with a diagonal, or with a pattern gap */
      auto const starts_diagonal = [&](int64_t const row, int64_t const column) -> bool {
        return (row == 0) or ((column == 0) and at_read_start and (row <= lane.free_overhang));
      };
      auto const starts_pattern_gap = [&](int64_t const row, int64_t const column) -> bool {
        return (row == 0) and not ((column == 0) and at_read_start);
      };
      /* the best way into (row, column) by `move`, from the given candidates */
      auto const relax = [&](std::size_t const target, int64_t const cost, Move const from) -> void {
        if (cost < cell_costs_[target])
          {
            cell_costs_[target] = cost;
            cell_moves_[target] = static_cast<unsigned char>(from);
          }
      };
      auto const extend = [&](int64_t const row, int64_t const column,
                              Move const into, int64_t const step) -> void {
        /* into: the move entering (row, column); its predecessor cell */
        auto const previous_row = (into == Move::read_gap) ? row : row - 1;
        auto const previous_column = (into == Move::pattern_gap) ? column : column - 1;
        auto const target = cell(row, column, into);
        auto const opening = (into == Move::diagonal) ? int64_t{0} : opening_weight;
        for (auto const from : {Move::diagonal, Move::read_gap, Move::pattern_gap})
          {
            auto const previous = cell_costs_[cell(previous_row, previous_column, from)];
            if (previous >= unreachable) { continue; }
            auto const reopens = ((into != Move::diagonal) and (from != into)) ? opening : int64_t{0};
            relax(target, previous + step + reopens, from);
          }
        auto const may_start = (into == Move::diagonal) ? starts_diagonal(previous_row, previous_column)
          : (into == Move::pattern_gap) ? starts_pattern_gap(previous_row, previous_column)
          : false;
        if (may_start) { relax(target, step + opening, Move::start); }
      };

      for (int64_t row = 1; row <= length; ++row)
        {
          for (int64_t column = 0; column <= width; ++column)
            {
              if (column > 0)
                {
                  auto const code = codes_[static_cast<std::size_t>(window_start + column - 1)];
                  auto const equal = bit_is_set(lane_set_.mask(code, candidate.lane),
                                                static_cast<int>(row - 1));
                  extend(row, column, Move::diagonal, equal ? int64_t{0} : difference_weight);
                  extend(row, column, Move::read_gap, difference_weight);
                }
              extend(row, column, Move::pattern_gap, difference_weight);
            }
        }

      /* the end: the lowest cost, then the longest pattern prefix */
      auto best_cost = unreachable;
      auto best_row = int64_t{0};
      auto best_move = Move::diagonal;
      auto const shortest = at_read_end ? length - lane.free_overhang : length;
      for (auto row = shortest; row <= length; ++row)
        {
          auto const diagonal = cell_costs_[cell(row, width, Move::diagonal)];
          if (diagonal <= best_cost)
            {
              best_cost = diagonal;
              best_row = row;
              best_move = Move::diagonal;
            }
        }
      if (not at_read_end)
        {
          auto const pattern_gap = cell_costs_[cell(length, width, Move::pattern_gap)];
          if (pattern_gap < best_cost)
            {
              best_cost = pattern_gap;
              best_row = length;
              best_move = Move::pattern_gap;
            }
        }
      if (best_cost >= unreachable) { return false; }

      /* trace back to the start */
      occurrence.columns.clear();
      auto row = best_row;
      auto column = width;
      auto move = best_move;
      while (true)
        {
          auto const from = static_cast<Move>(cell_moves_[cell(row, column, move)]);
          switch (move)
            {
            case Move::diagonal:
              occurrence.columns.push_back('M');
              --row;
              --column;
              break;
            case Move::read_gap:
              occurrence.columns.push_back('D');
              --column;
              break;
            case Move::pattern_gap:
              occurrence.columns.push_back('I');
              --row;
              break;
            case Move::start:
            default:
              assert(false);
              return false;
            }
          if (from == Move::start) { break; }
          move = from;
        }
      std::reverse(occurrence.columns.begin(), occurrence.columns.end());

      occurrence.lane = candidate.lane;
      occurrence.read_start = window_start + column;
      occurrence.read_end = candidate.end;
      occurrence.pattern_start = static_cast<int>(row);
      occurrence.pattern_end = static_cast<int>(best_row);
      occurrence.matches = 0;
      occurrence.mismatches = 0;
      occurrence.gap_columns = 0;
      occurrence.gap_openings = 0;
      auto read_position = occurrence.read_start;
      auto pattern_position = row;
      auto previous = 'M';
      for (auto const operation : occurrence.columns)
        {
          if (operation == 'M')
            {
              auto const code = codes_[static_cast<std::size_t>(read_position)];
              if (bit_is_set(lane_set_.mask(code, candidate.lane),
                             static_cast<int>(pattern_position)))
                {
                  ++occurrence.matches;
                }
              else
                {
                  ++occurrence.mismatches;
                }
              ++read_position;
              ++pattern_position;
            }
          else
            {
              ++occurrence.gap_columns;
              if (operation != previous) { ++occurrence.gap_openings; }
              if (operation == 'D') { ++read_position; }
              else { ++pattern_position; }
            }
          previous = operation;
        }
      assert(read_position == occurrence.read_end + 1);
      assert(pattern_position == occurrence.pattern_end);
      assert(((occurrence.mismatches + occurrence.gap_columns) * difference_weight)
             + occurrence.gap_openings == best_cost);
      return true;
    }


    /* --target_cov bounds truncation, not differences: it is compared with
       the span of the oligo inside the read, deleted bases included, so
       that --target_cov 1 excludes truncated occurrences and nothing else.
       (Elsewhere in vsearch, and in the tcov field, the aligned fraction
       counts letter pairs only, and a deletion lowers it.) */
    auto Scanner::passes_limits(Occurrence const & occurrence) const noexcept -> bool
    {
      auto const & limits = lane_set_.limits();
      auto const length = lane_set_.lanes()[occurrence.lane].pattern.size();
      auto const pairs = occurrence.matches + occurrence.mismatches;
      auto const span = occurrence.pattern_end - occurrence.pattern_start;
      return (pairs > 0) and
        (occurrence.mismatches + occurrence.gap_columns <= limits.max_diffs) and
        (occurrence.gap_openings <= limits.max_gap_openings) and
        (span >= limits.target_cov * static_cast<double>(length));
    }


    auto Scanner::verify(Candidate const & candidate, Occurrence & occurrence) -> bool
    {
      auto const aligned = (lane_set_.limits().model == Model::substitutions_only) ?
        align_substitutions(candidate, occurrence) :
        align_edit_distance(candidate, occurrence);
      return aligned and passes_limits(occurrence);
    }


    auto Scanner::scan(View<char> const read, std::vector<Occurrence> & found) -> void
    {
      found.clear();
      candidates_.clear();
      if (read.empty()) { return; }

      codes_.resize(read.size());
      std::transform(read.cbegin(), read.cend(), codes_.begin(),
                     [](char const nucleotide) -> unsigned char {
                       return maps::four_bit::map(nucleotide);
                     });
      std::fill(open_candidate_.begin(), open_candidate_.end(), no_candidate);

      if (lane_set_.limits().model == Model::substitutions_only)
        {
          scan_substitutions();
        }
      else
        {
          scan_edit_distance();
        }

      Occurrence occurrence;
      for (auto const & candidate : candidates_)
        {
          if (not verify(candidate, occurrence)) { continue; }
          found.push_back(occurrence);
        }

      std::sort(found.begin(), found.end(),
                [](Occurrence const & lhs, Occurrence const & rhs) -> bool {
                  return std::tie(lhs.read_start, lhs.read_end, lhs.lane)
                    < std::tie(rhs.read_start, rhs.read_end, rhs.lane);
                });
    }

  }  // namespace oligo
}  // namespace vsearch
