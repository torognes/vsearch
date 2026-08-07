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

#include "utils/view.hpp"
#include "vsearch.hpp"
#include "core/msa.hpp"
#include "core/db.hpp"
#include "core/fasta.hpp"
#include "utils/ascii_case.hpp"  // to_upper
#include "utils/cigar.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/span.hpp"
#include "utils/reverse_complement.hpp"
#include <array>
#include <algorithm>  // std::fill, std::fill_n, std::max, std::max_element
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE
#include <iterator> // std::distance, std::next
#include <numeric> // std::accumulate
#include <vector>

#ifndef NDEBUG
#include <limits>
#endif


/* Compute multiple sequence alignment (msa), profile, and consensus
   sequence of clustered sequences */

using prof_type = std::uint64_t;

/* One profile column per alignment position: the six abundance counters below,
   in storage order. The four nucleotide counters come first and in the order
   sym_nt_4bit encodes them (1 << A_counter == 'A' and so on), which is what
   lets the consensus read the winner's symbol straight out of its position. */
constexpr auto A_counter = 0;
constexpr auto C_counter = 1;
constexpr auto G_counter = 2;
constexpr auto U_counter = 3;  // note: T converted to U?
constexpr auto N_counter = 4;
constexpr auto gap_counter = 5;
constexpr auto profsize = 6;  // the six counters above


namespace {
/* The counters of alignment position `index`. The profile is a flat vector of
   profsize-wide columns, and every reader used to spell that stride by hand;
   the subspan also carries the bounds assertions an open-coded
   profsize * index + counter cannot. */
auto column(std::vector<prof_type> & profile, int const index) -> Span<prof_type> {
  return make_span(profile).subspan(static_cast<std::size_t>(index)
                                    * static_cast<std::size_t>(profsize),
                                    static_cast<std::size_t>(profsize));
}

auto column(std::vector<prof_type> const & profile, int const index) -> View<prof_type> {
  return make_view(profile).subspan(static_cast<std::size_t>(index)
                                    * static_cast<std::size_t>(profsize),
                                    static_cast<std::size_t>(profsize));
}


auto update_profile(char const nucleotide,
                    int const position_in_alignment,
                    prof_type const abundance,
                    std::vector<prof_type>& profile) -> void {
  auto const counters = column(profile, position_in_alignment);

  // refactoring: eliminate unused cases? No, T and U are merged, same as IUPAC and N
  switch (to_upper(nucleotide))
    {
    case 'A':
      counters[A_counter] += abundance;
      break;
    case 'C':
      counters[C_counter] += abundance;
      break;
    case 'G':
      counters[G_counter] += abundance;
      break;
    case 'T':
    case 'U':
      counters[U_counter] += abundance;
      break;
    case 'R':
    case 'Y':
    case 'S':
    case 'W':
    case 'K':
    case 'M':
    case 'B':
    case 'D':
    case 'H':
    case 'V':
    case 'N':
      counters[N_counter] += abundance;
      break;
    case '-':
      counters[gap_counter] += abundance;
      break;
    default:
      break;
    }
}


auto update_msa(char const nucleotide, int &position_in_alignment,
                std::vector<char>& alignment) -> void {
  alignment[static_cast<std::vector<char>::size_type>(position_in_alignment)] = nucleotide;
  ++position_in_alignment;
}
}  // anonymous namespace


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

}  // end of anonymous namespace


namespace {
auto find_max_insertions_per_position(View<struct msa_target_s> const targets,
                                      int const centroid_len) -> std::vector<int> {
  std::vector<int> max_insertions(static_cast<std::vector<int>::size_type>(centroid_len + 1));

  // the centroid (targets.front()) carries no cigar string, hence drop(1)
  for (auto const & a_msa_target : targets.drop(1)) {
    auto position = 0LL;
    auto const cigar_pairs = parse_cigar_string(a_msa_target.cigar);

    for (auto const & a_pair: cigar_pairs) {
      auto const operation = a_pair.first;
      auto const runlength = a_pair.second;
      switch (operation) {
      case Operation::match:
      case Operation::insertion:
        position += runlength;
        break;

      case Operation::deletion:
        assert(runlength <= std::numeric_limits<int>::max());
        max_insertions[static_cast<std::vector<int>::size_type>(position)] = std::max(static_cast<int>(runlength), max_insertions[static_cast<std::vector<int>::size_type>(position)]);
        break;
      }
    }
  }
  return max_insertions;
}


auto find_total_alignment_length(std::vector<int> const & max_insertions) -> int64_t {
  /* Sum in 64-bit: the total alignment length is centroid length plus the max
     insertions at every position, which can exceed INT_MAX for a large cluster
     even though each individual sequence is bounded. A 32-bit accumulator here
     wraps to a wrong (often negative) length that then under-sizes the profile
     and alignment buffers while the alignment walk uses the true extent -- a
     heap overflow. With a 64-bit result the buffers are sized correctly; a true
     length beyond int capacity fails the allocation cleanly before any walk. */
  auto const centroid_len = static_cast<int64_t>(max_insertions.size() - 1);
  return std::accumulate(max_insertions.begin(), max_insertions.end(), centroid_len);
}


auto find_longest_target_on_reverse_strand(View<struct msa_target_s> const targets,
                                           struct Database const & db) -> int64_t {
  int64_t longest_reversed = 0;
  for (auto const & target : targets)
    {
      if (target.strand == 0) { continue; }
      auto const len = static_cast<int64_t>(db.getsequencelen(static_cast<uint64_t>(target.seqno)));
      longest_reversed = std::max(len, longest_reversed);
    }
  return longest_reversed;
}


auto allocate_buffer_for_reverse_strand_target(View<struct msa_target_s> const targets,
                                               std::vector<char> & rc_buffer_v,
                                               struct Database const & db) -> Span<char> {
  /* Find longest target sequence on reverse strand and allocate buffer */
  auto const longest_reversed = find_longest_target_on_reverse_strand(targets, db);
  if (longest_reversed > 0)
    {
      rc_buffer_v.resize(static_cast<std::vector<char>::size_type>(longest_reversed + 1));
      return make_span(rc_buffer_v);
    }
  return {};
}


auto blank_line_before_each_msa(std::FILE * fp_msaout) -> void {
  if (fp_msaout == nullptr) { return ; }
  fprint(fp_msaout, '\n');
}


auto print_header_and_sequence(std::FILE * fp_msaout, char const * header_prefix,
                               int const target_seqno,
                               std::vector<char> const & aln_v,
                               struct Database const & db,
                               struct Parameters const & parameters) -> void {
  // header_prefix == "*" or "", resulting in ">*header" or ">header"
  if (fp_msaout == nullptr) { return ; }

  fasta_print_general(fp_msaout,
                      header_prefix,
                      make_view(aln_v).first(aln_v.size() - 1),
                      db.header_view(static_cast<uint64_t>(target_seqno)),
                      db.getabundance(static_cast<uint64_t>(target_seqno)),
                      0, -1.0, -1, -1, nullptr, 0.0, 0, parameters);
}


auto reverse_complement_target_if_need_be(int const strand, Span<char> const rc_buffer,
                                          View<char> const target_seq) -> View<char> {
  if (strand == 0) { return target_seq; }
  reverse_complement(rc_buffer.first(target_seq.size() + 1), target_seq);
  return View<char>{rc_buffer.first(target_seq.size())};
}


auto process_and_print_centroid(Span<char> const rc_buffer,
                                View<struct msa_target_s> const targets,
                                std::vector<int> const &max_insertions,
                                std::vector<prof_type> &profile,
                                std::vector<char> &aln_v,
                                std::FILE * fp_msaout,
                                struct Database const & db,
                                struct Parameters const & parameters) -> void {
  auto const centroid_len = static_cast<int>(max_insertions.size() - 1);
  auto const & target = targets.front();
  auto const target_seqno = target.seqno;
  auto const target_seq = reverse_complement_target_if_need_be(target.strand, rc_buffer,
                                                               db.sequence_view(static_cast<uint64_t>(target_seqno)));
  prof_type const target_abundance = parameters.opt_sizein ? db.getabundance(static_cast<uint64_t>(target_seqno)) : 1;
  auto position_in_alignment = 0;

  for (auto i = 0; i < centroid_len; ++i)
    {
      for (auto j = 0; j < max_insertions[static_cast<std::vector<int>::size_type>(i)]; ++j)
        {
          update_profile('-', position_in_alignment, target_abundance, profile);
          update_msa('-', position_in_alignment, aln_v);
        }
      update_profile(target_seq[static_cast<std::size_t>(i)], position_in_alignment, target_abundance, profile);
      update_msa(target_seq[static_cast<std::size_t>(i)], position_in_alignment, aln_v);
    }

  // insert
  for (auto j = 0; j < max_insertions[static_cast<std::vector<int>::size_type>(centroid_len)]; ++j)
    {
      update_profile('-', position_in_alignment, target_abundance, profile);
      update_msa('-', position_in_alignment, aln_v);
    }

  /* end of sequence string */
  aln_v[static_cast<std::vector<char>::size_type>(position_in_alignment)] = '\0';

  /* print header & sequence */
  print_header_and_sequence(fp_msaout, "*", target_seqno, aln_v, db, parameters);
}


auto insert_gaps_in_alignment_and_profile(bool const is_inserted,
                                          int const max_insertions_at_position,
                                          int & position_in_alignment,
                                          prof_type const target_abundance,
                                          std::vector<prof_type> & profile,
                                          std::vector<char> & aln_v) -> void {
  if (is_inserted) { return ; }
  for (auto i = 0; i < max_insertions_at_position; ++i) {
    update_profile('-', position_in_alignment, target_abundance, profile);
    update_msa('-', position_in_alignment, aln_v);
  }
}


auto compute_and_print_msa(View<struct msa_target_s> const targets,
                           std::vector<int> const &max_insertions,
                           std::vector<prof_type> &profile,
                           std::vector<char> &aln_v,
                           std::FILE * fp_msaout,
                           struct Database const & db,
                           struct Parameters const & parameters) -> void {

  blank_line_before_each_msa(fp_msaout);

  /* Find longest target sequence on reverse strand and allocate buffer */
  std::vector<char> rc_buffer_v;
  auto const rc_buffer = allocate_buffer_for_reverse_strand_target(targets, rc_buffer_v, db);

  // ------------------------------------------------------- deal with centroid
  process_and_print_centroid(rc_buffer, targets, max_insertions,
                             profile, aln_v, fp_msaout, db, parameters);

  // --------------------------------- deal with other sequences in the cluster
  for (auto const & target : targets.drop(1))
    {
      auto const target_seqno = target.seqno;
      auto const target_seq = reverse_complement_target_if_need_be(target.strand, rc_buffer,
                                                                   db.sequence_view(static_cast<uint64_t>(target_seqno)));
      prof_type const target_abundance = parameters.opt_sizein ? db.getabundance(static_cast<uint64_t>(target_seqno)) : 1;
      int position_in_alignment = 0;

      auto is_inserted = false;
      auto qpos = 0;
      auto tpos = 0;

      auto const * const cigar_end = target.cigar.end();
      auto const * position_in_cigar = target.cigar.begin();
      while (position_in_cigar < cigar_end)
        {
          // Consume digits (if any), return the position of the
          // first char (M, D, or I), store it, move cursor to the next byte.
          // Operations: match (M), insertion (I), or deletion (D)
          auto ** next_operation = &position_in_cigar;
          auto const runlength = find_runlength_of_leftmost_operation(position_in_cigar, next_operation);
          auto const operation = **next_operation;
          position_in_cigar = std::next(position_in_cigar);

          switch (operation) {
          case 'D':
            for (auto j = 0; j < runlength; ++j)
              {
                update_profile(target_seq[static_cast<std::size_t>(tpos)], position_in_alignment, target_abundance, profile);
                update_msa(target_seq[static_cast<std::size_t>(tpos)], position_in_alignment, aln_v);
                ++tpos;
              }
            for (auto j = runlength; j < max_insertions[static_cast<std::vector<int>::size_type>(qpos)]; ++j)
              {
                update_profile('-', position_in_alignment, target_abundance, profile);
                update_msa('-', position_in_alignment, aln_v);
              }
            is_inserted = true;
            break;
          case 'M':
            for (auto j = 0; j < runlength; ++j)
              {
                insert_gaps_in_alignment_and_profile(is_inserted, max_insertions[static_cast<std::vector<int>::size_type>(qpos)],
                                                     position_in_alignment, target_abundance,
                                                     profile, aln_v);
                update_profile(target_seq[static_cast<std::size_t>(tpos)], position_in_alignment, target_abundance, profile);
                update_msa(target_seq[static_cast<std::size_t>(tpos)], position_in_alignment, aln_v);
                ++tpos;
                ++qpos;
                is_inserted = false;
              }
            break;
          case 'I':
            for (auto j = 0; j < runlength; ++j)
              {
                insert_gaps_in_alignment_and_profile(is_inserted, max_insertions[static_cast<std::vector<int>::size_type>(qpos)],
                                                     position_in_alignment, target_abundance,
                                                     profile, aln_v);
                update_profile('-', position_in_alignment, target_abundance, profile);
                update_msa('-', position_in_alignment, aln_v);
                ++qpos;
                is_inserted = false;
              }
              break;
          default:
            break;
          }
        }

      insert_gaps_in_alignment_and_profile(is_inserted, max_insertions[static_cast<std::vector<int>::size_type>(qpos)],
                                           position_in_alignment, target_abundance,
                                           profile, aln_v);

      /* end of sequence string */
      aln_v[static_cast<std::vector<char>::size_type>(position_in_alignment)] = '\0';

      /* print header & sequence */
      print_header_and_sequence(fp_msaout, "", target_seqno, aln_v, db, parameters);
    }
}


auto compute_and_print_consensus(std::vector<int> const &max_insertions,
                                 std::vector<char> &aln_v,
                                 std::vector<char> &cons_v,
                                 std::vector<prof_type> const &profile,
                                 std::FILE * fp_msaout,
                                 struct Parameters const & parameters) -> void {
  static constexpr std::array<char, 16> sym_nt_4bit = {{'-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}};
  static constexpr char index_of_N = 15;  // 15th char in sym_nt_4bit[] (=> 'N')

  auto const alignment_length = static_cast<int>(aln_v.size() - 1);
  int conslen = 0;

  /* Censor part of the consensus sequence outside the centroid sequence */
  auto const left_censored = max_insertions.front();
  auto const right_censored = max_insertions.back();
  assert(left_censored >= 0);
  assert(right_censored >= 0);
  /* aln_v ends in a '\0' past the alignment, so the censoring and the
     consensus below run over its first alignment_length characters */
  auto const alignment = make_span(aln_v).first(static_cast<std::size_t>(alignment_length));
  std::fill_n(alignment.begin(), left_censored, '+');
  auto const right_censored_tail = alignment.last(static_cast<std::size_t>(right_censored));
  std::fill(right_censored_tail.begin(), right_censored_tail.end(), '+');

  for (auto i = left_censored; i < alignment_length - right_censored; ++i)
    {
      auto const counters = column(profile, i);

      /* find most common symbol of A, C, G and T. Strictly-greater kept the
         first of a tie, which is what max_element returns; a column with no
         nucleotide at all leaves best_sym at 0, so the count is tested. */
      auto const * const most_common =
        std::max_element(counters.begin(), std::next(counters.begin(), U_counter + 1));
      char best_sym = 0;
      auto best_count = *most_common;
      if (best_count > 0)
        {
          // 1, 2, 4, or 8 -- A_counter..U_counter are sym_nt_4bit's bit positions
          best_sym = static_cast<char>(
            1U << static_cast<unsigned int>(std::distance(counters.begin(), most_common)));
        }

      /* if no A, C, G, or T, check if there are any N's */
      auto const N_count = counters[N_counter];
      if ((best_count == 0) and (N_count > 0))
        {
          best_count = N_count;
          best_sym = index_of_N; // N
        }

      /* compare to the number of gap symbols */
      auto const gap_count = counters[gap_counter];
      if (best_count >= gap_count)
        {
          auto const index = static_cast<unsigned char>(best_sym);
          auto const sym = sym_nt_4bit[index];  // A, C, G, T, or N
          alignment[static_cast<std::size_t>(i)] = sym;
          cons_v[static_cast<std::vector<char>::size_type>(conslen)] = sym;
          ++conslen;
        }
      else
        {
          alignment[static_cast<std::size_t>(i)] = '-';
        }
    }

  aln_v.back() = '\0';
  cons_v[static_cast<std::vector<char>::size_type>(conslen)] = '\0';
  cons_v.resize(static_cast<std::vector<char>::size_type>(conslen + 1));

  if (fp_msaout != nullptr)
    {
      fasta_print(fp_msaout, "consensus", aln_v.data(), static_cast<uint64_t>(alignment_length), parameters);
    }
}


auto print_consensus_sequence(std::FILE *fp_consout, std::vector<char> const & cons_v,
                              int64_t const totalabundance, int const target_count,
                              int const cluster,
                              int const centroid_seqno,
                              struct Database const & db,
                              struct Parameters const & parameters) -> void {
  if (fp_consout == nullptr) { return ; }
  fasta_print_general(fp_consout,
                      "centroid=",
                      make_view(cons_v).first(cons_v.size() - 1),  // exclude the '\0' terminator slot
                      db.header_view(static_cast<uint64_t>(centroid_seqno)),
                      static_cast<uint64_t>(totalabundance),
                      cluster + 1,
                      -1.0,
                      target_count,
                      parameters.opt_clusterout_id ? cluster : -1,
                      nullptr, 0.0,
                      0,
                      parameters);
}


auto print_alignment_profile(std::FILE *fp_profile, std::vector<char> &aln_v,
                             std::vector<prof_type> const &profile,
                             int64_t const totalabundance, int const target_count,
                             int const cluster,
                             int const centroid_seqno,
                             struct Database const & db,
                             struct Parameters const & parameters) -> void {
  if (fp_profile == nullptr) { return ; }

  // Note: gaps before Ns in profile output
  // 0 = A, 1 = C, 2 = G, 3 = T, 4 = N, 5 = '-' (gap)
  static const std::array<int, 6> symbol_indexes = {0, 1, 2, 3, 5, 4};
  fasta_print_general(fp_profile,
                      "centroid=",
                      View<char>{},  // the profile output carries no centroid sequence
                      db.header_view(static_cast<uint64_t>(centroid_seqno)),
                      static_cast<uint64_t>(totalabundance),
                      cluster + 1,
                      -1.0,
                      target_count,
                      parameters.opt_clusterout_id ? cluster : -1,
                      nullptr, 0.0,
                      0,
                      parameters);

  aln_v.pop_back(); // remove last element ('\0')
  auto counter = 0;
  for (auto const nucleotide: aln_v) {
    fprint_integer(fp_profile, counter);
    fprint(fp_profile, '\t');
    fprint(fp_profile, nucleotide);
      // A, C, G and T, then gap '-', then N
      for (auto const symbol_index : symbol_indexes) {
        fprint(fp_profile, '\t');
        fprint_integer(fp_profile, column(profile, counter)[static_cast<std::size_t>(symbol_index)]);
      }
      fprint(fp_profile, '\n');
      ++counter;
    }
  fprint(fp_profile, '\n');
}
}  // anonymous namespace


auto msa(std::FILE * fp_msaout, std::FILE * fp_consout, std::FILE * fp_profile,
         int const cluster,
         View<struct msa_target_s> const targets,
         int64_t const totalabundance,
         struct Database const & db,
         struct Parameters const & parameters) -> void
{
  assert(not targets.empty());  // a cluster always holds at least its centroid
  assert(targets.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()));
  auto const target_count = static_cast<int>(targets.size());  // reported as 'seqs=' below
  int const centroid_seqno = targets.front().seqno;
  auto const centroid_length = static_cast<int>(db.getsequencelen(static_cast<uint64_t>(centroid_seqno)));

  /* find max insertions in front of each position in the centroid sequence */
  auto const max_insertions = find_max_insertions_per_position(targets, centroid_length);
  auto const alignment_length = find_total_alignment_length(max_insertions);

  /* allocate memory for profile (for consensus) and aligned seq */
  /* one profsize-wide column per alignment position, flattened; read through
     column() above, which is what gives the call sites the shape the brief
     below asks for. The brief stands for the storage itself: a
     vector<array<prof_type, profsize>> would make the stride a type-level fact
     and let the compiler see the column width, which a subspan cannot. */
  std::vector<prof_type> profile(static_cast<unsigned long>(profsize) * static_cast<unsigned long>(alignment_length));  // C++20 refactoring: std::vector<std::array<prof_type, profsize>>(alnlen);
  std::vector<char> aln_v(static_cast<std::vector<char>::size_type>(alignment_length + 1));
  std::vector<char> cons_v(static_cast<std::vector<char>::size_type>(alignment_length + 1));

  /* msaout: multiple sequence alignment ... */
  compute_and_print_msa(targets, max_insertions,
                        profile, aln_v,
                        fp_msaout, db, parameters);

  /* msaout: ... and consensus sequence at the end */
  compute_and_print_consensus(max_insertions,
                              aln_v,
                              cons_v,
                              profile,
                              fp_msaout,
                              parameters);

  /* consout: consensus sequence (dedicated input) */
  print_consensus_sequence(fp_consout, cons_v,
                           totalabundance, target_count,
                           cluster,
                           centroid_seqno, db, parameters);

  /* profile: multiple sequence alignment profile (dedicated input) */
  print_alignment_profile(fp_profile, aln_v,
                          profile,
                          totalabundance, target_count,
                          cluster,
                          centroid_seqno, db, parameters);
}
