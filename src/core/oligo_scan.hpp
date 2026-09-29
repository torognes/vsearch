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

#pragma once

#include "utils/view.hpp"  // View
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <limits>  // std::numeric_limits
#include <string>
#include <vector>


/* Approximate search of short oligos (primers, tags, barcodes) in long
   sequences: every occurrence of every oligo within a number of differences,
   on one or both strands, including occurrences truncated by a read end.

   The scan is Myers's bit-vector algorithm (1999), semi-global: the oligo is
   aligned end to end, the read is free at both ends. One 64-bit word holds
   the state of one oligo, which bounds oligos to 64 nucleotides. Oligos are
   scanned as "lanes" (one lane per oligo and strand, the minus strand being
   the reverse complement of the oligo scanned along the read), and the loop
   over lanes sits inside the loop over read positions: independent lanes
   give the CPU instruction-level parallelism, and the loop auto-vectorises
   where the instruction set has 64-bit vector compares.

   The scan only yields end positions and costs. Each occurrence is then
   aligned again by a small dynamic programming over a window of the read,
   which recovers its start, its columns, and the fewest gap openings among
   the alignments of lowest cost.

   Nucleotides are compared by overlap of their IUPAC base sets, as every
   vsearch aligner does (maps::four_bit::is_equivalent): an N in the read
   matches everything, a Y in the read matches a C, a T or a Y in the oligo.
   With --n_mismatch, a position involving an N is a mismatch. */

namespace vsearch {
  namespace oligo {

    /* one 64-bit word per lane */
    constexpr std::size_t max_oligo_length = 64;

    enum struct Strand : unsigned char { plus, minus };

    enum struct Strands : unsigned char { plus_only, both };

    /* edit_distance: substitutions, insertions and deletions;
       substitutions_only: --maxgaps 0, mismatches alone (Hamming distance) */
    enum struct Model : unsigned char { edit_distance, substitutions_only };

    /* how a position facing an N is counted (--n_mismatch) */
    enum struct NPositions : unsigned char { match, mismatch };


    /* What an occurrence must satisfy to be reported. */
    struct Limits
    {
      int64_t max_diffs = 2;             /* --maxdiffs: mismatches + gap columns */
      int64_t max_gap_openings = std::numeric_limits<int64_t>::max();  /* --maxgaps */
      double target_cov = 0.0;           /* --target_cov: oligo span in the read over its length */
      Model model = Model::edit_distance;
    };


    /* One searched pattern: an oligo, or its reverse complement. */
    struct Lane
    {
      std::size_t oligo = 0;       /* index of the oligo in the database */
      Strand strand = Strand::plus;
      std::string pattern;         /* as scanned along the read, 4-bit codes */
      int free_overhang = 0;       /* pattern bases allowed off a read end */
    };


    /* The lanes of a set of oligos, and their match masks. Built once, then
       shared read-only by every thread. */
    class LaneSet
    {
    public:
      LaneSet(std::vector<View<char>> const & oligos,
              Strands strands,
              NPositions n_positions,
              Limits const & limits);

      auto lanes() const noexcept -> std::vector<Lane> const & { return lanes_; }
      auto limits() const noexcept -> Limits const & { return limits_; }

      /* bit i is set when pattern position i matches read symbol code */
      auto mask(unsigned char const code, std::size_t const lane) const noexcept -> uint64_t
      {
        return masks_[(static_cast<std::size_t>(code) * lanes_.size()) + lane];
      }

      auto masks() const noexcept -> std::vector<uint64_t> const & { return masks_; }

    private:
      std::vector<Lane> lanes_;
      std::vector<uint64_t> masks_;  /* [code][lane], 16 codes */
      Limits limits_;
    };


    /* One reported occurrence, in the lane's orientation, i.e. read
       coordinates on the plus strand and pattern coordinates in the lane's
       pattern (the reverse complement of the oligo on a minus lane). */
    struct Occurrence
    {
      std::size_t lane = 0;
      int64_t read_start = 0;     /* first read base aligned, 0-based */
      int64_t read_end = 0;       /* last read base aligned, 0-based, inclusive */
      int pattern_start = 0;      /* first pattern base aligned, 0-based */
      int pattern_end = 0;        /* one past the last pattern base aligned */
      /* one letter per column, as in vsearch's cigar strings: M (two
         nucleotides), D (a read base facing a gap), I (a pattern base facing
         a gap) */
      std::string columns;
      int matches = 0;
      int mismatches = 0;
      int gap_columns = 0;
      int gap_openings = 0;
    };


    /* The per-thread scanner: its state is reused from read to read. */
    class Scanner
    {
    public:
      explicit Scanner(LaneSet const & lane_set);

      /* Every occurrence in the read that passes the limits, ordered by read
         position (read_start, then read_end, then lane). */
      auto scan(View<char> read, std::vector<Occurrence> & found) -> void;

    private:
      struct Candidate
      {
        std::size_t lane;
        int64_t end;       /* last read position of the occurrence */
        int64_t cost;
        int rows;          /* pattern bases up to the end (< length: truncated) */
      };

      auto note(std::size_t lane, int64_t end, int64_t cost, int rows) -> void;
      auto scan_edit_distance() -> void;
      auto scan_substitutions() -> void;
      auto note_truncated_ends_edit_distance() -> void;
      auto note_truncated_ends_substitutions() -> void;
      auto verify(Candidate const & candidate, Occurrence & occurrence) -> bool;
      auto align_edit_distance(Candidate const & candidate, Occurrence & occurrence) -> bool;
      auto align_substitutions(Candidate const & candidate, Occurrence & occurrence) -> bool;
      auto passes_limits(Occurrence const & occurrence) const noexcept -> bool;

      LaneSet const & lane_set_;
      int64_t threshold_ = 0;      /* max_diffs, capped at the longest possible cost */
      std::vector<unsigned char> codes_;  /* the read, as 4-bit codes */

      /* bit-vector state, one element per lane */
      std::vector<uint64_t> positive_;
      std::vector<uint64_t> negative_;
      std::vector<int64_t> score_;
      std::vector<uint64_t> high_bit_;
      /* substitutions_only: (threshold + 1) words per lane, [cost][lane] */
      std::vector<uint64_t> states_;

      /* collapsing of neighbouring end positions, per lane */
      std::vector<Candidate> candidates_;
      std::vector<int64_t> open_candidate_;  /* index in candidates_, or -1 */
      std::vector<int64_t> last_end_;

      /* alignment scratch */
      std::vector<int64_t> cell_costs_;
      std::vector<unsigned char> cell_moves_;
    };

  }  // namespace oligo
}  // namespace vsearch
