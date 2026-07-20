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

/* Internal interface shared between the fastq_mergepairs command
   (commands/fastq_mergepairs.cpp) and the paired-end merge engine
   (core/mergepairs.cpp): the per-pair working struct, the merge
   status/reason enums, the cooperative-abort control surface, and the two
   engine entry points the command drives directly (precompute_qual, process).
   Not part of the public library API (that is core/mergepairs.hpp). */

#include "core/mergepairs.hpp"  // QualityTables (shared with the merge engine)
#include <atomic>  // std::memory_order
#include <cstdint>  // int64_t
#include <vector>

enum class MergeAbortReason : std::uint8_t { quality_below_qmin, quality_above_qmax, more_fwd_than_rev };


/* reasons for not merging:
   - undefined
   - ok
   - input seq too short (after truncation)
   - input seq too long
   - too many Ns in input
   - overlap too short
   - too many differences (maxdiffs)
   - too high percentage of differences (maxdiffpct)
   - staggered
   - indels in overlap region
   - potential repeats in overlap region / multiple overlaps
   - merged sequence too short
   - merged sequence too long
   - expected error too high
   - alignment score too low, insignificant, potential indel
   - too few kmers on same diag found
*/

enum struct Reason : char {
  undefined,
  ok,
  minlen,
  maxlen,
  maxns,
  minovlen,
  maxdiffs,
  maxdiffpct,
  staggered,
  indel,
  repeat,
  minmergelen,
  maxmergelen,
  maxee,
  minscore,
  nokmers,
};

enum struct State: char {
  empty,
  filled,
  inprogress,
  processed,
};

struct merge_data_s
{
  std::vector<char> fwd_header;
  std::vector<char> rev_header;
  std::vector<char> fwd_sequence;
  std::vector<char> rev_sequence;
  std::vector<char> fwd_quality;
  std::vector<char> rev_quality;
  int64_t header_alloc = 0;
  int64_t seq_alloc = 0;
  int64_t fwd_length = 0;
  int64_t rev_length = 0;
  int64_t fwd_trunc = 0;
  int64_t rev_trunc = 0;
  int64_t fwd_abundance = 1;
  int64_t rev_abundance = 1;
  int64_t pair_no = 0;
  std::vector<char> merged_sequence;
  std::vector<char> merged_quality_v;
  int64_t merged_length = 0;
  int64_t merged_seq_alloc = 0;
  double ee_merged = 0;
  double ee_fwd = 0;
  double ee_rev = 0;
  int64_t fwd_errors = 0;
  int64_t rev_errors = 0;
  int64_t offset = 0;
  bool merged = false;
  Reason reason = Reason::undefined;
  State state = State::empty;
};

using merge_data_t = struct merge_data_s;


/* Cooperative-abort control. The state lives in core/mergepairs.cpp; the
   merge engine (get_qual/process) and the CLI worker loop both signal/poll it. */
auto request_merge_abort(MergeAbortReason reason, int value) -> void;
auto merge_aborted(std::memory_order order = std::memory_order_relaxed) -> bool;
auto merge_abort_reset() -> void;
auto report_merge_abort(struct Parameters const & parameters) -> void;

/* Merge engine entry points driven by both the CLI command and the library API. */
auto precompute_qual(struct Parameters const & parameters) -> QualityTables;
auto process(merge_data_t & a_read_pair,
             struct kh_handle_s & kmerhash,
             QualityTables const & tables,
             struct Parameters const & parameters) -> void;
