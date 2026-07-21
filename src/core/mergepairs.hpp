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

#include "utils/view.hpp"  // View<char> for read-only sequence/quality inputs
#include <array>
#include <string>

/* === Library API for embedding paired-end merging === */

/* A hard input error that made merge() reject the pair, distinct from an
   ordinary "did not overlap well enough" non-merge: a FASTQ quality symbol
   outside [fastq_qmin, fastq_qmax]. The CLI treats this as fatal; a library
   caller can inspect MergeResult::error and choose to skip the pair or abort its
   own run. (more_fwd_than_rev, the CLI reader's fwd/rev count mismatch, cannot
   occur on the single-pair path and so is not represented here.) */
enum struct MergeError {
  none,                 /* successful merge, or an ordinary non-merge */
  quality_below_qmin,
  quality_above_qmax,
};

/* Result of merging a single read pair, returned by MergePairs::merge().

   sequence and quality own their storage (std::string), so there is nothing
   for the caller to release: the buffers are freed automatically when the
   MergeResult goes out of scope (RAII). Both strings are empty when the merge
   fails (merged == false). quality holds ASCII-encoded quality symbols;
   sequence holds the merged DNA. Neither string carries a trailing NUL of its
   own, but .c_str() null-terminates as usual.

   A non-merge is ordinary (poor overlap, too many differences, ...) unless
   error != MergeError::none, which flags a hard input error (an out-of-range
   FASTQ quality value; error_value carries the offending value). */
struct MergeResult {
  bool merged = false;         /* true if merge succeeded */
  int merged_length = 0;       /* length of merged sequence */
  std::string sequence;        /* merged DNA (empty on failure) */
  std::string quality;         /* merged quality string, ASCII (empty on failure) */
  double ee_merged = 0.0;      /* expected errors in merged sequence */
  double ee_fwd = 0.0;         /* expected errors from forward read */
  double ee_rev = 0.0;         /* expected errors from reverse read */
  int fwd_errors = 0;          /* mismatches attributed to forward read */
  int rev_errors = 0;          /* mismatches attributed to reverse read */
  int overlap_length = 0;      /* length of overlap region */
  MergeError error = MergeError::none;  /* hard input error, if any (see MergeError) */
  int error_value = 0;         /* the offending quality value when error != none */
};

/* Number of ASCII quality symbols indexable in the score tables. */
constexpr auto n_quality_symbols = 128U;

/* Precomputed per-quality-symbol score tables built by a MergePairs session
   from the run's Parameters. Held privately by MergePairs; also produced
   directly by precompute_qual() for the fastq_mergepairs command (see
   core/mergepairs_internal.hpp). Once built it is read-only. Replaces the
   former file-scope globals in mergepairs.cpp. */
struct QualityTables {
  std::array<std::array<char,   n_quality_symbols>, n_quality_symbols> merge_qual_same {{}};
  std::array<std::array<char,   n_quality_symbols>, n_quality_symbols> merge_qual_diff {{}};
  std::array<std::array<double, n_quality_symbols>, n_quality_symbols> match_score     {{}};
  std::array<std::array<double, n_quality_symbols>, n_quality_symbols> mism_score      {{}};
  std::array<double, n_quality_symbols> q2p {{}};
};

/* One read of a pair, passed to MergePairs::merge().
   sequence and quality are read-only views over caller-owned buffers and MUST
   have the same length (one quality symbol per base). The views need not be
   null-terminated; only sequence.size() bases are read. */
struct MergeInput {
  View<char> sequence;
  View<char> quality;
};

/* A paired-end merging session.

   Construct one from a configured Parameters — this builds the quality score
   lookup tables once — then call merge() for each read pair. merge() is const
   and keeps no per-call mutable state, so a single MergePairs instance may be
   shared across threads.

   Relevant opt_* overrides (set on Parameters before vsearch_session_begin,
   then pass the same Parameters to the constructor and to merge()):
     opt_fastq_minovlen   — minimum overlap length (default 10)
     opt_fastq_maxdiffs   — max mismatches in overlap (default 10)
     opt_fastq_maxdiffpct — max mismatch % in overlap (default 100.0)
     opt_fastq_maxee      — max expected errors (default unlimited)
     opt_fastq_minlen     — min merged length (default 1)
     opt_fastq_maxlen     — max merged length (default unlimited)
     opt_fastq_maxns      — max Ns allowed (default unlimited)
     opt_fastq_ascii      — quality ASCII offset (default 33) */
class MergePairs {
public:
  explicit MergePairs(struct Parameters const & parameters);

  /* Merge a single forward/reverse read pair.
     parameters: the configured Parameters (same one passed to the
       constructor / vsearch_session_begin); supplies the merge tunables.
     fwd/rev: the forward and reverse reads (matching sequence + quality
       views; see MergeInput).
     Returns a MergeResult. On success result.merged is true and
     result.sequence / result.quality hold the merged read, owned by the
     result. On failure result.merged is false and both strings are empty; if
     the failure was a hard input error (a FASTQ quality value outside
     [fastq_qmin, fastq_qmax]) result.error says which and result.error_value
     carries the offending value, letting the caller distinguish it from an
     ordinary non-merge. Does not throw on an out-of-range quality.

     Thread-safe: a const MergePairs may be shared across threads. */
  auto merge(struct Parameters const & parameters,
             MergeInput const & fwd,
             MergeInput const & rev) const -> MergeResult;

private:
  QualityTables tables_;
};
