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

auto fastq_mergepairs(struct Parameters const & parameters) -> void;

/* === Library API for embedding paired-end merging === */

/* Result of merging a single read pair. */
struct merge_result_s {
  bool merged;                 /* true if merge succeeded */
  int merged_length;           /* length of merged sequence */
  char merged_sequence[10000]; /* merged DNA sequence (null-terminated) */
  char merged_quality[10000];  /* merged quality string (null-terminated) */
  double ee_merged;            /* expected errors in merged sequence */
  double ee_fwd;               /* expected errors from forward read */
  double ee_rev;               /* expected errors from reverse read */
  int fwd_errors;              /* mismatches attributed to forward read */
  int rev_errors;              /* mismatches attributed to reverse read */
  int overlap_length;          /* length of overlap region */
};

/* Initialize the quality score lookup table.
   Must be called once before mergepairs_single().
   Requires: opt_fastq_ascii, opt_fastq_qmin, opt_fastq_qmax set
   (vsearch_init_defaults provides correct values). */
auto mergepairs_init() -> void;

/* Merge a single forward/reverse read pair.
   fwd_seq/rev_seq: null-terminated DNA sequences.
   fwd_qual/rev_qual: null-terminated quality strings (ASCII-encoded).
   fwd_len/rev_len: sequence lengths.
   fwd_header/rev_header: null-terminated read headers.
   result: output struct populated on return.
   Returns 0 on success (merged), -1 on failure (not merged).

   Thread-safe: each call is independent (no shared mutable state).

   Relevant opt_* overrides (after vsearch_init_defaults):
     opt_fastq_minovlen   — minimum overlap length (default 10)
     opt_fastq_maxdiffs   — max mismatches in overlap (default 10)
     opt_fastq_maxdiffpct — max mismatch % in overlap (default 100.0)
     opt_fastq_maxee      — max expected errors (default unlimited)
     opt_fastq_minlen     — min merged length (default 1)
     opt_fastq_maxlen     — max merged length (default unlimited)
     opt_fastq_maxns      — max Ns allowed (default unlimited)
     opt_fastq_ascii      — quality ASCII offset (default 33) */
auto mergepairs_single(const char * fwd_seq,
                        const char * fwd_qual,
                        int fwd_len,
                        const char * rev_seq,
                        const char * rev_qual,
                        int rev_len,
                        const char * fwd_header,
                        const char * rev_header,
                        struct merge_result_s * result) -> int;
