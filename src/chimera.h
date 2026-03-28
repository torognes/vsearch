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

constexpr auto maxparents = 20; /* max, could be fewer */

auto chimera(struct Parameters const & parameters) -> void;

/* === Library API for embedding chimera detection === */

/* Result of chimera detection for a single query.
   Fields match vsearch's --uchimeout 18-column format.
   For non-chimeric results (flag='N'), only query_label and flag are
   populated; all other fields are zero.
   Labels may be silently truncated to 1023 characters. */
struct chimera_result_s {
  double score;               /* h-score */
  char query_label[1024];     /* query header (may truncate) */
  char parent_a_label[1024];  /* parent A header (or empty if none) */
  char parent_b_label[1024];  /* parent B header (or empty if none) */
  char closest_parent_label[1024]; /* closest parent header */
  double id_query_model;      /* query-to-model identity % */
  double id_query_a;          /* query-to-parentA identity % */
  double id_query_b;          /* query-to-parentB identity % */
  double id_a_b;              /* parentA-to-parentB identity % */
  double id_query_top;        /* query-to-closest-parent identity % */
  int left_yes, left_no, left_abstain;
  int right_yes, right_no, right_abstain;
  double divergence;
  char flag;                  /* 'Y', 'N', or '?' */
};

struct chimera_info_s;

/* Allocate an opaque chimera_info_s on the heap.
   Callers that cannot see the full struct definition use this
   instead of stack/member allocation. Free with chimera_info_free(). */
auto chimera_info_alloc() -> struct chimera_info_s *;

/* Free a chimera_info_s allocated by chimera_info_alloc().
   Does NOT call chimera_detect_cleanup — call that first if init was called. */
auto chimera_info_free(struct chimera_info_s * ci) -> void;

/* Initialize chimera detection subsystem for library use.
   Must call before chimera_detect_single.
   Sets search-shaping globals (opt_maxaccepts, opt_maxrejects, tophits, etc.)
   that chimera() normally sets but the library API path bypasses.
   Requires: global opt_* scoring/penalty variables already set,
   global db + dbindex loaded and indexed.
   ci is per-thread working state — each thread must have its own instance.
   Do NOT share ci across threads. */
auto chimera_detect_init(struct chimera_info_s * ci) -> void;

/* Detect chimera for a single query sequence.
   Supports both uchime_ref and uchime_denovo modes (based on opt_chimeras_denovo).
   ci: per-thread working state (from chimera_detect_init). NOT thread-safe if shared.
   query_seq: null-terminated query sequence (DNA, uppercase).
   query_head: null-terminated query header.
   query_len: length of query sequence.
   query_size: abundance (1 for uchime_ref, actual count for uchime_denovo).
   result: output struct populated on return.
   Returns 0 on success. */
auto chimera_detect_single(struct chimera_info_s * ci,
                           const char * query_seq,
                           const char * query_head,
                           int query_len,
                           int query_size,
                           struct chimera_result_s * result) -> int;

/* Clean up per-thread chimera working state.
   Frees all resources allocated by chimera_detect_init (SIMD aligners,
   unique k-mer finders, minheaps, CIGAR strings). */
auto chimera_detect_cleanup(struct chimera_info_s * ci) -> void;
