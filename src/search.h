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

auto usearch_global(struct Parameters const & parameters, char * cmdline, char * progheader) -> void;

/* === Library API for embedding global search === */

/* Result of a single search hit. */
struct search_result_s {
  int target;                  /* database sequence index */
  char target_label[1024];     /* target header (may truncate) */
  double id;                   /* percent identity (method per opt_iddef) */
  int matches;                 /* matching columns */
  int mismatches;              /* mismatching columns */
  int gaps;                    /* gap columns */
  int alignment_length;        /* total alignment length */
  int query_length;            /* query sequence length */
  int target_length;           /* target sequence length */
  bool accepted;               /* true if passed identity threshold */
};

struct searchinfo_s;

/* Allocate/free opaque per-thread search state. */
auto search_info_alloc() -> struct searchinfo_s *;
auto search_info_free(struct searchinfo_s * si) -> void;

/* Initialize search state for library use.
   Requires: global opt_* set, database loaded and indexed.
   One instance per thread — do NOT share across threads. */
auto search_init(struct searchinfo_s * si) -> void;

/* Search for a single query against the global database.
   Plus-strand only (library API does not support reverse complement).
   results: caller-allocated array of at least max_results elements.
   result_count: number of results populated on return.
   Results are ordered by identity (descending). */
auto search_single(struct searchinfo_s * si,
                    const char * query_seq,
                    const char * query_head,
                    int query_len,
                    int query_size,
                    struct search_result_s * results,
                    int max_results,
                    int * result_count) -> void;

/* Clean up per-thread search state.
   Call before search_info_free(). */
auto search_cleanup(struct searchinfo_s * si) -> void;
