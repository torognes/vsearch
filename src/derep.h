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

auto derep(struct Parameters const & parameters, char * input_filename, bool use_header) -> void;

/* === Library API for in-memory dereplication === */

/* Result for one unique sequence, populated by derep_get_results(). */
struct derep_result_s {
  const char * header;    /* representative header (owned by session — valid until cleanup) */
  const char * sequence;  /* normalized sequence: uppercase DNA, U→T (owned by session) */
  uint64_t abundance;     /* total abundance (sum of all identical sequences) */
  uint64_t seqlen;        /* sequence length */
  int count;              /* number of input sequences that collapsed into this */
};

/* Opaque session state. */
struct derep_session_s;

auto derep_session_alloc() -> struct derep_session_s *;
auto derep_session_free(struct derep_session_s * ds) -> void;

/* Initialize a dereplication session.
   Does NOT require a loaded database — sequences are added via derep_add_sequence.
   Call vsearch_init_defaults() + vsearch_apply_defaults_fixups() before this. */
auto derep_session_init(struct derep_session_s * ds) -> void;

/* Add a single sequence to the dereplication session.
   Sequences are normalized (uppercase, U→T) internally.
   abundance: the abundance of this input sequence (typically 1). */
auto derep_add_sequence(struct derep_session_s * ds,
                        const char * header,
                        const char * sequence,
                        int seqlen,
                        int64_t abundance) -> void;

/* Finalize and retrieve results, sorted by abundance (descending).
   results: caller-provided array of at least max_results elements.
   max_results: maximum number of unique sequences to return.
   result_count: set to actual number of results returned.
   Results point into session-owned memory (valid until cleanup). */
auto derep_get_results(struct derep_session_s * ds,
                       struct derep_result_s * results,
                       int max_results,
                       int * result_count) -> void;

/* Clean up session resources. Call before derep_session_free(). */
auto derep_session_cleanup(struct derep_session_s * ds) -> void;
