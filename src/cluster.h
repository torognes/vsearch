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

auto cluster_smallmem(char * cmdline, char * progheader) -> void;
auto cluster_fast(char * cmdline, char * progheader) -> void;
auto cluster_size(char * cmdline, char * progheader) -> void;
auto cluster_unoise(char * cmdline, char * progheader) -> void;

/* === Library API for embedding clustering === */

/* Result of assigning a single sequence to a cluster. */
struct cluster_result_s {
  bool is_centroid;            /* true if this sequence started a new cluster */
  int cluster_id;              /* cluster number (0-based) */
  int centroid_seqno;          /* database seqno of the cluster centroid */
  char centroid_label[1024];   /* centroid header (may truncate) */
  double identity;             /* identity to centroid (100.0 if is_centroid) */
  char cigar[4096];            /* CIGAR alignment string (empty if centroid) */
  bool cigar_truncated;        /* true if cigar was truncated to fit buffer */
};

/* Opaque session state for incremental clustering. */
struct cluster_session_s;

/* Allocate/free a clustering session. */
auto cluster_session_alloc() -> struct cluster_session_s *;
auto cluster_session_free(struct cluster_session_s * cs) -> void;

/* Initialize a clustering session.
   Requires: global opt_* set (including opt_id for identity threshold),
   database loaded, masked, and dbindex_prepare() called with bitmap=1.
   Do NOT call dbindex_addallsequences — centroids are indexed
   incrementally as they are discovered.

   Database must be pre-sorted by length (cluster_fast) or
   abundance (cluster_size) before loading. */
auto cluster_session_init(struct cluster_session_s * cs) -> void;

/* Assign a single database sequence to a cluster.
   Must be called sequentially (seqno 0, 1, 2, ...).
   Sequences matching an existing centroid (above opt_id) are assigned
   to that cluster. Unmatched sequences become new centroids.
   Single-threaded only. */
auto cluster_assign_single(struct cluster_session_s * cs,
                            int seqno,
                            struct cluster_result_s * result) -> void;

/* Clean up clustering session resources.
   Call before cluster_session_free(). */
auto cluster_session_cleanup(struct cluster_session_s * cs) -> void;
