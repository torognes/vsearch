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
#include <cstdio>  // std::FILE
#include <cstdint>  // int64_t


struct Database;
struct hit;

/* The query header and the query sequence (plus its reverse complement, when
   the minus strand was searched) are passed as views: the caller already knows
   both lengths, and a bare pointer would force each leaf to recover them with
   std::strlen. A view carrying the length also removes the separate qseqlen
   argument wherever the sequence itself is passed. */

/* Three shapes of hit argument, one per kind of writer, where the spelling used
   to be `struct hit const *` for all three:

   - View<struct hit>, for the writers that report a whole query at once. The
     view is the hits to report, already cut to --maxhits by the caller; an
     empty one is a query with no hit, which is not the same thing as a null
     pointer and now cannot be confused with one.
   - struct hit const &, for the per-hit writers the caller only ever reaches
     with a hit in hand.
   - struct hit const *, for the three per-hit writers that a caller does reach
     with no hit at all, to emit the query's "no hit" record under
     --output_no_hits. Null is a value here, not an oversight. */

auto results_show_alnout(std::FILE * output_handle,
                         View<struct hit> hits,
                         View<char> query_head,
                         View<char> qsequence,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void;

auto results_show_lcaout(std::FILE * output_handle,
                         View<struct hit> hits,
                         View<char> query_head,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void;

/* hit == nullptr: the query matched nothing */
auto results_show_blast6out_one(std::FILE * output_handle,
                                struct hit const * hit,
                                View<char> query_head,
                                int64_t qseqlen,
                                struct Database const & db) -> void;

/* hit == nullptr: the query matched nothing */
auto results_show_uc_one(std::FILE * output_handle,
                         struct hit const * hit,
                         View<char> query_head,
                         int64_t qseqlen,
                         int clusterno,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void;

/* hit == nullptr: the query matched nothing */
auto results_show_userout_one(std::FILE * output_handle,
                              struct hit const * hit,
                              View<char> query_head,
                              View<char> qsequence,
                              View<char> qsequence_rc,
                              struct Database const & db,
                              struct Parameters const & parameters) -> void;

auto results_show_fastapairs_one(std::FILE * output_handle,
                                 struct hit const & hit,
                                 View<char> query_head,
                                 View<char> qsequence,
                                 View<char> qsequence_rc,
                                 struct Database const & db,
                                 struct Parameters const & parameters) -> void;

auto results_show_qsegout_one(std::FILE * output_handle,
                              struct hit const & hit,
                              View<char> query_head,
                              View<char> qsequence,
                              View<char> qsequence_rc,
                              struct Parameters const & parameters) -> void;

auto results_show_tsegout_one(std::FILE * output_handle,
                              struct hit const & hit,
                              struct Database const & db,
                              struct Parameters const & parameters) -> void;

auto results_show_samheader(std::FILE * output_handle,
                            char const * dbname,
                            struct Database const & db,
                            struct Parameters const & parameters) -> void;

auto results_show_samout(std::FILE * output_handle,
                         View<struct hit> hits,
                         View<char> query_head,
                         View<char> qsequence,
                         View<char> qsequence_rc,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void;
