/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "maps.h"
#include <cstdio>  // std::FILE, std::size_t
#include <cstdint>  // int64_t
#include <vector>


struct Scoring {
  int64_t match = 0;
  int64_t mismatch = 0;
  int64_t gap_open_query_interior = 0;
};


class LinearMemoryAligner
{
private:
  char op = '\0';
  int64_t op_run = 0;
  int64_t cigar_alloc = 0;
  int64_t cigar_length = 0;
  std::vector<char> cigar_string;

  char * a_seq = nullptr;
  char * b_seq = nullptr;

  std::vector<int64_t> scorematrix;

  int64_t q = 0;  // general gap opening penalty (same as gap open query interior)
  int64_t r = 0;  // general gap extension penalty (same as gap extension query interior)

  /* gap penalties for open/extension query/target left/interior/right */
  int64_t go_q_l = 0;
  int64_t go_t_l = 0;
  int64_t go_q_i = 0;
  int64_t go_t_i = 0;
  int64_t go_q_r = 0;
  int64_t go_t_r = 0;
  int64_t ge_q_l = 0;
  int64_t ge_t_l = 0;
  int64_t ge_q_i = 0;
  int64_t ge_t_i = 0;
  int64_t ge_q_r = 0;
  int64_t ge_t_r = 0;

  std::size_t vector_alloc = 0;

  std::vector<int64_t> HH;
  std::vector<int64_t> EE;
  std::vector<int64_t> XX;
  std::vector<int64_t> YY;

  // initializers
  auto scorematrix_create(struct Scoring const & scoring) -> void;

  auto cigar_reset() -> void;

  auto cigar_flush() -> void;

  auto cigar_add(char _op, int64_t run) -> void;

  auto subst_score(int64_t lhs_pos, int64_t rhs_pos) -> int64_t;

  auto diff(int64_t a_start,
            int64_t b_start,
            int64_t a_len,
            int64_t b_len,
            bool gap_b_left,  /* gap open left of b      */
            bool gap_b_right, /* gap open right of b     */
            bool a_left,      /* includes left end of a  */
            bool a_right,     /* includes right end of a */
            bool b_left,      /* includes left end of b  */
            bool b_right) -> void;    /* includes right end of b */

  auto alloc_vectors(std::size_t size) -> void;

public:

  LinearMemoryAligner(struct Scoring const & scoring);

  ~LinearMemoryAligner();

  auto set_parameters(int64_t _gap_open_query_left,
                      int64_t _gap_open_target_left,
                      int64_t _gap_open_query_interior,
                      int64_t _gap_open_target_interior,
                      int64_t _gap_open_query_right,
                      int64_t _gap_open_target_right,
                      int64_t _gap_extension_query_left,
                      int64_t _gap_extension_target_left,
                      int64_t _gap_extension_query_interior,
                      int64_t _gap_extension_target_interior,
                      int64_t _gap_extension_query_right,
                      int64_t _gap_extension_target_right) -> void;

  auto align(char * _a_seq,
             char * _b_seq,
             int64_t a_len,
             int64_t b_len) -> char *;

  auto alignstats(char * cigar,
                  char * a_seq,
                  char * b_seq,
                  int64_t * nwscore,
                  int64_t * nwalignmentlength,
                  int64_t * nwmatches,
                  int64_t * nwmismatches,
                  int64_t * nwgaps) -> void;

};
