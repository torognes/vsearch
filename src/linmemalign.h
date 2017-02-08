/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

class LinearMemoryAligner
{
  char op;
  int64_t op_run;
  int64_t cigar_alloc;
  int64_t cigar_length;
  char * cigar_string;
  
  char * a_seq;
  char * b_seq;

  int64_t * scorematrix;

  int64_t q;
  int64_t r;

  /* gap penalties for open/extension query/target left/interior/right */
  int64_t go_q_l;
  int64_t go_t_l;
  int64_t go_q_i;
  int64_t go_t_i;
  int64_t go_q_r;
  int64_t go_t_r;
  int64_t ge_q_l;
  int64_t ge_t_l;
  int64_t ge_q_i;
  int64_t ge_t_i;
  int64_t ge_q_r;
  int64_t ge_t_r;

  size_t vector_alloc;

  int64_t * HH;
  int64_t * EE;
  int64_t * XX;
  int64_t * YY;
  
  void cigar_reset();

  void cigar_flush();

  void cigar_add(char _op, int64_t run);

  inline int64_t subst_score(int64_t x, int64_t y)
  {
    /* return substitution score for replacing symbol at position x in a
       with symbol at position y in b */
    return scorematrix[chrmap_4bit[(int)(b_seq[y])] * 16 + 
                       chrmap_4bit[(int)(a_seq[x])]];
  }

  void diff(int64_t a_start,
            int64_t b_start,
            int64_t a_len,
            int64_t b_len,
            bool gap_b_left,  /* gap open left of b      */
            bool gap_b_right, /* gap open right of b     */
            bool a_left,      /* includes left end of a  */
            bool a_right,     /* includes right end of a */
            bool b_left,      /* includes left end of b  */
            bool b_right);    /* includes right end of b */
  
  void alloc_vectors(size_t N);

  void show_matrix();

public:

  LinearMemoryAligner();

  ~LinearMemoryAligner();
  
  int64_t * scorematrix_create(int64_t match, int64_t mismatch);
  
  void set_parameters(int64_t * _scorematrix,
                      int64_t _gap_open_query_left,
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
                      int64_t _gap_extension_target_right);
  
  char * align(char * _a_seq,
               char * _b_seq,
               int64_t M,
               int64_t N);

  void alignstats(char * cigar,
                  char * a_seq,
                  char * b_seq,
                  int64_t * nwscore,
                  int64_t * nwalignmentlength,
                  int64_t * nwmatches,
                  int64_t * nwmismatches,
                  int64_t * nwgaps);

};
