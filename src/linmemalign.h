/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
  long op_run;
  long cigar_alloc;
  long cigar_length;
  char * cigar_string;
  
  char * a_seq;
  char * b_seq;

  long * scorematrix;

  long q;
  long r;

  /* gap penalties for open/extension query/target left/interior/right */
  long go_q_l;
  long go_t_l;
  long go_q_i;
  long go_t_i;
  long go_q_r;
  long go_t_r;
  long ge_q_l;
  long ge_t_l;
  long ge_q_i;
  long ge_t_i;
  long ge_q_r;
  long ge_t_r;

  size_t vector_alloc;

  long * HH;
  long * EE;
  long * XX;
  long * YY;
  
  void cigar_reset();

  void cigar_flush();

  void cigar_add(char _op, long run);

  inline long subst_score(long x, long y)
  {
    /* return substitution score for replacing symbol at position x in a
       with symbol at position y in b */
    return scorematrix[chrmap_4bit[(int)(b_seq[y])] * 16 + 
                       chrmap_4bit[(int)(a_seq[x])]];
  }

  void diff(long a_start,
            long b_start,
            long a_len,
            long b_len,
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
  
  long * scorematrix_create(long match, long mismatch);
  
  void set_parameters(long * _scorematrix,
                      long _gap_open_query_left,
                      long _gap_open_target_left,
                      long _gap_open_query_interior,
                      long _gap_open_target_interior,
                      long _gap_open_query_right,
                      long _gap_open_target_right,
                      long _gap_extension_query_left,
                      long _gap_extension_target_left,
                      long _gap_extension_query_interior,
                      long _gap_extension_target_interior,
                      long _gap_extension_query_right,
                      long _gap_extension_target_right);
  
  char * align(char * _a_seq,
               char * _b_seq,
               long M,
               long N);

  void alignstats(char * cigar,
                  char * a_seq,
                  char * b_seq,
                  long * nwscore,
                  long * nwalignmentlength,
                  long * nwmatches,
                  long * nwmismatches,
                  long * nwgaps);

};
