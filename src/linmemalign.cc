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

#include "vsearch.h"

/*

  Compute the optimal global alignment of two sequences
  in linear space using the divide and conquer method.
  
  These functions are based on the following articles:
  - Hirschberg (1975) Comm ACM 18:341-343
  - Myers & Miller (1988) CABIOS 4:11-17
  
  The method has been adapted for the use of different
  gap penalties for query/target/left/interior/right gaps.

  scorematrix consists of 16x16 int64_t integers
  
  Sequences and alignment matrix:
  A/a/i/query/q/downwards/vertical/top/bottom
  B/b/j/target/t/rightwards/horizontal/left/right

  f corresponds to score ending with gap in A/query
  EE corresponds to score ending with gap in B/target

*/

LinearMemoryAligner::LinearMemoryAligner()
{
  scorematrix = 0;

  cigar_alloc = 0;
  cigar_string = 0;

  vector_alloc = 0;
  HH = 0;
  EE = 0;
  XX = 0;
  YY = 0;
}

LinearMemoryAligner::~LinearMemoryAligner()
{
  if (cigar_string)
    xfree(cigar_string);
  if (HH)
    xfree(HH);
  if (EE)
    xfree(EE);
  if (XX)
    xfree(XX);
  if (YY)
    xfree(YY);
}

int64_t * LinearMemoryAligner::scorematrix_create(int64_t match, int64_t mismatch)
{
  int64_t * newscorematrix = (int64_t*) xmalloc(16*16*sizeof(int64_t));

  for(int i=0; i<16; i++)
    for(int j=0; j<16; j++)
      {
        int64_t value;
        if ((i==0) || (j==0) || (i>4) || (j>4))
          value = 0;
        else if (i==j)
          value = match;
        else
          value = mismatch;
        newscorematrix[16*i+j] = value;
      }
  return newscorematrix;
}

void LinearMemoryAligner::alloc_vectors(size_t x)
{
  if (vector_alloc < x)
    {
      vector_alloc = x;

      if (HH)
        xfree(HH);
      if (EE)
        xfree(EE);
      if (XX)
        xfree(XX);
      if (YY)
        xfree(YY);

      HH = (int64_t*) xmalloc(vector_alloc * (sizeof(int64_t)));
      EE = (int64_t*) xmalloc(vector_alloc * (sizeof(int64_t)));
      XX = (int64_t*) xmalloc(vector_alloc * (sizeof(int64_t)));
      YY = (int64_t*) xmalloc(vector_alloc * (sizeof(int64_t)));
    }
}

void LinearMemoryAligner::cigar_reset()
{
  if (cigar_alloc < 1)
    {
      cigar_alloc = 64;
      cigar_string = (char*) xrealloc(cigar_string, cigar_alloc);
    }
  cigar_string[0] = 0;
  cigar_length = 0;
  op = 0;
  op_run = 0;
}

void LinearMemoryAligner::cigar_flush()
{
  if (op_run > 0)
  {
    while (1)
    {
      /* try writing string until enough memory has been allocated */

      int64_t rest = cigar_alloc - cigar_length;
      int n;
      if (op_run > 1)
        n = snprintf(cigar_string + cigar_length,
                     rest,
                     "%" PRId64 "%c", op_run, op);
      else
        n = snprintf(cigar_string + cigar_length,
                     rest,
                     "%c", op);
      if (n < 0)
        {
          fatal("snprintf returned a negative number.\n");
        }
      else if (n >= rest)
        {
          cigar_alloc += MAX(n - rest + 1, 64);
          cigar_string = (char*) xrealloc(cigar_string, cigar_alloc);
        }
      else
        {
          cigar_length += n;
          break;
        }
    }
  }
}

void LinearMemoryAligner::cigar_add(char _op, int64_t run)
  {
    if (op == _op)
      op_run += run;
    else
      {
        cigar_flush();
        op = _op;
        op_run = run;
      }
  }

void LinearMemoryAligner::show_matrix()
{
  for(int i=0; i<16; i++)
    {
      printf("%2d:", i);
      for(int j=0; j<16; j++)
        printf(" %2" PRId64, scorematrix[16*i+j]);
      printf("\n");
    }
}

void LinearMemoryAligner::diff(int64_t a_start,
                               int64_t b_start,
                               int64_t a_len,
                               int64_t b_len,
                               bool gap_b_left,  /* gap open left of b      */
                               bool gap_b_right, /* gap open right of b     */
                               bool a_left,      /* includes left end of a  */
                               bool a_right,     /* includes right end of a */
                               bool b_left,      /* includes left end of b  */
                               bool b_right)     /* includes right end of b */
{
  if (b_len == 0)
    {
      /* B and possibly A is empty */
      if (a_len > 0)
        {
	  // Delete a_len from A
	  // AAA
	  // ---

          cigar_add('D', a_len);
        }
    }
  else if (a_len == 0)
    {
      /* A is empty, B is not */

      // Delete b_len from B
      // ---
      // BBB
      
      cigar_add('I', b_len);
    }
  else if (a_len == 1)
    {
      /*
        Convert 1 symbol from A to b_len symbols from B
        b_len >= 1
      */


      int64_t MaxScore;
      int64_t best;

      int64_t Score = 0;

      /* First possibility */

      // Delete 1 from A, Insert b_len from B
      // A----
      // -BBBB
      
      /* gap penalty for gap in B of length 1 */
      
      if (! gap_b_left)
        Score -= b_left ? go_t_l : go_t_i;
      
      Score -= b_left ? ge_t_l : ge_t_i;
      
      /* gap penalty for gap in A of length b_len */
      
      Score -= a_right ? go_q_r + b_len * ge_q_r : go_q_i + b_len * ge_q_i;

      MaxScore = Score;
      best = -1;


      /* Second possibility */
      
      // Insert b_len from B, Delete 1 from A
      // ----A
      // BBBB-
      
      /* gap penalty for gap in A of length b_len */
      
      Score -= a_left ? go_q_l + b_len * ge_q_l : go_q_i + b_len * ge_q_i;
      
      /* gap penalty for gap in B of length 1 */
      
      if (! gap_b_right)
        Score -= b_right ? go_t_r : go_t_i;
      
      Score -= b_right ? ge_t_r : ge_t_i;
      
      if (Score > MaxScore)
        {
          MaxScore = Score;
          best = b_len;
        }


      /* Third possibility */

      for (int64_t j = 0; j < b_len; j++)
	{
	  // Insert zero or more from B, replace 1, insert rest of B
	  // -A--
	  // BBBB

          Score = 0;

          if (j > 0)
            Score -= a_left ? go_q_l + j * ge_q_l : go_q_i + j * ge_q_i;
          
	  Score += subst_score(a_start, b_start + j);
          
          if (j < b_len - 1)
            Score -= a_right ? 
              go_q_r + (b_len-1-j) * ge_q_r : 
              go_q_i + (b_len-1-j) * ge_q_i;
          
	  if (Score > MaxScore)
	    {
	      MaxScore = Score;
	      best = j;
	    }
	}

      if (best == -1)
	{
	  cigar_add('D', 1);
          cigar_add('I', b_len);
	}
      else if (best == b_len)
	{
	  cigar_add('I', b_len);
          cigar_add('D', 1);
	}
      else
	{
	  if (best > 0)
	    cigar_add('I', best);
          cigar_add('M', 1);
	  if (best < b_len - 1)
	    cigar_add('I', b_len - 1 - best);
	}
    }
  else
    {
      /* a_len >= 2, b_len >= 1 */

      int64_t I = a_len / 2;
      int64_t i, j;

      // Compute HH & EE in forward phase
      // Upper part

      /* initialize HH and EE for values corresponding to 
         empty seq A vs B of j symbols,
         i.e. a gap of length j in A                 */

      HH[0] = 0;
      EE[0] = 0;

      for (j = 1; j <= b_len; j++)
	{
	  HH[j] = - (a_left ? go_q_l + j * ge_q_l : go_q_i + j * ge_q_i);
	  EE[j] = LONG_MIN;
	}

      /* compute matrix */

      for (i = 1; i <= I; i++)
	{
	  int64_t p = HH[0];
          
          int64_t h = - (b_left ?
                      (gap_b_left ? 0 : go_t_l) + i * ge_t_l :
                      (gap_b_left ? 0 : go_t_i) + i * ge_t_i);
          
          HH[0] = h;
	  int64_t f = LONG_MIN;
          
	  for (j = 1; j <= b_len; j++)
	    {
              f = MAX(f, h - go_q_i) - ge_q_i;
              if (b_right && (j==b_len))
                EE[j] = MAX(EE[j], HH[j] - go_t_r) - ge_t_r;
              else
                EE[j] = MAX(EE[j], HH[j] - go_t_i) - ge_t_i;

              h = p + subst_score(a_start + i - 1, b_start + j - 1);
              
	      if (f > h)
		h = f;
	      if (EE[j] > h)
		h = EE[j];
	      p = HH[j];
	      HH[j] = h;
	    }
	}

      EE[0] = HH[0];

      // Compute XX & YY in reverse phase
      // Lower part

      /* initialize XX and YY */

      XX[0] = 0;
      YY[0] = 0;

      for (j = 1; j <= b_len; j++)
	{
	  XX[j] = - (a_right ? go_q_r + j * ge_q_r : go_q_i + j * ge_q_i);
	  YY[j] = LONG_MIN;
	}
      
      /* compute matrix */

      for (i = 1; i <= a_len - I; i++)
	{
	  int64_t p = XX[0];

          int64_t h = - (b_right ?
                      (gap_b_right ? 0 : go_t_r) + i * ge_t_r :
                      (gap_b_right ? 0 : go_t_i) + i * ge_t_i);
          XX[0] = h;
	  int64_t f = LONG_MIN;

	  for (j = 1; j <= b_len; j++)
	    {
              f = MAX(f, h - go_q_i) - ge_q_i;
              if (b_left && (j==b_len))
                YY[j] = MAX(YY[j], XX[j] - go_t_l) - ge_t_l;
              else
                YY[j] = MAX(YY[j], XX[j] - go_t_i) - ge_t_i;

              h = p + subst_score(a_start + a_len - i, b_start + b_len - j);

	      if (f > h)
		h = f;
	      if (YY[j] > h)
		h = YY[j];
	      p = XX[j];
	      XX[j] = h;
	    }
	}

      YY[0] = XX[0];


      /* find maximum score along division line */

      int64_t MaxScore0 = LONG_MIN;
      int64_t best0 = -1;

      /* solutions with diagonal at break */
      
      for (j=0; j <= b_len; j++)
	{
	  int64_t Score = HH[j] + XX[b_len - j];

	  if (Score > MaxScore0)
	    {
	      MaxScore0 = Score;
	      best0 = j;
	    }
	}

      int64_t MaxScore1 = LONG_MIN;
      int64_t best1 = -1;

      /* solutions that end with a gap in b from both ends at break */

      for (j=0; j <= b_len; j++)
	{
          int64_t g;
          if (b_left && (j==0))
            g = go_t_l;
          else if (b_right && (j==b_len))
            g = go_t_r;
          else
            g = go_t_i;

	  int64_t Score = EE[j] + YY[b_len - j] + g;

	  if (Score > MaxScore1)
	    {
	      MaxScore1 = Score;
	      best1 = j;
	    }
	}
      
      int64_t P;
      int64_t best;

      if (MaxScore0 > MaxScore1)
        {
          P = 0;
          best = best0;
        }
      else if (MaxScore1 > MaxScore0)
        {
          P = 1;
          best = best1;
        }
      else
        {
          if (best0 <= best1)
            {
              P = 0;
              best = best0;
            }
          else
            {
              P = 1;
              best = best1;
            }
        }

      /* recursively compute upper left and lower right parts */

      if (P == 0)
	{
	  diff(a_start,               b_start,     
               I,                     best,      
               gap_b_left,            false,
               a_left,                false,
               b_left,                b_right && (best == b_len));

	  diff(a_start + I,           b_start + best,
               a_len - I,             b_len - best,
               false,                 gap_b_right,
               false,                 a_right,
               b_left && (best == 0), b_right);
        }
      else if (P == 1)
	{
	  diff(a_start,               b_start,      
               I - 1,                 best,       
               gap_b_left,            true,
               a_left,                false,
               b_left,                b_right && (best == b_len));

	  cigar_add('D', 2);

	  diff(a_start + I + 1,       b_start + best, 
               a_len - I - 1,         b_len - best, 
               true,                  gap_b_right,
               false,                 a_right,
               b_left && (best == 0), b_right);
	}
    }
}

void LinearMemoryAligner::set_parameters(int64_t * _scorematrix,
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
                                         int64_t _gap_extension_target_right)
{
  scorematrix = _scorematrix;

  /* a = query/q   b = t/target */

  go_q_l = _gap_open_query_left;
  go_t_l = _gap_open_target_left;
  go_q_i = _gap_open_query_interior;
  go_t_i = _gap_open_target_interior;
  go_q_r = _gap_open_query_right;
  go_t_r = _gap_open_target_right;
  ge_q_l = _gap_extension_query_left;
  ge_t_l = _gap_extension_target_left;
  ge_q_i = _gap_extension_query_interior;
  ge_t_i = _gap_extension_target_interior;
  ge_q_r = _gap_extension_query_right;
  ge_t_r = _gap_extension_target_right;

  q = _gap_open_query_interior;
  r = _gap_extension_query_interior;
}
  


char * LinearMemoryAligner::align(char * _a_seq,
                                  char * _b_seq,
                                  int64_t a_len,
                                  int64_t b_len)
{
  /* copy parameters */
  a_seq = _a_seq;
  b_seq = _b_seq;

  /* init cigar operations */
  cigar_reset();

  /* allocate enough memory for vectors */
  alloc_vectors(b_len+1);

  /* perform alignment */
  diff(0, 0, a_len, b_len, 0, 0, true, true, true, true);

  /* ensure entire cigar has been written */
  cigar_flush();

  /* return cigar */
  return cigar_string;
}

void LinearMemoryAligner::alignstats(char * cigar,
                                     char * _a_seq,
                                     char * _b_seq,
                                     int64_t * _nwscore,
                                     int64_t * _nwalignmentlength,
                                     int64_t * _nwmatches,
                                     int64_t * _nwmismatches,
                                     int64_t * _nwgaps)
{
  a_seq = _a_seq;
  b_seq = _b_seq;

  int64_t nwscore = 0;
  int64_t nwalignmentlength = 0;
  int64_t nwmatches = 0;
  int64_t nwmismatches = 0;
  int64_t nwgaps = 0;

  int64_t a_pos = 0;
  int64_t b_pos = 0;

  char * p = cigar;

  int64_t g;

  while (*p)
    {
      int64_t run = 1;
      int scanlength = 0;
      sscanf(p, "%" PRId64 "%n", &run, &scanlength);
      p += scanlength;
      switch (*p++)
        {
        case 'M':
          nwalignmentlength += run;
          for(int64_t k=0; k<run; k++)
            {
              nwscore += subst_score(a_pos, b_pos);

              if (chrmap_4bit[(int)(a_seq[a_pos])] ==
                  chrmap_4bit[(int)(b_seq[b_pos])])
                nwmatches++;
              else
                nwmismatches++;
              
              a_pos++;
              b_pos++;
            }
          break;
              
        case 'I':
          if ((a_pos == 0) && (b_pos == 0))
            g = go_q_l + run * ge_q_l;
          else if (*p == 0)
            g = go_q_r + run * ge_q_r;
          else
            g = go_q_i + run * ge_q_i;
          nwscore -= g;
          nwgaps++;
          nwalignmentlength += run;
          b_pos += run;
          break;
              
        case 'D':
          if ((a_pos == 0) && (b_pos == 0))
            g = go_t_l + run * ge_t_l;
          else if (*p == 0)
            g = go_t_r + run * ge_t_r;
          else
            g = go_t_i + run * ge_t_i;
          nwscore -= g;
          nwgaps++;
          nwalignmentlength += run;
          a_pos += run;
          break;
        }
    }

  *_nwscore = nwscore;
  *_nwalignmentlength = nwalignmentlength;
  *_nwmatches = nwmatches;
  *_nwmismatches = nwmismatches;
  *_nwgaps = nwgaps;
}

