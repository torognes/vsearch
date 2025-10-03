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

#include "vsearch.h"
#include "linmemalign.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include <algorithm>  // std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t
#include <cstdio>  // std::printf, std::size_t, std::snprintf, std::sscanf
#include <limits>
// #include <vector>


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

constexpr auto minimal_length = int64_t{64};


LinearMemoryAligner::LinearMemoryAligner(struct Scoring const & scoring)
    : go_q_l(scoring.gap_open_query_left),
      go_t_l(scoring.gap_open_target_left),
      go_q_i(scoring.gap_open_query_interior),
      go_t_i(scoring.gap_open_target_interior),
      go_q_r(scoring.gap_open_query_right),
      go_t_r(scoring.gap_open_target_right),
      ge_q_l(scoring.gap_extension_query_left),
      ge_t_l(scoring.gap_extension_target_left),
      ge_q_i(scoring.gap_extension_query_interior),
      ge_t_i(scoring.gap_extension_target_interior),
      ge_q_r(scoring.gap_extension_query_right),
      ge_t_r(scoring.gap_extension_target_right)
{
  scorematrix_fill(scoring);
}


/*
  Expected score matrix (if option N is mismatch):

     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
     -  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
0  - M  X  X  0  X  0  0  0  X  0  0  0  0  0  0  X
1  A X  M  X  0  X  0  0  0  X  0  0  0  0  0  0  X
2  C X  X  M  0  X  0  0  0  X  0  0  0  0  0  0  X
3  M 0  0  0  M  0  0  0  0  0  0  0  0  0  0  0  X
4  G X  X  X  0  M  0  0  0  X  0  0  0  0  0  0  X
5  R 0  0  0  0  0  M  0  0  0  0  0  0  0  0  0  X
6  S 0  0  0  0  0  0  M  0  0  0  0  0  0  0  0  X
7  V 0  0  0  0  0  0  0  M  0  0  0  0  0  0  0  X
8  T X  X  X  0  X  0  0  0  M  0  0  0  0  0  0  X
9  W 0  0  0  0  0  0  0  0  0  M  0  0  0  0  0  X
10 Y 0  0  0  0  0  0  0  0  0  0  M  0  0  0  0  X
11 H 0  0  0  0  0  0  0  0  0  0  0  M  0  0  0  X
12 K 0  0  0  0  0  0  0  0  0  0  0  0  M  0  0  X
13 D 0  0  0  0  0  0  0  0  0  0  0  0  0  M  0  X
14 B 0  0  0  0  0  0  0  0  0  0  0  0  0  0  M  X
15 N X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  M

  M = match, X = mismatch

  Map from ascii to 4-bit nucleotide code
  -:  0
  A:  1
  C:  2
  M:  3
  G:  4
  R:  5
  S:  6
  V:  7
  T:  8
  W:  9
  Y: 10
  H: 11
  K: 12
  D: 13
  B: 14
  N: 15
*/

auto LinearMemoryAligner::scorematrix_fill(struct Scoring const & scoring) -> void {
 std::vector<char> const nucleotides = {'-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

  // fill-in the score matrix
  for (auto const row_nuc : nucleotides) {
    auto const row = map_4bit(row_nuc);
    for (auto const column_nuc : nucleotides) {
      auto const column = map_4bit(column_nuc);
      if (is_ambiguous_4bit(row) or is_ambiguous_4bit(column)) {
        continue;  // then score is 0; already zero-initialized
      }
      if (row == column) { // diagonal
        scorematrix[row][column] = scoring.match;
      }
      else {
        scorematrix[row][column] = scoring.mismatch;
      }
    }
  }

  // if alignment with N is set to be a mismatch
  if (opt_n_mismatch) {
    // last column
    for (auto & row : scorematrix) {
      row.back() = scoring.mismatch;
    }
    // last row
    auto & last_row = scorematrix.back();
    std::fill(last_row.begin(), last_row.end(), scoring.mismatch);
  }
}


auto LinearMemoryAligner::alloc_vectors(std::size_t const size) -> void {
  if (vector_alloc >= size) { return; }
  vector_alloc = size;
  HH.resize(vector_alloc);
  EE.resize(vector_alloc);
  XX.resize(vector_alloc);
  YY.resize(vector_alloc);
}


auto LinearMemoryAligner::cigar_reset() -> void
{
  if (cigar_alloc < 1)
    {
      cigar_alloc = minimal_length;
      cigar_string.resize(cigar_alloc);
    }
  cigar_string[0] = '\0';
  cigar_length = 0;
  op = '\0';
  op_run = 0;
}


auto LinearMemoryAligner::cigar_flush() -> void
{
  if (op_run <= 0) { return; }
  while (true)
    {
      /* try writing string until enough memory has been allocated */

      auto const rest = cigar_alloc - cigar_length;
      auto n = 0;
      if (op_run > 1)
        {
          n = std::snprintf(&cigar_string[cigar_length],
                       rest,
                       "%" PRId64 "%c", op_run, op);
        }
      else
        {
          n = std::snprintf(&cigar_string[cigar_length],
                       rest,
                       "%c", op);
        }
      if (n < 0)
        {
          fatal("snprintf returned a negative number.\n");
        }
      else if (n >= rest)
        {
          cigar_alloc += std::max(n - rest + 1, minimal_length);
          cigar_string.resize(cigar_alloc);
        }
      else
        {
          cigar_length += n;
          break;
        }
    }
}


auto LinearMemoryAligner::subst_score(char const lhs, char const rhs) -> int64_t
{
  /* return substitution score for replacing char lhs (sequence a),
     with char rhs (sequence b) */
  return scorematrix[map_4bit(rhs)][map_4bit(lhs)];
}


auto LinearMemoryAligner::cigar_add(char _op, int64_t run) -> void
{
  if (op == _op)
    {
      op_run += run;
    }
  else
    {
      cigar_flush();
      op = _op;
      op_run = run;
    }
}


auto LinearMemoryAligner::diff(int64_t a_start,
                               int64_t b_start,
                               int64_t a_len,
                               int64_t b_len,
                               bool gap_b_left,  /* gap open left of b      */
                               bool gap_b_right, /* gap open right of b     */
                               bool a_left,      /* includes left end of a  */
                               bool a_right,     /* includes right end of a */
                               bool b_left,      /* includes left end of b  */
                               bool b_right) -> void  /* includes right end of b */
{
  static constexpr auto int64_min = std::numeric_limits<int64_t>::min();
  // auto span_A = Span{std::next(a_seq, a_start), a_len};
  // auto span_B = Span{std::next(b_seq, b_start), b_len};

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


      int64_t MaxScore = 0;
      int64_t best = 0;

      int64_t Score = 0;

      /* First possibility */

      // Delete 1 from A, Insert b_len from B
      // A----
      // -BBBB

      /* gap penalty for gap in B of length 1 */

      if (not gap_b_left)
        {
          Score -= b_left ? go_t_l : go_t_i;
        }

      Score -= b_left ? ge_t_l : ge_t_i;

      /* gap penalty for gap in A of length b_len */

      Score -= a_right ? go_q_r + (b_len * ge_q_r) : go_q_i + (b_len * ge_q_i);

      MaxScore = Score;
      best = -1;


      /* Second possibility */

      // Insert b_len from B, Delete 1 from A
      // ----A
      // BBBB-

      /* gap penalty for gap in A of length b_len */

      Score -= a_left ? go_q_l + (b_len * ge_q_l) : go_q_i + (b_len * ge_q_i);

      /* gap penalty for gap in B of length 1 */

      if (not gap_b_right)
        {
          Score -= b_right ? go_t_r : go_t_i;
        }

      Score -= b_right ? ge_t_r : ge_t_i;

      if (Score > MaxScore)
        {
          MaxScore = Score;
          best = b_len;
        }


      /* Third possibility */

      for (int64_t i = 0; i < b_len; i++)
        {
          // Insert zero or more from B, replace 1, insert rest of B
          // -A--
          // BBBB

          Score = 0;

          if (i > 0)
            {
              Score -= a_left ? go_q_l + (i * ge_q_l) : go_q_i + (i * ge_q_i);
            }

          Score += subst_score(a_seq[a_start], b_seq[b_start + i]);

          if (i < b_len - 1)
            {
              Score -= a_right ?
                go_q_r + ((b_len - 1 - i) * ge_q_r) :
                go_q_i + ((b_len - 1 - i) * ge_q_i);
            }

          if (Score > MaxScore)
            {
              MaxScore = Score;
              best = i;
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
            {
              cigar_add('I', best);
            }
          cigar_add('M', 1);
          if (best < b_len - 1)
            {
              cigar_add('I', b_len - 1 - best);
            }
        }
    }
  else
    {
      /* a_len >= 2, b_len >= 1 */

      int64_t const I = a_len / 2;  // rename: median?

      // Compute HH & EE in forward phase
      // Upper part

      /* initialize HH and EE for values corresponding to
         empty seq A vs B of i symbols,
         i.e. a gap of length i in A                 */

      HH[0] = 0;
      EE[0] = 0;

      for (int64_t i = 1; i <= b_len; i++)
        {
          HH[i] = - (a_left ? go_q_l + (i * ge_q_l) : go_q_i + (i * ge_q_i));
          EE[i] = int64_min;
        }

      /* compute matrix */

      for (int64_t i = 1; i <= I; i++)
        {
          auto p = HH[0];

          int64_t h = - (b_left ?
                         (gap_b_left ? 0 : go_t_l) + (i * ge_t_l) :
                         (gap_b_left ? 0 : go_t_i) + (i * ge_t_i));

          HH[0] = h;
          auto f = int64_min;

          for (int64_t j = 1; j <= b_len; j++)
            {
              f = std::max(f, h - go_q_i) - ge_q_i;
              if (b_right and (j == b_len))
                {
                  EE[j] = std::max(EE[j], HH[j] - go_t_r) - ge_t_r;
                }
              else
                {
                  EE[j] = std::max(EE[j], HH[j] - go_t_i) - ge_t_i;
                }

              h = p + subst_score(a_seq[a_start + i - 1], b_seq[b_start + j - 1]);

              h = std::max(f, h);
              h = std::max(EE[j], h);
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

      for (int64_t i = 1; i <= b_len; i++)
        {
          XX[i] = - (a_right ? go_q_r + (i * ge_q_r) : go_q_i + (i * ge_q_i));
          YY[i] = int64_min;
        }

      /* compute matrix */

      for (int64_t i = 1; i <= a_len - I; i++)
        {
          auto p = XX[0];

          int64_t h = - (b_right ?
                         (gap_b_right ? 0 : go_t_r) + (i * ge_t_r) :
                         (gap_b_right ? 0 : go_t_i) + (i * ge_t_i));
          XX[0] = h;
          auto f = int64_min;

          for (int64_t j = 1; j <= b_len; j++)
            {
              f = std::max(f, h - go_q_i) - ge_q_i;
              if (b_left and (j == b_len))
                {
                  YY[j] = std::max(YY[j], XX[j] - go_t_l) - ge_t_l;
                }
              else
                {
                  YY[j] = std::max(YY[j], XX[j] - go_t_i) - ge_t_i;
                }

              h = p + subst_score(a_seq[a_start + a_len - i], b_seq[b_start + b_len - j]);

              h = std::max(f, h);
              h = std::max(YY[j], h);
              p = XX[j];
              XX[j] = h;
            }
        }

      YY[0] = XX[0];


      /* find maximum score along division line */

      auto MaxScore0 = int64_min;
      int64_t best0 = -1;

      /* solutions with diagonal at break */

      for (int64_t i = 0; i <= b_len; i++)
        {
          auto const Score = HH[i] + XX[b_len - i];

          if (Score > MaxScore0)
            {
              MaxScore0 = Score;
              best0 = i;
            }
        }

      auto MaxScore1 = int64_min;
      int64_t best1 = -1;

      /* solutions that end with a gap in b from both ends at break */

      for (int64_t i = 0; i <= b_len; i++)
        {
          int64_t g = 0;
          if (b_left and (i == 0))
            {
              g = go_t_l;
            }
          else if (b_right and (i == b_len))
            {
              g = go_t_r;
            }
          else
            {
              g = go_t_i;
            }

          auto const Score = EE[i] + YY[b_len - i] + g;

          if (Score > MaxScore1)
            {
              MaxScore1 = Score;
              best1 = i;
            }
        }

      int64_t P = 0;  // rename: is_parted?? convert to bool?
      int64_t best = 0;

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
               b_left,                b_right and (best == b_len));

          diff(a_start + I,           b_start + best,
               a_len - I,             b_len - best,
               false,                 gap_b_right,
               false,                 a_right,
               b_left and (best == 0), b_right);
        }
      else if (P == 1)
        {
          diff(a_start,               b_start,
               I - 1,                 best,
               gap_b_left,            true,
               a_left,                false,
               b_left,                b_right and (best == b_len));

          cigar_add('D', 2);

          diff(a_start + I + 1,       b_start + best,
               a_len - I - 1,         b_len - best,
               true,                  gap_b_right,
               false,                 a_right,
               b_left and (best == 0), b_right);
        }
    }
}


auto LinearMemoryAligner::align(char * _a_seq,
                                char * _b_seq,
                                int64_t a_len,
                                int64_t b_len) -> char *
{
  /* copy parameters */
  a_seq = _a_seq;
  b_seq = _b_seq;

  /* init cigar operations */
  cigar_reset();

  /* allocate enough memory for vectors */
  alloc_vectors(b_len + 1);

  /* perform alignment */
  diff(0, 0, a_len, b_len, false, false, true, true, true, true);

  /* ensure entire cigar has been written */
  cigar_flush();

  /* return cigar */
  return cigar_string.data();
}


auto LinearMemoryAligner::alignstats(char * cigar,
                                     char * _a_seq,
                                     char * _b_seq,
                                     int64_t * _nwscore,
                                     int64_t * _nwalignmentlength,
                                     int64_t * _nwmatches,
                                     int64_t * _nwmismatches,
                                     int64_t * _nwgaps) -> void
{
  static constexpr auto is_N = 15;  // 4-bit code for 'N' or 'n'
  a_seq = _a_seq;
  b_seq = _b_seq;

  int64_t nwscore = 0;
  int64_t nwalignmentlength = 0;
  int64_t nwmatches = 0;
  int64_t nwmismatches = 0;
  int64_t nwgaps = 0;

  int64_t a_pos = 0;
  int64_t b_pos = 0;

  auto * p = cigar;

  int64_t g = 0;

  while (*p != '\0')
    {
      int64_t runlength = 1;
      auto scanlength = 0;
      std::sscanf(p, "%" PRId64 "%n", &runlength, &scanlength);
      p += scanlength;
      switch (*p++)
        {
        case 'M':
          nwalignmentlength += runlength;
          for (int64_t k = 0; k < runlength; k++)
            {
              auto const a_nuc = a_seq[a_pos];
              auto const b_nuc = b_seq[b_pos];
              nwscore += subst_score(a_nuc, b_nuc);

              if (opt_n_mismatch and ((map_4bit(a_nuc) == is_N) or
                                     (map_4bit(b_nuc) == is_N)))
                {
                  ++nwmismatches;
                }
              else if ((map_4bit(a_nuc) &
                        map_4bit(b_nuc)) != 0U)
                {
                  ++nwmatches;
                }
              else
                {
                  ++nwmismatches;
                }

              ++a_pos;
              ++b_pos;
            }
          break;

        case 'I':
          if ((a_pos == 0) and (b_pos == 0))
            {
              g = go_q_l + (runlength * ge_q_l);
            }
          else if (*p == '\0')  // last operation?
            {
              g = go_q_r + (runlength * ge_q_r);
            }
          else
            {
              g = go_q_i + (runlength * ge_q_i);
            }
          nwscore -= g;
          ++nwgaps;
          nwalignmentlength += runlength;
          b_pos += runlength;
          break;

        case 'D':
          if ((a_pos == 0) and (b_pos == 0))
            {
              g = go_t_l + (runlength * ge_t_l);
            }
          else if (*p == '\0')  // last operation?
            {
              g = go_t_r + (runlength * ge_t_r);
            }
          else
            {
              g = go_t_i + (runlength * ge_t_i);
            }
          nwscore -= g;
          ++nwgaps;
          nwalignmentlength += runlength;
          a_pos += runlength;
          break;
        }  // end of switch
    }  // end of cigar parsing

  *_nwscore = nwscore;
  *_nwalignmentlength = nwalignmentlength;
  *_nwmatches = nwmatches;
  *_nwmismatches = nwmismatches;
  *_nwgaps = nwgaps;
}


// TODO: include guards span.hpp, linmemalign.h  *DONE*
//       scorematrix as vector of vectors? fix scorematrix_create? *DONE*
//       pass nucleotides to subst_score(char const lhs, char const rhs) *DONE*
//       struct scoring as class private member (rename struct members everywhere),
//       inject struct Span in diff(),
//       pass a Pair<nucleotides>
//       design a struct Pair<sequences>?

// struct Pair {
//   Span seq_A;
//   Span seq_B;
// };
