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

#include "vsearch.hpp"
#include "core/linmemalign.hpp"
#include "utils/cigar.hpp"  // find_runlength_of_leftmost_operation
#include "utils/decimal_digits.hpp"  // decimal::Buffer, decimal::to_decimal
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/score_4bit.hpp"  // vsearch::score_4bit, SubstitutionScores, nucleotide_codes_4bit
#include "utils/view.hpp"  // View<char>
#include <algorithm>  // std::copy, std::max
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // int64_t
#include <iterator>  // std::next
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


auto scoring_from_options(struct Parameters const & parameters) -> struct Scoring
{
  struct Scoring scoring;
  scoring.match = parameters.opt_match;
  scoring.mismatch = parameters.opt_mismatch;
  scoring.gap_open_query_left = parameters.opt_gap_open_query_left;
  scoring.gap_open_target_left = parameters.opt_gap_open_target_left;
  scoring.gap_open_query_interior = parameters.opt_gap_open_query_interior;
  scoring.gap_open_target_interior = parameters.opt_gap_open_target_interior;
  scoring.gap_open_query_right = parameters.opt_gap_open_query_right;
  scoring.gap_open_target_right = parameters.opt_gap_open_target_right;
  scoring.gap_extension_query_left = parameters.opt_gap_extension_query_left;
  scoring.gap_extension_target_left = parameters.opt_gap_extension_target_left;
  scoring.gap_extension_query_interior = parameters.opt_gap_extension_query_interior;
  scoring.gap_extension_target_interior = parameters.opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = parameters.opt_gap_extension_query_right;
  scoring.gap_extension_target_right = parameters.opt_gap_extension_target_right;
  scoring.n_mismatch = parameters.opt_n_mismatch;
  return scoring;
}


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
      ge_t_r(scoring.gap_extension_target_right),
      n_mismatch(scoring.n_mismatch)
{
  scorematrix_fill(scoring);
}


/*
  Score matrix built below (shown with option N as mismatch):

     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
     -  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
0  - 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
1  A 0  M  X  0  X  0  0  0  X  0  0  0  0  0  0  X
2  C 0  X  M  0  X  0  0  0  X  0  0  0  0  0  0  X
3  M 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
4  G 0  X  X  0  M  0  0  0  X  0  0  0  0  0  0  X
5  R 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
6  S 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
7  V 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
8  T 0  X  X  0  X  0  0  0  M  0  0  0  0  0  0  X
9  W 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
10 Y 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
11 H 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
12 K 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
13 D 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
14 B 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  X
15 N X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X

  M = match, X = mismatch, 0 = neither

  Only the unambiguous nucleotides A, C, G and T/U can match or
  mismatch: any column involving code 0 (no nucleotide) or an ambiguity
  code scores zero, its own diagonal included -- R aligned to R is not a
  match. --n_mismatch overrides that for N alone: every column with an N
  becomes a mismatch, N against N included. Without --n_mismatch the
  last row and column above are all zeros like the other ambiguous ones.
  The SIMD aligner builds the same matrix at 16-bit width
  (search16_init); pairs move between the aligners on size and
  representability grounds alone, so the two must agree cell for cell --
  which they do by construction, both taking their cells from
  vsearch::score_4bit (utils/score_4bit.hpp).
  Checked against raw alignment scores: with match = +2, an identical
  5010-nt pair scores 10020, the same pair with five R/R (or N/N)
  columns scores 10010, and with --n_mismatch five N/N columns score
  9990 (five mismatches at -4).

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
  // fill-in the score matrix
  static_assert(matrix_size == vsearch::nucleotide_codes_4bit.size(),
                "one matrix row and column per 4-bit nucleotide code");
  vsearch::SubstitutionScores<int64_t> const scores {scoring.match, scoring.mismatch,
                                                     scoring.n_mismatch,};
  for (auto const row : vsearch::nucleotide_codes_4bit) {
    for (auto const column : vsearch::nucleotide_codes_4bit) {
      scorematrix[(matrix_size * std::size_t{row}) + column] =
          vsearch::score_4bit(row, column, scores);
    }
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
      cigar_string.resize(static_cast<std::size_t>(cigar_alloc));
    }
  cigar_string[0] = '\0';
  cigar_length = 0;
  op = '\0';
  op_run = 0;
}


auto LinearMemoryAligner::cigar_flush() -> void
{
  if (op_run <= 0) { return; }

  /* A run of length 1 is written as the operation alone, as "%c" did; a longer
     one is preceded by its length, as "%" PRId64 "%c" did. The digits come from
     decimal_digits.hpp, so their width is known before anything is written --
     which is what replaces the grow-and-retry loop this used to be: snprintf
     had to be *called* to discover how much room the number needed. */
  decimal::Buffer digits {};
  auto const run = (op_run > 1) ? decimal::to_decimal(digits, op_run) : View<char>{};
  auto const width = static_cast<int64_t>(run.size()) + 1;  /* + the operation */

  /* One character beyond the run, for the terminator: cigar_length counts the
     characters written, and cigar_reset()/the readers expect cigar_string to
     stay NUL-terminated. */
  if (cigar_length + width + 1 > cigar_alloc)
    {
      cigar_alloc += std::max(cigar_length + width + 1 - cigar_alloc, minimal_length);
      cigar_string.resize(static_cast<std::size_t>(cigar_alloc));
    }

  auto cursor = std::next(cigar_string.begin(), static_cast<std::ptrdiff_t>(cigar_length));
  cursor = std::copy(run.cbegin(), run.cend(), cursor);
  *cursor = op;
  ++cursor;
  *cursor = '\0';
  cigar_length += width;
}


auto LinearMemoryAligner::subst_score(char const lhs, char const rhs) -> int64_t
{
  /* return substitution score for replacing char lhs (sequence a),
     with char rhs (sequence b) */
  return scorematrix[(matrix_size * std::size_t{map_4bit(rhs)}) + map_4bit(lhs)];
}


auto LinearMemoryAligner::cigar_add(char const _op, int64_t const run) -> void
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


auto LinearMemoryAligner::diff(int64_t const a_start,
                               int64_t const b_start,
                               int64_t const a_len,
                               int64_t const b_len,
                               bool const gap_b_left,  /* gap open left of b      */
                               bool const gap_b_right, /* gap open right of b     */
                               bool const a_left,      /* includes left end of a  */
                               bool const a_right,     /* includes right end of a */
                               bool const b_left,      /* includes left end of b  */
                               bool const b_right) -> void  /* includes right end of b */
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

          Score += subst_score(a_seq[static_cast<std::size_t>(a_start)],
                               b_seq[static_cast<std::size_t>(b_start + i)]);

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

      /* Substitution scores for the two fill nests below, with the two
         lookups subst_score performs hoisted to where each is invariant.
         map_4bit is a cross-TU call vsearch cannot inline (no LTO), and
         subst_score runs once per cell here, so the nests fetch the map
         once and translate the A-side nucleotide once per row instead --
         the same hoist reverse_complement uses. subst_score stays as the
         off-the-hot-path spelling; its cell is
         scorematrix[16 * b_code + a_code], which is what the nests index. */
      auto const * const nucleotide_codes = chrmap_4bit();

      // Compute HH & EE in forward phase
      // Upper part

      /* initialize HH and EE for values corresponding to
         empty seq A vs B of i symbols,
         i.e. a gap of length i in A                 */

      HH[0] = 0;
      EE[0] = 0;

      for (int64_t i = 1; i <= b_len; i++)
        {
          HH[static_cast<std::size_t>(i)] = - (a_left ? go_q_l + (i * ge_q_l) : go_q_i + (i * ge_q_i));
          EE[static_cast<std::size_t>(i)] = int64_min;
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

          auto const a_code = nucleotide_codes[static_cast<unsigned char>(
              a_seq[static_cast<std::size_t>(a_start + i - 1)])];

          for (int64_t j = 1; j <= b_len; j++)
            {
              auto const jdx = static_cast<std::size_t>(j);
              f = std::max(f, h - go_q_i) - ge_q_i;
              if (b_right and (j == b_len))
                {
                  EE[jdx] = std::max(EE[jdx], HH[jdx] - go_t_r) - ge_t_r;
                }
              else
                {
                  EE[jdx] = std::max(EE[jdx], HH[jdx] - go_t_i) - ge_t_i;
                }

              auto const b_code = nucleotide_codes[static_cast<unsigned char>(
                  b_seq[static_cast<std::size_t>(b_start + j - 1)])];
              h = p + scorematrix[(matrix_size * std::size_t{b_code}) + a_code];

              h = std::max(f, h);
              h = std::max(EE[jdx], h);
              p = HH[jdx];
              HH[jdx] = h;
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
          XX[static_cast<std::size_t>(i)] = - (a_right ? go_q_r + (i * ge_q_r) : go_q_i + (i * ge_q_i));
          YY[static_cast<std::size_t>(i)] = int64_min;
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

          auto const a_code = nucleotide_codes[static_cast<unsigned char>(
              a_seq[static_cast<std::size_t>(a_start + a_len - i)])];

          for (int64_t j = 1; j <= b_len; j++)
            {
              auto const jdx = static_cast<std::size_t>(j);
              f = std::max(f, h - go_q_i) - ge_q_i;
              if (b_left and (j == b_len))
                {
                  YY[jdx] = std::max(YY[jdx], XX[jdx] - go_t_l) - ge_t_l;
                }
              else
                {
                  YY[jdx] = std::max(YY[jdx], XX[jdx] - go_t_i) - ge_t_i;
                }

              auto const b_code = nucleotide_codes[static_cast<unsigned char>(
                  b_seq[static_cast<std::size_t>(b_start + b_len - j)])];
              h = p + scorematrix[(matrix_size * std::size_t{b_code}) + a_code];

              h = std::max(f, h);
              h = std::max(YY[jdx], h);
              p = XX[jdx];
              XX[jdx] = h;
            }
        }

      YY[0] = XX[0];


      /* find maximum score along division line */

      auto MaxScore0 = int64_min;
      int64_t best0 = -1;

      /* solutions with diagonal at break */

      for (int64_t i = 0; i <= b_len; i++)
        {
          auto const Score = HH[static_cast<std::size_t>(i)] + XX[static_cast<std::size_t>(b_len - i)];

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

          auto const Score = EE[static_cast<std::size_t>(i)] + YY[static_cast<std::size_t>(b_len - i)] + g;

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


auto LinearMemoryAligner::align(View<char> const a_sequence,
                                View<char> const b_sequence) -> char *
{
  /* copy parameters */
  a_seq = a_sequence;
  b_seq = b_sequence;
  auto const a_len = static_cast<int64_t>(a_sequence.size());
  auto const b_len = static_cast<int64_t>(b_sequence.size());

  /* init cigar operations */
  cigar_reset();

  /* allocate enough memory for vectors */
  alloc_vectors(static_cast<std::size_t>(b_len + 1));

  /* perform alignment */
  diff(0, 0, a_len, b_len, false, false, true, true, true, true);

  /* ensure entire cigar has been written */
  cigar_flush();

  /* return cigar */
  return cigar_string.data();
}


auto LinearMemoryAligner::alignstats(char const * cigar,
                                     View<char> const a_sequence,
                                     View<char> const b_sequence) -> AlignStats
{
  static constexpr auto is_N = 15;  // 4-bit code for 'N' or 'n'
  a_seq = a_sequence;
  b_seq = b_sequence;

  int64_t nwscore = 0;
  int64_t nwalignmentlength = 0;
  int64_t nwmatches = 0;
  int64_t nwmismatches = 0;
  int64_t nwgaps = 0;

  int64_t a_pos = 0;
  int64_t b_pos = 0;

  auto const * p = cigar;

  int64_t g = 0;

  while (*p != '\0')
    {
      char const * next_operation = nullptr;
      auto const runlength = find_runlength_of_leftmost_operation(p, next_operation);
      p = next_operation;
      auto const operation = *p;
      p = std::next(p);
      switch (operation)
        {
        case 'M':
          nwalignmentlength += runlength;
          for (int64_t k = 0; k < runlength; k++)
            {
              auto const a_nuc = a_seq[static_cast<std::size_t>(a_pos)];
              auto const b_nuc = b_seq[static_cast<std::size_t>(b_pos)];
              nwscore += subst_score(a_nuc, b_nuc);

              if (n_mismatch and ((map_4bit(a_nuc) == is_N) or
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
        default:
          break;
        }  // end of switch
    }  // end of cigar parsing

  /* not brace-initialized: AlignStats has default member initializers, which
     in C++11 make it a non-aggregate (relaxed only in C++14) */
  AlignStats stats;
  stats.score = nwscore;
  stats.alignmentlength = nwalignmentlength;
  stats.matches = nwmatches;
  stats.mismatches = nwmismatches;
  stats.gaps = nwgaps;
  return stats;
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
