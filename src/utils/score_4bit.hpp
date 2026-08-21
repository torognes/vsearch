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

#include "utils/maps.hpp"  // is_ambiguous_4bit
#include <array>
#include <cassert>


/* The one substitution rule both aligners build their score matrices from:
   search16_init (core/align_simd.cpp) at 16-bit width, and
   LinearMemoryAligner::scorematrix_fill (core/linmemalign.cpp) at 64-bit
   width. Pairs move between those aligners on size and representability
   grounds alone, so the two matrices must agree cell for cell; both taking
   their cells from score_4bit() makes drift impossible. The full matrix the
   rule produces is drawn above scorematrix_fill.

   In namespace vsearch for the same reason as warn(): a new header costs
   nothing to be born namespaced (TBD_20260725_vsearch_namespace.md). */
namespace vsearch
{

  /* The sixteen 4-bit nucleotide codes (utils/maps.hpp): 1, 2, 4 and 8 are
     A, C, G and T/U, 0 maps no accepted nucleotide, and the rest are the
     ambiguity combinations up to N = 15. Code n sits at index n, so a fill
     iterating this list can use its loop variable as a cell index too --
     and since the elements are born unsigned char, no call site needs a
     cast to feed them to score_4bit(). The values cannot be pinned by a
     static_assert, because std::array::operator[] is not constexpr at
     C++11 (libc++ enforces that); what each fill does pin is its matrix
     dimension against this list's size(), which is constexpr. */
  constexpr std::array<unsigned char, 16> nucleotide_codes_4bit =
    {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}};


  /* One aligner's substitution configuration, bundled so the two same-typed
     scores cannot be transposed at a call site. */
  template <typename Cell>
  struct SubstitutionScores
  {
    Cell match;
    Cell mismatch;
    bool n_mismatch;  // any column with an N is a mismatch (--n_mismatch)
  };


  /* Score one aligned column of 4-bit nucleotide codes:

     - an unambiguous nucleotide (A, C, G, T/U) against itself: match,
     - two different unambiguous nucleotides: mismatch,
     - any column involving code 0 or an ambiguity code: zero, its own
       diagonal included -- R against R is not a match,
     - under n_mismatch, any column involving N: mismatch, N against N
       included; this rule wins over the zero rule.

     The two code parameters are interchangeable on purpose: the rule is
     symmetric, so transposing them at a call site cannot change a result. */
  template <typename Cell>
  auto score_4bit(unsigned char const lhs_code,
                  unsigned char const rhs_code,
                  SubstitutionScores<Cell> const & scores) noexcept -> Cell
  {
    constexpr unsigned char code_of_N = 15;
    assert(lhs_code <= code_of_N);  // 4-bit codes index a 16-entry table
    assert(rhs_code <= code_of_N);
    if (scores.n_mismatch and ((lhs_code == code_of_N) or (rhs_code == code_of_N)))
      {
        return scores.mismatch;
      }
    if (is_ambiguous_4bit(lhs_code) or is_ambiguous_4bit(rhs_code))
      {
        return Cell{0};
      }
    return (lhs_code == rhs_code) ? scores.match : scores.mismatch;
  }

}  // namespace vsearch
