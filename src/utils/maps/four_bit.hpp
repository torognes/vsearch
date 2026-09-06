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

#include "utils/maps/to_uchar.hpp"
#include <array>
#include <cassert>  // assert
#include <cstddef>


/* The 4-bit nucleotide codes and the three relations vsearch builds on
   them. Both tables live here because is_ambiguous() answers a question
   about a code this map produces, not about an ascii byte.

   The table is a function-local `static constexpr`: constant-initialized
   into .rodata, so it carries no thread-safe-static guard, and vague-linked
   through the inline accessor, so the translation units that use it share
   one copy. See utils/maps/to_uchar.hpp for why it is shaped this way.

   C++17 refactoring: the three blocks below become
   `namespace vsearch::maps::four_bit`. */
namespace vsearch
{
  namespace maps
  {
    namespace four_bit
    {

      constexpr std::size_t table_size = 256;

      /* the sixteen codes 0 to 15 the map produces; N is 15 */
      constexpr std::size_t code_count = 16;


      /* The raw table, for the few loops that hand a whole map to
         std::transform or to a SIMD gather. Reaching for it is the
         exception: it puts the caller back in charge of the index cast,
         which is what map() exists to prevent. */
      inline auto get_map() noexcept
        -> std::array<unsigned char, table_size> const & {
        static constexpr std::array<unsigned char, table_size> table =
          {{
            /*
              Map from ascii to 4-bit nucleotide code

              Aa:  1    0001
              Bb: 14    1110   (not A) ex: 'B' & 'A' == 0000 while 'B' & anyother != 0000
              Cc:  2    0010
              Dd: 13    1101   (not C)
              Gg:  4    0100
              Hh: 11    1011   (not G)
              Kk: 12    1100   (G or T)
              Mm:  3    0011   (A or C)
              Nn: 15    1111   ex: 'N' & any != 0000
              Rr:  5    0101   (A or G)
              Ss:  6    0110   (S or G)
              Tt:  8    1000
              Uu:  8    1000
              Vv:  7    0111   (not T)
              Ww:  9    1001   (A or T)
              Yy: 10    1010   (C or T) ex: 'Y' & 'C' or 'T' == 0000
              Others: 0

              @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
              P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
            */

            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  1, 14,  2, 13,  0,  0,  4, 11,  0,  0, 12,  0,  3, 15,  0,
            0,  0,  5,  6,  8,  8,  7,  9,  0, 10,  0,  0,  0,  0,  0,  0,
            0,  1, 14,  2, 13,  0,  0,  4, 11,  0,  0, 12,  0,  3, 15,  0,
            0,  0,  5,  6,  8,  8,  7,  9,  0, 10,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          },};
        return table;
      }


      /* THE entry point -- the 4-bit code of one nucleotide.
         Every call site uses it, so the char -> unsigned char cast
         happens in to_uchar() and nowhere else in the tree. */
      inline auto map(char const nucleotide) noexcept -> unsigned char {
        return get_map()[to_uchar(nucleotide)];
      }


      /* Is this 4-bit code an ambiguity code? Takes a code, not an ascii
         byte, so callers that already mapped do not map twice. */
      inline auto is_ambiguous(unsigned char const code) noexcept -> bool {
        static constexpr std::array<bool, code_count> ambiguous =
          {{
            true,
            false,  // Aa
            false,  // Cc
            true,
            false,  // Gg
            true,
            true,
            true,
            false,  // TtUu
            true,
            true,
            true,
            true,
            true,
            true,
            true,
          },};
        assert(code < code_count);
        return ambiguous[code];
      }


      /* Can these two nucleotides stand for a common base? The codes are
         bit sets, so a non-empty intersection is the answer. */
      inline auto is_equivalent(char const lhs, char const rhs) noexcept
        -> bool {
        auto const lhs_code = map(lhs);
        auto const rhs_code = map(rhs);
        return ((lhs_code & rhs_code) != 0);
      }


      /* Do these two nucleotides have the very same code? Stricter than
         is_equivalent(): R and A are equivalent, but not the same. */
      inline auto is_same(char const lhs, char const rhs) noexcept
        -> bool {
        return map(lhs) == map(rhs);
      }

    }  // namespace four_bit
  }  // namespace maps
}  // namespace vsearch
