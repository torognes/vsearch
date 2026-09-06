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
#include <cstddef>


/* One character map, its lookup and nothing else -- a translation unit
   includes the maps it names and pays for no other. The table is a
   function-local `static constexpr`: constant-initialized into .rodata, so
   it carries no thread-safe-static guard, and vague-linked through the
   inline accessor, so the 35 translation units that use the maps share one
   copy. See utils/maps/to_uchar.hpp for why it is shaped this way.

   C++17 refactoring: the three blocks below become
   `namespace vsearch::maps::complement`. */
namespace vsearch
{
  namespace maps
  {
    namespace complement
    {

      constexpr std::size_t table_size = 256;


      /* The raw table, for the few loops that hand a whole map to
         std::transform or to a SIMD gather. Reaching for it is the
         exception: it puts the caller back in charge of the index cast,
         which is what map() exists to prevent. */
      inline auto get_map() noexcept
        -> std::array<unsigned char, table_size> const & {
        static constexpr std::array<unsigned char, table_size> table =
          {{
            /*

              Map from ascii to ascii, complementary nucleotide

              @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
              P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
            */

            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

            'N','T','V','G','H','N','N','C','D','N','N','M','N','K','N','N',
            'N','N','Y','S','A','A','B','W','N','R','N','N','N','N','N','N',
            'N','t','v','g','h','N','N','c','d','N','N','m','N','k','n','N',
            'N','N','y','s','a','a','b','w','N','r','N','N','N','N','N','N',

            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
          },};
        return table;
      }


      /* THE entry point -- complementary nucleotide, case preserved.
         Every call site uses it, so the char -> unsigned char cast
         happens in to_uchar() and nowhere else in the tree. */
      inline auto map(char const nucleotide) noexcept -> char {
        return static_cast<char>(get_map()[to_uchar(nucleotide)]);
      }

    }  // namespace complement
  }  // namespace maps
}  // namespace vsearch
