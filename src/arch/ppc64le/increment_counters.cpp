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

#include "arch/increment_counters.hpp"
#include "arch/intrinsics.hpp"
#include "vsearch.hpp"
#include <cstring>  // std::memcpy
#include <iterator>  // std::next


// ppc64le backend: AltiVec/VSX intrinsics (altivec.h, via arch/intrinsics.hpp). Single
// plain-named variant (no runtime dispatch off x86).
auto increment_counters_from_bitmap(Span<count_t> const counters,
                                    View<unsigned char> const bitmap) -> void
{
  __vector unsigned char const shuffle_pattern =
    { 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
  __vector unsigned char const bit_selectors =
    { 0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f,
      0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f };
  __vector unsigned char const all_ones =
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff };

  auto const * bits = reinterpret_cast<unsigned short const *>(bitmap.data());
  auto * counter_vector = reinterpret_cast<__vector signed short *>(counters.data());
  auto const rounds = (counters.size() + counters_per_round - 1) / counters_per_round;

  for (auto round = std::size_t{0}; round < rounds; round++)
    {
      __vector unsigned char bit_word;

      std::memcpy(&bit_word, bits, bytes_per_round);
      bits = std::next(bits);
      __vector unsigned char const spread = vec_perm(bit_word, bit_word, shuffle_pattern);
      __vector unsigned char const selected = vec_or(spread, bit_selectors);
      __vector __bool char const mask = vec_cmpeq(selected, all_ones);
      /* vec_unpack* widen the boolean mask to __vector __bool short; the two
         casts below only reinterpret those bits as signed counters, so they
         are reinterpret_cast and not static_cast -- AltiVec vector types of
         different element type have no conversion for static_cast to perform,
         and GCC rejects it outright */
      auto const mask_low = reinterpret_cast<__vector signed short>(vec_unpackl(mask));
      auto const mask_high = reinterpret_cast<__vector signed short>(vec_unpackh(mask));
      *counter_vector = vec_subs(*counter_vector, mask_low);
      counter_vector = std::next(counter_vector);
      *counter_vector = vec_subs(*counter_vector, mask_high);
      counter_vector = std::next(counter_vector);
    }
}
