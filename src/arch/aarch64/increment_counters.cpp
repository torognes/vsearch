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
#include <iterator>  // std::next


// aarch64 backend: NEON intrinsics (arm_neon.h, via arch/intrinsics.hpp). Single
// plain-named variant (no runtime dispatch off x86).
auto increment_counters_from_bitmap(count_t * counters,
                                    unsigned char const * bitmap,
                                    unsigned int const totalbits) -> void
{
  uint8x16_t const bit_selectors =
    { 0x01, 0x01, 0x02, 0x02, 0x04, 0x04, 0x08, 0x08,
      0x10, 0x10, 0x20, 0x20, 0x40, 0x40, 0x80, 0x80 };

  auto const * bits = reinterpret_cast<unsigned short const *>(bitmap);
  auto * counter_vector = reinterpret_cast<int16x8_t *>(counters);
  auto const rounds = (totalbits + counters_per_round - 1) / counters_per_round;

  for (auto round = 0U; round < rounds; round++)
    {
      // load and duplicate short
      auto const bit_word = vdupq_n_u16(*bits);
      bits = std::next(bits);

      // cast to bytes
      auto const bit_bytes = vreinterpretq_u8_u16(bit_word);

      // bit test with mask giving 0x00 or 0xff
      auto const mask = vtstq_u8(bit_bytes, bit_selectors);

      // transpose to duplicate even bytes
      auto const mask_even = vtrn1q_u8(mask, mask);

      // transpose to duplicate odd bytes
      auto const mask_odd = vtrn2q_u8(mask, mask);

      // cast to signed 0x0000 or 0xffff
      auto const mask_low = vreinterpretq_s16_u8(mask_even);

      // cast to signed 0x0000 or 0xffff
      auto const mask_high = vreinterpretq_s16_u8(mask_odd);

      // subtract signed 0 or -1 (i.e add 0 or 1) with saturation to counter
      *counter_vector = vqsubq_s16(*counter_vector, mask_low);
      counter_vector = std::next(counter_vector);

      // subtract signed 0 or 1 (i.e. add 0 or 1) with saturation to counter
      *counter_vector = vqsubq_s16(*counter_vector, mask_high);
      counter_vector = std::next(counter_vector);
    }
}
