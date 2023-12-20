/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <cstdint>  // int32_t


/* This file contains code dependent on special cpu features. */
/* The file may be compiled several times with different cpu options. */

#ifdef __aarch64__

void increment_counters_from_bitmap(count_t * counters,
                                    unsigned char * bitmap,
                                    unsigned int totalbits)
{
  const uint8x16_t c1 =
    { 0x01, 0x01, 0x02, 0x02, 0x04, 0x04, 0x08, 0x08,
      0x10, 0x10, 0x20, 0x20, 0x40, 0x40, 0x80, 0x80 };

  unsigned short * p = (unsigned short *)(bitmap);
  int16x8_t * q = (int16x8_t *)(counters);
  const auto r = (totalbits + 15) / 16;

  for(auto j = 0U; j < r; j++)
    {
      // load and duplicate short
      uint16x8_t r0 = vdupq_n_u16(*p);
      ++p;

      // cast to bytes
      uint8x16_t r1 = vreinterpretq_u8_u16(r0);

      // bit test with mask giving 0x00 or 0xff
      uint8x16_t r2 = vtstq_u8(r1, c1);

      // transpose to duplicate even bytes
      uint8x16_t r3 = vtrn1q_u8(r2, r2);

      // transpose to duplicate odd bytes
      uint8x16_t r4 = vtrn2q_u8(r2, r2);

      // cast to signed 0x0000 or 0xffff
      int16x8_t r5 = vreinterpretq_s16_u8(r3);

      // cast to signed 0x0000 or 0xffff
      int16x8_t r6 = vreinterpretq_s16_u8(r4);

      // subtract signed 0 or -1 (i.e add 0 or 1) with saturation to counter
      *q = vqsubq_s16(*q, r5);
      ++q;

      // subtract signed 0 or 1 (i.e. add 0 or 1) with saturation to counter
      *q = vqsubq_s16(*q, r6);
      ++q;
    }
}

#elif defined __PPC__

void increment_counters_from_bitmap(count_t * counters,
                                    unsigned char * bitmap,
                                    unsigned int totalbits)
{
  const __vector unsigned char c1 =
    { 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
  const __vector unsigned char c2 =
    { 0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f,
      0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f };
  const __vector unsigned char c3 =
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff };

  unsigned short * p = (unsigned short *)(bitmap);
  __vector signed short * q = (__vector signed short *) (counters);
  const auto r = (totalbits + 15) / 16;

  for(auto j = 0U; j < r; j++)
    {
      __vector unsigned char r0;

      memcpy(&r0, p, 2);
      ++p;
      __vector unsigned char r1 = vec_perm(r0, r0, c1);
      __vector unsigned char r2 = vec_or(r1, c2);
      __vector __bool char r3 = vec_cmpeq(r2, c3);
      __vector signed short r4 = (__vector signed short) vec_unpackl(r3);
      __vector signed short r5 = (__vector signed short) vec_unpackh(r3);
      *q = vec_subs(*q, r4);
      ++q;
      *q = vec_subs(*q, r5);
      ++q;
    }
}

#elif __x86_64__

#include <emmintrin.h>

#ifdef SSSE3
void increment_counters_from_bitmap_ssse3(count_t * counters,
                                          unsigned char * bitmap,
                                          unsigned int totalbits)
#else
void increment_counters_from_bitmap_sse2(count_t * counters,
                                         unsigned char * bitmap,
                                         unsigned int totalbits)
#endif
{
  /*
    Increment selected elements in an array of 16 bit counters.
    The counters to increment are indicated by 1's in the bitmap.

    We read 16 bytes from the bitmap, but use only two bytes (16 bits).
    Convert these 16 bits into 16 bytes with either 0x00 or 0xFF.
    Extend these to 16 words (32 bytes) with either 0x0000 or 0xFFFF.
    Use these values to increment 16 words in an array by subtraction.

    See article below for some hints:
    http://stackoverflow.com/questions/21622212/
    how-to-perform-the-inverse-of-mm256-movemask-epi8-vpmovmskb

    Because the efficient PSHUFB instruction is a SSSE3 instruction
    lacking in many AMD cpus, we provide slightly slower alternative
    SSE2 code.
  */

  // 0xffffffff -> 1111'1111'1111'1111'1111'1111'1111'1111 (32 bits)
  static constexpr auto all_ones = static_cast<int32_t>(0xffffffff);
  // 0x7fbfdfef -> 0111'1111'1011'1111'1101'1111'1110'1111 (32 bits)
  static constexpr auto mask1 = static_cast<int32_t>(0x7fbfdfef);
  // 0xf7fbfdfe -> 1111'0111'1111'1011'1111'1101'1111'1110 (32 bits)
  static constexpr auto mask2 = static_cast<int32_t>(0xf7fbfdfe);

#ifdef SSSE3
  const auto c1 = _mm_set_epi32(0x01010101, 0x01010101, 0x00000000, 0x00000000);
#endif
  const auto c2 = _mm_set_epi32(mask1, mask2, mask1, mask2);
  const auto c3 = _mm_set_epi32(all_ones, all_ones, all_ones, all_ones);

  auto * p = (unsigned short *)(bitmap);
  auto * q = (__m128i *)(counters);
  const auto r = (totalbits + 15) / 16;

  for(auto j = 0U; j < r; j++)
    {
      const auto xmm0 = _mm_loadu_si128((__m128i*)p++);
#ifdef SSSE3
      const auto xmm1 = _mm_shuffle_epi8(xmm0, c1);
#else
      const auto xmm6 = _mm_unpacklo_epi8(xmm0, xmm0);
      const auto xmm7 = _mm_unpacklo_epi16(xmm6, xmm6);
      const auto xmm1 = _mm_unpacklo_epi32(xmm7, xmm7);
#endif
      const auto xmm2 = _mm_or_si128(xmm1, c2);
      const auto xmm3 = _mm_cmpeq_epi8(xmm2, c3);
      const auto xmm4 = _mm_unpacklo_epi8(xmm3, xmm3);
      const auto xmm5 = _mm_unpackhi_epi8(xmm3, xmm3);
      *q = _mm_subs_epi16(*q, xmm4);
      ++q;
      *q = _mm_subs_epi16(*q, xmm5);
      ++q;
    }
}

#else

#error Unknown architecture

#endif
