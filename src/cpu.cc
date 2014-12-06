/*
  Copyright (C) 2014 Torbjorn Rognes

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "vsearch.h"

/* This file contains code dependent on special cpu features (e.g. ssse3) */
/* The file will be compiled several times with different cpu options */

#ifdef SSSE3
void increment_counters_from_bitmap_ssse3(unsigned short * counters,
                                          unsigned char * bitmap,
                                          unsigned int totalbits)
#else
void increment_counters_from_bitmap_sse2(unsigned short * counters,
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

#ifdef SSSE3
  const __m128i c1 =
    _mm_set_epi32(0x01010101, 0x01010101, 0x00000000, 0x00000000);
#endif

  const __m128i c2 = 
    _mm_set_epi32(0x7fbfdfef, 0xf7fbfdfe, 0x7fbfdfef, 0xf7fbfdfe);

  const __m128i c3 =
    _mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);

  unsigned short * p = (unsigned short *)(bitmap);
  __m128i * q = (__m128i *)(counters);
  int r = (totalbits + 15) / 16;
   
  for(int j=0; j<r; j++)
    {
      __m128i xmm0, xmm1, xmm2, xmm3, xmm4, xmm5;
      xmm0 = _mm_loadu_si128((__m128i*)p++);
#ifdef SSSE3
      xmm1 = _mm_shuffle_epi8(xmm0, c1);
#else
      __m128i xmm6, xmm7;
      xmm6 = _mm_unpacklo_epi8(xmm0, xmm0);
      xmm7 = _mm_unpacklo_epi16(xmm6, xmm6);
      xmm1 = _mm_unpacklo_epi32(xmm7, xmm7);
#endif
      xmm2 = _mm_or_si128(xmm1, c2);
      xmm3 = _mm_cmpeq_epi8(xmm2, c3);
      xmm4 = _mm_unpacklo_epi8(xmm3, xmm3);
      xmm5 = _mm_unpackhi_epi8(xmm3, xmm3);
      *q = _mm_subs_epi16(*q, xmm4);
      q++;
      *q = _mm_subs_epi16(*q, xmm5);
      q++;
    }
}
