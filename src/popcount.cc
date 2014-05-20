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

/*                                                                                       
   Unable to get the Mac gcc compiler v 4.2.1 to produce the real                           
   popcnt instruction. Therefore resorting to assembly code.                             
*/

#define popcnt_asm(x,y)                                         \
  __asm__ __volatile__ ("popcnt %1,%0" : "=r"(y) : "r"(x));

unsigned long popcount(unsigned long x)
{
  unsigned long y;
  popcnt_asm(x,y);
  return y;
}

void pprint(__m128i x)
{
  unsigned char * p = (unsigned char *) & x;

  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
  printf("%02x", *p++);
}

void pshow(char * name, __m128i x)
{
  printf("%s: ", name);
  pprint(x);
  printf("\n");
}


unsigned long popcount_128(__m128i x)
{

  //  pshow("x", x);

  __m128i mask1 = _mm_set_epi8(0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 
                               0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55);
  
  __m128i mask2 = _mm_set_epi8(0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 
                               0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33);
  
  __m128i mask4 = _mm_set_epi8(0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 
                               0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f);
  
  __m128i zero = _mm_setzero_si128();

  /* add together 2 bits: 0+1, 2+3, 3+4, ... 126+127 */

  __m128i a = _mm_srli_epi64(x, 1);
  __m128i b = _mm_and_si128(x, mask1);
  __m128i c = _mm_and_si128(a, mask1);
  __m128i d = _mm_add_epi64(b, c);

  //  pshow("d", d);

  /* add together 4 bits: (0+1)+(2+3), ... (124+125)+(126+127) */

  __m128i e = _mm_srli_epi64(d, 2);
  __m128i f = _mm_and_si128(d, mask2);
  __m128i g = _mm_and_si128(e, mask2);
  __m128i h = _mm_add_epi64(f, g);

  //  pshow("h", h);

  /* add together 8 bits: (0..3)+(4..7), ... (120..123)+(124..127) */

  __m128i i = _mm_srli_epi64(h, 4);
  __m128i j = _mm_add_epi64(h, i);
  __m128i k = _mm_and_si128(j, mask4);

  //  pshow("k", k);

  /* add together 8 bytes: (0..63) and (64..127) */

  __m128i l = _mm_sad_epu8(k, zero);

  //  pshow("l", l);

  /* add together 64-bit values into final 128 bit value */

  __m128i m = _mm_srli_si128(l, 8);
  __m128i n = _mm_add_epi64(m, l);

  //  pshow("n", n);

  /* return low 64 bits - return value is always in range 0 to 128 */

  unsigned long zz = (unsigned long) _mm_movepi64_pi64(n);

  //  printf("z: %lu\n", zz);

  return zz;
}
