/*
 Copyright (C) 2015 Jakob Frielingsdorf

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

 Contact: Jakob Frielingsdorf <jfrielingsdorf@gmail.com>
 */

#include "align_simd_dprofile.h"

#include "align_simd_helper.h"

#include "score_matrix.h"

/*
  Using 16-bit signed values, from -32768 to +32767.
  match: positive
  mismatch: negative
  gap penalties: positive (open, extend, query/target, left/interior/right)
  optimal global alignment (NW)
  maximize score
*/

void dprofile_fill16_aa(CELL * dprofile_word,
                        CELL * score_matrix_word,
                        BYTE * dseq)
{
#if 0
  dumpscorematrix(score_matrix_word);

  for (int j=0; j<CDEPTH; j++)
    {
      for(int z=0; z<CHANNELS; z++)
        fprintf(stderr, " [%c]", sym_aa_5bit[dseq[j*CHANNELS+z]]);
      fprintf(stderr, "\n");
    }
#endif
  __m128i xmm[CHANNELS];
  __m128i xmm_t[CHANNELS];

  for (int j = 0; j<CDEPTH; j++)
    {
//      union { TODO seems to be slower than the loop
//          __m128i v;
//          int16_t a[CHANNELS];
//      } d;
//
//      __m128i tmp = _mm_loadu_si128( (__m128i *) (dseq + (j * CHANNELS)) );
//      tmp = _mm_unpacklo_epi8(tmp, _mm_setzero_si128());
//      _mm_store_si128( &d.v, _mm_slli_epi16( tmp, 4 ) ); // 5 with amino acids

      int d[CHANNELS];
      for (int z = 0; z<CHANNELS; z++)
        d[z] = dseq[j*CHANNELS+z]<<5;

      for (int i = 0; i<32; i += 8)
        {
          for (int x = 0; x<CHANNELS; ++x)
            {
              xmm[x] = _mm_load_si128((__m128i *)(score_matrix_word+d[x]+i));
            }

          for (int x = 0; x<CHANNELS; x += 2)
            {
              xmm_t[x+0] = _mm_unpacklo_epi16(xmm[x+0], xmm[x+1]);
              xmm_t[x+1] = _mm_unpackhi_epi16(xmm[x+0], xmm[x+1]);
            }

          for (int x = 0; x<CHANNELS; x += 4)
            {
              xmm[x+0] = _mm_unpacklo_epi32(xmm_t[x+0], xmm_t[x+2]);
              xmm[x+1] = _mm_unpackhi_epi32(xmm_t[x+0], xmm_t[x+2]);
              xmm[x+2] = _mm_unpacklo_epi32(xmm_t[x+1], xmm_t[x+3]);
              xmm[x+3] = _mm_unpackhi_epi32(xmm_t[x+1], xmm_t[x+3]);
            }

          for (int x = 0; x<(CHANNELS/2); x++)
            {
              xmm_t[(x*2)+0] = _mm_unpacklo_epi64(xmm[x+0], xmm[x+4]);
              xmm_t[(x*2)+1] = _mm_unpackhi_epi64(xmm[x+0], xmm[x+4]);
            }

          for (int x = 0; x<CHANNELS; x++)
            {
              _mm_store_si128(((__m128i *)(dprofile_word)+CDEPTH*(i+x)+j), xmm_t[x]);
            }
        }
    }
#if 0
  dprofile_dump16(dprofile_word);
#endif
}

void dprofile_fill16(CELL * dprofile_word,
                     CELL * score_matrix_word,
                     BYTE * dseq)
{
#if 0
  dumpscorematrix(score_matrix_word);

  for (int j=0; j<CDEPTH; j++)
    {
      for(int z=0; z<CHANNELS; z++)
        fprintf(stderr, " [%c]", sym_nt_4bit[dseq[j*CHANNELS+z]]);
      fprintf(stderr, "\n");
    }
#endif

  __m128i xmm[CHANNELS];
  __m128i xmm_t[CHANNELS];

  for (int j = 0; j<CDEPTH; j++)
    {
//      union { TODO seems to be slower than the loop
//          __m128i v;
//          int16_t a[CHANNELS];
//      } d;
//
//      __m128i tmp = _mm_loadu_si128( (__m128i *) (dseq + (j * CHANNELS)) );
//      tmp = _mm_unpacklo_epi8(tmp, _mm_setzero_si128());
//      _mm_store_si128( &d.v, _mm_slli_epi16( tmp, 4 ) ); // 5 with amino acids

      int d[CHANNELS];
      for (int z = 0; z<CHANNELS; z++)
        d[z] = dseq[j*CHANNELS+z]<<4;

      for (int i = 0; i<16; i += 8) // TODO 16 as a constant makes it faster
        {
          for (int x = 0; x<CHANNELS; ++x)
            {
              xmm[x] = _mm_load_si128((__m128i *)(score_matrix_word+d[x]+i));
            }

          for (int x = 0; x<CHANNELS; x += 2)
            {
              xmm_t[x+0] = _mm_unpacklo_epi16(xmm[x+0], xmm[x+1]);
              xmm_t[x+1] = _mm_unpackhi_epi16(xmm[x+0], xmm[x+1]);
            }

          for (int x = 0; x<CHANNELS; x += 4)
            {
              xmm[x+0] = _mm_unpacklo_epi32(xmm_t[x+0], xmm_t[x+2]);
              xmm[x+1] = _mm_unpackhi_epi32(xmm_t[x+0], xmm_t[x+2]);
              xmm[x+2] = _mm_unpacklo_epi32(xmm_t[x+1], xmm_t[x+3]);
              xmm[x+3] = _mm_unpackhi_epi32(xmm_t[x+1], xmm_t[x+3]);
            }

          for (int x = 0; x<(CHANNELS/2); x++)
            {
              xmm_t[(x*2)+0] = _mm_unpacklo_epi64(xmm[x+0], xmm[x+4]);
              xmm_t[(x*2)+1] = _mm_unpackhi_epi64(xmm[x+0], xmm[x+4]);
            }

          for (int x = 0; x<CHANNELS; x++)
            {
              _mm_store_si128(((__m128i *)(dprofile_word)+CDEPTH*(i+x)+j), xmm_t[x]);
              // CDEPTH * CHANNELS * (i + x) + CHANNELS * j
            }
        }
    }

//  __m128i xmm0,  xmm1,  xmm2,  xmm3,  xmm4,  xmm5,  xmm6,  xmm7;
//  __m128i xmm8,  xmm9,  xmm10, xmm11, xmm12, xmm13, xmm14, xmm15;
//  __m128i xmm16, xmm17, xmm18, xmm19, xmm20, xmm21, xmm22, xmm23;
//  __m128i xmm24, xmm25, xmm26, xmm27, xmm28, xmm29, xmm30, xmm31;
//
//  /* does not require ssse3 */
//  /* approx 4*(5*8+2*40)=480 instructions */
//
//  for (int j=0; j<CDEPTH; j++)
//  {
//    int d[CHANNELS];
//    for(int z=0; z<CHANNELS; z++)
//      d[z] = dseq[j*CHANNELS+z] << 4;
//
//    for(int i=0; i<16; i += 8)
//    {
//      xmm0  = _mm_load_si128((__m128i*)(score_matrix_word + d[0] + i));
//      xmm1  = _mm_load_si128((__m128i*)(score_matrix_word + d[1] + i));
//      xmm2  = _mm_load_si128((__m128i*)(score_matrix_word + d[2] + i));
//      xmm3  = _mm_load_si128((__m128i*)(score_matrix_word + d[3] + i));
//      xmm4  = _mm_load_si128((__m128i*)(score_matrix_word + d[4] + i));
//      xmm5  = _mm_load_si128((__m128i*)(score_matrix_word + d[5] + i));
//      xmm6  = _mm_load_si128((__m128i*)(score_matrix_word + d[6] + i));
//      xmm7  = _mm_load_si128((__m128i*)(score_matrix_word + d[7] + i));
//
//      xmm8  = _mm_unpacklo_epi16(xmm0,  xmm1);
//      xmm9  = _mm_unpackhi_epi16(xmm0,  xmm1);
//      xmm10 = _mm_unpacklo_epi16(xmm2,  xmm3);
//      xmm11 = _mm_unpackhi_epi16(xmm2,  xmm3);
//      xmm12 = _mm_unpacklo_epi16(xmm4,  xmm5);
//      xmm13 = _mm_unpackhi_epi16(xmm4,  xmm5);
//      xmm14 = _mm_unpacklo_epi16(xmm6,  xmm7);
//      xmm15 = _mm_unpackhi_epi16(xmm6,  xmm7);
//
//      xmm16 = _mm_unpacklo_epi32(xmm8,  xmm10);
//      xmm17 = _mm_unpackhi_epi32(xmm8,  xmm10);
//      xmm18 = _mm_unpacklo_epi32(xmm12, xmm14);
//      xmm19 = _mm_unpackhi_epi32(xmm12, xmm14);
//      xmm20 = _mm_unpacklo_epi32(xmm9,  xmm11);
//      xmm21 = _mm_unpackhi_epi32(xmm9,  xmm11);
//      xmm22 = _mm_unpacklo_epi32(xmm13, xmm15);
//      xmm23 = _mm_unpackhi_epi32(xmm13, xmm15);
//
//      xmm24 = _mm_unpacklo_epi64(xmm16, xmm18);
//      xmm25 = _mm_unpackhi_epi64(xmm16, xmm18);
//      xmm26 = _mm_unpacklo_epi64(xmm17, xmm19);
//      xmm27 = _mm_unpackhi_epi64(xmm17, xmm19);
//      xmm28 = _mm_unpacklo_epi64(xmm20, xmm22);
//      xmm29 = _mm_unpackhi_epi64(xmm20, xmm22);
//      xmm30 = _mm_unpacklo_epi64(xmm21, xmm23);
//      xmm31 = _mm_unpackhi_epi64(xmm21, xmm23);
//
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+0) + CHANNELS*j), xmm24);
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+1) + CHANNELS*j), xmm25);
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+2) + CHANNELS*j), xmm26);
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+3) + CHANNELS*j), xmm27);
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+4) + CHANNELS*j), xmm28);
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+5) + CHANNELS*j), xmm29);
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+6) + CHANNELS*j), xmm30);
//      _mm_store_si128((__m128i*)(dprofile_word +
//                                 CDEPTH*CHANNELS*(i+7) + CHANNELS*j), xmm31);
//    }
//  }
#if 0
  dprofile_dump16(dprofile_word);
#endif
}
