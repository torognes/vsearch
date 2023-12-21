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
#include <limits>


/*
  Using 16-bit signed values, from -32768 to +32767.
  match: positive
  mismatch: negative
  gap penalties: positive (open, extend, query/target, left/interior/right)
  optimal global alignment (NW)
  maximize score
*/

constexpr auto CHANNELS = 8;
constexpr auto CDEPTH = 4;

/*
  Due to memory usage, limit the product of the length of the sequences.
  If the product of the query length and any target sequence length
  is above the limit, the alignment will not be computed and a score
  of SHRT_MAX will be returned as the score.
  If an overflow occurs during alignment computation, a score of
  SHRT_MAX will also be returned.

  The limit is set to 5 000 * 5 000 = 25 000 000. This will allocate up to
  200 MB per thread. It will align pairs of sequences less than 5000 nt long
  using the SIMD implementation, larger alignments will be performed with
  the linear memory aligner.
*/

constexpr auto MAXSEQLENPRODUCT = 25000000LL;

static int64_t scorematrix[16][16];

/*
  The macros below usually operate on 128-bit vectors of 8 signed
  short 16-bit integers. Additions and subtractions should be
  saturated.  The shift operation should shift left by 2 bytes (one
  short int) and shift in zeros. The v_mask_gt operation should
  compare two vectors of signed shorts and return a 16-bit bitmask
  with pairs of 2 bits set for each element greater in the first than
  in the second argument.
*/

#ifdef __PPC__

typedef __vector signed short VECTOR_SHORT;

const __vector unsigned char perm_merge_long_low =
  {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
   0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17};

const __vector unsigned char perm_merge_long_high =
  {0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
   0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

#define v_init(a,b,c,d,e,f,g,h) (const VECTOR_SHORT){a,b,c,d,e,f,g,h}
#define v_load(a) vec_ld(0, (VECTOR_SHORT *)(a))
#define v_store(a, b) vec_st((__vector unsigned char)(b), 0,    \
                             (__vector unsigned char *)(a))
#define v_add(a, b) vec_adds((a), (b))
#define v_sub(a, b) vec_subs((a), (b))
#define v_sub_unsigned(a, b) ((VECTOR_SHORT)                            \
                              vec_subs((__vector unsigned short) (a),   \
                                       (__vector unsigned short) (b)))
#define v_max(a, b) vec_max((a), (b))
#define v_min(a, b) vec_min((a), (b))
#define v_dup(a) vec_splat((VECTOR_SHORT){(short)(a), 0, 0, 0, 0, 0, 0, 0}, 0);
#define v_zero vec_splat_s16(0)
#define v_and(a, b) vec_and((a), (b))
#define v_xor(a, b) vec_xor((a), (b))
#define v_shift_left(a) vec_sld((a), v_zero, 2)

#elif defined __aarch64__

typedef int16x8_t VECTOR_SHORT;

const uint16x8_t neon_mask =
  {0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000};

#define v_init(a,b,c,d,e,f,g,h) (const VECTOR_SHORT){a,b,c,d,e,f,g,h}
#define v_load(a) vld1q_s16((const int16_t *)(a))
#define v_store(a, b) vst1q_s16((int16_t *)(a), (b))
#define v_merge_lo_16(a, b) vzip1q_s16((a),(b))
#define v_merge_hi_16(a, b) vzip2q_s16((a),(b))
#define v_merge_lo_32(a, b) vreinterpretq_s16_s32(vzip1q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)))
#define v_merge_hi_32(a, b) vreinterpretq_s16_s32(vzip2q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)))
#define v_merge_lo_64(a, b) vreinterpretq_s16_s64(vcombine_s64(vget_low_s64(vreinterpretq_s64_s16(a)), vget_low_s64(vreinterpretq_s64_s16(b))))
#define v_merge_hi_64(a, b) vreinterpretq_s16_s64(vcombine_s64(vget_high_s64(vreinterpretq_s64_s16(a)), vget_high_s64(vreinterpretq_s64_s16(b))))
#define v_add(a, b) vqaddq_s16((a), (b))
#define v_sub(a, b) vqsubq_s16((a), (b))
#define v_sub_unsigned(a, b) vreinterpretq_s16_u16(vqsubq_u16(vreinterpretq_u16_s16(a), vreinterpretq_u16_s16(b)))
#define v_max(a, b) vmaxq_s16((a), (b))
#define v_min(a, b) vminq_s16((a), (b))
#define v_dup(a) vdupq_n_s16(a)
#define v_zero v_dup(0)
#define v_and(a, b) vandq_s16((a), (b))
#define v_xor(a, b) veorq_s16((a), (b))
#define v_shift_left(a) vextq_s16((v_zero), (a), 7)
#define v_mask_gt(a, b) vaddvq_u16(vandq_u16((vcgtq_s16((a), (b))), neon_mask))

#elif __x86_64__

typedef __m128i VECTOR_SHORT;

#define v_init(a,b,c,d,e,f,g,h) _mm_set_epi16(h,g,f,e,d,c,b,a)
#define v_load(a) _mm_load_si128((VECTOR_SHORT *)(a))
#define v_store(a, b) _mm_store_si128((VECTOR_SHORT *)(a), (b))
#define v_merge_lo_16(a, b) _mm_unpacklo_epi16((a),(b))
#define v_merge_hi_16(a, b) _mm_unpackhi_epi16((a),(b))
#define v_merge_lo_32(a, b) _mm_unpacklo_epi32((a),(b))
#define v_merge_hi_32(a, b) _mm_unpackhi_epi32((a),(b))
#define v_merge_lo_64(a, b) _mm_unpacklo_epi64((a),(b))
#define v_merge_hi_64(a, b) _mm_unpackhi_epi64((a),(b))
#define v_add(a, b) _mm_adds_epi16((a), (b))
#define v_sub(a, b) _mm_subs_epi16((a), (b))
#define v_sub_unsigned(a, b) _mm_subs_epu16((a), (b))
#define v_max(a, b) _mm_max_epi16((a), (b))
#define v_min(a, b) _mm_min_epi16((a), (b))
#define v_dup(a) _mm_set1_epi16(a)
#define v_zero v_dup(0)
#define v_and(a, b) _mm_and_si128((a), (b))
#define v_xor(a, b) _mm_xor_si128((a), (b))
#define v_shift_left(a) _mm_slli_si128((a), 2)
#define v_mask_gt(a, b) _mm_movemask_epi8(_mm_cmpgt_epi16((a), (b)))

#else

#error Unknown Architecture

#endif

struct s16info_s
{
  VECTOR_SHORT matrix[32];
  VECTOR_SHORT * hearray;
  VECTOR_SHORT * dprofile;
  VECTOR_SHORT ** qtable;
  unsigned short * dir;
  char * qseq;
  uint64_t diralloc;

  char * cigar;
  char * cigarend;
  int64_t cigaralloc;
  int opcount;
  char op;

  int qlen;
  int maxdlen;
  CELL penalty_gap_open_query_left;
  CELL penalty_gap_open_target_left;
  CELL penalty_gap_open_query_interior;
  CELL penalty_gap_open_target_interior;
  CELL penalty_gap_open_query_right;
  CELL penalty_gap_open_target_right;
  CELL penalty_gap_extension_query_left;
  CELL penalty_gap_extension_target_left;
  CELL penalty_gap_extension_query_interior;
  CELL penalty_gap_extension_target_interior;
  CELL penalty_gap_extension_query_right;
  CELL penalty_gap_extension_target_right;
};

void _mm_print(VECTOR_SHORT x)
{
  auto * y = (unsigned short*)&x;
  for (int i=0; i<8; i++)
    {
      printf("%s%6d", (i>0?" ":""), y[7-i]);
    }
}

void _mm_print2(VECTOR_SHORT x)
{
  auto * y = (signed short*)&x;
  for (int i=0; i<8; i++)
    {
      printf("%s%2d", (i>0?" ":""), y[7-i]);
    }
}

void dprofile_dump16(CELL * dprofile)
{
  char * s = sym_nt_4bit;
  printf("\ndprofile:\n");
  for(int i=0; i<16; i++)
    {
      printf("%c: ",s[i]);
      for(int k=0; k<CDEPTH; k++)
        {
          printf("[");
          for(int j=0; j<CHANNELS; j++)
            {
              printf(" %3d", dprofile[CHANNELS*CDEPTH*i + CHANNELS*k + j]);
            }
          printf("]");
        }
      printf("\n");
    }
}

void dumpscorematrix(CELL * m)
{
  for(int i = 0; i < 16; i++)
    {
      printf("%2d %c", i, sym_nt_4bit[i]);
      for(int j = 0; j < 16; j++)
        {
          printf(" %2d", m[16 * i + j]);
        }
      printf("\n");
    }
}

void dprofile_fill16(CELL * dprofile_word,
                     CELL * score_matrix_word,
                     BYTE * dseq)
{
#if 0
  dumpscorematrix(score_matrix_word);

  for (int j = 0; j < CDEPTH; j++)
    {
      for(int z = 0; z < CHANNELS; z++)
        fprintf(stderr, " [%c]", sym_nt_4bit[dseq[j * CHANNELS + z]]);
      fprintf(stderr, "\n");
    }
#endif

  for (int j = 0; j < CDEPTH; j++)
    {
      int d[CHANNELS];
      for(int z = 0; z < CHANNELS; z++)
        {
          d[z] = dseq[j * CHANNELS + z] << 4;
        }

      for(int i = 0; i < 16; i += 8)
        {

#ifdef __PPC__
          __vector signed short reg0;
          __vector signed short reg1;
          __vector signed short reg2;
          __vector signed short reg3;
          __vector signed short reg4;
          __vector signed short reg5;
          __vector signed short reg6;
          __vector signed short reg7;

          __vector signed int reg8;
          __vector signed int reg9;
          __vector signed int reg10;
          __vector signed int reg11;
          __vector signed int reg12;
          __vector signed int reg13;
          __vector signed int reg14;
          __vector signed int reg15;

          __vector signed long long reg16;
          __vector signed long long reg17;
          __vector signed long long reg18;
          __vector signed long long reg19;
          __vector signed long long reg20;
          __vector signed long long reg21;
          __vector signed long long reg22;
          __vector signed long long reg23;

          __vector signed long long reg24;
          __vector signed long long reg25;
          __vector signed long long reg26;
          __vector signed long long reg27;
          __vector signed long long reg28;
          __vector signed long long reg29;
          __vector signed long long reg30;
          __vector signed long long reg31;
#else
          VECTOR_SHORT reg0;
          VECTOR_SHORT reg1;
          VECTOR_SHORT reg2;
          VECTOR_SHORT reg3;
          VECTOR_SHORT reg4;
          VECTOR_SHORT reg5;
          VECTOR_SHORT reg6;
          VECTOR_SHORT reg7;
          VECTOR_SHORT reg8;
          VECTOR_SHORT reg9;
          VECTOR_SHORT reg10;
          VECTOR_SHORT reg11;
          VECTOR_SHORT reg12;
          VECTOR_SHORT reg13;
          VECTOR_SHORT reg14;
          VECTOR_SHORT reg15;
          VECTOR_SHORT reg16;
          VECTOR_SHORT reg17;
          VECTOR_SHORT reg18;
          VECTOR_SHORT reg19;
          VECTOR_SHORT reg20;
          VECTOR_SHORT reg21;
          VECTOR_SHORT reg22;
          VECTOR_SHORT reg23;
          VECTOR_SHORT reg24;
          VECTOR_SHORT reg25;
          VECTOR_SHORT reg26;
          VECTOR_SHORT reg27;
          VECTOR_SHORT reg28;
          VECTOR_SHORT reg29;
          VECTOR_SHORT reg30;
          VECTOR_SHORT reg31;
#endif

          reg0 = v_load(score_matrix_word + d[0] + i);
          reg1 = v_load(score_matrix_word + d[1] + i);
          reg2 = v_load(score_matrix_word + d[2] + i);
          reg3 = v_load(score_matrix_word + d[3] + i);
          reg4 = v_load(score_matrix_word + d[4] + i);
          reg5 = v_load(score_matrix_word + d[5] + i);
          reg6 = v_load(score_matrix_word + d[6] + i);
          reg7 = v_load(score_matrix_word + d[7] + i);

#ifdef __PPC__
          reg8  = (__vector signed int) vec_mergeh(reg0, reg1);
          reg9  = (__vector signed int) vec_mergel(reg0, reg1);
          reg10 = (__vector signed int) vec_mergeh(reg2, reg3);
          reg11 = (__vector signed int) vec_mergel(reg2, reg3);
          reg12 = (__vector signed int) vec_mergeh(reg4, reg5);
          reg13 = (__vector signed int) vec_mergel(reg4, reg5);
          reg14 = (__vector signed int) vec_mergeh(reg6, reg7);
          reg15 = (__vector signed int) vec_mergel(reg6, reg7);

          reg16 = (__vector signed long long) vec_mergeh(reg8,  reg10);
          reg17 = (__vector signed long long) vec_mergel(reg8,  reg10);
          reg18 = (__vector signed long long) vec_mergeh(reg12, reg14);
          reg19 = (__vector signed long long) vec_mergel(reg12, reg14);
          reg20 = (__vector signed long long) vec_mergeh(reg9,  reg11);
          reg21 = (__vector signed long long) vec_mergel(reg9,  reg11);
          reg22 = (__vector signed long long) vec_mergeh(reg13, reg15);
          reg23 = (__vector signed long long) vec_mergel(reg13, reg15);

          reg24 = (__vector signed long long) vec_perm
            (reg16, reg18, perm_merge_long_low);
          reg25 = (__vector signed long long) vec_perm
            (reg16, reg18, perm_merge_long_high);
          reg26 = (__vector signed long long) vec_perm
            (reg17, reg19, perm_merge_long_low);
          reg27 = (__vector signed long long) vec_perm
            (reg17, reg19, perm_merge_long_high);
          reg28 = (__vector signed long long) vec_perm
            (reg20, reg22, perm_merge_long_low);
          reg29 = (__vector signed long long) vec_perm
            (reg20, reg22, perm_merge_long_high);
          reg30 = (__vector signed long long) vec_perm
            (reg21, reg23, perm_merge_long_low);
          reg31 = (__vector signed long long) vec_perm
            (reg21, reg23, perm_merge_long_high);
#else
          reg8  = v_merge_lo_16(reg0,  reg1);
          reg9  = v_merge_hi_16(reg0,  reg1);
          reg10 = v_merge_lo_16(reg2,  reg3);
          reg11 = v_merge_hi_16(reg2,  reg3);
          reg12 = v_merge_lo_16(reg4,  reg5);
          reg13 = v_merge_hi_16(reg4,  reg5);
          reg14 = v_merge_lo_16(reg6,  reg7);
          reg15 = v_merge_hi_16(reg6,  reg7);

          reg16 = v_merge_lo_32(reg8,  reg10);
          reg17 = v_merge_hi_32(reg8,  reg10);
          reg18 = v_merge_lo_32(reg12, reg14);
          reg19 = v_merge_hi_32(reg12, reg14);
          reg20 = v_merge_lo_32(reg9,  reg11);
          reg21 = v_merge_hi_32(reg9,  reg11);
          reg22 = v_merge_lo_32(reg13, reg15);
          reg23 = v_merge_hi_32(reg13, reg15);

          reg24 = v_merge_lo_64(reg16, reg18);
          reg25 = v_merge_hi_64(reg16, reg18);
          reg26 = v_merge_lo_64(reg17, reg19);
          reg27 = v_merge_hi_64(reg17, reg19);
          reg28 = v_merge_lo_64(reg20, reg22);
          reg29 = v_merge_hi_64(reg20, reg22);
          reg30 = v_merge_lo_64(reg21, reg23);
          reg31 = v_merge_hi_64(reg21, reg23);
#endif

          v_store(dprofile_word + CDEPTH*CHANNELS*(i+0) + CHANNELS*j, reg24);
          v_store(dprofile_word + CDEPTH*CHANNELS*(i+1) + CHANNELS*j, reg25);
          v_store(dprofile_word + CDEPTH*CHANNELS*(i+2) + CHANNELS*j, reg26);
          v_store(dprofile_word + CDEPTH*CHANNELS*(i+3) + CHANNELS*j, reg27);
          v_store(dprofile_word + CDEPTH*CHANNELS*(i+4) + CHANNELS*j, reg28);
          v_store(dprofile_word + CDEPTH*CHANNELS*(i+5) + CHANNELS*j, reg29);
          v_store(dprofile_word + CDEPTH*CHANNELS*(i+6) + CHANNELS*j, reg30);
          v_store(dprofile_word + CDEPTH*CHANNELS*(i+7) + CHANNELS*j, reg31);
        }
    }
#if 0
  dprofile_dump16(dprofile_word);
#endif
}

/*
  The direction bits are set as follows:
  in DIR[0..1] if F>H initially (must go up) (4th pri)
  in DIR[2..3] if E>max(H,F) (must go left) (3rd pri)
  in DIR[4..5] if new F>H (must extend up) (2nd pri)
  in DIR[6..7] if new E>H (must extend left) (1st pri)
  no bits set: go diagonally
*/

/*
  On PPC the fifth parameter is a vector for the result in the lower 64 bits.
  On x86_64 the fifth parameter is the address to write the result to.
*/

#ifdef __PPC__

/* Handle differences between GNU and IBM compilers */

#ifdef __IBMCPP__
#define VECTORBYTEPERMUTE vec_bperm
#else
#define VECTORBYTEPERMUTE vec_vbpermq
#endif

/* The VSX vec_bperm instruction puts the 16 selected bits of the first
   source into bits 48-63 of the destination. */

const __vector unsigned char perm  = { 120, 112, 104,  96,  88,  80,  72,  64,
  56,  48,  40,  32,  24,  16,   8,   0 };

#define ALIGNCORE(H, N, F, V, RES, QR_q, R_q, QR_t, R_t, H_MIN, H_MAX)  \
  {                                                                     \
    __vector unsigned short W, X, Y, Z;                                 \
    __vector unsigned int WX, YZ;                                       \
    __vector short VV;                                                  \
    VV = v_load(&V);                                                    \
    H = v_add(H, VV);                                                   \
    W = (__vector unsigned short) VECTORBYTEPERMUTE                     \
      ((__vector unsigned char) vec_cmpgt(F, H), perm);                 \
    H = v_max(H, F);                                                    \
    X = (__vector unsigned short) VECTORBYTEPERMUTE                     \
      ((__vector unsigned char) vec_cmpgt(E, H), perm);                 \
    H = v_max(H, E);                                                    \
    H_MIN = v_min(H_MIN, H);                                            \
    H_MAX = v_max(H_MAX, H);                                            \
    N = H;                                                              \
    HF = v_sub(H, QR_t);                                                \
    F = v_sub(F, R_t);                                                  \
    Y = (__vector unsigned short) VECTORBYTEPERMUTE                     \
      ((__vector unsigned char) vec_cmpgt(F, HF), perm);                \
    F = v_max(F, HF);                                                   \
    HE = v_sub(H, QR_q);                                                \
    E = v_sub(E, R_q);                                                  \
    Z = (__vector unsigned short) VECTORBYTEPERMUTE                     \
      ((__vector unsigned char) vec_cmpgt(E, HE), perm);                \
    E = v_max(E, HE);                                                   \
    WX = (__vector unsigned int) vec_mergel(W, X);                      \
    YZ = (__vector unsigned int) vec_mergel(Y, Z);                      \
    RES = (__vector unsigned long long) vec_mergeh(WX, YZ);             \
  }

#else

/* x86_64 & aarch64 */

#define ALIGNCORE(H, N, F, V, PATH, QR_q, R_q, QR_t, R_t, H_MIN, H_MAX) \
  H = v_add(H, V);                                                      \
  *((PATH)+0) = v_mask_gt(F, H);                                        \
  (H) = v_max(H, F);                                                    \
  *((PATH)+1) = v_mask_gt(E, H);                                        \
  (H) = v_max(H, E);                                                    \
  (H_MIN) = v_min(H_MIN, H);                                            \
  (H_MAX) = v_max(H_MAX, H);                                            \
  (N) = H;                                                              \
  HF = v_sub(H, QR_t);                                                  \
  (F) = v_sub(F, R_t);                                                  \
  *((PATH)+2) = v_mask_gt(F, HF);                                       \
  (F) = v_max(F, HF);                                                   \
  HE = v_sub(H, QR_q);                                                  \
  E = v_sub(E, R_q);                                                    \
  *((PATH)+3) = v_mask_gt(E, HE);                                       \
  E = v_max(E, HE);

#endif

void aligncolumns_first(VECTOR_SHORT * Sm,
                        VECTOR_SHORT * hep,
                        VECTOR_SHORT ** qp,
                        VECTOR_SHORT QR_q_i,
                        VECTOR_SHORT R_q_i,
                        VECTOR_SHORT QR_q_r,
                        VECTOR_SHORT R_q_r,
                        VECTOR_SHORT QR_t_0,
                        VECTOR_SHORT R_t_0,
                        VECTOR_SHORT QR_t_1,
                        VECTOR_SHORT R_t_1,
                        VECTOR_SHORT QR_t_2,
                        VECTOR_SHORT R_t_2,
                        VECTOR_SHORT QR_t_3,
                        VECTOR_SHORT R_t_3,
                        VECTOR_SHORT h0,
                        VECTOR_SHORT h1,
                        VECTOR_SHORT h2,
                        VECTOR_SHORT h3,
                        VECTOR_SHORT f0,
                        VECTOR_SHORT f1,
                        VECTOR_SHORT f2,
                        VECTOR_SHORT f3,
                        VECTOR_SHORT * _h_min,
                        VECTOR_SHORT * _h_max,
                        VECTOR_SHORT Mm,
                        VECTOR_SHORT M_QR_t_left,
                        VECTOR_SHORT M_R_t_left,
                        VECTOR_SHORT M_QR_q_interior,
                        VECTOR_SHORT M_QR_q_right,
                        int64_t ql,
                        unsigned short * dir)
{

  VECTOR_SHORT h4;
  VECTOR_SHORT h5;
  VECTOR_SHORT h6;
  VECTOR_SHORT h7;
  VECTOR_SHORT h8;
  VECTOR_SHORT E;
  VECTOR_SHORT HE;
  VECTOR_SHORT HF;
  VECTOR_SHORT * vp;

  VECTOR_SHORT h_min = v_zero;
  VECTOR_SHORT h_max = v_zero;

#ifdef __PPC__
  __vector unsigned long long RES1;
  __vector unsigned long long RES2;
  __vector unsigned long long RES;
#endif

  int64_t i;

  f0 = v_sub(f0, QR_t_0);
  f1 = v_sub(f1, QR_t_1);
  f2 = v_sub(f2, QR_t_2);
  f3 = v_sub(f3, QR_t_3);


  for(i=0; i < ql - 1; i++)
    {
      vp = qp[i+0];

      h4 = hep[2*i+0];

      E  = hep[2*i+1];

      /*
        Initialize selected h and e values for next/this round.
        First zero those cells where a new sequence starts
        by using an unsigned saturated subtraction of a huge value to
        set it to zero.
        Then use signed subtraction to obtain the correct value.
      */

      h4 = v_sub_unsigned(h4, Mm);
      h4 = v_sub(h4, M_QR_t_left);

      E  = v_sub_unsigned(E, Mm);
      E  = v_sub(E, M_QR_t_left);
      E  = v_sub(E, M_QR_q_interior);

      M_QR_t_left = v_add(M_QR_t_left, M_R_t_left);

#ifdef __PPC__
      ALIGNCORE(h0, h5, f0, vp[0], RES1,
                QR_q_i, R_q_i, QR_t_0, R_t_0, h_min, h_max);
      ALIGNCORE(h1, h6, f1, vp[1], RES2,
                QR_q_i, R_q_i, QR_t_1, R_t_1, h_min, h_max);
      RES = vec_perm(RES1, RES2, perm_merge_long_low);
      v_store((dir + 16*i + 0), RES);
      ALIGNCORE(h2, h7, f2, vp[2], RES1,
                QR_q_i, R_q_i, QR_t_2, R_t_2, h_min, h_max);
      ALIGNCORE(h3, h8, f3, vp[3], RES2,
                QR_q_i, R_q_i, QR_t_3, R_t_3, h_min, h_max);
      RES = vec_perm(RES1, RES2, perm_merge_long_low);
      v_store((dir + 16*i + 8), RES);
#else
      ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+0,
                QR_q_i, R_q_i, QR_t_0, R_t_0, h_min, h_max);
      ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+4,
                QR_q_i, R_q_i, QR_t_1, R_t_1, h_min, h_max);
      ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+8,
                QR_q_i, R_q_i, QR_t_2, R_t_2, h_min, h_max);
      ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
                QR_q_i, R_q_i, QR_t_3, R_t_3, h_min, h_max);
#endif

      hep[2*i+0] = h8;
      hep[2*i+1] = E;

      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  /* the final round - using query gap penalties for right end */

  vp = qp[i+0];

  E  = hep[2*i+1];

  E  = v_sub_unsigned(E, Mm);
  E  = v_sub(E, M_QR_t_left);
  E  = v_sub(E, M_QR_q_right);


#ifdef __PPC__
  ALIGNCORE(h0, h5, f0, vp[0], RES1,
            QR_q_r, R_q_r, QR_t_0, R_t_0, h_min, h_max);
  ALIGNCORE(h1, h6, f1, vp[1], RES2,
            QR_q_r, R_q_r, QR_t_1, R_t_1, h_min, h_max);
  RES = vec_perm(RES1, RES2, perm_merge_long_low);
  v_store((dir + 16*i + 0), RES);
  ALIGNCORE(h2, h7, f2, vp[2], RES1,
            QR_q_r, R_q_r, QR_t_2, R_t_2, h_min, h_max);
  ALIGNCORE(h3, h8, f3, vp[3], RES2,
            QR_q_r, R_q_r, QR_t_3, R_t_3, h_min, h_max);
  RES = vec_perm(RES1, RES2, perm_merge_long_low);
  v_store((dir + 16*i + 8), RES);
#else
  ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0,
            QR_q_r, R_q_r, QR_t_0, R_t_0, h_min, h_max);
  ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4,
            QR_q_r, R_q_r, QR_t_1, R_t_1, h_min, h_max);
  ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8,
            QR_q_r, R_q_r, QR_t_2, R_t_2, h_min, h_max);
  ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
            QR_q_r, R_q_r, QR_t_3, R_t_3, h_min, h_max);
#endif


  hep[2*i+0] = h8;
  hep[2*i+1] = E;

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;

  *_h_min = h_min;
  *_h_max = h_max;
}

void aligncolumns_rest(VECTOR_SHORT * Sm,
                       VECTOR_SHORT * hep,
                       VECTOR_SHORT ** qp,
                       VECTOR_SHORT QR_q_i,
                       VECTOR_SHORT R_q_i,
                       VECTOR_SHORT QR_q_r,
                       VECTOR_SHORT R_q_r,
                       VECTOR_SHORT QR_t_0,
                       VECTOR_SHORT R_t_0,
                       VECTOR_SHORT QR_t_1,
                       VECTOR_SHORT R_t_1,
                       VECTOR_SHORT QR_t_2,
                       VECTOR_SHORT R_t_2,
                       VECTOR_SHORT QR_t_3,
                       VECTOR_SHORT R_t_3,
                       VECTOR_SHORT h0,
                       VECTOR_SHORT h1,
                       VECTOR_SHORT h2,
                       VECTOR_SHORT h3,
                       VECTOR_SHORT f0,
                       VECTOR_SHORT f1,
                       VECTOR_SHORT f2,
                       VECTOR_SHORT f3,
                       VECTOR_SHORT * _h_min,
                       VECTOR_SHORT * _h_max,
                       int64_t ql,
                       unsigned short * dir)
{
  VECTOR_SHORT h4;
  VECTOR_SHORT h5;
  VECTOR_SHORT h6;
  VECTOR_SHORT h7;
  VECTOR_SHORT h8;
  VECTOR_SHORT E;
  VECTOR_SHORT HE;
  VECTOR_SHORT HF;
  VECTOR_SHORT * vp;

  VECTOR_SHORT h_min = v_zero;
  VECTOR_SHORT h_max = v_zero;

#ifdef __PPC__
  __vector unsigned long long RES1;
  __vector unsigned long long RES2;
  __vector unsigned long long RES;
#endif

  int64_t i;

  f0 = v_sub(f0, QR_t_0);
  f1 = v_sub(f1, QR_t_1);
  f2 = v_sub(f2, QR_t_2);
  f3 = v_sub(f3, QR_t_3);

  for(i=0; i < ql - 1; i++)
    {
      vp = qp[i+0];

      h4 = hep[2*i+0];

      E  = hep[2*i+1];

#ifdef __PPC__
      ALIGNCORE(h0, h5, f0, vp[0], RES1,
                QR_q_i, R_q_i, QR_t_0, R_t_0, h_min, h_max);
      ALIGNCORE(h1, h6, f1, vp[1], RES2,
                QR_q_i, R_q_i, QR_t_1, R_t_1, h_min, h_max);
      RES = vec_perm(RES1, RES2, perm_merge_long_low);
      v_store((dir + 16*i + 0), RES);
      ALIGNCORE(h2, h7, f2, vp[2], RES1,
                QR_q_i, R_q_i, QR_t_2, R_t_2, h_min, h_max);
      ALIGNCORE(h3, h8, f3, vp[3], RES2,
                QR_q_i, R_q_i, QR_t_3, R_t_3, h_min, h_max);
      RES = vec_perm(RES1, RES2, perm_merge_long_low);
      v_store((dir + 16*i + 8), RES);
#else
      ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0,
                QR_q_i, R_q_i, QR_t_0, R_t_0, h_min, h_max);
      ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4,
                QR_q_i, R_q_i, QR_t_1, R_t_1, h_min, h_max);
      ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8,
                QR_q_i, R_q_i, QR_t_2, R_t_2, h_min, h_max);
      ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
                QR_q_i, R_q_i, QR_t_3, R_t_3, h_min, h_max);
#endif

      hep[2*i+0] = h8;
      hep[2*i+1] = E;

      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  /* the final round - using query gap penalties for right end */

  vp = qp[i+0];

  E  = hep[2*i+1];

#ifdef __PPC__
  ALIGNCORE(h0, h5, f0, vp[0], RES1,
            QR_q_r, R_q_r, QR_t_0, R_t_0, h_min, h_max);
  ALIGNCORE(h1, h6, f1, vp[1], RES2,
            QR_q_r, R_q_r, QR_t_1, R_t_1, h_min, h_max);
  RES = vec_perm(RES1, RES2, perm_merge_long_low);
  v_store((dir + 16*i + 0), RES);
  ALIGNCORE(h2, h7, f2, vp[2], RES1,
            QR_q_r, R_q_r, QR_t_2, R_t_2, h_min, h_max);
  ALIGNCORE(h3, h8, f3, vp[3], RES2,
            QR_q_r, R_q_r, QR_t_3, R_t_3, h_min, h_max);
  RES = vec_perm(RES1, RES2, perm_merge_long_low);
  v_store((dir + 16*i + 8), RES);
#else
  ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0,
            QR_q_r, R_q_r, QR_t_0, R_t_0, h_min, h_max);
  ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4,
            QR_q_r, R_q_r, QR_t_1, R_t_1, h_min, h_max);
  ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8,
            QR_q_r, R_q_r, QR_t_2, R_t_2, h_min, h_max);
  ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
            QR_q_r, R_q_r, QR_t_3, R_t_3, h_min, h_max);
#endif

  hep[2*i+0] = h8;
  hep[2*i+1] = E;

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;

  *_h_min = h_min;
  *_h_max = h_max;
}

inline void pushop(s16info_s * s, char newop)
{
  if (newop == s->op)
    {
      s->opcount++;
    }
  else
    {
      *--s->cigarend = s->op;
      if (s->opcount > 1)
        {
          const int size = 11;
          char buf[size];
          int len = snprintf(buf, size, "%d", s->opcount);
          s->cigarend -= len;
          memcpy(s->cigarend, buf, len);
        }
      s->op = newop;
      s->opcount = 1;
    }
}

inline void finishop(s16info_s * s)
{
  if (s->op && s->opcount)
    {
      *--s->cigarend = s->op;
      if (s->opcount > 1)
        {
          const int size = 11;
          char buf[size];
          int len = snprintf(buf, size, "%d", s->opcount);
          s->cigarend -= len;
          memcpy(s->cigarend, buf, len);
        }
      s->op = 0;
      s->opcount = 0;
    }
}

void backtrack16(s16info_s * s,
                 char * dseq,
                 uint64_t dlen,
                 uint64_t offset,
                 uint64_t channel,
                 unsigned short * paligned,
                 unsigned short * pmatches,
                 unsigned short * pmismatches,
                 unsigned short * pgaps)
{
  unsigned short * dirbuffer = s->dir;
  uint64_t dirbuffersize = s->qlen * s->maxdlen * 4;
  uint64_t qlen = s->qlen;
  char * qseq = s->qseq;

  uint64_t maskup      = 3ULL << (2*channel+ 0);
  uint64_t maskleft    = 3ULL << (2*channel+16);
  uint64_t maskextup   = 3ULL << (2*channel+32);
  uint64_t maskextleft = 3ULL << (2*channel+48);

#if 0

  printf("Dumping backtracking array\n");

  for(uint64_t i=0; i<qlen; i++)
    {
      for(uint64_t j=0; j<dlen; j++)
        {
          uint64_t d = *((uint64_t *) (dirbuffer +
                                       (offset + 16*s->qlen*(j/4) +
                                        16*i + 4*(j&3)) % dirbuffersize));
          if (d & maskup)
            {
              if (d & maskleft)
                printf("+");
              else
                printf("^");
            }
          else if (d & maskleft)
            {
              printf("<");
            }
          else
            {
              printf("\\");
            }
        }
      printf("\n");
    }

  printf("Dumping gap extension array\n");

  for(uint64_t i=0; i<qlen; i++)
    {
      for(uint64_t j=0; j<dlen; j++)
        {
          uint64_t d = *((uint64_t *) (dirbuffer +
                                       (offset + 16*s->qlen*(j/4) +
                                        16*i + 4*(j&3)) % dirbuffersize));
          if (d & maskextup)
            {
              if (d & maskextleft)
                printf("+");
              else
                printf("^");
            }
          else if (d & maskextleft)
            {
              printf("<");
            }
          else
            {
              printf("\\");
            }
        }
      printf("\n");
    }

#endif

  unsigned short aligned = 0;
  unsigned short matches = 0;
  unsigned short mismatches = 0;
  unsigned short gaps = 0;

  int64_t i = qlen - 1;
  int64_t j = dlen - 1;

  s->cigarend = s->cigar + s->qlen + s->maxdlen + 1;
  s->op = 0;
  s->opcount = 1;

  while ((i>=0) && (j>=0))
    {
      aligned++;

      uint64_t d = *((uint64_t *) (dirbuffer +
                                   (offset + 16*s->qlen*(j/4) +
                                    16*i + 4*(j&3)) % dirbuffersize));

      if ((s->op == 'I') && (d & maskextleft))
        {
          j--;
          pushop(s, 'I');
        }
      else if ((s->op == 'D') && (d & maskextup))
        {
          i--;
          pushop(s, 'D');
        }
      else if (d & maskleft)
        {
          if (s->op != 'I')
            {
              gaps++;
            }
          j--;
          pushop(s, 'I');
        }
      else if (d & maskup)
        {
          if (s->op != 'D')
            {
              gaps++;
            }
          i--;
          pushop(s, 'D');
        }
      else
        {
          if (chrmap_4bit[(int)(qseq[i])] & chrmap_4bit[(int)(dseq[j])])
            {
              matches++;
            }
          else
            {
              mismatches++;
            }
          i--;
          j--;
          pushop(s, 'M');
        }
    }

  while(i>=0)
    {
      aligned++;
      if (s->op != 'D')
        {
          gaps++;
        }
      i--;
      pushop(s, 'D');
    }

  while(j>=0)
    {
      aligned++;
      if (s->op != 'I')
        {
          gaps++;
        }
      j--;
      pushop(s, 'I');
    }

  finishop(s);

  /* move cigar to beginning of allocated memory area */
  int cigarlen = s->cigar + s->qlen + s->maxdlen - s->cigarend;
  memmove(s->cigar, s->cigarend, cigarlen + 1);

  * paligned = aligned;
  * pmatches = matches;
  * pmismatches = mismatches;
  * pgaps = gaps;
}

struct s16info_s * search16_init(CELL score_match,
                                 CELL score_mismatch,
                                 CELL penalty_gap_open_query_left,
                                 CELL penalty_gap_open_target_left,
                                 CELL penalty_gap_open_query_interior,
                                 CELL penalty_gap_open_target_interior,
                                 CELL penalty_gap_open_query_right,
                                 CELL penalty_gap_open_target_right,
                                 CELL penalty_gap_extension_query_left,
                                 CELL penalty_gap_extension_target_left,
                                 CELL penalty_gap_extension_query_interior,
                                 CELL penalty_gap_extension_target_interior,
                                 CELL penalty_gap_extension_query_right,
                                 CELL penalty_gap_extension_target_right)
{
  (void) score_match;
  (void) score_mismatch;

  /* prepare alloc of qtable, dprofile, hearray, dir */
  auto * s = (struct s16info_s *)
    xmalloc(sizeof(struct s16info_s));

  s->dprofile = (VECTOR_SHORT *) xmalloc(2*4*8*16);
  s->qlen = 0;
  s->qseq = nullptr;
  s->maxdlen = 0;
  s->dir = nullptr;
  s->diralloc = 0;
  s->hearray = nullptr;
  s->qtable = nullptr;
  s->cigar = nullptr;
  s->cigarend = nullptr;
  s->cigaralloc = 0;

  for(int i=0; i<16; i++)
    {
      for(int j=0; j<16; j++)
        {
          CELL value;
          if (ambiguous_4bit[i] or ambiguous_4bit[j])
            {
              value = 0;
            }
          else if (i == j)
            {
              value = opt_match;
            }
          else
            {
              value = opt_mismatch;
            }
          ((CELL*)(&s->matrix))[16*i+j] = value;
          scorematrix[i][j] = value;
        }
    }


  s->penalty_gap_open_query_left =
    penalty_gap_open_query_left;
  s->penalty_gap_open_query_interior =
    penalty_gap_open_query_interior;
  s->penalty_gap_open_query_right =
    penalty_gap_open_query_right;

  s->penalty_gap_open_target_left =
    penalty_gap_open_target_left;
  s->penalty_gap_open_target_interior =
    penalty_gap_open_target_interior;
  s->penalty_gap_open_target_right =
    penalty_gap_open_target_right;

  s->penalty_gap_extension_query_left =
    penalty_gap_extension_query_left;
  s->penalty_gap_extension_query_interior =
    penalty_gap_extension_query_interior;
  s->penalty_gap_extension_query_right =
    penalty_gap_extension_query_right;

  s->penalty_gap_extension_target_left =
    penalty_gap_extension_target_left;
  s->penalty_gap_extension_target_interior =
    penalty_gap_extension_target_interior;
  s->penalty_gap_extension_target_right =
    penalty_gap_extension_target_right;

  return s;
}

void search16_exit(s16info_s * s)
{
  /* free mem for dprofile, hearray, dir, qtable */
  if (s->dir)
    {
      xfree(s->dir);
    }
  if (s->hearray)
    {
      xfree(s->hearray);
    }
  if (s->dprofile)
    {
      xfree(s->dprofile);
    }
  if (s->qtable)
    {
      xfree(s->qtable);
    }
  if (s->cigar)
    {
      xfree(s->cigar);
    }
  xfree(s);
}

void search16_qprep(s16info_s * s, char * qseq, int qlen)
{
  s->qlen = qlen;
  s->qseq = qseq;

  if (s->hearray)
    {
      xfree(s->hearray);
    }
  s->hearray = (VECTOR_SHORT *) xmalloc(2 * s->qlen * sizeof(VECTOR_SHORT));
  memset(s->hearray, 0, 2 * s->qlen * sizeof(VECTOR_SHORT));

  if (s->qtable)
    {
      xfree(s->qtable);
    }
  s->qtable = (VECTOR_SHORT **) xmalloc(s->qlen * sizeof(VECTOR_SHORT*));

  for(int i = 0; i < qlen; i++)
    {
      s->qtable[i] = s->dprofile + 4 * chrmap_4bit[(int)(qseq[i])];
    }
}

void search16(s16info_s * s,
              unsigned int sequences,
              unsigned int * seqnos,
              CELL * pscores,
              unsigned short * paligned,
              unsigned short * pmatches,
              unsigned short * pmismatches,
              unsigned short * pgaps,
              char ** pcigar)
{
  CELL ** q_start = (CELL**) s->qtable;
  CELL * dprofile = (CELL*) s->dprofile;
  CELL * hearray = (CELL*) s->hearray;
  uint64_t qlen = s->qlen;

  if (qlen == 0)
    {
      for (unsigned int cand_id = 0; cand_id < sequences; cand_id++)
        {
          unsigned int seqno = seqnos[cand_id];
          int64_t length = db_getsequencelen(seqno);

          paligned[cand_id] = length;
          pmatches[cand_id] = 0;
          pmismatches[cand_id] = 0;
          pgaps[cand_id] = length;

          if (length == 0)
            {
              pscores[cand_id] = 0;
            }
          else
            {
              pscores[cand_id] =
                MAX(- s->penalty_gap_open_target_left -
                    length * s->penalty_gap_extension_target_left,
                    - s->penalty_gap_open_target_right -
                    length * s->penalty_gap_extension_target_right);
            }

          char * cigar = nullptr;
          if (length > 0)
            {
              int ret = xsprintf(&cigar, "%ldI", length);
              if ((ret < 2) or !cigar)
                {
                  fatal("Unable to allocate enough memory.");
                }
            }
          else
            {
              cigar = (char *) xmalloc(1);
              cigar[0] = 0;
            }
          pcigar[cand_id] = cigar;
        }
      return;
    }

  /* find longest target sequence and reallocate direction buffer */
  uint64_t maxdlen = 0;
  for(int64_t i = 0; i < sequences; i++)
    {
      uint64_t dlen = db_getsequencelen(seqnos[i]);
      /* skip the very long sequences */
      if ((int64_t)(s->qlen) * dlen <= MAXSEQLENPRODUCT)
        {
          if (dlen > maxdlen)
            {
              maxdlen = dlen;
            }
        }
    }
  maxdlen = 4 * ((maxdlen + 3) / 4);
  s->maxdlen = maxdlen;
  uint64_t dirbuffersize = s->qlen * s->maxdlen * 4;

  if (dirbuffersize > s->diralloc)
    {
      s->diralloc = dirbuffersize;
      if (s->dir)
        {
          xfree(s->dir);
        }
      s->dir = (unsigned short*) xmalloc(dirbuffersize *
                                         sizeof(unsigned short));
    }

  unsigned short * dirbuffer = s->dir;

  if (s->qlen + s->maxdlen + 1 > s->cigaralloc)
    {
      s->cigaralloc = s->qlen + s->maxdlen + 1;
      if (s->cigar)
        {
          xfree(s->cigar);
        }
      s->cigar = (char *) xmalloc(s->cigaralloc);
    }

  VECTOR_SHORT M;
  VECTOR_SHORT T0;

  VECTOR_SHORT M_QR_target_left;
  VECTOR_SHORT M_R_target_left;
  VECTOR_SHORT M_QR_query_interior;
  VECTOR_SHORT M_QR_query_right;

  VECTOR_SHORT R_query_left;
  VECTOR_SHORT QR_query_interior;
  VECTOR_SHORT R_query_interior;
  VECTOR_SHORT QR_query_right;
  VECTOR_SHORT R_query_right;
  VECTOR_SHORT QR_target_left;
  VECTOR_SHORT R_target_left;
  VECTOR_SHORT QR_target_interior;
  VECTOR_SHORT R_target_interior;
  VECTOR_SHORT QR_target_right;
  VECTOR_SHORT R_target_right;
  VECTOR_SHORT QR_target[4];
  VECTOR_SHORT R_target[4];

  VECTOR_SHORT *hep, **qp;

  BYTE * d_begin[CHANNELS];
  BYTE * d_end[CHANNELS];
  uint64_t d_offset[CHANNELS];
  BYTE * d_address[CHANNELS];
  uint64_t d_length[CHANNELS];
  int64_t seq_id[CHANNELS];
  bool overflow[CHANNELS];

  VECTOR_SHORT dseqalloc[CDEPTH];
  VECTOR_SHORT S[4];

  BYTE * dseq = (BYTE*) & dseqalloc;
  BYTE zero = 0;

  uint64_t next_id = 0;
  uint64_t done = 0;

  T0 = v_init(-1, 0, 0, 0, 0, 0, 0, 0);

  R_query_left = v_dup(s->penalty_gap_extension_query_left);

  QR_query_interior = v_dup((s->penalty_gap_open_query_interior +
                             s->penalty_gap_extension_query_interior));
  R_query_interior  = v_dup(s->penalty_gap_extension_query_interior);

  QR_query_right  = v_dup((s->penalty_gap_open_query_right +
                           s->penalty_gap_extension_query_right));
  R_query_right  = v_dup(s->penalty_gap_extension_query_right);

  QR_target_left  = v_dup((s->penalty_gap_open_target_left +
                           s->penalty_gap_extension_target_left));
  R_target_left  = v_dup(s->penalty_gap_extension_target_left);

  QR_target_interior = v_dup((s->penalty_gap_open_target_interior +
                              s->penalty_gap_extension_target_interior));
  R_target_interior = v_dup(s->penalty_gap_extension_target_interior);

  QR_target_right  = v_dup((s->penalty_gap_open_target_right +
                            s->penalty_gap_extension_target_right));
  R_target_right  = v_dup(s->penalty_gap_extension_target_right);

  hep = (VECTOR_SHORT*) hearray;
  qp = (VECTOR_SHORT**) q_start;

  for (int c=0; c<CHANNELS; c++)
    {
      d_begin[c] = &zero;
      d_end[c] = d_begin[c];
      d_address[c] = nullptr;
      d_offset[c] = 0;
      d_length[c] = 0;
      seq_id[c] = -1;
      overflow[c] = false;
    }

  short gap_penalty_max = 0;

  gap_penalty_max = MAX(gap_penalty_max,
                        s->penalty_gap_open_query_left +
                        s->penalty_gap_extension_query_left);
  gap_penalty_max = MAX(gap_penalty_max,
                        s->penalty_gap_open_query_interior +
                        s->penalty_gap_extension_query_interior);
  gap_penalty_max = MAX(gap_penalty_max,
                        s->penalty_gap_open_query_right +
                        s->penalty_gap_extension_query_right);
  gap_penalty_max = MAX(gap_penalty_max,
                        s->penalty_gap_open_target_left +
                        s->penalty_gap_extension_target_left);
  gap_penalty_max = MAX(gap_penalty_max,
                        s->penalty_gap_open_target_interior +
                        s->penalty_gap_extension_target_interior);
  gap_penalty_max = MAX(gap_penalty_max,
                        s->penalty_gap_open_target_right +
                        s->penalty_gap_extension_target_right);

  short score_min = std::numeric_limits<short>::min() + gap_penalty_max;
  short score_max = std::numeric_limits<short>::max();

  for(int i=0; i<4; i++)
    {
      S[i] = v_zero;
      dseqalloc[i] = v_zero;
    }

  VECTOR_SHORT H0 = v_zero;
  VECTOR_SHORT H1 = v_zero;
  VECTOR_SHORT H2 = v_zero;
  VECTOR_SHORT H3 = v_zero;

  VECTOR_SHORT F0 = v_zero;
  VECTOR_SHORT F1 = v_zero;
  VECTOR_SHORT F2 = v_zero;
  VECTOR_SHORT F3 = v_zero;

  bool easy = false;

  unsigned short * dir = dirbuffer;

  while(true)
    {
      if (easy)
        {
          /* fill all channels with symbols from the database sequences */

          for(int c = 0; c < CHANNELS; c++)
            {
              for(int j = 0; j < CDEPTH; j++)
                {
                  if (d_begin[c] < d_end[c])
                    {
                      dseq[CHANNELS * j + c] = chrmap_4bit[*(d_begin[c]++)];
                    }
                  else
                    {
                      dseq[CHANNELS * j + c] = 0;
                    }
                }
              if (d_begin[c] == d_end[c])
                {
                  easy = false;
                }
            }

          dprofile_fill16(dprofile, (CELL*) s->matrix, dseq);

          /* create vectors of gap penalties for target depending on whether
             any of the database sequences ended in these four columns */

          if (easy)
            {
              for(unsigned int j = 0; j < CDEPTH; j++)
                {
                  QR_target[j] = QR_target_interior;
                  R_target[j]  = R_target_interior;
                }
            }
          else
            {
              /* one or more sequences ended */

              VECTOR_SHORT QR_diff = v_sub(QR_target_right,
                                           QR_target_interior);
              VECTOR_SHORT R_diff  = v_sub(R_target_right,
                                           R_target_interior);
              for(unsigned int j = 0; j < CDEPTH; j++)
                {
                  VECTOR_SHORT MM = v_zero;
                  VECTOR_SHORT TT = T0;
                  for(int c = 0; c < CHANNELS; c++)
                    {
                      if ((d_begin[c] == d_end[c]) &&
                          (j >= ((d_length[c]+3) % 4)))
                        {
                          MM = v_xor(MM, TT);
                        }
                      TT = v_shift_left(TT);
                    }
                  QR_target[j] = v_add(QR_target_interior,
                                       v_and(QR_diff, MM));
                  R_target[j]  = v_add(R_target_interior,
                                       v_and(R_diff, MM));
                }
            }

          VECTOR_SHORT h_min;
          VECTOR_SHORT h_max;

          aligncolumns_rest(S, hep, qp,
                            QR_query_interior, R_query_interior,
                            QR_query_right, R_query_right,
                            QR_target[0], R_target[0],
                            QR_target[1], R_target[1],
                            QR_target[2], R_target[2],
                            QR_target[3], R_target[3],
                            H0, H1, H2, H3,
                            F0, F1, F2, F3,
                            & h_min, & h_max,
                            qlen, dir);

          VECTOR_SHORT h_min_vector;
          VECTOR_SHORT h_max_vector;
          v_store(& h_min_vector, h_min);
          v_store(& h_max_vector, h_max);
          for(int c = 0; c < CHANNELS; c++)
            {
              if (! overflow[c])
                {
                  signed short h_min_c = ((signed short *)(& h_min_vector))[c];
                  signed short h_max_c = ((signed short *)(& h_max_vector))[c];
                  if ((h_min_c <= score_min) or
                      (h_max_c >= score_max))
                    {
                      overflow[c] = true;
                    }
                }
            }
        }
      else
        {
          /* One or more sequences ended in the previous block.
             We have to switch over to a new sequence           */

          easy = true;

          M = v_zero;

          VECTOR_SHORT T = T0;
          for (int c = 0; c < CHANNELS; c++)
            {
              if (d_begin[c] < d_end[c])
                {
                  /* this channel has more sequence */

                  for(int j = 0; j < CDEPTH; j++)
                    {
                      if (d_begin[c] < d_end[c])
                        {
                          dseq[CHANNELS * j + c] = chrmap_4bit[*(d_begin[c]++)];
                        }
                      else
                        {
                          dseq[CHANNELS * j + c] = 0;
                        }
                    }
                  if (d_begin[c] == d_end[c])
                    {
                      easy = false;
                    }
                }
              else
                {
                  /* sequence in channel c ended. change of sequence */

                  M = v_xor(M, T);

                  int64_t cand_id = seq_id[c];

                  if (cand_id >= 0)
                    {
                      /* save score */

                      char * dbseq = (char*) d_address[c];
                      int64_t dbseqlen = d_length[c];
                      int64_t z = (dbseqlen+3) % 4;
                      int64_t score = ((CELL*)S)[z*CHANNELS+c];

                      if (overflow[c])
                        {
                          pscores[cand_id] = std::numeric_limits<short>::max();
                          paligned[cand_id] = 0;
                          pmatches[cand_id] = 0;
                          pmismatches[cand_id] = 0;
                          pgaps[cand_id] = 0;
                          pcigar[cand_id] = xstrdup("");
                        }
                      else
                        {
                          pscores[cand_id] = score;
                          backtrack16(s, dbseq, dbseqlen, d_offset[c], c,
                                      paligned + cand_id,
                                      pmatches + cand_id,
                                      pmismatches + cand_id,
                                      pgaps + cand_id);
                          pcigar[cand_id] =
                            (char *) xmalloc(strlen(s->cigar)+1);
                          strcpy(pcigar[cand_id], s->cigar);
                        }

                      done++;
                    }

                  /* get next sequence of reasonable length */

                  int64_t length = 0;

                  while ((length == 0) && (next_id < sequences))
                    {
                      cand_id = next_id++;
                      length = db_getsequencelen(seqnos[cand_id]);
                      if ((length == 0) or (s->qlen * length > MAXSEQLENPRODUCT))
                        {
                          pscores[cand_id] = std::numeric_limits<short>::max();
                          paligned[cand_id] = 0;
                          pmatches[cand_id] = 0;
                          pmismatches[cand_id] = 0;
                          pgaps[cand_id] = 0;
                          pcigar[cand_id] = xstrdup("");
                          length = 0;
                          done++;
                        }
                    }

                  if (length > 0)
                    {
                      seq_id[c] = cand_id;
                      char * address = db_getsequence(seqnos[cand_id]);
                      d_address[c] = (BYTE*) address;
                      d_length[c] = length;
                      d_begin[c] = (unsigned char*) address;
                      d_end[c] = (unsigned char*) address + length;
                      d_offset[c] = dir - dirbuffer;
                      overflow[c] = false;

                      ((CELL*)&H0)[c] = 0;
                      ((CELL*)&H1)[c] = - s->penalty_gap_open_query_left
                        - 1*s->penalty_gap_extension_query_left;
                      ((CELL*)&H2)[c] = - s->penalty_gap_open_query_left
                        - 2*s->penalty_gap_extension_query_left;
                      ((CELL*)&H3)[c] = - s->penalty_gap_open_query_left
                        - 3*s->penalty_gap_extension_query_left;

                      ((CELL*)&F0)[c] = - s->penalty_gap_open_query_left
                        - 1*s->penalty_gap_extension_query_left;
                      ((CELL*)&F1)[c] = - s->penalty_gap_open_query_left
                        - 2*s->penalty_gap_extension_query_left;
                      ((CELL*)&F2)[c] = - s->penalty_gap_open_query_left
                        - 3*s->penalty_gap_extension_query_left;
                      ((CELL*)&F3)[c] = - s->penalty_gap_open_query_left
                        - 4*s->penalty_gap_extension_query_left;

                      /* fill channel */

                      for(int j = 0; j < CDEPTH; j++)
                        {
                          if (d_begin[c] < d_end[c])
                            {
                              dseq[CHANNELS * j + c] =
                                chrmap_4bit[*(d_begin[c]++)];
                            }
                          else
                            {
                              dseq[CHANNELS * j + c] = 0;
                            }
                        }
                      if (d_begin[c] == d_end[c])
                        {
                          easy = false;
                        }
                    }
                  else
                    {
                      /* no more sequences, empty channel */

                      seq_id[c] = -1;
                      d_address[c] = nullptr;
                      d_begin[c] = &zero;
                      d_end[c] = d_begin[c];
                      d_length[c] = 0;
                      d_offset[c] = 0;
                      for (int j = 0; j < CDEPTH; j++)
                        {
                          dseq[CHANNELS * j + c] = 0;
                        }
                    }
                }
              T = v_shift_left(T);
            }

          if (done == sequences)
            {
              break;
            }

          /* make masked versions of QR and R for gaps in target */

          M_QR_target_left = v_and(M, QR_target_left);
          M_R_target_left = v_and(M, R_target_left);

          /* make masked versions of QR for gaps in query at target left end */

          M_QR_query_interior = v_and(M, QR_query_interior);
          M_QR_query_right = v_and(M, QR_query_right);

          dprofile_fill16(dprofile, (CELL*) s->matrix, dseq);

          /* create vectors of gap penalties for target depending on whether
             any of the database sequences ended in these four columns */

          if (easy)
            {
              for(unsigned int j = 0; j < CDEPTH; j++)
                {
                  QR_target[j] = QR_target_interior;
                  R_target[j]  = R_target_interior;
                }
            }
          else
            {
              /* one or more sequences ended */

              VECTOR_SHORT QR_diff = v_sub(QR_target_right,
                                           QR_target_interior);
              VECTOR_SHORT R_diff  = v_sub(R_target_right,
                                           R_target_interior);
              for(unsigned int j = 0; j < CDEPTH; j++)
                {
                  VECTOR_SHORT MM = v_zero;
                  VECTOR_SHORT TT = T0;
                  for(int c = 0; c < CHANNELS; c++)
                    {
                      if ((d_begin[c] == d_end[c]) &&
                          (j >= ((d_length[c] + 3) % 4)))
                        {
                          MM = v_xor(MM, TT);
                        }
                      TT = v_shift_left(TT);
                    }
                  QR_target[j] = v_add(QR_target_interior,
                                       v_and(QR_diff, MM));
                  R_target[j]  = v_add(R_target_interior,
                                       v_and(R_diff, MM));
                }
            }

          VECTOR_SHORT h_min;
          VECTOR_SHORT h_max;

          aligncolumns_first(S, hep, qp,
                             QR_query_interior, R_query_interior,
                             QR_query_right, R_query_right,
                             QR_target[0], R_target[0],
                             QR_target[1], R_target[1],
                             QR_target[2], R_target[2],
                             QR_target[3], R_target[3],
                             H0, H1, H2, H3,
                             F0, F1, F2, F3,
                             & h_min, & h_max,
                             M,
                             M_QR_target_left, M_R_target_left,
                             M_QR_query_interior,
                             M_QR_query_right,
                             qlen, dir);

          VECTOR_SHORT h_min_vector;
          VECTOR_SHORT h_max_vector;
          v_store(& h_min_vector, h_min);
          v_store(& h_max_vector, h_max);
          for(int c = 0; c < CHANNELS; c++)
            {
              if (! overflow[c])
                {
                  signed short h_min_c = ((signed short *)(& h_min_vector))[c];
                  signed short h_max_c = ((signed short *)(& h_max_vector))[c];
                  if ((h_min_c <= score_min) or
                      (h_max_c >= score_max))
                    {
                      overflow[c] = true;
                    }
                }
            }
        }

      H0 = v_sub(H3, R_query_left);
      H1 = v_sub(H0, R_query_left);
      H2 = v_sub(H1, R_query_left);
      H3 = v_sub(H2, R_query_left);

      F0 = v_sub(F3, R_query_left);
      F1 = v_sub(F0, R_query_left);
      F2 = v_sub(F1, R_query_left);
      F3 = v_sub(F2, R_query_left);

      dir += 4 * 4 * s->qlen;

      if (dir >= dirbuffer + dirbuffersize)
        {
          dir -= dirbuffersize;
        }
    }
}
