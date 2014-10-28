/*
    Copyright (C) 2012-2014 Torbjorn Rognes & Frederic Mahe

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

#define CHANNELS 8
#define CDEPTH 4

static long scorematrix[16][16];

/*
  Using 16-bit signed values, from -32768 to +32767.
  match: positive
  mismatch: negative
  gap penalties: positive (open, extend, query/target, left/interior/right)
  optimal global alignment (NW)
  maximize score
*/

void _mm_print(__m128i x)
{
  unsigned short * y = (unsigned short*)&x;
  for (int i=0; i<8; i++)
    printf("%s%6d", (i>0?" ":""), y[7-i]);
}

void _mm_print2(__m128i x)
{
  signed short * y = (signed short*)&x;
  for (int i=0; i<8; i++)
    printf("%s%2d", (i>0?" ":""), y[7-i]);
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
        printf(" %3d", dprofile[CHANNELS*CDEPTH*i + CHANNELS*k + j]);
      printf("]");
    }
    printf("\n");
  }
}

void dumpscorematrix(CELL * m)
{
  for(int i=0; i<16; i++)
    {
      printf("%2d %c", i, sym_nt_4bit[i]);
      for(int j=0; j<16; j++)
        printf(" %2d", m[16*i+j]);
      printf("\n");
    }
}


void dprofile_fill16(CELL * dprofile_word,
                     CELL * score_matrix_word,
                     BYTE * dseq)
{
  __m128i xmm0,  xmm1,  xmm2,  xmm3,  xmm4,  xmm5,  xmm6,  xmm7;
  __m128i xmm8,  xmm9,  xmm10, xmm11, xmm12, xmm13, xmm14, xmm15;
  __m128i xmm16, xmm17, xmm18, xmm19, xmm20, xmm21, xmm22, xmm23;
  __m128i xmm24, xmm25, xmm26, xmm27, xmm28, xmm29, xmm30, xmm31;

  /* does not require ssse3 */
  /* approx 4*(5*8+2*40)=480 instructions */
  
#if 0
  dumpscorematrix(score_matrix_word);
  
  for (int j=0; j<CDEPTH; j++)
    {
      for(int z=0; z<CHANNELS; z++)
        fprintf(stderr, " [%c]", sym_nt_4bit[dseq[j*CHANNELS+z]]);
      fprintf(stderr, "\n");
    }
#endif

  for (int j=0; j<CDEPTH; j++)
  {
    int d[CHANNELS];
    for(int z=0; z<CHANNELS; z++)
      d[z] = dseq[j*CHANNELS+z] << 4;
      
    for(int i=0; i<16; i += 8)
    {
      xmm0  = _mm_load_si128((__m128i*)(score_matrix_word + d[0] + i));
      xmm1  = _mm_load_si128((__m128i*)(score_matrix_word + d[1] + i));
      xmm2  = _mm_load_si128((__m128i*)(score_matrix_word + d[2] + i));
      xmm3  = _mm_load_si128((__m128i*)(score_matrix_word + d[3] + i));
      xmm4  = _mm_load_si128((__m128i*)(score_matrix_word + d[4] + i));
      xmm5  = _mm_load_si128((__m128i*)(score_matrix_word + d[5] + i));
      xmm6  = _mm_load_si128((__m128i*)(score_matrix_word + d[6] + i));
      xmm7  = _mm_load_si128((__m128i*)(score_matrix_word + d[7] + i));
      
      xmm8  = _mm_unpacklo_epi16(xmm0,  xmm1);
      xmm9  = _mm_unpackhi_epi16(xmm0,  xmm1);
      xmm10 = _mm_unpacklo_epi16(xmm2,  xmm3);
      xmm11 = _mm_unpackhi_epi16(xmm2,  xmm3);
      xmm12 = _mm_unpacklo_epi16(xmm4,  xmm5);
      xmm13 = _mm_unpackhi_epi16(xmm4,  xmm5);
      xmm14 = _mm_unpacklo_epi16(xmm6,  xmm7);
      xmm15 = _mm_unpackhi_epi16(xmm6,  xmm7);
      
      xmm16 = _mm_unpacklo_epi32(xmm8,  xmm10);
      xmm17 = _mm_unpackhi_epi32(xmm8,  xmm10);
      xmm18 = _mm_unpacklo_epi32(xmm12, xmm14);
      xmm19 = _mm_unpackhi_epi32(xmm12, xmm14);
      xmm20 = _mm_unpacklo_epi32(xmm9,  xmm11);
      xmm21 = _mm_unpackhi_epi32(xmm9,  xmm11);
      xmm22 = _mm_unpacklo_epi32(xmm13, xmm15);
      xmm23 = _mm_unpackhi_epi32(xmm13, xmm15);
      
      xmm24 = _mm_unpacklo_epi64(xmm16, xmm18);
      xmm25 = _mm_unpackhi_epi64(xmm16, xmm18);
      xmm26 = _mm_unpacklo_epi64(xmm17, xmm19);
      xmm27 = _mm_unpackhi_epi64(xmm17, xmm19);
      xmm28 = _mm_unpacklo_epi64(xmm20, xmm22);
      xmm29 = _mm_unpackhi_epi64(xmm20, xmm22);
      xmm30 = _mm_unpacklo_epi64(xmm21, xmm23);
      xmm31 = _mm_unpackhi_epi64(xmm21, xmm23);
      
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+0) + CHANNELS*j), xmm24);
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+1) + CHANNELS*j), xmm25);
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+2) + CHANNELS*j), xmm26);
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+3) + CHANNELS*j), xmm27);
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+4) + CHANNELS*j), xmm28);
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+5) + CHANNELS*j), xmm29);
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+6) + CHANNELS*j), xmm30);
      _mm_store_si128((__m128i*)(dprofile_word +
                                 CDEPTH*CHANNELS*(i+7) + CHANNELS*j), xmm31);
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

#define ALIGNCORE(H, N, F, V, PATH, QR_q, R_q, QR_t, R_t)       \
  H = _mm_adds_epi16(H, V);                                     \
  *(PATH+0) = _mm_movemask_epi8(_mm_cmpgt_epi16(F, H));         \
  H = _mm_max_epi16(H, F);                                      \
  *(PATH+1) = _mm_movemask_epi8(_mm_cmpgt_epi16(E, H));         \
  H = _mm_max_epi16(H, E);                                      \
  N = H;                                                        \
  HF = _mm_subs_epi16(H, QR_t);                                 \
  F = _mm_subs_epi16(F, R_t);                                   \
  *(PATH+2) = _mm_movemask_epi8(_mm_cmpgt_epi16(F, HF));        \
  F = _mm_max_epi16(F, HF);                                     \
  HE = _mm_subs_epi16(H, QR_q);                                 \
  E = _mm_subs_epi16(E, R_q);                                   \
  *(PATH+3) = _mm_movemask_epi8(_mm_cmpgt_epi16(E, HE));        \
  E = _mm_max_epi16(E, HE);


void aligncolumns_first(__m128i * Sm,
                        __m128i * hep,
                        __m128i ** qp,
                        __m128i QR_q_i,
                        __m128i R_q_i,
                        __m128i QR_q_r,
                        __m128i R_q_r,
                        __m128i QR_t_0,
                        __m128i R_t_0,
                        __m128i QR_t_1,
                        __m128i R_t_1,
                        __m128i QR_t_2,
                        __m128i R_t_2,
                        __m128i QR_t_3,
                        __m128i R_t_3,
                        __m128i h0,
                        __m128i h1,
                        __m128i h2,
                        __m128i h3,
                        __m128i E0,
                        __m128i Mm,
                        __m128i M_QR_t_left,
                        __m128i M_R_t_left,
                        long ql,
                        unsigned short * dir)
{
  __m128i h4, h5, h6, h7, h8, f0, f1, f2, f3, E, HE, HF;
  __m128i * vp;
  long i;

  f0 = _mm_set1_epi16(SHRT_MIN);
  f1 = f0;
  f2 = f0;
  f3 = f0;

  for(i=0; i < ql - 1; i++)
    {
      vp = qp[i+0];

      h4 = hep[2*i+0];
      h4 = _mm_subs_epu16(h4, Mm);
      h4 = _mm_subs_epi16(h4, M_QR_t_left);
      M_QR_t_left = _mm_adds_epi16(M_QR_t_left, M_R_t_left);

      E  = hep[2*i+1];
      E  = _mm_subs_epu16(E, Mm);
      E  = _mm_adds_epi16(E, E0);

      ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0, QR_q_i, R_q_i, QR_t_0, R_t_0);
      ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4, QR_q_i, R_q_i, QR_t_1, R_t_1);
      ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8, QR_q_i, R_q_i, QR_t_2, R_t_2);
      ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12, QR_q_i, R_q_i, QR_t_3, R_t_3);
      
      hep[2*i+0] = h8;
      hep[2*i+1] = E;

      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  /* the final round - using alternative query gap penalties */
  
  vp = qp[i+0];

  E  = hep[2*i+1];
  E  = _mm_subs_epu16(E, Mm);
  E  = _mm_adds_epi16(E, E0);
  
  ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0, QR_q_r, R_q_r, QR_t_0, R_t_0);
  ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4, QR_q_r, R_q_r, QR_t_1, R_t_1);
  ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8, QR_q_r, R_q_r, QR_t_2, R_t_2);
  ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12, QR_q_r, R_q_r, QR_t_3, R_t_3);

  hep[2*i+0] = h8;
  hep[2*i+1] = E;

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}

void aligncolumns_rest(__m128i * Sm,
                       __m128i * hep,
                       __m128i ** qp,
                       __m128i QR_q_i,
                       __m128i R_q_i,
                       __m128i QR_q_r,
                       __m128i R_q_r,
                       __m128i QR_t_0,
                       __m128i R_t_0,
                       __m128i QR_t_1,
                       __m128i R_t_1,
                       __m128i QR_t_2,
                       __m128i R_t_2,
                       __m128i QR_t_3,
                       __m128i R_t_3,
                       __m128i h0,
                       __m128i h1,
                       __m128i h2,
                       __m128i h3,
                       long ql,
                       unsigned short * dir)
{
  __m128i h4, h5, h6, h7, h8, f0, f1, f2, f3, E, HE, HF;
  __m128i * vp;
  long i;

  f0 = _mm_set1_epi16(SHRT_MIN);
  f1 = f0;
  f2 = f0;
  f3 = f0;

  for(i=0; i < ql - 1; i++)
    {
      vp = qp[i+0];

      h4 = hep[2*i+0];

      E  = hep[2*i+1];

      ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0, QR_q_i, R_q_i, QR_t_0, R_t_0);
      ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4, QR_q_i, R_q_i, QR_t_1, R_t_1);
      ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8, QR_q_i, R_q_i, QR_t_2, R_t_2);
      ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12, QR_q_i, R_q_i, QR_t_3, R_t_3);
      
      hep[2*i+0] = h8;
      hep[2*i+1] = E;

      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  /* the final round - using alternative query gap penalties */
  
  vp = qp[i+0];

  E  = hep[2*i+1];
  
  ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0, QR_q_r, R_q_r, QR_t_0, R_t_0);
  ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4, QR_q_r, R_q_r, QR_t_1, R_t_1);
  ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8, QR_q_r, R_q_r, QR_t_2, R_t_2);
  ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12, QR_q_r, R_q_r, QR_t_3, R_t_3);

  hep[2*i+0] = h8;
  hep[2*i+1] = E;

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}

struct s16info_s
{
  __m128i matrix[32];
  __m128i * hearray;
  __m128i * dprofile;
  __m128i ** qtable;
  unsigned short * dir;
  char * qseq;

  char * cigar;
  char * cigarend;
  int opcount;
  char op;

  int qlen;
  int maxqlen;
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

inline void pushop(s16info_s * s, char newop)
{
  if (newop == s->op)
    s->opcount++;
  else
    {
      *--s->cigarend = s->op;
      if (s->opcount > 1)
        {
          char buf[11];
          int len = sprintf(buf, "%d", s->opcount);
          s->cigarend -= len;
          strncpy(s->cigarend, buf, len);
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
          char buf[11];
          int len = sprintf(buf, "%d", s->opcount);
          s->cigarend -= len;
          strncpy(s->cigarend, buf, len);
        }
      s->op = 0;
      s->opcount = 0;
    }
}

void backtrack16(s16info_s * s,
                 char * dseq,
                 unsigned long dlen,
                 unsigned long offset,
                 unsigned long channel,
                 unsigned short * paligned,
                 unsigned short * pmatches,
                 unsigned short * pmismatches,
                 unsigned short * pgaps)
{
  unsigned short * dirbuffer = s->dir;
  unsigned long dirbuffersize = s->maxqlen * s->maxdlen * 4;
  unsigned long qlen = s->qlen;
  char * qseq = s->qseq;

  unsigned long maskup      = 3UL << (2*channel+ 0);
  unsigned long maskleft    = 3UL << (2*channel+16);
  unsigned long maskextup   = 3UL << (2*channel+32);
  unsigned long maskextleft = 3UL << (2*channel+48);

#if 0

  printf("Dumping backtracking array\n");

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long d = *((unsigned long *) (dirbuffer + 
                                             (offset + 16*s->maxqlen*(j/4) + 
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

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long d = *((unsigned long *) (dirbuffer + 
                                             (offset + 16*s->maxqlen*(j/4) + 
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

  long i = qlen - 1;
  long j = dlen - 1;

  s->cigarend = s->cigar + s->maxqlen + s->maxdlen + 1;
  s->op = 0;
  s->opcount = 1;

  while ((i>=0) && (j>=0))
  {
    aligned++;

    unsigned long d = *((unsigned long *) (dirbuffer + 
                                           (offset + 16*s->maxqlen*(j/4) + 
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
        gaps++;
      j--;
      pushop(s, 'I');
    }
    else if (d & maskup)
    {
      if (s->op != 'D')
        gaps++;
      i--;
      pushop(s, 'D');
    }
    else
    {
      if (chrmap_4bit[(int)(qseq[i])] == chrmap_4bit[(int)(dseq[j])])
        matches++;
      else
        mismatches++;
      i--;
      j--;
      pushop(s, 'M');
    }
  }
  
  while(i>=0)
    {
      aligned++;
      if (s->op != 'D')
        gaps++;
      i--;
      pushop(s, 'D');
    }

  while(j>=0)
    {
      aligned++;
      if (s->op != 'I')
        gaps++;
      j--;
      pushop(s, 'I');
    }

  finishop(s);

  /* move cigar to beginning of allocated memory area */
  int cigarlen = s->cigar + s->maxqlen + s->maxdlen - s->cigarend;
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
  /* prepare alloc of qtable, dprofile, hearray, dir */
  struct s16info_s * s = (struct s16info_s *)
    xmalloc(sizeof(struct s16info_s));

  s->maxdlen = 4 * ((db_getlongestsequence() + 3) / 4);
  s->dprofile = (__m128i *) xmalloc(2*4*8*16);
  s->qlen = 0;
  s->qseq = 0;
  s->maxqlen = 0;
  s->dir = 0;
  s->hearray = 0;
  s->qtable = 0;
  s->cigar = 0;
  s->cigarend = 0;

  for(int i=0; i<16; i++)
    for(int j=0; j<16; j++)
      {
        CELL value;
        if (i==j)
          value = match_score;
        else if ((i==0) || (j==0) || (i>4) || (j>4))
          value = 0;
        else
          value = mismatch_score;
        ((CELL*)(&s->matrix))[16*i+j] = value;
      }
  
  for(int i=0; i<16; i++)
    for(int j=0; j<16; j++)
      {
        CELL value;
        if ((i==0) || (j==0) || (i>4) || (j>4))
          value = 0;
        else if (i==j)
          value = match_score;
        else
          value = mismatch_score;
        scorematrix[i][j] = value;
      }
  
  s->penalty_gap_open_query_left = penalty_gap_open_query_left;
  s->penalty_gap_open_query_interior = penalty_gap_open_query_interior;
  s->penalty_gap_open_query_right = penalty_gap_open_query_right;

  s->penalty_gap_open_target_left = penalty_gap_open_target_left;
  s->penalty_gap_open_target_interior = penalty_gap_open_target_interior;
  s->penalty_gap_open_target_right = penalty_gap_open_target_right;

  s->penalty_gap_extension_query_left = penalty_gap_extension_query_left;
  s->penalty_gap_extension_query_interior = penalty_gap_extension_query_interior;
  s->penalty_gap_extension_query_right = penalty_gap_extension_query_right;

  s->penalty_gap_extension_target_left = penalty_gap_extension_target_left;
  s->penalty_gap_extension_target_interior = penalty_gap_extension_target_interior;
  s->penalty_gap_extension_target_right = penalty_gap_extension_target_right;

  return s;
}

void search16_exit(s16info_s * s)
{
  /* free mem for dprofile, hearray, dir, qtable */
  if (s->dir)
    free(s->dir);
  if (s->hearray)
    free(s->hearray);
  if (s->dprofile)
    free(s->dprofile);
  if (s->qtable)
    free(s->qtable);
  if (s->cigar)
    free(s->cigar);
  free(s);
}

void search16_qprep(s16info_s * s, char * qseq, int qlen)
{
  s->qlen = qlen;
  s->qseq = qseq;

  if (s->qlen > s->maxqlen)
    {
      if (s->maxqlen == 0)
        s->maxqlen = 256;
      while (s->qlen > s->maxqlen)
        s->maxqlen <<= 1;

      if (s->dir)
        free(s->dir);
      s->dir = (unsigned short*) xmalloc(s->maxqlen * s->maxdlen *
                                         4 * sizeof(unsigned short));

      if (s->hearray)
        free(s->hearray);
      s->hearray = (__m128i *) xmalloc(2 * s->maxqlen * sizeof(__m128i));
      memset(s->hearray, 0, 2 * s->maxqlen * sizeof(__m128i));

      if (s->qtable)
        free(s->qtable);
      s->qtable = (__m128i **) xmalloc(s->maxqlen * sizeof(__m128i*));

      if (s->cigar)
        free(s->cigar);
      s->cigar = (char *) xmalloc(s->maxqlen + s->maxdlen + 1);
    }

  for(int i = 0; i < qlen; i++)
    s->qtable[i] = s->dprofile + 4 * chrmap_4bit[(int)(qseq[i])];
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
  unsigned long qlen = s->qlen;
  unsigned long dirbuffersize = s->maxqlen * s->maxdlen * 4;
  unsigned short * dirbuffer = s->dir;
  
  __m128i T, M, T0, E0;

  __m128i M_QR_target_left, M_R_target_left;

  __m128i QR_query_left, R_query_left;
  __m128i QR_query_interior, R_query_interior;
  __m128i QR_query_right, R_query_right;
  __m128i QR_target_left, R_target_left;
  __m128i QR_target_interior, R_target_interior;
  __m128i QR_target_right, R_target_right;
  __m128i QR_target[4], R_target[4];
  
  __m128i *hep, **qp;

  BYTE * d_begin[CHANNELS];
  BYTE * d_end[CHANNELS];
  unsigned long d_offset[CHANNELS];
  BYTE * d_address[CHANNELS];
  unsigned long d_length[CHANNELS];
  long seq_id[CHANNELS];
  
  __m128i dseqalloc[CDEPTH];
  __m128i S[4];

  BYTE * dseq = (BYTE*) & dseqalloc;
  BYTE zero = 0;

  unsigned long next_id = 0;
  unsigned long done = 0;
  
  T0 = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 0xffff);

  QR_query_left = _mm_set1_epi16(s->penalty_gap_open_query_left + 
                                 s->penalty_gap_extension_query_left);
  R_query_left = _mm_set1_epi16(s->penalty_gap_extension_query_left);
  
  QR_query_interior = _mm_set1_epi16(s->penalty_gap_open_query_interior + 
                                     s->penalty_gap_extension_query_interior);
  R_query_interior  = _mm_set1_epi16(s->penalty_gap_extension_query_interior);

  QR_query_right  = _mm_set1_epi16(s->penalty_gap_open_query_right + 
                                   s->penalty_gap_extension_query_right);
  R_query_right  = _mm_set1_epi16(s->penalty_gap_extension_query_right);
  
  QR_target_left  = _mm_set1_epi16(s->penalty_gap_open_target_left + 
                                   s->penalty_gap_extension_target_left);
  R_target_left  = _mm_set1_epi16(s->penalty_gap_extension_target_left);
  
  QR_target_interior = _mm_set1_epi16(s->penalty_gap_open_target_interior + 
                                     s->penalty_gap_extension_target_interior);
  R_target_interior = _mm_set1_epi16(s->penalty_gap_extension_target_interior);
  
  QR_target_right  = _mm_set1_epi16(s->penalty_gap_open_target_right + 
                                   s->penalty_gap_extension_target_right);
  R_target_right  = _mm_set1_epi16(s->penalty_gap_extension_target_right);
  
  hep = (__m128i*) hearray;
  qp = (__m128i**) q_start;

  for (int c=0; c<CHANNELS; c++)
    {
      d_begin[c] = &zero;
      d_end[c] = d_begin[c];
      d_address[c] = 0;
      d_offset[c] = 0;
      d_length[c] = 0;
      seq_id[c] = -1;
    }
  
  for(int i=0; i<4; i++)
    {
      S[i] = _mm_setzero_si128();
      dseqalloc[i] = _mm_setzero_si128();
    }
  
  __m128i VECTOR_SHRT_MIN = _mm_set1_epi16(SHRT_MIN);
  __m128i H0 = _mm_setzero_si128();
  __m128i H1 = _mm_setzero_si128();
  __m128i H2 = _mm_setzero_si128();
  __m128i H3 = _mm_setzero_si128();
  
  int easy = 0;

  unsigned short * dir = dirbuffer;

  
  while(1)
  {
    if (easy)
    {
      /* fill all channels with symbols from the database sequences */

      for(int c=0; c<CHANNELS; c++)
      {
        for(int j=0; j<CDEPTH; j++)
        {
          if (d_begin[c] < d_end[c])
            dseq[CHANNELS*j+c] = chrmap_4bit[*(d_begin[c]++)];
          else
            dseq[CHANNELS*j+c] = 0;
        }
        if (d_begin[c] == d_end[c])
          easy = 0;
      }

      dprofile_fill16(dprofile, (CELL*) s->matrix, dseq);

      /* create vectors of gap penalties for target depending on whether
         any of the database sequences ended in these four columns */

      if (easy)
        {
          for(unsigned int j=0; j<CDEPTH; j++)
            {
              QR_target[j] = QR_target_interior;
              R_target[j]  = R_target_interior;
            }
        }
      else
        {
          /* one or more sequences ended */
          __m128i QR_diff = _mm_subs_epi16(QR_target_right,
                                           QR_target_interior);
          __m128i R_diff  = _mm_subs_epi16(R_target_right,
                                           R_target_interior);
          for(unsigned int j=0; j<CDEPTH; j++)
            {
              __m128i M = _mm_setzero_si128();
              __m128i T = T0;
              for(int c=0; c<CHANNELS; c++)
                {
                  if ((d_begin[c] == d_end[c]) &&
                      (j >= ((d_length[c]+3) % 4)))
                    {
                      M = _mm_xor_si128(M, T);
                    }
                  T = _mm_slli_si128(T, 2);
                }
              QR_target[j] = _mm_adds_epi16(QR_target_interior, 
                                           _mm_and_si128(QR_diff, M));
              R_target[j]  = _mm_adds_epi16(R_target_interior,
                                           _mm_and_si128(R_diff, M));
            }
        }
      
      aligncolumns_rest(S, hep, qp,
                        QR_query_interior, R_query_interior, 
                        QR_query_right, R_query_right, 
                        QR_target[0], R_target[0],
                        QR_target[1], R_target[1],
                        QR_target[2], R_target[2],
                        QR_target[3], R_target[3],
                        H0, H1, H2, H3,
                        qlen, dir);
    }
    else
    {
      /* One or more sequences ended in the previous block.
         We have to switch over to a new sequence           */

      easy = 1;
      
      M = _mm_setzero_si128();
      T = T0;
      for (int c=0; c<CHANNELS; c++)
      {
        if (d_begin[c] < d_end[c])
        {
          /* this channel has more sequence */

          for(int j=0; j<CDEPTH; j++)
          {
            if (d_begin[c] < d_end[c])
              dseq[CHANNELS*j+c] = chrmap_4bit[*(d_begin[c]++)];
            else
              dseq[CHANNELS*j+c] = 0;
          }
          if (d_begin[c] == d_end[c])
            easy = 0;
        }
        else
        {
          /* sequence in channel c ended. change of sequence */

          M = _mm_xor_si128(M, T);

          long cand_id = seq_id[c];
          
          if (cand_id >= 0)
          {
            /* save score */

            char * dbseq = (char*) d_address[c];
            long dbseqlen = d_length[c];
            long z = (dbseqlen+3) % 4;
            long score = ((CELL*)S)[z*CHANNELS+c];

            pscores[cand_id] = score;
            
            if ((score > SHRT_MIN) && (score < SHRT_MAX))
              {
                backtrack16(s, dbseq, dbseqlen, d_offset[c], c,
                            paligned + cand_id,
                            pmatches + cand_id,
                            pmismatches + cand_id,
                            pgaps + cand_id);
                pcigar[cand_id] = (char *) xmalloc(strlen(s->cigar)+1);
                strcpy(pcigar[cand_id], s->cigar);
              }
            else
              {
                paligned[cand_id] = 0;
                pmatches[cand_id] = 0;
                pmismatches[cand_id] = 0;
                pgaps[cand_id] = 0;
                pcigar[cand_id] = strdup("");
              }

            done++;
          }

          if (next_id < sequences)
          {

            long seqno;
            char* address;
            long length = 0;

            /* get next sequence with length>0 */
            
            while ((length==0) && (next_id < sequences))
              {
                seq_id[c] = next_id;
                seqno = seqnos[next_id];

                address = db_getsequence(seqno);
                length = db_getsequencelen(seqno);

                next_id++;
              }

            d_address[c] = (BYTE*) address;
            d_length[c] = length;
            d_begin[c] = (unsigned char*) address;
            d_end[c] = (unsigned char*) address + length;
            d_offset[c] = dir - dirbuffer;
            
            ((CELL*)&H0)[c] = 0;
            ((CELL*)&H1)[c] = - s->penalty_gap_open_query_left
              - 1*s->penalty_gap_extension_query_left;
            ((CELL*)&H2)[c] = - s->penalty_gap_open_query_left
              - 2*s->penalty_gap_extension_query_left;
            ((CELL*)&H3)[c] = - s->penalty_gap_open_query_left
              - 3*s->penalty_gap_extension_query_left;
            
            /* fill channel */

            for(int j=0; j<CDEPTH; j++)
            {
              if (d_begin[c] < d_end[c])
                dseq[CHANNELS*j+c] = chrmap_4bit[*(d_begin[c]++)];
              else
                dseq[CHANNELS*j+c] = 0;
            }
            if (d_begin[c] == d_end[c])
              easy = 0;
          }
          else
          {
            /* no more sequences, empty channel */

            seq_id[c] = -1;
            d_address[c] = 0;
            d_begin[c] = &zero;
            d_end[c] = d_begin[c];
            d_length[c] = 0;
            d_offset[c] = 0;
            for (int j=0; j<CDEPTH; j++)
              dseq[CHANNELS*j+c] = 0;
          }
        }

        T = _mm_slli_si128(T, 2);
      }

      if (done == sequences)
        break;
          
      /* make masked versions of QR, R and E0 */

      M_QR_target_left = _mm_and_si128(M, QR_target_left);
      M_R_target_left = _mm_and_si128(M, R_target_left);
      E0 = _mm_and_si128(M, VECTOR_SHRT_MIN);
      
      dprofile_fill16(dprofile, (CELL*) s->matrix, dseq);
      
      /* create vectors of gap penalties for target depending on whether
         any of the database sequences ended in these four columns */

      if (easy)
        {
          for(unsigned int j=0; j<CDEPTH; j++)
            {
              QR_target[j] = QR_target_interior;
              R_target[j]  = R_target_interior;
            }
        }
      else
        {
          /* one or more sequences ended */
          __m128i QR_diff = _mm_subs_epi16(QR_target_right,
                                           QR_target_interior);
          __m128i R_diff  = _mm_subs_epi16(R_target_right,
                                           R_target_interior);
          for(unsigned int j=0; j<CDEPTH; j++)
            {
              __m128i M = _mm_setzero_si128();
              __m128i T = T0;
              for(int c=0; c<CHANNELS; c++)
                {
                  if ((d_begin[c] == d_end[c]) &&
                      (j >= ((d_length[c]+3) % 4)))
                    {
                      M = _mm_xor_si128(M, T);
                    }
                  T = _mm_slli_si128(T, 2);
                }
              QR_target[j] = _mm_adds_epi16(QR_target_interior, 
                                           _mm_and_si128(QR_diff, M));
              R_target[j]  = _mm_adds_epi16(R_target_interior,
                                           _mm_and_si128(R_diff, M));
            }
        }
      
      aligncolumns_first(S, hep, qp, 
                         QR_query_interior, R_query_interior, 
                         QR_query_right, R_query_right, 
                         QR_target[0], R_target[0],
                         QR_target[1], R_target[1],
                         QR_target[2], R_target[2],
                         QR_target[3], R_target[3],
                         H0, H1, H2, H3,
                         E0,
                         M,
                         M_QR_target_left, M_R_target_left,
                         qlen, dir);
      
    }
    
    H0 = _mm_subs_epi16(H3, R_query_left);
    H1 = _mm_subs_epi16(H0, R_query_left);
    H2 = _mm_subs_epi16(H1, R_query_left);
    H3 = _mm_subs_epi16(H2, R_query_left);

    dir += 4 * 4 * s->maxqlen;
    
    if (dir >= dirbuffer + dirbuffersize)
      dir -= dirbuffersize;
  }
}
