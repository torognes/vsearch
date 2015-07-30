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

#include "align_simd.h"

#include <string.h>

#include "align_simd_dprofile.h"
#include "align_simd_helper.h"
#include "score_matrix.h"

#include "db.h"
#include "maps.h"
#include "util.h"

#define MAXSEQLENPRODUCT 25000000
//#define MAXSEQLENPRODUCT 160000

static void (*dprofile_fill16_func)(CELL *, CELL *, BYTE *);
unsigned int * chrmap;

struct s16info_s * search16_init(int score_match,
                                 int score_mismatch,
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
  ScoreMatrix::instance.init(score_match, score_mismatch, MATRIX_MODE_NUC);

  return search16_init_2(penalty_gap_open_query_left,
                       penalty_gap_open_target_left,
                       penalty_gap_open_query_interior,
                       penalty_gap_open_target_interior,
                       penalty_gap_open_query_right,
                       penalty_gap_open_target_right,
                       penalty_gap_extension_query_left,
                       penalty_gap_extension_target_left,
                       penalty_gap_extension_query_interior,
                       penalty_gap_extension_target_interior,
                       penalty_gap_extension_query_right,
                       penalty_gap_extension_target_right);
}

struct s16info_s * search16_init_2(CELL penalty_gap_open_query_left,
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
  struct s16info_s * s = (struct s16info_s *)xmalloc(sizeof(struct s16info_s));

  s->dprofile = (__m128i *)xmalloc(sizeof(CELL)*CDEPTH*CHANNELS*ScoreMatrix::instance.get_dimension());
  s->qlen = 0;
  s->qseq = 0;
  s->maxdlen = 0;
  s->dir = 0;
  s->diralloc = 0;
  s->hearray = 0;
  s->qtable = 0;
  s->cigar = 0;
  s->cigarend = 0;
  s->cigaralloc = 0;

  s->matrix = (__m128i *)ScoreMatrix::instance.score_matrix_16;

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

  if(ScoreMatrix::instance.is_nucleotide_mode())
    {
      dprofile_fill16_func = &dprofile_fill16;
      chrmap = chrmap_4bit;
    }
  else
    {
      dprofile_fill16_func = &dprofile_fill16_aa;
      chrmap = chrmap_aa_5bit;
    }

  return s;
}

/*
  The direction bits are set as follows:
  in DIR[0..1] if F>H initially (must go up) (4th pri)
  in DIR[2..3] if E>max(H,F) (must go left) (3rd pri)
  in DIR[4..5] if new F>H (must extend up) (2nd pri)
  in DIR[6..7] if new E>H (must extend left) (1st pri)
  no bits set: go diagonally
*/

#define ALIGNCORE(H, N, F, V, PATH, QR_q, R_q, QR_t, R_t, H_MIN, H_MAX) \
  H = _mm_adds_epi16(H, V);                                             \
  *(PATH+0) = _mm_movemask_epi8(_mm_cmpgt_epi16(F, H));                 \
  H = _mm_max_epi16(H, F);                                              \
  *(PATH+1) = _mm_movemask_epi8(_mm_cmpgt_epi16(E, H));                 \
  H = _mm_max_epi16(H, E);                                              \
  H_MIN = _mm_min_epi16(H_MIN, H);                                      \
  H_MAX = _mm_max_epi16(H_MAX, H);                                      \
  N = H;                                                                \
  HF = _mm_subs_epi16(H, QR_t);                                         \
  F = _mm_subs_epi16(F, R_t);                                           \
  *(PATH+2) = _mm_movemask_epi8(_mm_cmpgt_epi16(F, HF));                \
  F = _mm_max_epi16(F, HF);                                             \
  HE = _mm_subs_epi16(H, QR_q);                                         \
  E = _mm_subs_epi16(E, R_q);                                           \
  *(PATH+3) = _mm_movemask_epi8(_mm_cmpgt_epi16(E, HE));                \
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
                        __m128i f0,
                        __m128i f1,
                        __m128i f2,
                        __m128i f3,
                        __m128i * _h_min,
                        __m128i * _h_max,
                        __m128i Mm,
                        __m128i M_QR_t_left,
                        __m128i M_R_t_left,
                        __m128i M_QR_q_interior,
                        __m128i M_QR_q_right,
                        long ql,
                        unsigned short * dir)
{
  __m128i h4, h5, h6, h7, h8, E, HE, HF;
  __m128i * vp;
  __m128i h_min = _mm_setzero_si128();
  __m128i h_max = _mm_setzero_si128();
  long i;

  f0 = _mm_subs_epi16(f0, QR_t_0);
  f1 = _mm_subs_epi16(f1, QR_t_1);
  f2 = _mm_subs_epi16(f2, QR_t_2);
  f3 = _mm_subs_epi16(f3, QR_t_3);

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

      h4 = _mm_subs_epu16(h4, Mm);
      h4 = _mm_subs_epi16(h4, M_QR_t_left);

      E  = _mm_subs_epu16(E, Mm);
      E  = _mm_subs_epi16(E, M_QR_t_left);
      E  = _mm_subs_epi16(E, M_QR_q_interior);

      M_QR_t_left = _mm_adds_epi16(M_QR_t_left, M_R_t_left);

      ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0,
                QR_q_i, R_q_i, QR_t_0, R_t_0, h_min, h_max);
      ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4,
                QR_q_i, R_q_i, QR_t_1, R_t_1, h_min, h_max);
      ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8,
                QR_q_i, R_q_i, QR_t_2, R_t_2, h_min, h_max);
      ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
                QR_q_i, R_q_i, QR_t_3, R_t_3, h_min, h_max);

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

  E  = _mm_subs_epu16(E, Mm);
  E  = _mm_subs_epi16(E, M_QR_t_left);
  E  = _mm_subs_epi16(E, M_QR_q_right);

  ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0,
            QR_q_r, R_q_r, QR_t_0, R_t_0, h_min, h_max);
  ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4,
            QR_q_r, R_q_r, QR_t_1, R_t_1, h_min, h_max);
  ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8,
            QR_q_r, R_q_r, QR_t_2, R_t_2, h_min, h_max);
  ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
            QR_q_r, R_q_r, QR_t_3, R_t_3, h_min, h_max);

  hep[2*i+0] = h8;
  hep[2*i+1] = E;

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;

  *_h_min = h_min;
  *_h_max = h_max;
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
                       __m128i f0,
                       __m128i f1,
                       __m128i f2,
                       __m128i f3,
                       __m128i * _h_min,
                       __m128i * _h_max,
                       long ql,
                       unsigned short * dir)
{
  __m128i h4, h5, h6, h7, h8, E, HE, HF;
  __m128i * vp;
  __m128i h_min = _mm_setzero_si128();
  __m128i h_max = _mm_setzero_si128();
  long i;

  f0 = _mm_subs_epi16(f0, QR_t_0);
  f1 = _mm_subs_epi16(f1, QR_t_1);
  f2 = _mm_subs_epi16(f2, QR_t_2);
  f3 = _mm_subs_epi16(f3, QR_t_3);

  for(i=0; i < ql - 1; i++)
    {
      vp = qp[i+0];

      h4 = hep[2*i+0];

      E  = hep[2*i+1];

      ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0,
                QR_q_i, R_q_i, QR_t_0, R_t_0, h_min, h_max);
      ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4,
                QR_q_i, R_q_i, QR_t_1, R_t_1, h_min, h_max);
      ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8,
                QR_q_i, R_q_i, QR_t_2, R_t_2, h_min, h_max);
      ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
                QR_q_i, R_q_i, QR_t_3, R_t_3, h_min, h_max);

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

  ALIGNCORE(h0, h5, f0, vp[0], dir+16*i+ 0,
            QR_q_r, R_q_r, QR_t_0, R_t_0, h_min, h_max);
  ALIGNCORE(h1, h6, f1, vp[1], dir+16*i+ 4,
            QR_q_r, R_q_r, QR_t_1, R_t_1, h_min, h_max);
  ALIGNCORE(h2, h7, f2, vp[2], dir+16*i+ 8,
            QR_q_r, R_q_r, QR_t_2, R_t_2, h_min, h_max);
  ALIGNCORE(h3, h8, f3, vp[3], dir+16*i+12,
            QR_q_r, R_q_r, QR_t_3, R_t_3, h_min, h_max);

  hep[2*i+0] = h8;
  hep[2*i+1] = E;

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;

  *_h_min = h_min;
  *_h_max = h_max;
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

  if (s->hearray)
    free(s->hearray);
  s->hearray = (__m128i *) xmalloc(2 * s->qlen * sizeof(__m128i));
  memset(s->hearray, 0, 2 * s->qlen * sizeof(__m128i));

  if (s->qtable)
    free(s->qtable);
  s->qtable = (__m128i **) xmalloc(s->qlen * sizeof(__m128i*));

  for(int i = 0; i < qlen; i++)
    s->qtable[i] = s->dprofile + CDEPTH * chrmap[(int)(qseq[i])];
}

static void check_for_overflows(bool overflow[CHANNELS],
                                __m128i h_min, __m128i h_max,
                                short score_min, short score_max)
{
  for (int c = 0; c<CHANNELS; c++)
    {
      if (!overflow[c])
        {
          signed short h_min_array[CHANNELS];
          signed short h_max_array[CHANNELS];
          _mm_storeu_si128((__m128i *)h_min_array, h_min);
          _mm_storeu_si128((__m128i *)h_max_array, h_max);
          signed short h_min_c = h_min_array[c];
          signed short h_max_c = h_max_array[c];
          if ((h_min_c<=score_min)||(h_max_c>=score_max))
            overflow[c] = true;
        }
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
  unsigned long qlen = s->qlen;

  /* find longest target sequence and reallocate direction buffer */
  unsigned long maxdlen = 0;
  for(long i = 0; i < sequences; i++)
    {
      unsigned long dlen = db_getsequencelen(seqnos[i]);
      /* skip the very long sequences */
      if ((long)(s->qlen) * dlen <= MAXSEQLENPRODUCT)
        {
          if (dlen > maxdlen)
            maxdlen = dlen;
        }
    }
  maxdlen = CDEPTH * ((maxdlen + 3) / CDEPTH);
  s->maxdlen = maxdlen;
  unsigned long dirbuffersize = s->qlen * s->maxdlen * CDEPTH;

  if (dirbuffersize > s->diralloc)
    {
      s->diralloc = dirbuffersize;
      if (s->dir)
        free(s->dir);
      s->dir = (unsigned short*) xmalloc(dirbuffersize * sizeof(unsigned short));
    }

  unsigned short * dirbuffer = s->dir;

  if (s->qlen + s->maxdlen + 1 > s->cigaralloc)
    {
      s->cigaralloc = s->qlen + s->maxdlen + 1;
      if (s->cigar)
        free(s->cigar);
      s->cigar = (char *) xmalloc(s->cigaralloc);
    }

  __m128i T, M, T0;

  __m128i M_QR_target_left, M_R_target_left;
  __m128i M_QR_query_interior;
  __m128i M_QR_query_right;

  __m128i R_query_left;
  __m128i QR_query_interior, R_query_interior;
  __m128i QR_query_right, R_query_right;
  __m128i QR_target_left, R_target_left;
  __m128i QR_target_interior, R_target_interior;
  __m128i QR_target_right, R_target_right;
  __m128i QR_target[CDEPTH], R_target[CDEPTH];

  __m128i *hep, **qp;

  BYTE * d_begin[CHANNELS];
  BYTE * d_end[CHANNELS];
  unsigned long d_offset[CHANNELS];
  BYTE * d_address[CHANNELS];
  unsigned long d_length[CHANNELS];
  long seq_id[CHANNELS];
  bool overflow[CHANNELS];

  __m128i dseqalloc[CDEPTH];
  __m128i S[CDEPTH];

  BYTE * dseq = (BYTE*) & dseqalloc;
  BYTE zero = 0;

  unsigned long next_id = 0;
  unsigned long done = 0;

  T0 = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 0xffff);

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

  short score_min = SHRT_MIN + gap_penalty_max;
  short score_max = SHRT_MAX;

  for(int i=0; i<CDEPTH; i++)
    {
      S[i] = _mm_setzero_si128();
      dseqalloc[i] = _mm_setzero_si128();
    }

  __m128i H0 = _mm_setzero_si128();
  __m128i H1 = _mm_setzero_si128();
  __m128i H2 = _mm_setzero_si128();
  __m128i H3 = _mm_setzero_si128();

  __m128i F0 = _mm_setzero_si128();
  __m128i F1 = _mm_setzero_si128();
  __m128i F2 = _mm_setzero_si128();
  __m128i F3 = _mm_setzero_si128();

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
            dseq[CHANNELS*j+c] = chrmap[*(d_begin[c]++)];
          else
            dseq[CHANNELS*j+c] = 0;
        }
        if (d_begin[c] == d_end[c])
          easy = 0;
      }

      dprofile_fill16_func(dprofile, (CELL*) s->matrix, dseq);

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
                      (j >= ((d_length[c]+3) % CDEPTH)))
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

      __m128i h_min, h_max;

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

      check_for_overflows(overflow, h_min, h_max, score_min, score_max);
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
              dseq[CHANNELS*j+c] = chrmap[*(d_begin[c]++)];
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
            long z = (dbseqlen+3) % CDEPTH;
            long score = ((CELL*)S)[z*CHANNELS+c];

            if (overflow[c])
              {
                pscores[cand_id] = SHRT_MAX;
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
                pcigar[cand_id] = (char *) xmalloc(strlen(s->cigar)+1);
                strcpy(pcigar[cand_id], s->cigar);
              }

            done++;
          }

          /* get next sequence of reasonable length */

          long length = 0;

          while ((length == 0) && (next_id < sequences))
          {
            cand_id = next_id++;
            length = db_getsequencelen(seqnos[cand_id]);
            if ((length==0) || (s->qlen * length > MAXSEQLENPRODUCT))
              {
                pscores[cand_id] = SHRT_MAX;
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

              for(int j=0; j<CDEPTH; j++)
                {
                  if (d_begin[c] < d_end[c])
                    dseq[CHANNELS*j+c] = chrmap[*(d_begin[c]++)];
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

      /* make masked versions of QR and R for gaps in target */

      M_QR_target_left = _mm_and_si128(M, QR_target_left);
      M_R_target_left = _mm_and_si128(M, R_target_left);

      /* make masked versions of QR for gaps in query at target left end */

      M_QR_query_interior = _mm_and_si128(M, QR_query_interior);
      M_QR_query_right = _mm_and_si128(M, QR_query_right);

      dprofile_fill16_func(dprofile, (CELL*) s->matrix, dseq);

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
                      (j >= ((d_length[c]+3) % CDEPTH)))
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

      __m128i h_min, h_max;

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

      check_for_overflows(overflow, h_min, h_max, score_min, score_max);
    }

    H0 = _mm_subs_epi16(H3, R_query_left);
    H1 = _mm_subs_epi16(H0, R_query_left);
    H2 = _mm_subs_epi16(H1, R_query_left);
    H3 = _mm_subs_epi16(H2, R_query_left);

    F0 = _mm_subs_epi16(F3, R_query_left);
    F1 = _mm_subs_epi16(F0, R_query_left);
    F2 = _mm_subs_epi16(F1, R_query_left);
    F3 = _mm_subs_epi16(F2, R_query_left);

    dir += 4 * 4 * s->qlen;

    if (dir >= dirbuffer + dirbuffersize)
      dir -= dirbuffersize;
  }
}

