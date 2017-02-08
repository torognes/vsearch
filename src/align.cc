/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

struct nwinfo_s
{
  int64_t dir_alloc;
  int64_t hearray_alloc;
  char * dir;
  int64_t * hearray;
};

static const char maskup      =  1;
static const char maskleft    =  2;
static const char maskextup   =  4;
static const char maskextleft =  8;

inline void pushop(char newop, char ** cigarendp, char * op, int * count)
{
  if (newop == *op)
    (*count)++;
  else
  {
    *--*cigarendp = *op;
    if (*count > 1)
    {
      char buf[25];
      int len = sprintf(buf, "%d", *count);
      *cigarendp -= len;
      strncpy(*cigarendp, buf, (size_t)len);
    }
    *op = newop;
    *count = 1;
  }
}

inline void finishop(char ** cigarendp, char * op, int * count)
{
  if ((op) && (count))
  {
    *--*cigarendp = *op;
    if (*count > 1)
    {
      char buf[25];
      int len = sprintf(buf, "%d", *count);
      *cigarendp -= len;
      strncpy(*cigarendp, buf, (size_t)len);
    }
    *op = 0;
    *count = 0;
  }
}

/*

  Needleman-Wunsch aligner

  finds a global alignment with maximum score
  positive score for matches; negative score mismatches
  gap penalties are positive, but counts negatively

  alignment priority when backtracking (from lower right corner):
  1. left/insert/e (gap in query sequence (qseq))
  2. up/delete/f (gap in database sequence (dseq))
  3. align/diag/h (match/mismatch)
  
  qseq: the reference/query/upper/vertical/from sequence
  dseq: the sample/database/lower/horisontal/to sequence

  default (interior) scores:
  match: +2
  mismatch: -4
  gap open: 20
  gap extend: 2

  input

  dseq: pointer to start of database sequence
  dend: pointer after database sequence
  qseq: pointer to start of query sequence
  qend: pointer after database sequence
  score_matrix: 16x16 matrix of longs with scores for aligning two symbols
  gapopen: positive number indicating penalty for opening a gap of length zero
  gapextend: positive number indicating penalty for extending a gap

  output

  nwscore: the global alignment score
  nwdiff: number of non-identical nucleotides in one optimal global alignment
  nwalignmentlength: the length of one optimal alignment
  nwalignment: cigar string with one optimal alignment

*/

struct nwinfo_s * nw_init()
{
  struct nwinfo_s * nw = (struct nwinfo_s *) xmalloc(sizeof(struct nwinfo_s));
  nw->dir = 0;
  nw->dir_alloc = 0;
  nw->hearray = 0;
  nw->hearray_alloc = 0;
  return nw;
}

void nw_exit(struct nwinfo_s * nw)
{
  if (nw->dir)
    xfree(nw->dir);
  if (nw->hearray)
    xfree(nw->hearray);
  xfree(nw);
}

inline int nt_identical(char a, char b)
{
  if (chrmap_4bit[(int)a] == chrmap_4bit[(int)b])
    return 1;
  else
    return 0;
}

inline int64_t getscore(int64_t * score_matrix, char a, char b)
{
  return score_matrix[(chrmap_4bit[(int)a]<<4) + chrmap_4bit[(int)b]];
}

void nw_align(char * dseq,
              char * dend,
              char * qseq,
              char * qend,
              int64_t * score_matrix,
              int64_t gapopen_q_left,
              int64_t gapopen_q_interior,
              int64_t gapopen_q_right,
              int64_t gapopen_t_left,
              int64_t gapopen_t_interior,
              int64_t gapopen_t_right,
              int64_t gapextend_q_left,
              int64_t gapextend_q_interior,
              int64_t gapextend_q_right,
              int64_t gapextend_t_left,
              int64_t gapextend_t_interior,
              int64_t gapextend_t_right,
              int64_t * nwscore,
              int64_t * nwdiff,
              int64_t * nwgaps,
              int64_t * nwindels,
              int64_t * nwalignmentlength,
              char ** nwalignment,
              int64_t queryno,
              int64_t dbseqno,
              struct nwinfo_s * nw)
{

  int64_t h, n, e, f, h_e, h_f;
  int64_t *hep;
  
  int64_t qlen = qend - qseq;
  int64_t dlen = dend - dseq;

  if (qlen * dlen > nw->dir_alloc)
    {
      nw->dir_alloc = qlen * dlen;
      nw->dir = (char *) xrealloc(nw->dir, (size_t)nw->dir_alloc);
    }

  int64_t need = 2 * qlen * (int64_t) sizeof(int64_t);
  if (need > nw->hearray_alloc)
    {
      nw->hearray_alloc = need;
      nw->hearray = (int64_t *) xrealloc(nw->hearray, (size_t)nw->hearray_alloc);
    }

  memset(nw->dir, 0, (size_t)(qlen*dlen));

  int64_t i, j;

  for(i=0; i<qlen; i++)
  {
    nw->hearray[2*i] = -gapopen_t_left - (i+1) * gapextend_t_left;
    if (i < qlen-1)
      nw->hearray[2*i+1] =
        - gapopen_t_left - (i+1) * gapextend_t_left
        - gapopen_q_interior - gapextend_q_interior;
    else
      nw->hearray[2*i+1] =
        - gapopen_t_left - (i+1) * gapextend_t_left
        - gapopen_q_right - gapextend_q_right;
  }

  for(j=0; j<dlen; j++)
  {
    hep = nw->hearray;

    if (j == 0)
      h = 0;
    else
      h = - gapopen_q_left - j * gapextend_q_left;

    if (j < dlen-1)
      f = - gapopen_q_left - (j+1) * gapextend_q_left
        - gapopen_t_interior - gapextend_t_interior;
    else
      f = - gapopen_q_left - (j+1) * gapextend_q_left
        - gapopen_t_right - gapextend_t_right;

    for(i=0; i<qlen; i++)
    {
      char * d = nw->dir + qlen*j+i;
      
      n = *hep;
      e = *(hep+1);
      h += getscore(score_matrix, dseq[j], qseq[i]);
      
      if (f > h)
        {  
          h = f;
          *d |= maskup;
        }
      
      if (e > h)
        {
          h = e;
          *d |= maskleft;
        }

      *hep = h;

      if (i < qlen-1)
        {
          h_e = h - gapopen_q_interior - gapextend_q_interior;
          e -= gapextend_q_interior;
        }
      else
        {
          h_e = h - gapopen_q_right - gapextend_q_right;
          e -= gapextend_q_right;
        }
      
      if (j < dlen-1)
        {
          h_f = h - gapopen_t_interior - gapextend_t_interior;
          f -= gapextend_t_interior;
        }
      else
        {
          h_f = h - gapopen_t_right - gapextend_t_right;
          f -= gapextend_t_right;
        }

      if (f > h_f)
        *d |= maskextup;
      else
        f = h_f;

      if (e > h_e)
        *d |= maskextleft;
      else
        e = h_e;

      *(hep+1) = e;
      h = n;
      hep += 2;
    }
  }
  
  int64_t dist = nw->hearray[2*qlen-2];
  
  /* backtrack: count differences and save alignment in cigar string */

  int64_t score = 0;
  int64_t alength = 0;
  int64_t matches = 0;
  int64_t gaps = 0;
  int64_t indels = 0;

  char * cigar = (char *) xmalloc((size_t)(qlen + dlen + 1));
  char * cigarend = cigar+qlen+dlen+1;

  char op = 0;
  int count = 0;
  *(--cigarend) = 0;

  i = qlen;
  j = dlen;

  while ((i>0) && (j>0))
  {
    int64_t gapopen_q   = (i < qlen) ? gapopen_q_interior   : gapopen_q_right;
    int64_t gapextend_q = (i < qlen) ? gapextend_q_interior : gapextend_q_right;
    int64_t gapopen_t   = (j < dlen) ? gapopen_t_interior   : gapopen_t_right;
    int64_t gapextend_t = (j < dlen) ? gapextend_t_interior : gapextend_t_right;
  
    int d = nw->dir[qlen*(j-1)+(i-1)];

    alength++;

    if ((op == 'I') && (d & maskextleft))
    {
      score -= gapextend_q;
      indels++;
      j--;
      pushop('I', &cigarend, &op, &count);
    }
    else if ((op == 'D') && (d & maskextup))
    {
      score -= gapextend_t;
      indels++;
      i--;
      pushop('D', &cigarend, &op, &count);
    }
    else if (d & maskleft)
    {
      score -= gapextend_q;
      indels++;
      if (op != 'I')
        {
          score -= gapopen_q;
          gaps++;
        }
      j--;
      pushop('I', &cigarend, &op, &count);
    }
    else if (d & maskup)
    {
      score -= gapextend_t;
      indels++;
      if (op != 'D')
        {
          score -= gapopen_t;
          gaps++;
        }
      i--;
      pushop('D', &cigarend, &op, &count);
    }
    else
    {
      score += getscore(score_matrix, dseq[j-1], qseq[i-1]);
      if (nt_identical(dseq[j-1], qseq[i-1]))
        matches++;
      i--;
      j--;
      pushop('M', &cigarend, &op, &count);
    }
  }
  
  while(i>0)
  {
    alength++;
    score -= gapextend_t_left;
    indels++;
    if (op != 'D')
      {
        score -= gapopen_t_left;
        gaps++;
      }
    i--;
    pushop('D', &cigarend, &op, &count);
  }

  while(j>0)
  {
    alength++;
    score -= gapextend_q_left;
    indels++;
    if (op != 'I')
      {
        score -= gapopen_q_left;
        gaps++;
      }
    j--;
    pushop('I', &cigarend, &op, &count);
  }

  finishop(&cigarend, &op, &count);

  /* move and reallocate cigar */

  int64_t cigarlength = cigar+qlen+dlen-cigarend;
  memmove(cigar, cigarend, (size_t)(cigarlength+1));
  cigar = (char*) xrealloc(cigar, (size_t)(cigarlength+1));

  * nwscore = dist;
  * nwdiff = alength - matches;
  * nwalignmentlength = alength;
  * nwalignment = cigar;
  * nwgaps = gaps;
  * nwindels = indels;

#if 1
  if (score != dist)
  {
    fprintf(stderr, "WARNING: Error with query no %" PRId64 " and db sequence no %" PRId64 ":\n", queryno, dbseqno);
    fprintf(stderr, "Initial and recomputed alignment score disagreement: %" PRId64 " %" PRId64 "\n", dist, score);
    fprintf(stderr, "Alignment: %s\n", cigar);

    if (opt_log)
      {
        fprintf(fp_log, "WARNING: Error with query no %" PRId64 " and db sequence no %" PRId64 ":\n", queryno, dbseqno);
        fprintf(fp_log, "Initial and recomputed alignment score disagreement: %" PRId64 " %" PRId64 "\n", dist, score);
        fprintf(fp_log, "Alignment: %s\n", cigar);
        fprintf(fp_log, "\n");
      }
  }
#endif
}
