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
      strncpy(*cigarendp, buf, len);
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
      strncpy(*cigarendp, buf, len);
    }
    *op = 0;
    *count = 0;
  }
}

const unsigned char maskup      = 1;
const unsigned char maskleft    = 2;
const unsigned char maskextup   = 4;
const unsigned char maskextleft = 8;

/*

  Needleman/Wunsch/Sellers aligner

  finds a global alignment with minimum cost
  there should be positive costs/penalties for gaps and for mismatches
  matches should have zero cost (0)

  alignment priority when backtracking (from lower right corner):
  1. left/insert/e (gap in query sequence (qseq))
  2. align/diag/h (match/mismatch)
  3. up/delete/f (gap in database sequence (dseq))
  
  qseq: the reference/query/upper/vertical/from sequence
  dseq: the sample/database/lower/horisontal/to sequence

  typical costs:
  match: 0
  mismatch: 3
  gapopen: 4
  gapextend: 3

  input

  dseq: pointer to start of database sequence
  dend: pointer after database sequence
  qseq: pointer to start of query sequence
  qend: pointer after database sequence
  score_matrix: 32x32 matrix of longs with scores for aligning two symbols
  gapopen: positive number indicating penalty for opening a gap of length zero
  gapextend: positive number indicating penalty for extending a gap

  output

  nwscore: the global alignment score
  nwdiff: number of non-identical nucleotides in one optimal global alignment
  nwalignmentlength: the length of one optimal alignment
  nwalignment: cigar string with one optimal alignment

*/

long dir_alloc = 0;
unsigned long hearray_alloc = 0;

unsigned char * dir;
unsigned long * hearray;

void nw_init()
{
  dir = 0;
  dir_alloc = 0;
  hearray = 0;
  hearray_alloc = 0;
}

void nw_exit()
{
  if (dir)
    free(dir);
  if (hearray)
    free(hearray);
}


void nw_align(char * dseq,
              char * dend,
              char * qseq,
              char * qend,
              long * score_matrix,
              unsigned long gapopen,
              unsigned long gapextend,
              unsigned long * nwscore,
              unsigned long * nwdiff,
              unsigned long * nwgaps,
              unsigned long * nwindels,
              unsigned long * nwalignmentlength,
              char ** nwalignment,
              unsigned long queryno,
              unsigned long dbseqno)
{

  long h, n, e, f;
  long unsigned *hep;

  long qlen = qend - qseq;
  long dlen = dend - dseq;

  if (qlen * dlen > dir_alloc)
    {
      dir_alloc = qlen * dlen;
      dir = (unsigned char *) xrealloc(dir, dir_alloc);
    }

  if (2 * qlen * sizeof(long) > hearray_alloc)
    {
      hearray_alloc = 2 * qlen * sizeof(long);
      hearray = (unsigned long *) xrealloc(hearray, hearray_alloc);
    }

  memset(dir, 0, qlen*dlen);

  long i, j;

  for(i=0; i<qlen; i++)
  {
    hearray[2*i]   = 1 * gapopen + (i+1) * gapextend; // H (N)
    hearray[2*i+1] = 2 * gapopen + (i+2) * gapextend; // E
  }

  for(j=0; j<dlen; j++)
  {
    hep = hearray;
    f = 2 * gapopen + (j+2) * gapextend;
    h = (j == 0) ? 0 : (gapopen + j * gapextend);
    
    for(i=0; i<qlen; i++)
    {
      unsigned char * d = dir + qlen*j+i;
      
      n = *hep;
      e = *(hep+1);
      h += score_matrix[(dseq[j]<<5) + qseq[i]];
      
      *d |= (f < h ? maskup : 0);
      h = MIN(h, f);
      h = MIN(h, e);
      *d |= (e == h ? maskleft : 0);

      *hep = h;
      
      h += gapopen + gapextend;
      e += gapextend;
      f += gapextend;
      
      *d |= (f < h ? maskextup : 0);
      *d |= (e < h ? maskextleft : 0);
      f = MIN(h,f);
      e = MIN(h,e);
      
      *(hep+1) = e;
      h = n;
      hep += 2;
    }
  }
  
  long dist = hearray[2*qlen-2];
  
  /* backtrack: count differences and save alignment in cigar string */

  long score = 0;
  long alength = 0;
  long matches = 0;
  long gaps = 0;
  long indels = 0;

  char * cigar = (char *) xmalloc(qlen + dlen + 1);
  char * cigarend = cigar+qlen+dlen+1;

  char op = 0;
  int count = 0;
  *(--cigarend) = 0;

  i = qlen;
  j = dlen;

  while ((i>0) && (j>0))
  {
    int d = dir[qlen*(j-1)+(i-1)];

    alength++;

    if ((op == 'I') && (d & maskextleft))
    {
      score += gapextend;
      indels++;
      j--;
      pushop('I', &cigarend, &op, &count);
    }
    else if ((op == 'D') && (d & maskextup))
    {
      score += gapextend;
      indels++;
      i--;
      pushop('D', &cigarend, &op, &count);
    }
    else if (d & maskleft)
    {
      score += gapextend;
      indels++;
      if (op != 'I')
        {
          score += gapopen;
          gaps++;
        }
      j--;
      pushop('I', &cigarend, &op, &count);
    }
    else if (d & maskup)
    {
      score += gapextend;
      indels++;
      if (op != 'D')
        {
          score += gapopen;
          gaps++;
        }
      i--;
      pushop('D', &cigarend, &op, &count);
    }
    else
    {
      score += score_matrix[(dseq[j-1] << 5) + qseq[i-1]];
      if (qseq[i-1] == dseq[j-1])
        matches++;
      i--;
      j--;
      pushop('M', &cigarend, &op, &count);
    }
  }
  
  while(i>0)
  {
    alength++;
    score += gapextend;
    indels++;
    if (op != 'D')
      {
        score += gapopen;
        gaps++;
      }
    i--;
    pushop('D', &cigarend, &op, &count);
  }
  
  while(j>0)
  {
    alength++;
    score += gapextend;
    indels++;
    if (op != 'I')
      {
        score += gapopen;
        gaps++;
      }
    j--;
    pushop('I', &cigarend, &op, &count);
  }

  finishop(&cigarend, &op, &count);

  /* move and reallocate cigar */

  long cigarlength = cigar+qlen+dlen-cigarend;
  memmove(cigar, cigarend, cigarlength+1);
  cigar = (char*) xrealloc(cigar, cigarlength+1);

  * nwscore = dist;
  * nwdiff = alength - matches;
  * nwalignmentlength = alength;
  * nwalignment = cigar;
  * nwgaps = gaps;
  * nwindels = indels;

  if (score != dist)
  {
    fprintf(stderr, "Error with query no %lu and db sequence no %lu:\n", queryno, dbseqno);
    fprintf(stderr, "Initial and recomputed alignment score disagreement: %lu %lu\n", dist, score);
    fprintf(stderr, "Alignment: %s\n", cigar);
  }
}
