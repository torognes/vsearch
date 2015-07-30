/*
    Copyright (C) 2012-2015 Torbjorn Rognes & Frederic Mahe

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

#include "align_simd.h"

#include <string.h>

#include "align_simd_helper.h"

/*
  Using 16-bit signed values, from -32768 to +32767.
  match: positive
  mismatch: negative
  gap penalties: positive (open, extend, query/target, left/interior/right)
  optimal global alignment (NW)
  maximize score
*/

//#define DEBUG

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

/*
  The direction bits are set as follows:
  in DIR[0..1] if F>H initially (must go up) (4th pri)
  in DIR[2..3] if E>max(H,F) (must go left) (3rd pri)
  in DIR[4..5] if new F>H (must extend up) (2nd pri)
  in DIR[6..7] if new E>H (must extend left) (1st pri)
  no bits set: go diagonally
*/

static inline void pushop(s16info_s * s, char newop)
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

static inline void finishop(s16info_s * s)
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
  unsigned long dirbuffersize = s->qlen * s->maxdlen * 4;
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

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long d = *((unsigned long *) (dirbuffer + 
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

  long i = qlen - 1;
  long j = dlen - 1;

  s->cigarend = s->cigar + s->qlen + s->maxdlen + 1;
  s->op = 0;
  s->opcount = 1;

  while ((i>=0) && (j>=0))
  {
    aligned++;

    unsigned long d = *((unsigned long *) (dirbuffer +  // TODO check this part regarding amino acid sequences
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
      if (chrmap[(int)(qseq[i])] == chrmap[(int)(dseq[j])])
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
  int cigarlen = s->cigar + s->qlen + s->maxdlen - s->cigarend;
  memmove(s->cigar, s->cigarend, cigarlen + 1);
  
  * paligned = aligned;
  * pmatches = matches;
  * pmismatches = mismatches;
  * pgaps = gaps;
}
