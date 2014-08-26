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

static int k;

unsigned int * kmerhash;
unsigned int * kmerindex;
unsigned int kmerhashsize;
unsigned int kmerindexsize;

static char sym_nt[] = "ACGT";

inline int dbindex_getkmermatchcount(int kmer)
{
  return kmerhash[kmer+1] - kmerhash[kmer];
}

inline int dbindex_getkmermatch(int kmer, int matchno)
{
  return kmerindex[kmerhash[kmer]+matchno];
}

void fprint_kmer(FILE * f, unsigned int kk, unsigned long kmer)
{
  unsigned long x = kmer;
  for(unsigned int i=0; i<kk; i++)
    fprintf(f, "%c", sym_nt[(x >> (2*(kk-i-1))) & 3]);
}

unsigned int extract_sequence_kmer(char * seq, unsigned int kk, unsigned int pos)
{
  unsigned int kmer = 0;
  unsigned int mask = (1<<(2*kk)) - 1;
  char * s = seq + pos;
  char * e = s + kk;

  while (s < e)
    {
      kmer <<= 2;
      kmer |= chrmap_2bit[(int)(*s++)];
    }

  kmer &= mask;

  return kmer;
}

void dbindex_build()
{
  k = wordlength;

  show_rusage();

  unsigned int seqcount = db_getsequencecount();
  
  kmerhashsize = 1 << (2*k);

  unsigned int * kmercount = (unsigned int *) xmalloc(kmerhashsize * sizeof(unsigned int));
  memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  struct uhandle_s * uh = unique_init();

  /* status of kmers for each sequence: 0 (zero), 1 (one) or 2 (two or more) */
  unsigned char * kmerstatus = (unsigned char *) xmalloc(kmerhashsize * sizeof(unsigned char));
  
  /* first scan, just count occurences */
  
  progress_init("Counting unique k-mers", seqcount);

  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      char * sequence;
      long seqlen;
      
      db_getsequenceandlength(seqno, & sequence, & seqlen);

      if (wordlength > 9)
        {
	  unsigned int uniquecount;
	  unsigned int * uniquelist;
	  unique_count(uh, k, seqlen, sequence, & uniquecount, & uniquelist);
          for(unsigned int i=0; i<uniquecount; i++)
	    kmercount[uniquelist[i]]++;
        }
      else
        {
          memset(kmerstatus, 0, kmerhashsize * sizeof(unsigned char));

          unsigned int kmer = 0;
          unsigned int mask = kmerhashsize - 1;
          unsigned int i = 0;
          
          char * s = sequence + i;
          char * e1 = sequence + k - 1;
          char * e2 = sequence + seqlen;
          if (e2 < e1)
            e1 = e2;
          
          while (s < e1)
            {
              kmer <<= 2;
	      kmer |= chrmap_2bit[(int)(*s++)];
            }
          
          while (s < e2)
            {
              kmer <<= 2;
	      kmer |= chrmap_2bit[(int)(*s++)];
              kmer &= mask;
              
              if (kmerstatus[kmer] == 0)
                {
                  /* first occurence */
                  kmerstatus[kmer] = 1;
                  kmercount[kmer]++;
                }
              else if (kmerstatus[kmer] == 1)
                {
                  /* second occurence - not unique */
                  kmerstatus[kmer] = 2;
                  /* correct the count of unique kmers */
                  kmercount[kmer]--;
                }
            }
        }
      progress_update(seqno);
    }
 
  progress_done();

  /* hash setup */

  /* convert hash counts to position in index */

  kmerhash = (unsigned int *) xmalloc((kmerhashsize+1) * sizeof(unsigned int));

  unsigned int sum = 0;

  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      kmerhash[i] = sum;
      sum += kmercount[i];
    }

  kmerindexsize = sum;
  kmerhash[kmerhashsize] = sum;

  fprintf(stderr, "Unique %u-mers: %u\n", k, kmerindexsize);

  kmerindex = (unsigned int *) xmalloc(kmerindexsize * sizeof(unsigned int));

  show_rusage();

  progress_init("Creating index of unique k-mers", seqcount);

  /* second scan, fill in actual index of the unique kmers */

  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      char * sequence;
      long seqlen;
      
      db_getsequenceandlength(seqno, & sequence, & seqlen);

      if (wordlength > 9)
        {
	  unsigned int uniquecount;
	  unsigned int * uniquelist;
	  unique_count(uh, k, seqlen, sequence, & uniquecount, & uniquelist);
          for(unsigned int i=0; i<uniquecount; i++)
	    kmerindex[kmerhash[uniquelist[i]]++] = seqno;
        }
      else
        {

          /* count kmers */
          
          memset(kmerstatus, 0, kmerhashsize * sizeof(unsigned char));
          
          unsigned int kmer = 0;
          unsigned int mask = kmerhashsize - 1;
          unsigned int i = 0;
          
          char * s = sequence + i;
          char * e1 = sequence + k - 1;
          char * e2 = sequence + seqlen;
          if (e2 < e1)
            e1 = e2;
          
          while (s < e1)
            {
              kmer <<= 2;
	      kmer |= chrmap_2bit[(int)(*s++)];
            }
          
          while (s < e2)
            {
              kmer <<= 2;
	      kmer |= chrmap_2bit[(int)(*s++)];
              kmer &= mask;
              
              if (kmerstatus[kmer] == 0)
                {
                  /* first occurence */
                  if (kmercount[kmer] > 0)
                    {
                      /* list not already full */
                      kmerstatus[kmer] = 1;
                      kmerindex[kmerhash[kmer]++] = seqno;
                      kmercount[kmer]--;
                    }
                  else
                    {
                      /* already full - cannot be unique */
                      kmerstatus[kmer] = 2;
                    }
                }
              else if (kmerstatus[kmer] == 1)
                {
                  /* second occurence */
                  /* not unique after all - adjust counts */
                  kmerstatus[kmer] = 2;
                  kmerhash[kmer]--;
                  kmercount[kmer]++;
                }
            }
        }
      progress_update(seqno);
    }
  
  progress_done();

  free(kmercount);
  free(kmerstatus);

  /* reset kmerhash pointers (move up) */

  unsigned int temp = 0;
  unsigned int next;

  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      next = kmerhash[i];
      kmerhash[i] = temp;
      temp = next;
    }

  show_rusage();
  unique_exit(uh);
}

void dbindex_free()
{
  free(kmerhash);
  free(kmerindex);
}
