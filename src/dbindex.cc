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

int k;

unsigned int * kmerhash;
unsigned int * kmerindex;
unsigned int kmerhashsize;
unsigned int kmerindexsize;

inline int dbindex_getkmermatchcount(int kmer)
{
  return kmerhash[kmer+1] - kmerhash[kmer];
}

inline int dbindex_getkmermatch(int kmer, int matchno)
{
  return kmerindex[kmerhash[kmer]+matchno];
}

void fprint_kmer(FILE * f, unsigned int k, unsigned long kmer)
{
  unsigned long x = kmer;
  char sym[] = "ACGT";
  for(unsigned int i=0; i<k; i++)
    fprintf(f, "%c", sym[(x >> (2*(k-i-1))) & 3]);
}

unsigned int extract_sequence_kmer(char * seq, unsigned int k, unsigned int pos)
{
  unsigned int kmer = 0;
  unsigned int mask = (1<<(2*k)) - 1;
  char * s = seq + pos;
  char * e = s + k;

  while (s < e)
    {
      kmer <<= 2;
      kmer |= *s++;
    }

  kmer &= mask;

  return kmer;
}

void dbindex_build()
{
  k = wordlength;

  show_rusage();

  fprintf(stderr, "Performing first database pass - counting unique %u-mers: ", k);

  unsigned int seqcount = db_getsequencecount();
  
  kmerhashsize = 1 << (2*k);

  unsigned int * kmercount = (unsigned int *) xmalloc(kmerhashsize * sizeof(unsigned int));
  memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  count_kmers_init();

  /* status of kmers for each sequence: 0 (zero), 1 (one) or 2 (two or more) */
  unsigned char * kmerstatus = (unsigned char *) xmalloc(kmerhashsize * sizeof(unsigned char));
  
  /* first scan, just count occurences */
  
  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      char * sequence;
      long seqlen;
      
      db_getsequenceandlength(seqno, & sequence, & seqlen);

      if (wordlength > 9)
        {
          count_kmers(k, sequence, seqlen);
          for(unsigned int i=0; i< count_kmers_gethashsize(); i++)
            if (kmercounthash[i].count == 1)
              kmercount[kmercounthash[i].kmer]++;
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
              kmer |= *s++;
            }
          
          while (s < e2)
            {
              kmer <<= 2;
              kmer |= *s++;
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
    }
  

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

  fprintf(stderr,"%u\n", kmerindexsize);

  kmerindex = (unsigned int *) xmalloc(kmerindexsize * sizeof(unsigned int));

  show_rusage();

  fprintf(stderr,"Performing second database pass - filling index\n");

  /* second scan, fill in actual index of the unique kmers */

  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      char * sequence;
      long seqlen;
      
      db_getsequenceandlength(seqno, & sequence, & seqlen);

      if (wordlength > 9)
        {
          count_kmers(k, sequence, seqlen);
          for(unsigned int i=0; i<count_kmers_gethashsize(); i++)
            if (kmercounthash[i].count == 1)
              kmerindex[kmerhash[kmercounthash[i].kmer]++] = seqno;
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
              kmer |= *s++;
            }
          
          while (s < e2)
            {
              kmer <<= 2;
              kmer |= *s++;
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
    }
  
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

#if 0
  
  /* print kmer hash and index */

  for(unsigned int i=0; i<kmerhashsize; i++)
    {
      unsigned int kmer = i;
      for(unsigned int j=0; j<k; j++)
        putchar(sym_nt[1+((kmer >> 2*(k-1-j)) & 3)]);
      fprintf(stderr,":");
      for(unsigned int j = kmerhash[i]; j < kmerhash[i+1]; j++)
        fprintf(stderr," %d", kmerindex[j]);
      fprintf(stderr,"\n");
    }

#endif

  show_rusage();
  count_kmers_exit();
}

void dbindex_free()
{
  free(kmerhash);
  free(kmerindex);
}
