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

/* Find the unique kmers or words in a given sequence. Currently the
   definintion of "unique kmers" are the kmers that appear exactly once;
   kmers that appear twice or more often are ignored.
   Maybe the definition should be all kmers that appear in the input
   sequence, but that each of those kmers shoudl appear only
   once in the resulting list.
   Experiments indicate that the current (first) definition gives
   both the highest sensitivity and precision. */

#define HASH CityHash64

struct bucket_s
{
  unsigned int kmer;
  unsigned int count;
};

struct uhandle_s
{
  struct bucket_s * hash;
  unsigned int * list;
  unsigned int hash_mask;
  int size;
  int alloc;

  unsigned long bitmap_size;
  unsigned long * bitmap;
};

struct uhandle_s * unique_init()
{
  uhandle_s * uh = (struct uhandle_s *) xmalloc(sizeof(struct uhandle_s));

  uh->alloc = 2048;
  uh->size = 0;
  uh->hash_mask = uh->alloc - 1;
  uh->hash = (struct bucket_s *) xmalloc(sizeof(struct bucket_s) * uh->alloc);
  uh->list = (unsigned int *) xmalloc(sizeof(unsigned int) * uh->alloc);

  uh->bitmap_size = 0;
  uh->bitmap = 0;

  return uh;
}

void unique_exit(struct uhandle_s * uh)
{
  if (uh->bitmap)
    free(uh->bitmap);
  if (uh->hash)
    free(uh->hash);
  if (uh->list)
    free(uh->list);
  free(uh);
}

int unique_compare(const void * a, const void * b)
{
  unsigned int * x = (unsigned int*) a;
  unsigned int * y = (unsigned int*) b;

  if (x<y)
    return -1;
  else
    if (x>y)
      return +1;
    else
      return 0;
}


void unique_count_bitmap(struct uhandle_s * uh, 
                         int k,
                         int seqlen,
                         char * seq,
                         unsigned int * listlen,
                         unsigned int * * list)
{
  /* if necessary, reallocate list of unique kmers */

  if (uh->alloc < seqlen)
    {
      while (uh->alloc < seqlen)
        uh->alloc *= 2;
      uh->list = (unsigned int *)
        xrealloc(uh->list, sizeof(unsigned int) * uh->alloc);
    }
  
  unsigned long size = 1UL << (k << 1UL);
  
  /* reallocate bitmap arrays if necessary */
  
  if (uh->bitmap_size < size)
    {
      uh->bitmap = (unsigned long *) xrealloc(uh->bitmap, size >> 2UL);
      uh->bitmap_size = size;
    }
  
  memset(uh->bitmap, 0, size >> 2UL);
  
  unsigned long bad = 0;
  unsigned long kmer = 0;
  unsigned long mask = size - 1UL;
  char * s = seq;
  char * e1 = s + k-1;
  char * e2 = s + seqlen;
  if (e2 < e1)
    e1 = e2;
  
  while (s < e1)
    {
      bad <<= 2UL;
      bad |= (chrmap_4bit[(int)(*s)] > 4);

      kmer <<= 2UL;
      kmer |= chrmap_2bit[(int)(*s++)];
    }
      
  int candidates = 0;

  while (s < e2)
    {
      bad <<= 2UL;
      bad |= (chrmap_4bit[(int)(*s)] > 4);
      bad &= mask;

      kmer <<= 2UL;
      kmer |= chrmap_2bit[(int)(*s++)];
      kmer &= mask;
      
      if (!bad)
        {
          unsigned long x = kmer >> 6UL;
          unsigned long y = 1UL << (kmer & 63UL);
              
          unsigned long uniq = uh->bitmap[2*x];
          unsigned long seen = uh->bitmap[2*x+1];
          
          uh->bitmap[2*x] = (uniq & ~y) | (y & ~seen);
          uh->bitmap[2*x+1] = seen | y;
          
          uh->list[candidates++] = kmer;
        }
    }
      
  int unique = 0;

  for (int i = 0; i < candidates; i++)
    {
      unsigned long kmer = uh->list[i];
      unsigned long x = kmer >> 6UL;
      unsigned long y = 1UL << (kmer & 63UL);
      
      if (uh->bitmap[2*x] & y)
        uh->list[unique++] = uh->list[i];
    }
  
  *listlen = unique;
  *list = uh->list;
}

void unique_count_hash(struct uhandle_s * uh, 
                       int k,
                       int seqlen,
                       char * seq,
                       unsigned int * listlen,
                       unsigned int * * list)
{
  unsigned long unique;

  /* if necessary, reallocate hash table and list of unique kmers */

  if (uh->alloc < 2*seqlen)
    {
      while (uh->alloc < 2*seqlen)
        uh->alloc *= 2;
      uh->hash = (struct bucket_s *)
        xrealloc(uh->hash, sizeof(struct bucket_s) * uh->alloc);
      uh->list = (unsigned int *)
        xrealloc(uh->list, sizeof(unsigned int) * uh->alloc);
    }
  
  /* hashtable variant */

  uh->size = 1;
  while (uh->size < 2*seqlen)
    uh->size *= 2;
  uh->hash_mask = uh->size - 1;
      
  memset(uh->hash, 0, sizeof(struct bucket_s) * uh->size);
      
  unsigned long bad = 0;
  unsigned long j;
  unsigned int kmer = 0;
  unsigned int mask = (1<<(2*k)) - 1;
  char * s = seq;
  char * e1 = s + k-1;
  char * e2 = s + seqlen;
  if (e2 < e1)
    e1 = e2;
      
  while (s < e1)
    {
      bad <<= 2UL;
      bad |= (chrmap_4bit[(int)(*s)] > 4);

      kmer <<= 2;
      kmer |= chrmap_2bit[(int)(*s++)];
    }
      
  while (s < e2)
    {
      bad <<= 2UL;
      bad |= (chrmap_4bit[(int)(*s)] > 4);
      bad &= mask;

      kmer <<= 2;
      kmer |= chrmap_2bit[(int)(*s++)];
      kmer &= mask;
          
      if (!bad)
        {
          /* find free appropriate bucket in hash */
          j = HASH((char*)&kmer, (k+3)/4) & uh->hash_mask;
          while((uh->hash[j].count) && (uh->hash[j].kmer != kmer))
            j = (j + 1) & uh->hash_mask;
              
          uh->hash[j].kmer = kmer;
          uh->hash[j].count++;
        }
    }
      
  unique = 0;
  for(int i=0; i<uh->size; i++)
    if (uh->hash[i].count == 1)
      uh->list[unique++] = uh->hash[i].kmer;
      
  *listlen = unique;
  *list = uh->list;
}

void unique_count(struct uhandle_s * uh, 
                  int k,
                  int seqlen,
                  char * seq,
                  unsigned int * listlen,
                  unsigned int * * list)
{
  if (k<10)
    unique_count_bitmap(uh, k, seqlen, seq, listlen, list);
  else
    unique_count_hash(uh, k, seqlen, seq, listlen, list);
}

int unique_count_shared(struct uhandle_s * uh,
                        int k,
                        int listlen,
                        unsigned int * list)
{
  /* counts how many of the kmers in list are present in the
     (already computed) hash or bitmap */
  
  int count = 0;
  if (k<10)
    {
      for(int i = 0; i<listlen; i++)
        {
          unsigned int kmer = list[i];
          unsigned long x = kmer >> 6UL;
          unsigned long y = 1UL << (kmer & 63UL);
          if (uh->bitmap[2*x] & y)
            count++;
        }
    }
  else
    {
      for(int i = 0; i<listlen; i++)
        {
          unsigned int kmer = list[i];
          unsigned long j = HASH((char*)&kmer, (k+3)/4) & uh->hash_mask;
          while((uh->hash[j].count) && (uh->hash[j].kmer != kmer))
            j = (j + 1) & uh->hash_mask;
          if (uh->hash[j].count == 1)
            count++;
        }
    }
  return count;
}
