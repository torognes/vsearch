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

#define HASH CityHash64

struct uhandle_s * unique_init()
{
  uhandle_s * uh = (struct uhandle_s *) xmalloc(sizeof(struct uhandle_s));
  uh->alloc = 512;
  uh->size = 0;
  uh->hash_mask = uh->alloc - 1;
  uh->hash = (struct bucket_s *) xmalloc(sizeof(struct bucket_s) * uh->alloc);
  uh->list = (unsigned int *) xmalloc(sizeof(unsigned int) * uh->alloc);
  return uh;
}

void unique_exit(struct uhandle_s * uh)
{
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


void unique_count(struct uhandle_s * uh, 
		  int k,
		  int seqlen,
		  char * seq,
		  unsigned int * listlen,
		  unsigned int * * list)
{
  if (uh->alloc < 2*seqlen)
    {
      while (uh->alloc < 2*seqlen)
	uh->alloc *= 2;
      uh->hash = (struct bucket_s *)
	xrealloc(uh->hash, sizeof(struct bucket_s) * uh->alloc);
      uh->list = (unsigned int *)
	xrealloc(uh->list, sizeof(unsigned int) * uh->alloc);
    }

  uh->size = 1;
  while (uh->size < 2*seqlen)
    uh->size *= 2;
  uh->hash_mask = uh->size - 1;

  memset(uh->hash, 0, sizeof(struct bucket_s) * uh->size);

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
      kmer <<= 2;
      kmer |= chrmap_2bit[(int)(*s++)];
    }

  while (s < e2)
    {
      kmer <<= 2;
      kmer |= chrmap_2bit[(int)(*s++)];
      kmer &= mask;

      /* find free appropriate bucket in hash */
      j = HASH((char*)&kmer, (k+3)/4) & uh->hash_mask;
      while((uh->hash[j].count) && (uh->hash[j].kmer != kmer))
        j = (j + 1) & uh->hash_mask;

      uh->hash[j].kmer = kmer;
      uh->hash[j].count++;
    }

  int unique = 0;
  for(int i=0; i<uh->size; i++)
    if (uh->hash[i].count == 1)
      uh->list[unique++] = uh->hash[i].kmer;

  //  qsort(uh->list, unique, sizeof(unsigned int), unique_compare);

  *listlen = unique;
  *list = uh->list;
}      

