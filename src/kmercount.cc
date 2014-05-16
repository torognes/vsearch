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

struct kmercountelem * kmercounthash;

unsigned int kmercounthashsize;

size_t kmercounthash_alloc;

unsigned int count_kmers_gethashsize()
{
  return kmercounthashsize;
}

void count_kmers_init()
{
  kmercounthash_alloc = 0;
  kmercounthashsize = 0;
  kmercounthash = 0;
}

void count_kmers_exit()
{
  free(kmercounthash);
  kmercounthash = 0;
}

void count_kmers(unsigned int k, char * seq, unsigned int seqlen)
{
  kmercounthashsize = 2 * seqlen;

  if (kmercounthashsize > kmercounthash_alloc)
    {
      /* resize memory fro hash if necessary */
      kmercounthash_alloc = kmercounthashsize;
      kmercounthash = (struct kmercountelem *) 
	xrealloc(kmercounthash, sizeof(struct kmercountelem) * kmercounthashsize);
    }

  memset(kmercounthash, 0, sizeof(struct kmercountelem) * kmercounthashsize);

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
      kmer |= *s++;
    }

  while (s < e2)
    {
      kmer <<= 2;
      kmer |= *s++;
      kmer &= mask;

      unsigned long j = hash_fnv_1a_64((unsigned char*)&kmer,
				       (k+3)/4) % kmercounthashsize;
      
      while((kmercounthash[j].count) && (kmercounthash[j].kmer != kmer))
	j = (j + 1) % kmercounthashsize;

      kmercounthash[j].kmer = kmer;
      kmercounthash[j].count++;
    }
}      

