/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <algorithm>  // std::min
#include <cstdint> // uint64_t
#include <cstring>  // std::memset


/*
  Find the unique kmers or words in a given sequence.
  Unique is now defined as all different words occuring at least once.
  Earlier it was defined as those words occuring exactly once, but
  that caused a problem when searching for sequences with many repeats.
*/

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

  uint64_t bitmap_size;
  uint64_t * bitmap;
};

auto unique_init() -> struct uhandle_s *
{
  auto * uh = (struct uhandle_s *) xmalloc(sizeof(struct uhandle_s));

  uh->alloc = 2048;
  uh->size = 0;
  uh->hash_mask = uh->alloc - 1;
  uh->hash = (struct bucket_s *) xmalloc(sizeof(struct bucket_s) * uh->alloc);
  uh->list = (unsigned int *) xmalloc(sizeof(unsigned int) * uh->alloc);

  uh->bitmap_size = 0;
  uh->bitmap = nullptr;

  return uh;
}

auto unique_exit(struct uhandle_s * unique_handle) -> void
{
  if (unique_handle->bitmap)
    {
      xfree(unique_handle->bitmap);
    }
  if (unique_handle->hash)
    {
      xfree(unique_handle->hash);
    }
  if (unique_handle->list)
    {
      xfree(unique_handle->list);
    }
  xfree(unique_handle);
}

auto unique_compare(const void * a, const void * b) -> int
{
  auto * x = (unsigned int *) a;
  auto * y = (unsigned int *) b;

  if (x < y)
    {
      return -1;
    }
  else
    if (x > y)
      {
        return +1;
      }
    else
      {
        return 0;
      }
}


auto unique_count_bitmap(struct uhandle_s * unique_handle,
                         int k,
                         int seqlen,
                         char * seq,
                         unsigned int * listlen,
                         unsigned int * * list,
                         int seqmask) -> void
{
  /* if necessary, reallocate list of unique kmers */

  if (unique_handle->alloc < seqlen)
    {
      while (unique_handle->alloc < seqlen)
        {
          unique_handle->alloc *= 2;
        }
      unique_handle->list = (unsigned int *)
        xrealloc(unique_handle->list, sizeof(unsigned int) * unique_handle->alloc);
    }

  uint64_t const size = 1ULL << (k << 1ULL);

  /* reallocate bitmap arrays if necessary */

  if (unique_handle->bitmap_size < size)
    {
      unique_handle->bitmap = (uint64_t *) xrealloc(unique_handle->bitmap, size >> 3ULL);
      unique_handle->bitmap_size = size;
    }

  memset(unique_handle->bitmap, 0, size >> 3ULL);

  uint64_t bad = 0;
  uint64_t kmer = 0;
  uint64_t const mask = size - 1ULL;
  char * s = seq;
  char * e1 = s + k - 1;
  char * e2 = s + seqlen;
  e1 = std::min(e2, e1);

  unsigned int * maskmap = (seqmask != MASK_NONE) ?
    chrmap_mask_lower : chrmap_mask_ambig;

  while (s < e1)
    {
      bad <<= 2ULL;
      bad |= maskmap[(int) (*s)];

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[(int) (*s++)];
    }

  int unique = 0;

  while (s < e2)
    {
      bad <<= 2ULL;
      bad |= maskmap[(int) (*s)];
      bad &= mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[(int) (*s++)];
      kmer &= mask;

      if (! bad)
        {
          uint64_t const x = kmer >> 6ULL;
          uint64_t const y = 1ULL << (kmer & 63ULL);
          if (! (unique_handle->bitmap[x] & y))
            {
              /* not seen before */
              unique_handle->list[unique++] = kmer;
              unique_handle->bitmap[x] |= y;
            }
        }
    }

  *listlen = unique;
  *list = unique_handle->list;
}

auto unique_count_hash(struct uhandle_s * uh,
                       int k,
                       int seqlen,
                       char * seq,
                       unsigned int * listlen,
                       unsigned int * * list,
                       int seqmask) -> void
{
  /* if necessary, reallocate hash table and list of unique kmers */

  if (uh->alloc < 2 * seqlen)
    {
      while (uh->alloc < 2 * seqlen)
        {
          uh->alloc *= 2;
        }
      uh->hash = (struct bucket_s *)
        xrealloc(uh->hash, sizeof(struct bucket_s) * uh->alloc);
      uh->list = (unsigned int *)
        xrealloc(uh->list, sizeof(unsigned int) * uh->alloc);
    }

  /* hashtable variant */

  uh->size = 1;
  while (uh->size < 2*seqlen)
    {
      uh->size *= 2;
    }
  uh->hash_mask = uh->size - 1;

  memset(uh->hash, 0, sizeof(struct bucket_s) * uh->size);

  uint64_t bad = 0;
  uint64_t j = 0;
  unsigned int kmer = 0;
  unsigned int const mask = (1ULL << (2ULL * k)) - 1ULL;
  char * s = seq;
  char * e1 = s + k - 1;
  char * e2 = s + seqlen;
  e1 = std::min(e2, e1);

  unsigned int * maskmap = (seqmask != MASK_NONE) ?
    chrmap_mask_lower : chrmap_mask_ambig;

  while (s < e1)
    {
      bad <<= 2ULL;
      bad |= maskmap[(int) (*s)];

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[(int) (*s++)];
    }

  uint64_t unique = 0;

  while (s < e2)
    {
      bad <<= 2ULL;
      bad |= maskmap[(int) (*s)];
      bad &= mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[(int) (*s++)];
      kmer &= mask;

      if (! bad)
        {
          /* find free appropriate bucket in hash */
          j = HASH((char *) &kmer, (k + 3) / 4) & uh->hash_mask;
          while((uh->hash[j].count) && (uh->hash[j].kmer != kmer))
            {
              j = (j + 1) & uh->hash_mask;
            }

          if (! (uh->hash[j].count))
            {
              /* not seen before */
              uh->list[unique++] = kmer;
              uh->hash[j].kmer = kmer;
              uh->hash[j].count = 1;
            }
        }
    }

  *listlen = unique;
  *list = uh->list;
}

auto unique_count(struct uhandle_s * unique_handle,
                  int wordlength,
                  int seqlen,
                  char * seq,
                  unsigned int * listlen,
                  unsigned int * * list,
                  int seqmask) -> void
{
  if (wordlength < 10)
    {
      unique_count_bitmap(unique_handle, wordlength, seqlen, seq, listlen, list, seqmask);
    }
  else
    {
      unique_count_hash(unique_handle, wordlength, seqlen, seq, listlen, list, seqmask);
    }
}

auto unique_count_shared(struct uhandle_s * unique_handle,
                        int wordlength,
                        int listlen,
                        unsigned int * list) -> int
{
  /* counts how many of the kmers in list are present in the
     (already computed) hash or bitmap */

  int count = 0;
  if (wordlength < 10)
    {
      for (int i = 0; i < listlen; i++)
        {
          unsigned int const kmer = list[i];
          uint64_t const x = kmer >> 6ULL;
          uint64_t const y = 1ULL << (kmer & 63ULL);
          if (unique_handle->bitmap[x] & y)
            {
              ++count;
            }
        }
    }
  else
    {
      for (int i = 0; i < listlen; i++)
        {
          unsigned int kmer = list[i];
          uint64_t j = HASH((char *) &kmer, (wordlength + 3) / 4) & unique_handle->hash_mask;
          while((unique_handle->hash[j].count) && (unique_handle->hash[j].kmer != kmer))
            {
              j = (j + 1) & unique_handle->hash_mask;
            }
          if (unique_handle->hash[j].count)
            {
              ++count;
            }
        }
    }
  return count;
}
