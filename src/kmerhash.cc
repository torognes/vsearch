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

#define HASH CityHash64

struct kh_bucket_s
{
  unsigned int kmer;
  unsigned int pos; /* 1-based position, 0 = empty */
};

struct kh_handle_s
{
  struct kh_bucket_s * hash;
  unsigned int hash_mask;
  int size;
  int alloc;
  int maxpos;
};

struct kh_handle_s * kh_init()
{
  kh_handle_s * kh =
    (struct kh_handle_s *) xmalloc(sizeof(struct kh_handle_s));

  kh->maxpos = 0;
  kh->alloc = 256;
  kh->size = 0;
  kh->hash_mask = kh->alloc - 1;
  kh->hash =
    (struct kh_bucket_s *) xmalloc(kh->alloc * sizeof(struct kh_bucket_s));

  return kh;
}

void kh_exit(struct kh_handle_s * kh)
{
  if (kh->hash)
    xfree(kh->hash);
  xfree(kh);
}

inline void kh_insert_kmer(struct kh_handle_s * kh,
                           int k,
                           unsigned int kmer,
                           unsigned int pos)
{
  /* find free bucket in hash */
  unsigned int j = HASH((char*)&kmer, (k+3)/4) & kh->hash_mask;
  while(kh->hash[j].pos)
    j = (j + 1) & kh->hash_mask;

  kh->hash[j].kmer = kmer;
  kh->hash[j].pos = pos;
}

void kh_insert_kmers(struct kh_handle_s * kh, int k, char * seq, int len)
{
  int kmers = 1 << (2 * k);
  unsigned int kmer_mask = kmers - 1;

  /* reallocate hash table if necessary */

  if (kh->alloc < 2 * len)
    {
      while (kh->alloc < 2 * len)
        kh->alloc *= 2;
      kh->hash = (struct kh_bucket_s *)
        xrealloc(kh->hash, kh->alloc * sizeof(struct kh_bucket_s));
    }

  kh->size = 1;
  while(kh->size < 2 * len)
    kh->size *= 2;
  kh->hash_mask = kh->size - 1;

  kh->maxpos = len;

  memset(kh->hash, 0, kh->size * sizeof(struct kh_bucket_s));

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;
  char * s = seq;

  unsigned int * maskmap = chrmap_mask_ambig;

  for (int pos = 0; pos < len; pos++)
    {
      int c = *s++;

      bad <<= 2ULL;
      bad |= maskmap[c];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[c];
      kmer &= kmer_mask;

      if (!bad)
        {
          /* 1-based pos of start of kmer */
          kh_insert_kmer(kh, k, kmer, pos - k + 1 + 1);
        }
    }
}

int kh_find_best_diagonal(struct kh_handle_s * kh, int k, char * seq, int len)
{
  int diag_counts[kh->maxpos];

  memset(diag_counts, 0, kh->maxpos * sizeof(int));

  int kmers = 1 << (2 * k);
  unsigned int kmer_mask = kmers - 1;

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;
  char * s = seq + len - 1;

  unsigned int * maskmap = chrmap_mask_ambig;

  for (int pos = 0; pos < len; pos++)
    {
      int c = *s--;

      bad <<= 2ULL;
      bad |= maskmap[c];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[chrmap_complement[c]];
      kmer &= kmer_mask;

      if (!bad)
        {
          /* find matching buckets in hash */
          unsigned int j = HASH((char*)&kmer, (k+3)/4) & kh->hash_mask;
          while(kh->hash[j].pos)
            {
              if (kh->hash[j].kmer == kmer)
                {
                  int fpos = kh->hash[j].pos - 1;
                  int diag = fpos - (pos - k + 1);
                  if (diag >= 0)
                    diag_counts[diag]++;
                }
              j = (j + 1) & kh->hash_mask;
            }
        }
    }

  int best_diag_count = -1;
  int best_diag = -1;
  int good_diags = 0;

  for(int d = 0; d < kh->maxpos - k + 1; d++)
    {
      int diag_len = kh->maxpos - d;
      int minmatch = MAX(1, diag_len - k + 1 - k * MAX(diag_len / 20, 0));
      int c = diag_counts[d];

      if (c >= minmatch)
        good_diags++;

      if (c > best_diag_count)
        {
          best_diag_count = c;
          best_diag = d;
        }
    }

  if (good_diags == 1)
    return best_diag;
  else
    return -1;
}

void kh_find_diagonals(struct kh_handle_s * kh,
                       int k,
                       char * seq,
                       int len,
                       int * diags)
{
  memset(diags, 0, (kh->maxpos+len) * sizeof(int));

  int kmers = 1 << (2 * k);
  unsigned int kmer_mask = kmers - 1;

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;
  char * s = seq + len - 1;

  for (int pos = 0; pos < len; pos++)
    {
      int c = *s--;

      bad <<= 2ULL;
      bad |= chrmap_mask_ambig[c];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[chrmap_complement[c]];
      kmer &= kmer_mask;

      if (!bad)
        {
          /* find matching buckets in hash */
          unsigned int j = HASH((char*)&kmer, (k+3)/4) & kh->hash_mask;
          while(kh->hash[j].pos)
            {
              if (kh->hash[j].kmer == kmer)
                {
                  int fpos = kh->hash[j].pos - 1;
                  int diag = len + fpos - (pos - k + 1);
                  if (diag >= 0)
                    diags[diag]++;
                }
              j = (j + 1) & kh->hash_mask;
            }
        }
    }
}
