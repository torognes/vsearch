/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "city.h"
#include "maps.h"
#include <algorithm>  // std::max
#include <cstring>  // std::memset
#include <vector>


using Hash = decltype(&CityHash64);
static Hash hash_function = CityHash64;

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


auto kh_init() -> struct kh_handle_s *
{
  auto * kmer_hash =
    (struct kh_handle_s *) xmalloc(sizeof(struct kh_handle_s));

  kmer_hash->maxpos = 0;
  kmer_hash->alloc = 256;
  kmer_hash->size = 0;
  kmer_hash->hash_mask = kmer_hash->alloc - 1;
  kmer_hash->hash =
    (struct kh_bucket_s *) xmalloc(kmer_hash->alloc * sizeof(struct kh_bucket_s));

  return kmer_hash;
}


auto kh_exit(struct kh_handle_s * kmer_hash) -> void
{
  if (kmer_hash->hash != nullptr)
    {
      xfree(kmer_hash->hash);
    }
  xfree(kmer_hash);
}


inline auto kh_insert_kmer(struct kh_handle_s * kmer_hash,
                           int k_offset,
                           unsigned int kmer,
                           unsigned int pos) -> void
{
  /* find free bucket in hash */
  auto bucket = hash_function((char *) &kmer, (k_offset + 3) / 4) & kmer_hash->hash_mask;
  while (kmer_hash->hash[bucket].pos != 0U)
    {
      bucket = (bucket + 1) & kmer_hash->hash_mask;
    }

  kmer_hash->hash[bucket].kmer = kmer;
  kmer_hash->hash[bucket].pos = pos;
}


auto kh_insert_kmers(struct kh_handle_s * kmer_hash, int k_offset, char * seq, int len) -> void
{
  int const kmers = 1U << (2U * k_offset);
  unsigned int const kmer_mask = kmers - 1;

  /* reallocate hash table if necessary */

  if (kmer_hash->alloc < 2 * len)
    {
      while (kmer_hash->alloc < 2 * len)
        {
          kmer_hash->alloc *= 2;
        }
      kmer_hash->hash = (struct kh_bucket_s *)
        xrealloc(kmer_hash->hash, kmer_hash->alloc * sizeof(struct kh_bucket_s));
    }

  kmer_hash->size = 1;
  while (kmer_hash->size < 2 * len)
    {
      kmer_hash->size *= 2;
    }
  kmer_hash->hash_mask = kmer_hash->size - 1;

  kmer_hash->maxpos = len;

  std::memset(kmer_hash->hash, 0, kmer_hash->size * sizeof(struct kh_bucket_s));

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;
  char * s = seq;

  unsigned int * maskmap = chrmap_mask_ambig;

  for (int pos = 0; pos < len; pos++)
    {
      int const c = *s++;

      bad <<= 2ULL;
      bad |= maskmap[c];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[c];
      kmer &= kmer_mask;

      if (bad == 0U)
        {
          /* 1-based pos of start of kmer */
          kh_insert_kmer(kmer_hash, k_offset, kmer, pos - k_offset + 1 + 1);
        }
    }
}


auto kh_find_best_diagonal(struct kh_handle_s * kmer_hash, int k_offset, char * seq, int len) -> int
{
  std::vector<int> diag_counts(kmer_hash->maxpos, 0);

  int const kmers = 1U << (2U * k_offset);
  unsigned int const kmer_mask = kmers - 1;

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;
  char * s = seq + len - 1;

  unsigned int * maskmap = chrmap_mask_ambig;

  for (int pos = 0; pos < len; pos++)
    {
      int const c = *s--;

      bad <<= 2ULL;
      bad |= maskmap[c];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[chrmap_complement[c]];
      kmer &= kmer_mask;

      if (bad == 0U)
        {
          /* find matching buckets in hash */
          unsigned int j = hash_function((char *) &kmer, (k_offset + 3) / 4) & kmer_hash->hash_mask;
          while (kmer_hash->hash[j].pos != 0U)
            {
              if (kmer_hash->hash[j].kmer == kmer)
                {
                  int const fpos = kmer_hash->hash[j].pos - 1;
                  int const diag = fpos - (pos - k_offset + 1);
                  if (diag >= 0)
                    {
                      ++diag_counts[diag];
                    }
                }
              j = (j + 1) & kmer_hash->hash_mask;
            }
        }
    }

  int best_diag_count = -1;
  int best_diag = -1;
  int good_diags = 0;

  for (int d = 0; d < kmer_hash->maxpos - k_offset + 1; d++)
    {
      int const diag_len = kmer_hash->maxpos - d;
      int const minmatch = std::max(1, diag_len - k_offset + 1 - (k_offset * std::max(diag_len / 20, 0)));
      int const c = diag_counts[d];

      if (c >= minmatch)
        {
          ++good_diags;
        }

      if (c > best_diag_count)
        {
          best_diag_count = c;
          best_diag = d;
        }
    }

  if (good_diags == 1)
    {
      return best_diag;
    }
  return -1;
}


auto kh_find_diagonals(struct kh_handle_s * kmer_hash,
                       int k_offset,
                       char * seq,
                       int len,
                       int * diags) -> void
{
  std::memset(diags, 0, (kmer_hash->maxpos+len) * sizeof(int));

  int const kmers = 1U << (2U * k_offset);
  unsigned int const kmer_mask = kmers - 1;

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;
  char * s = seq + len - 1;

  for (int pos = 0; pos < len; pos++)
    {
      int const c = *s--;

      bad <<= 2ULL;
      bad |= chrmap_mask_ambig[c];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= chrmap_2bit[chrmap_complement[c]];
      kmer &= kmer_mask;

      if (bad == 0U)
        {
          /* find matching buckets in hash */
          unsigned int j = hash_function((char *) &kmer, (k_offset + 3) / 4) & kmer_hash->hash_mask;
          while (kmer_hash->hash[j].pos != 0U)
            {
              if (kmer_hash->hash[j].kmer == kmer)
                {
                  int const fpos = kmer_hash->hash[j].pos - 1;
                  int const diag = len + fpos - (pos - k_offset + 1);
                  if (diag >= 0)
                    {
                      ++diags[diag];
                    }
                }
              j = (j + 1) & kmer_hash->hash_mask;
            }
        }
    }
}
