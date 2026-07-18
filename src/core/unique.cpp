/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "core/unique.hpp"
#include "vendored/city.h"
#include "core/mask.hpp"  // Masking
#include "utils/maps.hpp"  // chrmap_2bit, chrmap_mask_lower, chrmap_mask_ambig
#include <algorithm>  // std::min, std::fill_n
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t


/*
  Find the unique kmers or words in a given sequence.
  Unique is now defined as all different words occuring at least once.
  Earlier it was defined as those words occuring exactly once, but
  that caused a problem when searching for sequences with many repeats.
*/


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // refactoring:
  // replace with std::unordered_map (default hashing)
  // if performance are bad, see Victor_Ciura's Cpp Talk "So You Think You Can Hash"
  // then make a CityHash hasher object and use it with std::unordered_map
  using Hash = decltype(&CityHash64);
  constexpr Hash hash_function = CityHash64;

}  // end of anonymous namespace


auto Uniquer::count_bitmap(int const wordlength,
                           int const seqlen,
                           char const * seq,
                           unsigned int * listlen,
                           unsigned int const * * list,
                           Masking const seqmask) -> void
{
  /* if necessary, grow the list of unique kmers (at most seqlen entries) */

  if (list_.size() < static_cast<std::size_t>(seqlen))
    {
      list_.resize(static_cast<std::size_t>(seqlen));
    }

  uint64_t const size = 1ULL << (wordlength << 1ULL);

  /* (re)allocate and clear the bitmap: one bit per possible kmer, packed in
     64-bit words (size >> 6 of them). assign() reuses the buffer when the word
     count is unchanged (the common case: it depends only on wordlength). */

  bitmap_.assign(size >> 6ULL, 0);

  uint64_t bad = 0;
  uint64_t kmer = 0;
  uint64_t const mask = size - 1ULL;
  auto const * s = seq;
  auto const * e1 = s + wordlength - 1;
  auto const * e2 = s + seqlen;
  e1 = std::min(e2, e1);

  auto const * mask_map = (seqmask != Masking::none) ?
    chrmap_mask_lower() : chrmap_mask_ambig();
  auto const * two_bit_map = chrmap_2bit();

  auto * const bitmap = bitmap_.data();
  auto * const list_data = list_.data();

  while (s < e1)
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(*s)];

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(*s)];
      ++s;
    }

  auto unique = 0;

  while (s < e2)
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(*s)];
      bad &= mask;

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(*s)];
      ++s;
      kmer &= mask;

      if (bad == 0U)
        {
          uint64_t const x = kmer >> 6ULL;
          uint64_t const y = 1ULL << (kmer & 63ULL);
          if ((bitmap[x] & y) == 0U)
            {
              /* not seen before */
              list_data[unique] = static_cast<unsigned int>(kmer);
              ++unique;
              bitmap[x] |= y;
            }
        }
    }

  *listlen = static_cast<unsigned int>(unique);
  *list = list_data;
}


auto Uniquer::count_hash(int const wordlength,
                         int const seqlen,
                         char const * seq,
                         unsigned int * listlen,
                         unsigned int const * * list,
                         Masking const seqmask) -> void
{
  /* size the hash table and the list of unique kmers to the sequence. needed
     and size are 64-bit: the hash grows to 2 * sequence length, which exceeds
     INT_MAX for sequences above ~1.07 Gnt (the doubling would otherwise overflow
     int before reaching the target). */

  int64_t const needed = 2 * static_cast<int64_t>(seqlen);
  int64_t size = 1;
  while (size < needed)
    {
      size *= 2;
    }
  hash_mask_ = static_cast<unsigned int>(size - 1);

  /* the buffers only grow; the first 'table_size' buckets are cleared each call
     (the probe never looks past hash_mask_, so any leftover tail is untouched) */

  auto const table_size = static_cast<std::size_t>(size);
  if (hash_.size() < table_size)
    {
      hash_.resize(table_size);
    }
  if (list_.size() < table_size)
    {
      list_.resize(table_size);
    }
  std::fill_n(hash_.begin(), table_size, bucket{});

  uint64_t bad = 0;
  auto kmer = 0U;
  unsigned int const mask = static_cast<unsigned int>((1ULL << (2ULL * static_cast<unsigned int>(wordlength))) - 1ULL);
  auto const * s = seq;
  auto const * e1 = s + wordlength - 1;
  auto const * e2 = s + seqlen;
  e1 = std::min(e2, e1);

  auto const * mask_map = (seqmask != Masking::none) ?
    chrmap_mask_lower() : chrmap_mask_ambig();
  auto const * two_bit_map = chrmap_2bit();

  auto * const hash = hash_.data();
  auto * const list_data = list_.data();

  while (s < e1)
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(*s)];

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(*s)];
      ++s;
    }

  uint64_t unique = 0;

  while (s < e2)
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(*s)];
      bad &= mask;

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(*s)];
      ++s;
      kmer &= mask;

      if (bad == 0U)
        {
          /* find free appropriate bucket in hash */
          uint64_t j = hash_function(reinterpret_cast<char const *>(&kmer), static_cast<size_t>((wordlength + 3) / 4)) & hash_mask_;
          while ((hash[j].count != 0U) && (hash[j].kmer != kmer))
            {
              j = (j + 1) & hash_mask_;
            }

          if (hash[j].count == 0U)
            {
              /* not seen before */
              list_data[unique] = kmer;
              ++unique;
              hash[j].kmer = kmer;
              hash[j].count = 1;
            }
        }
    }

  *listlen = static_cast<unsigned int>(unique);
  *list = list_data;
}


auto Uniquer::count(int const wordlength,
                    int const seqlen,
                    char const * seq,
                    unsigned int * listlen,
                    unsigned int const * * list,
                    Masking const seqmask) -> void
{
  if (wordlength < 10)
    {
      count_bitmap(wordlength, seqlen, seq, listlen, list, seqmask);
    }
  else
    {
      count_hash(wordlength, seqlen, seq, listlen, list, seqmask);
    }
}


auto Uniquer::count_shared(int const wordlength,
                           int const listlen,
                           unsigned int const * list) const noexcept -> unsigned int
{
  /* counts how many of the kmers in list are present in the
     (already computed) hash or bitmap */

  auto count = 0U;
  if (wordlength < 10)
    {
      auto const * const bitmap = bitmap_.data();
      for (auto i = 0; i < listlen; i++)
        {
          auto const kmer = list[i];
          uint64_t const x = kmer >> 6ULL;
          uint64_t const y = 1ULL << (kmer & 63ULL);
          if ((bitmap[x] & y) != 0U)
            {
              ++count;
            }
        }
    }
  else
    {
      auto const * const hash = hash_.data();
      for (auto i = 0; i < listlen; i++)
        {
          auto kmer = list[i];
          uint64_t j = hash_function(reinterpret_cast<char const *>(&kmer), static_cast<size_t>((wordlength + 3) / 4)) & hash_mask_;
          while ((hash[j].count != 0U) && (hash[j].kmer != kmer))
            {
              j = (j + 1) & hash_mask_;
            }
          if (hash[j].count != 0U)
            {
              ++count;
            }
        }
    }
  return count;
}
