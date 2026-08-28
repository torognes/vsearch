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
#include "utils/cityhash.hpp"  // hash_packed_kmer
#include "core/mask.hpp"  // Masking
#include "utils/maps.hpp"  // chrmap_2bit, chrmap_mask_lower, chrmap_mask_ambig
#include "utils/hash_table_size.hpp"  // table_size_half
#include <algorithm>  // std::min, std::fill, std::fill_n
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // int64_t, uint64_t


/*
  Find the unique kmers or words in a given sequence.
  Unique is now defined as all different words occurring at least once.
  Earlier it was defined as those words occurring exactly once, but
  that caused a problem when searching for sequences with many repeats.
*/


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // refactoring: deliberately not std::unordered_map. These tables are
  // per-sequence scratch, sized up front so the load factor never exceeds 0.5
  // and no rehash is ever needed, and count() picks between a bitmap, an
  // epoch-stamp table and this one by wordlength -- none of which a node-based
  // container expresses. Reasoning: DONE_20260825_flatmap_helpers.md.
  using Hash = decltype(&hash_packed_kmer);
  constexpr Hash hash_function = hash_packed_kmer;

}  // end of anonymous namespace


auto Uniquer::count_bitmap(int const wordlength,
                           View<char> const seq,
                           Masking const seqmask) -> View<unsigned int>
{
  /* if necessary, grow the list of unique kmers (at most seq.size() entries) */

  if (list_.size() < seq.size())
    {
      list_.resize(seq.size());
    }

  uint64_t const size = 1ULL << (wordlength << 1ULL);

  /* (re)allocate and clear the bitmap: one bit per possible kmer, packed in
     64-bit words (size >> 6 of them). assign() reuses the buffer when the word
     count is unchanged (the common case: it depends only on wordlength). */

  bitmap_.assign(size >> 6ULL, 0);

  uint64_t bad = 0;
  uint64_t kmer = 0;
  uint64_t const mask = size - 1ULL;
  /* the wordlength - 1 leading bases only prime the rolling kmer; a sequence
     shorter than that primes it as far as it goes and yields no kmer */
  auto const primer_length = std::min(static_cast<std::size_t>(wordlength - 1), seq.size());

  auto const * mask_map = (seqmask != Masking::none) ?
    chrmap_mask_lower() : chrmap_mask_ambig();
  auto const * two_bit_map = chrmap_2bit();

  auto * const bitmap = bitmap_.data();
  auto * const list_data = list_.data();

  for (auto const nucleotide : seq.first(primer_length))
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(nucleotide)];

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(nucleotide)];
    }

  auto unique = 0;

  for (auto const nucleotide : seq.drop(primer_length))
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(nucleotide)];
      bad &= mask;

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(nucleotide)];
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

  return View<unsigned int>{list_data, static_cast<std::size_t>(unique)};
}


auto Uniquer::count_stamps(int const wordlength,
                           View<char> const seq,
                           Masking const seqmask) -> View<unsigned int>
{
  /* if necessary, grow the list of unique kmers (at most seq.size() entries) */

  if (list_.size() < seq.size())
    {
      list_.resize(seq.size());
    }

  uint64_t const size = 1ULL << (static_cast<uint64_t>(wordlength) << 1ULL);

  /* One stamp per possible kmer; a kmer was seen in THIS call iff its stamp
     equals the current epoch, so bumping the epoch invalidates the whole
     table at once and nothing is cleared between calls. The table is zeroed
     when (re)grown and when the 16-bit epoch wraps around, i.e. every 65535
     calls -- amortized to a fraction of one entry per kmer scanned. */
  if (stamp_of_kmer_.size() < size)
    {
      stamp_of_kmer_.assign(size, 0);
      epoch_ = 0;
    }
  ++epoch_;
  if (epoch_ == 0)
    {
      std::fill(stamp_of_kmer_.begin(), stamp_of_kmer_.end(), static_cast<uint16_t>(0));
      epoch_ = 1;
    }

  uint64_t bad = 0;
  uint64_t kmer = 0;
  uint64_t const mask = size - 1ULL;
  /* see count_bitmap: the leading wordlength - 1 bases only prime the kmer */
  auto const primer_length = std::min(static_cast<std::size_t>(wordlength - 1), seq.size());

  auto const * mask_map = (seqmask != Masking::none) ?
    chrmap_mask_lower() : chrmap_mask_ambig();
  auto const * two_bit_map = chrmap_2bit();

  auto * const stamps = stamp_of_kmer_.data();
  auto * const list_data = list_.data();
  auto const epoch = epoch_;

  for (auto const nucleotide : seq.first(primer_length))
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(nucleotide)];

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(nucleotide)];
    }

  auto unique = 0;

  for (auto const nucleotide : seq.drop(primer_length))
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(nucleotide)];
      bad &= mask;

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(nucleotide)];
      kmer &= mask;

      if ((bad == 0U) and (stamps[kmer] != epoch))
        {
          /* not seen before */
          stamps[kmer] = epoch;
          list_data[unique] = static_cast<unsigned int>(kmer);
          ++unique;
        }
    }

  return View<unsigned int>{list_data, static_cast<std::size_t>(unique)};
}


auto Uniquer::count_hash(int const wordlength,
                         View<char> const seq,
                         Masking const seqmask) -> View<unsigned int>
{
  /* size the hash table and the list of unique kmers to the sequence. The
     64-bit width is load-bearing: the hash grows to 2 * sequence length, which
     exceeds INT_MAX for sequences above ~1.07 Gnt (the doubling would otherwise
     overflow int before reaching the target). See utils/hash_table_size.hpp. */

  auto const size = vsearch::table_size_half(seq.size());
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
  auto const mask = static_cast<unsigned int>((1ULL << (2ULL * static_cast<unsigned int>(wordlength))) - 1ULL);
  /* see count_bitmap: the leading wordlength - 1 bases only prime the kmer */
  auto const primer_length = std::min(static_cast<std::size_t>(wordlength - 1), seq.size());

  auto const * mask_map = (seqmask != Masking::none) ?
    chrmap_mask_lower() : chrmap_mask_ambig();
  auto const * two_bit_map = chrmap_2bit();

  auto * const hash = hash_.data();
  auto * const list_data = list_.data();

  for (auto const nucleotide : seq.first(primer_length))
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(nucleotide)];

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(nucleotide)];
    }

  uint64_t unique = 0;

  for (auto const nucleotide : seq.drop(primer_length))
    {
      bad <<= 2ULL;
      bad |= mask_map[static_cast<unsigned char>(nucleotide)];
      bad &= mask;

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(nucleotide)];
      kmer &= mask;

      if (bad == 0U)
        {
          /* find free appropriate bucket in hash */
          uint64_t j = hash_function(kmer, wordlength) & hash_mask_;
          while (is_occupied(hash[j]) && (hash[j].kmer != kmer))
            {
              j = (j + 1) & hash_mask_;
            }

          if (not is_occupied(hash[j]))
            {
              /* not seen before */
              list_data[unique] = kmer;
              ++unique;
              hash[j].kmer = kmer;
              hash[j].count = 1;
            }
        }
    }

  return View<unsigned int>{list_data, static_cast<std::size_t>(unique)};
}


/* Dispatch by table size. Up to wordlength 9 the one-bit-per-kmer bitmap is
   at most 32 KiB and is recleared on every call. Up to wordlength 12 the
   stamp table is at most 32 MiB: too large to clear per call, but the epoch
   trick never clears it, and a direct lookup beats hashing every kmer
   occurrence (the CityHash call and probe were a fifth of --orient's
   runtime, whose default word length is 12). Beyond that the table would
   grow four-fold per unit of word length, so the hash takes over. */
constexpr auto max_bitmap_wordlength = 9;
constexpr auto max_stamps_wordlength = 12;

auto Uniquer::count(int const wordlength,
                    View<char> const seq,
                    Masking const seqmask) -> View<unsigned int>
{
  if (wordlength <= max_bitmap_wordlength)
    {
      return count_bitmap(wordlength, seq, seqmask);
    }
  if (wordlength <= max_stamps_wordlength)
    {
      return count_stamps(wordlength, seq, seqmask);
    }
  return count_hash(wordlength, seq, seqmask);
}


auto Uniquer::count_shared(int const wordlength,
                           View<unsigned int> const list) const noexcept -> unsigned int
{
  /* counts how many of the kmers in list are present in the (already
     computed) bitmap, stamp table or hash -- picked by the same word-length
     thresholds as count() above */

  auto count = 0U;
  if (wordlength <= max_bitmap_wordlength)
    {
      auto const * const bitmap = bitmap_.data();
      for (auto const kmer : list)
        {
          uint64_t const x = kmer >> 6ULL;
          uint64_t const y = 1ULL << (kmer & 63ULL);
          if ((bitmap[x] & y) != 0U)
            {
              ++count;
            }
        }
    }
  else if (wordlength <= max_stamps_wordlength)
    {
      auto const * const stamps = stamp_of_kmer_.data();
      for (auto const kmer : list)
        {
          if (stamps[kmer] == epoch_)
            {
              ++count;
            }
        }
    }
  else
    {
      auto const * const hash = hash_.data();
      for (auto kmer : list)
        {
          uint64_t j = hash_function(kmer, wordlength) & hash_mask_;
          while (is_occupied(hash[j]) && (hash[j].kmer != kmer))
            {
              j = (j + 1) & hash_mask_;
            }
          if (is_occupied(hash[j]))
            {
              ++count;
            }
        }
    }
  return count;
}
