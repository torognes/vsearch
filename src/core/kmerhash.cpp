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

#include "core/kmerhash.hpp"
#include "utils/cityhash.hpp"  // hash_packed_kmer
#include "utils/kmer_hash_struct.hpp"
#include "utils/maps.hpp"
#include "utils/hash_table_size.hpp"  // table_size_half
#include "utils/view.hpp"  // View<char>
#include <cstddef>
#include <cstdint>
#include <vector>


using Hash = decltype(&hash_packed_kmer);
static constexpr Hash hash_function = hash_packed_kmer;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto reset_buckets(std::vector<struct kh_bucket_s> & hash) -> void {
    auto const current_size = hash.size();
    hash.clear();
    hash.resize(current_size);
  }

}  // end of anonymous namespace


namespace {
/* the empty-bucket sentinel, and the reason positions are stored 1-based:
   position 0 is a real position, so it cannot double as "no entry here". */
inline auto is_occupied(struct kh_bucket_s const & entry) noexcept -> bool
{
  return entry.pos != 0U;
}


inline auto kh_insert_kmer(struct kh_handle_s & kmer_hash,
                           int const k_offset,
                           unsigned int const kmer,
                           unsigned int const pos) -> void
{
  /* find free bucket in hash */
  auto bucket = hash_function(kmer, k_offset) & kmer_hash.hash_mask;
  while (is_occupied(kmer_hash.hash[bucket]))
    {
      bucket = (bucket + 1) & kmer_hash.hash_mask;
    }

  kmer_hash.hash[bucket].kmer = kmer;
  kmer_hash.hash[bucket].pos = pos;
}
}  // anonymous namespace


auto kh_insert_kmers(struct kh_handle_s & kmer_hash, int const k_offset, View<char> const seq) -> void
{
  int const kmers = static_cast<int>(1U << (2U * static_cast<unsigned int>(k_offset)));
  auto const kmer_mask = static_cast<unsigned int>(kmers - 1);

  reset_buckets(kmer_hash.hash);

  /* reallocate hash table if necessary */

  int64_t const needed = 2 * static_cast<int64_t>(seq.size());
  auto const wanted = static_cast<int64_t>(vsearch::table_size_half(seq.size()));
  if (kmer_hash.alloc < needed)
    {
      /* 'alloc' starts at a power of two and only ever grows, so the smallest
         power of two that fits is the same answer the doubling loop reached */
      kmer_hash.alloc = wanted;
      kmer_hash.hash.resize(static_cast<std::size_t>(kmer_hash.alloc));
    }

  kmer_hash.size = wanted;
  kmer_hash.hash_mask = static_cast<unsigned int>(kmer_hash.size - 1);

  kmer_hash.maxpos = static_cast<int>(seq.size());


  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;

  auto const * two_bit_map = chrmap_2bit();
  auto const * mask_ambig_map = chrmap_mask_ambig();
  int pos = 0;
  for (auto const nucleotide : seq)
    {
      bad <<= 2ULL;
      bad |= mask_ambig_map[static_cast<unsigned char>(nucleotide)];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= two_bit_map[static_cast<unsigned char>(nucleotide)];
      kmer &= kmer_mask;

      if (bad == 0U)
        {
          /* 1-based pos of start of kmer */
          kh_insert_kmer(kmer_hash, k_offset, kmer, static_cast<unsigned int>(pos - k_offset + 1 + 1));
        }
      ++pos;
    }
}


auto kh_find_diagonals(struct kh_handle_s const & kmer_hash,
                       int const k_offset,
                       View<char> const seq,
                       std::vector<int> & diags) -> void
{

  int const kmers = static_cast<int>(1U << (2U * static_cast<unsigned int>(k_offset)));
  auto const kmer_mask = static_cast<unsigned int>(kmers - 1);

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;

  /* the diagonal is a signed offset into diags, so the length is taken
     as an int once here rather than at each of the call sites */
  int const len = static_cast<int>(seq.size());

  auto const * two_bit_map = chrmap_2bit();
  auto const * mask_ambig_map = chrmap_mask_ambig();
  auto const * complement_map = chrmap_complement();
  auto seq_cursor = seq.crbegin();
  for (int pos = 0; pos < len; pos++)
    {
      char const nucleotide = *seq_cursor;
      ++seq_cursor;

      bad <<= 2ULL;
      bad |= mask_ambig_map[static_cast<unsigned char>(nucleotide)];
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= two_bit_map[complement_map[static_cast<unsigned char>(nucleotide)]];
      kmer &= kmer_mask;

      if (bad == 0U)
        {
          /* find matching buckets in hash */
          auto j = static_cast<unsigned int>(hash_function(kmer, k_offset) & kmer_hash.hash_mask);
          while (is_occupied(kmer_hash.hash[j]))
            {
              if (kmer_hash.hash[j].kmer == kmer)
                {
                  int const fpos = static_cast<int>(kmer_hash.hash[j].pos) - 1;
                  int const diag = len + fpos - (pos - k_offset + 1);
                  if (diag >= 0)
                    {
                      ++diags[static_cast<std::size_t>(diag)];
                    }
                }
              j = (j + 1) & kmer_hash.hash_mask;
            }
        }
    }
}

