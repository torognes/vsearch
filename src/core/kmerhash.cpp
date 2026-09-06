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
#include "utils/kmer_hash_struct.hpp"
#include "utils/maps.hpp"
#include "utils/maps/mask_ambig.hpp"
#include "utils/maps/two_bit.hpp"
#include "utils/view.hpp"  // View<char>
#include <cassert>
#include <cstddef>
#include <limits>  // std::numeric_limits
#include <vector>

namespace maps = vsearch::maps;
namespace mask_ambig = vsearch::maps::mask_ambig;
namespace two_bit = vsearch::maps::two_bit;


/* A packed k-mer of length k occupies 2 * k bits, so the whole k-mer space is
   4^k values: at the k the merge core uses (5) that is 1024 of them, few
   enough that the k-mer can address an array directly. 'chain_head[kmer]'
   then names a forward position carrying that k-mer, and the other positions
   carrying it are threaded through 'chain_next', which is indexed by
   position. Answering "which forward positions carry this k-mer" costs one
   load and a chain walk in which every step is a match -- where an
   open-addressing table also had to hash the k-mer, probe linearly, and
   compare the k-mer of every other entry that had landed in the same run.

   Positions are stored biased by one, because position 0 is a real position
   and so cannot double as the "end of chain" sentinel. */


auto kh_insert_kmers(struct kh_handle_s & kmer_hash, int const k_offset, View<char> const seq) -> void
{
  assert(k_offset > 0);
  assert(k_offset <= kmer_hash_max_k);
  /* a position is stored in an unsigned int, biased by one */
  assert(seq.size() < std::numeric_limits<unsigned int>::max());

  int const kmers = static_cast<int>(1U << (2U * static_cast<unsigned int>(k_offset)));
  auto const kmer_mask = static_cast<unsigned int>(kmers - 1);

  /* the heads have to be cleared: a stale one would name a position of the
     previous read. The links do not, and clearing them would be wasted work
     -- a chain is only ever entered through its head, so every link followed
     at look-up time was written by this call. */
  kmer_hash.chain_head.assign(static_cast<std::size_t>(kmers), 0U);
  if (kmer_hash.chain_next.size() <= seq.size())
    {
      kmer_hash.chain_next.resize(seq.size() + 1, 0U);
    }

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;

  int pos = 0;
  for (auto const nucleotide : seq)
    {
      bad <<= 2ULL;
      bad |= mask_ambig::map(nucleotide);
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= two_bit::map(nucleotide);
      kmer &= kmer_mask;

      if (bad == 0U)
        {
          /* 1-based pos of start of kmer, biased by one */
          auto const entry = static_cast<unsigned int>(pos - k_offset + 1 + 1);
          kmer_hash.chain_next[entry] = kmer_hash.chain_head[kmer];
          kmer_hash.chain_head[kmer] = entry;
        }
      ++pos;
    }
}


auto kh_find_diagonals(struct kh_handle_s const & kmer_hash,
                       int const k_offset,
                       View<char> const seq,
                       std::vector<int> & diags) -> void
{
  assert(k_offset > 0);
  assert(k_offset <= kmer_hash_max_k);

  int const kmers = static_cast<int>(1U << (2U * static_cast<unsigned int>(k_offset)));
  auto const kmer_mask = static_cast<unsigned int>(kmers - 1);

  unsigned int bad = kmer_mask;
  unsigned int kmer = 0;

  /* the diagonal is a signed offset into diags, so the length is taken
     as an int once here rather than at each of the call sites */
  int const len = static_cast<int>(seq.size());

  auto const * complement_map = chrmap_complement();
  auto seq_cursor = seq.crbegin();
  for (int pos = 0; pos < len; pos++)
    {
      char const nucleotide = *seq_cursor;
      ++seq_cursor;

      bad <<= 2ULL;
      bad |= mask_ambig::map(nucleotide);
      bad &= kmer_mask;

      kmer <<= 2ULL;
      kmer |= two_bit::get_map()[complement_map[maps::to_uchar(nucleotide)]];
      kmer &= kmer_mask;

      if (bad == 0U)
        {
          /* walk the forward positions carrying exactly this k-mer. The
             diagonal counters are a bag, so the order in which the chain
             visits those positions cannot change the result. */
          int const base_diag = len - (pos - k_offset + 1);
          for (auto entry = kmer_hash.chain_head[kmer];
               entry != 0U;
               entry = kmer_hash.chain_next[entry])
            {
              /* 'entry - 1' undoes the bias, giving the 1-based start
                 position of the k-mer in the forward read */
              int const diag = base_diag + static_cast<int>(entry) - 1;
              if (diag >= 0)
                {
                  ++diags[static_cast<std::size_t>(diag)];
                }
            }
        }
    }
}

