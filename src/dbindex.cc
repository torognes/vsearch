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

#include "vsearch.h"
#include "bitmap.h"
#include "dbindex.h"
#include "unique.h"
#include "utils/progress.hpp"
#include <array>
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstring>  // std::memset
#include <iterator>  // std::next


constexpr unsigned int bitmap_threshold = 8;


auto Dbindex::getbitmap(unsigned int const kmer) const -> unsigned char *
{
  auto * a_bitmap_s = *std::next(kmerbitmap, kmer);
  if (a_bitmap_s != nullptr)
    {
      return a_bitmap_s->bitmap;
    }
  return nullptr;
}


auto Dbindex::getmatchcount(unsigned int const kmer) const -> unsigned int
{
  return *std::next(kmercount, kmer);
}


auto Dbindex::getmatchlist(unsigned int const kmer) const -> unsigned int *
{
  return std::next(kmerindex, static_cast<std::iterator_traits<unsigned int *>::difference_type>(*std::next(kmerhash, kmer)));
}


auto Dbindex::getmapping(unsigned int const index) const -> unsigned int
{
  return *std::next(map, index);
}


auto Dbindex::getcount() const -> unsigned int
{
  return count;
}


auto fprint_kmer(std::FILE * output_handle, unsigned int const kmer_length, uint64_t const kmer) -> void
{
  static constexpr std::array<char, 4> sym_nt_2bit = {{'A', 'C', 'G', 'T'}};
  for (auto i = 0U; i < kmer_length; ++i)
    {
      std::fprintf(output_handle, "%c", sym_nt_2bit[(kmer >> (2 * (kmer_length - i - 1))) & 3]);
    }
}

auto Dbindex::add_sequence(unsigned int seqno, int seqmask) -> void
{
#if 0
  std::printf("Adding seqno %d as index element no %d\n", seqno, count);
#endif

  unsigned int uniquecount = 0;
  unsigned int const * uniquelist = nullptr;
  unique_count(uhandle, static_cast<int>(wordlength),
               static_cast<int>(db_getsequencelen(seqno)), db_getsequence(seqno),
               &uniquecount, &uniquelist, seqmask);
  map[count] = seqno;
  for (auto i = 0U; i < uniquecount; i++)
    {
      auto const kmer = uniquelist[i];
      if (kmerbitmap[kmer] != nullptr)
        {
          ++kmercount[kmer];
          bitmap_set(kmerbitmap[kmer], count);
        }
      else
        {
          kmerindex[kmerhash[kmer] + kmercount[kmer]] = count;
          ++kmercount[kmer];
        }
    }
  ++count;
}


auto Dbindex::add_all_sequences(int seqmask, struct Parameters const & parameters) -> void
{
  unsigned int const seqcount = static_cast<unsigned int>(db_getsequencecount());
  Progress progress("Creating k-mer index", seqcount, parameters);
  for (auto seqno = 0U; seqno < seqcount ; seqno++)
    {
      add_sequence(seqno, seqmask);
      progress.update(seqno);
    }
}


auto Dbindex::prepare(int use_bitmap, int seqmask, struct Parameters const & parameters) -> void
{
  /* Release any state from a previous prepare first (mirrors db_init ->
     db_free), so a second prepare without an intervening clear() does
     not leak the earlier five buffers. clear() is a no-op on the
     first call (all members are null) (L2a). */
  clear();

  uhandle = unique_init();

  unsigned int const seqcount = static_cast<unsigned int>(db_getsequencecount());
  /* this is the FASTA-database path; the effective index word length is the
     configured one (a UDB database sets wordlength in udb_read instead). */
  wordlength = static_cast<unsigned int>(parameters.opt_wordlength);
  hashsize = 1U << (2 * wordlength);

  /* allocate memory for kmer count array */
  kmercount = static_cast<unsigned int *>(xmalloc(hashsize * sizeof(unsigned int)));
  std::memset(kmercount, 0, hashsize * sizeof(unsigned int));

  /* first scan, just count occurences */
  {
    Progress progress("Counting k-mers", seqcount, parameters);
    for (auto seqno = 0U; seqno < seqcount ; seqno++)
      {
        unsigned int uniquecount = 0;
        unsigned int const * uniquelist = nullptr;
        unique_count(uhandle, static_cast<int>(wordlength),
                     static_cast<int>(db_getsequencelen(seqno)), db_getsequence(seqno),
                     &uniquecount, &uniquelist, seqmask);
        for (auto i = 0U; i < uniquecount; i++)
          {
            ++kmercount[uniquelist[i]];
          }
        progress.update(seqno);
      }
  }

#if 0
  /* dump kmer counts */
  auto * f = fopen_output("kmercounts.txt");
  for (auto kmer = 0U; kmer < hashsize; kmer++)
    {
      fprint_kmer(f, 8, kmer);
      std::fprintf(f, "\t%d\t%d\n", kmer, kmercount[kmer]);
    }
  fclose_output(f);
#endif

  /* determine minimum kmer count for bitmap usage */
  unsigned int const bitmap_mincount = (use_bitmap != 0) ? (seqcount / bitmap_threshold) : (seqcount + 1);

  /* allocate and zero bitmap pointers */
  kmerbitmap = static_cast<struct bitmap_s **>(xmalloc(hashsize * sizeof(struct bitmap_s *)));
  std::memset(kmerbitmap, 0, hashsize * sizeof(struct bitmap_s *));

  /* hash / bitmap setup */
  /* convert hash counts to position in index */
  kmerhash = static_cast<uint64_t *>(xmalloc((hashsize + 1) * sizeof(uint64_t)));
  uint64_t sum = 0;
  for (auto i = 0U; i < hashsize; i++)
    {
      kmerhash[i] = sum;
      if (kmercount[i] >= bitmap_mincount)
        {
          kmerbitmap[i] = bitmap_init(seqcount + 127); // pad for xmm
          bitmap_reset_all(kmerbitmap[i]);
        }
      else
        {
          sum += kmercount[i];
        }
    }
  indexsize = sum;
  kmerhash[hashsize] = sum;

#if 0
  if (not parameters.opt_quiet)
    std::fprintf(stderr, "Unique %u-mers: %u\n", wordlength, indexsize);
#endif

  /* reset counts */
  std::memset(kmercount, 0, hashsize * sizeof(unsigned int));

  /* allocate space for actual data */
  kmerindex = static_cast<unsigned int *>(xmalloc(indexsize * sizeof(unsigned int)));

  /* allocate space for mapping from indexno to seqno */
  map = static_cast<unsigned int *>(xmalloc(seqcount * sizeof(unsigned int)));

  count = 0;

  // memory-intensive: the k-mer index has been allocated
}


auto Dbindex::clear() -> void
{
  /* Free and null every owned buffer so the routine is idempotent (a second
     call, or a call before any successful prepare, is a safe no-op) and a
     subsequent prepare() starts from a clean slate. xfree() fatals on
     a null pointer, so each free is guarded; the bitmap loop and unique handle
     are likewise guarded because they may not have been allocated yet (L2a). */
  if (kmerhash != nullptr) { xfree(kmerhash); kmerhash = nullptr; }
  if (kmerindex != nullptr) { xfree(kmerindex); kmerindex = nullptr; }
  if (kmercount != nullptr) { xfree(kmercount); kmercount = nullptr; }
  if (map != nullptr) { xfree(map); map = nullptr; }

  if (kmerbitmap != nullptr)
    {
      for (auto kmer = 0U; kmer < hashsize; kmer++)
        {
          if (kmerbitmap[kmer] != nullptr)
            {
              bitmap_free(kmerbitmap[kmer]);
            }
        }
      xfree(kmerbitmap);
      kmerbitmap = nullptr;
    }

  if (uhandle != nullptr)
    {
      unique_exit(uhandle);
      uhandle = nullptr;
    }

  hashsize = 0;
  indexsize = 0;
}


Dbindex::~Dbindex()
{
  clear();
}
