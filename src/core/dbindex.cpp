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

#include "core/mask.hpp"
#include "vsearch.hpp"
#include "core/bitmap.hpp"
#include "core/db.hpp"
#include "core/dbindex.hpp"
#include "core/unique.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include <algorithm>  // std::max
#include <array>
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE
#include <iterator>  // std::next
#include <vector>


constexpr unsigned int bitmap_threshold = 8;


auto bitmap_min_matches(unsigned int const seqcount) noexcept -> unsigned int
{
  /* A bitmap pays for itself only for a k-mer shared by a large share of the
     database, and never for a k-mer that does not occur at all, hence the floor
     of one match. Without it a database of fewer than bitmap_threshold
     sequences divides down to a threshold of zero, which every one of the
     4^wordlength slots meets, absent k-mers included: at --orient's default
     word length of 12 that is 16.7 M bitmaps (+0.5 GB, and a full SIMD scan
     per absent k-mer at search time) in place of the handful of k-mers the
     sequences really contain. The floor cannot change which k-mers get a
     bitmap: a k-mer that occurs has one match or more, so it falls on the same
     side of a threshold of zero and of one alike. */
  return std::max(1U, seqcount / bitmap_threshold);
}


auto Dbindex::getbitmap(unsigned int const kmer) const -> unsigned char const *
{
  auto const & a_bitmap = kmerbitmap[kmer];
  if (not a_bitmap.empty())
    {
      return a_bitmap.data();
    }
  return nullptr;
}


auto Dbindex::getmatchcount(unsigned int const kmer) const -> unsigned int
{
  return kmercount[kmer];
}


auto Dbindex::getmatchlist(unsigned int const kmer) const -> unsigned int const *
{
  return std::next(kmerindex.data(), static_cast<std::iterator_traits<unsigned int const *>::difference_type>(kmerhash[kmer]));
}


auto Dbindex::getmapping(unsigned int const index) const -> unsigned int
{
  return map[index];
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
      fprint(output_handle, sym_nt_2bit[(kmer >> (2 * (kmer_length - i - 1))) & 3]);
    }
}

auto Dbindex::add_sequence(unsigned int const seqno, Masking const seqmask, struct Database const & db) -> void
{
  auto const uniquelist = uhandle.count(static_cast<int>(wordlength),
                                        db.sequence_view(seqno), seqmask);
  map[count] = seqno;
  for (auto const kmer : uniquelist)
    {
      if (not kmerbitmap[kmer].empty())
        {
          ++kmercount[kmer];
          kmerbitmap[kmer].set(count);
        }
      else
        {
          kmerindex[kmerhash[kmer] + kmercount[kmer]] = count;
          ++kmercount[kmer];
        }
    }
  ++count;
}


auto Dbindex::add_all_sequences(Masking const seqmask, struct Database const & db, struct Parameters const & parameters) -> void
{
  auto const seqcount = static_cast<unsigned int>(db.getsequencecount());
  Progress progress("Creating k-mer index", seqcount, parameters);
  for (auto seqno = 0U; seqno < seqcount ; seqno++)
    {
      add_sequence(seqno, seqmask, db);
      progress.update(seqno);
    }
}


auto Dbindex::prepare(int const use_bitmap, Masking const seqmask, struct Database const & db, struct Parameters const & parameters) -> void
{
  /* Release any state from a previous prepare first (mirrors Database::init ->
     clear()), so a second prepare without an intervening clear() does
     not leak the earlier five buffers. clear() is a no-op on the
     first call (all members are null) (L2a). */
  clear();

  auto const seqcount = static_cast<unsigned int>(db.getsequencecount());
  /* this is the FASTA-database path; the effective index word length is the
     configured one (a UDB database sets wordlength in udb_read instead). */
  wordlength = static_cast<unsigned int>(parameters.opt_wordlength);
  hashsize = 1U << (2 * wordlength);

  /* allocate memory for kmer count array */
  kmercount.assign(hashsize, 0U);

  /* first scan, just count occurrences */
  {
    Progress progress("Counting k-mers", seqcount, parameters);
    for (auto seqno = 0U; seqno < seqcount ; seqno++)
      {
        auto const uniquelist = uhandle.count(static_cast<int>(wordlength),
                                              db.sequence_view(seqno), seqmask);
        for (auto const kmer : uniquelist)
          {
            ++kmercount[kmer];
          }
        progress.update(seqno);
      }
  }

  /* determine minimum kmer count for bitmap usage */
  unsigned int const bitmap_mincount = (use_bitmap != 0) ? bitmap_min_matches(seqcount) : (seqcount + 1);

  /* allocate empty (list-form) bitmap slots for every kmer */
  kmerbitmap = std::vector<Bitmap>(hashsize);

  /* hash / bitmap setup */
  /* convert hash counts to position in index */
  kmerhash.resize(hashsize + 1);
  uint64_t sum = 0;
  for (auto i = 0U; i < hashsize; i++)
    {
      kmerhash[i] = sum;
      if (kmercount[i] >= bitmap_mincount)
        {
          kmerbitmap[i] = Bitmap(seqcount + 127); // pad for xmm
        }
      else
        {
          sum += kmercount[i];
        }
    }
  indexsize = sum;
  kmerhash[hashsize] = sum;

  /* reset counts */
  kmercount.assign(hashsize, 0U);

  /* allocate space for actual data */
  kmerindex.resize(indexsize);

  /* allocate space for mapping from indexno to seqno */
  map.resize(seqcount);

  count = 0;

  // memory-intensive: the k-mer index has been allocated
}


auto Dbindex::clear() -> void
{
  /* Release every owned buffer so the routine is idempotent (a second call, or
     a call before any successful prepare, is a safe no-op) and a subsequent
     prepare() starts from a clean slate. Every buffer is now a RAII value
     member: the k-mer arrays and bitmaps are std::vector, and the unique-kmer
     finder is a Uniquer. clear() drops the contents and shrink_to_fit() returns
     the capacity to the allocator, so a rebuild does not carry the previous
     index's peak memory. */
  kmerhash.clear();   kmerhash.shrink_to_fit();
  kmerindex.clear();  kmerindex.shrink_to_fit();
  kmercount.clear();  kmercount.shrink_to_fit();
  map.clear();        map.shrink_to_fit();

  kmerbitmap.clear();
  kmerbitmap.shrink_to_fit();

  uhandle = Uniquer();

  hashsize = 0;
  indexsize = 0;
}


Dbindex::~Dbindex()
{
  clear();
}
