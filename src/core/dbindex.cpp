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
#include "utils/decimal_digits.hpp"  // decimal::to_text
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/warn.hpp"  // vsearch::warn
#include <algorithm>  // std::max
#include <array>
#include <cassert>  // assert
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE
#include <iterator>  // std::next
#include <string>  // std::string
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


auto Dbindex::bitmap_slots_reset(unsigned int const slots) -> void
{
  kmerbitmap_slot.assign(slots, 0U);
  bitmap_pool.clear();
  bitmap_pool.shrink_to_fit();
  bitmap_width = 0;
}


auto Dbindex::set_bitmap_width(unsigned int const sequences) -> void
{
  /* pad for xmm: the SIMD counter routines read whole registers, so they may
     touch up to 127 bits past the last sequence */
  static constexpr unsigned int simd_padding = 127;
  bitmap_width = sequences + simd_padding;
}


auto Dbindex::bitmap_create(unsigned int const kmer) -> Bitmap &
{
  assert(kmer < kmerbitmap_slot.size());
  assert(kmerbitmap_slot[kmer] == 0U);  /* at most one bitmap per k-mer */
  assert(bitmap_width != 0U);  /* set_bitmap_width comes first */
  bitmap_pool.emplace_back(bitmap_width);
  /* at most one bitmap per slot, so the 1-based index fits an unsigned int for
     every supported word length (4^15 slots at most) */
  assert(bitmap_pool.size() <= kmerbitmap_slot.size());
  kmerbitmap_slot[kmer] = static_cast<unsigned int>(bitmap_pool.size());
  return bitmap_pool.back();
}


auto Dbindex::has_bitmap(unsigned int const kmer) const -> bool
{
  assert(kmer < kmerbitmap_slot.size());
  return kmerbitmap_slot[kmer] != 0U;
}


auto Dbindex::bitmap_of(unsigned int const kmer) -> Bitmap &
{
  assert(has_bitmap(kmer));
  return bitmap_pool[kmerbitmap_slot[kmer] - 1];
}


auto Dbindex::bitmap_of(unsigned int const kmer) const -> Bitmap const &
{
  assert(has_bitmap(kmer));
  return bitmap_pool[kmerbitmap_slot[kmer] - 1];
}


auto Dbindex::getbitmap(unsigned int const kmer) const -> unsigned char const *
{
  auto const slot = kmerbitmap_slot[kmer];
  if (slot != 0U)
    {
      return bitmap_pool[slot - 1].data();
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
      if (has_bitmap(kmer))
        {
          ++kmercount[kmer];
          bitmap_of(kmer).set(count);
        }
      else
        {
          kmerindex[kmerhash[kmer] + kmercount[kmer]] = count;
          ++kmercount[kmer];
        }
    }
  ++count;
}


/* anonymous namespace: limit visibility and usage to this translation unit */
namespace {

  /* A sequence whose unique-k-mer list comes back empty is absent from the
     k-mer index, so it is never selected as a candidate target. The message
     says exactly that and no more: such a sequence can still be *reached*, by
     a query that samples no k-mers itself and is therefore compared against
     every database sequence (documented in option_qmask.md), so "can never be
     matched" would be an overclaim. Three situations produce an empty list --
     every k-mer overlaps a masked (lowercase) base, every k-mer contains a
     symbol outside ACGTU, or the sequence is shorter than the word length (see
     Uniquer::count(), core/unique.cpp) -- and telling them apart would cost a
     second pass over the whole database (a count() with Masking::none), so the
     warning reports the observable fact and lists the possible causes. When
     masking is off, lowercase cannot be one of them and is left out.

     One summary line, not one warning per sequence: a low-complexity or
     heavily masked dataset would otherwise emit millions of them, which is why
     the existing discard reports (core/db.cpp) are shaped the same way.

     Motivating case: torognes/vsearch#570, where the all-lowercase BOLD COI
     reference database was silently indexed as 2.2 M unhittable sequences and
     --sintax produced 13 days of empty classifications. */
  /* The two counts travel together and are both plain integers, so they are
     named in a struct rather than passed as two adjacent swappable arguments
     (the shape FilterCounts uses in core/filter.cpp; clang-tidy's
     readability-suspicious-call-argument flags the two-argument form). No
     default member initializers: they would make this a non-aggregate before
     C++14, and the call site brace-initializes it. */
  struct KmerlessCounts
  {
    uint64_t without_kmers;
    uint64_t sequences;
  };


  auto warn_sequences_without_kmers(KmerlessCounts const & counts,
                                    unsigned int const wordlength,
                                    Masking const seqmask) -> void
  {
    if (counts.without_kmers == 0)
      {
        return;
      }
    std::string message = decimal::to_text(counts.without_kmers) + " of "
      + decimal::to_text(counts.sequences)
      + " sequences yielded no k-mer for the index and cannot be selected as"
        " candidate targets.\n"
      + "Possible causes: all their k-mers ";
    if (seqmask != Masking::none)
      {
        message += "are masked (lowercase) or ";
      }
    message += "contain a symbol other than A, C, G, T or U, or the sequences"
      " are shorter than the word length (" + decimal::to_text(wordlength) + ").";
    vsearch::warn(message);
  }

}  // end of anonymous namespace


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


auto Dbindex::prepare(Masking const seqmask, struct Database const & db, struct Parameters const & parameters) -> void
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
  uint64_t sequences_without_kmers = 0;
  {
    Progress progress("Counting k-mers", seqcount, parameters);
    for (auto seqno = 0U; seqno < seqcount ; seqno++)
      {
        auto const uniquelist = uhandle.count(static_cast<int>(wordlength),
                                              db.sequence_view(seqno), seqmask);
        if (uniquelist.empty())
          {
            ++sequences_without_kmers;
          }
        for (auto const kmer : uniquelist)
          {
            ++kmercount[kmer];
          }
        progress.update(seqno);
      }
  }

  /* after the Progress object is destroyed, so the warning does not land in
     the middle of a progress line (warn() opens with a newline anyway) */
  warn_sequences_without_kmers(KmerlessCounts{sequences_without_kmers, seqcount},
                               wordlength, seqmask);

  /* determine minimum kmer count for bitmap usage */
  unsigned int const bitmap_mincount = bitmap_min_matches(seqcount);

  /* allocate empty (list-form) bitmap slots for every kmer */
  bitmap_slots_reset(hashsize);
  set_bitmap_width(seqcount);

  /* hash / bitmap setup */
  /* convert hash counts to position in index */
  /* filled by push_back into reserved space rather than resize + assign,
     because resize would value-initialise every entry (1 GB of zeros at word
     length 15) only for the loop below to overwrite all of them */
  kmerhash.clear();
  kmerhash.reserve(hashsize + 1);
  uint64_t sum = 0;
  for (auto i = 0U; i < hashsize; i++)
    {
      kmerhash.push_back(sum);
      if (kmercount[i] >= bitmap_mincount)
        {
          bitmap_create(i);
        }
      else
        {
          sum += kmercount[i];
        }
    }
  indexsize = sum;
  kmerhash.push_back(sum);
  /* one entry per slot plus the end marker. With the buffer filled rather than
     value-initialised, size() is the number of entries actually written, so
     this checks the loop covered every slot -- the guarantee resize() used to
     provide by zeroing. */
  assert(kmerhash.size() == hashsize + 1U);

  /* reset counts */
  /* left as one bulk assign rather than folded into the loop above: a single
     fill of 4 bytes per slot is cheaper than a checked store per slot, sharply
     so in a _GLIBCXX_DEBUG build (measured: folding it in cost 49 ms per index
     build at word length 12, and gained nothing in a release build) */
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

  kmerbitmap_slot.clear();  kmerbitmap_slot.shrink_to_fit();
  bitmap_pool.clear();      bitmap_pool.shrink_to_fit();
  bitmap_width = 0;

  uhandle = Uniquer();

  hashsize = 0;
  indexsize = 0;
}


Dbindex::~Dbindex()
{
  clear();
}
