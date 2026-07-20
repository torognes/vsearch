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

#include "utils/view.hpp"
#include "utils/span.hpp"
#include "vsearch.hpp"
#include "utils/progress.hpp"
#include "core/bitmap.hpp"
#include "core/db.hpp"
#include "core/dbhash.hpp"
#include "utils/seqcmp.hpp"
#include "utils/cityhash.hpp"
#include "utils/string_normalize.hpp"
#include <cstdint>  // int64_t, uint64_t
#include <cstring>  // std::memset
#include <vector>


auto Dbhash::open(uint64_t const maxelements) -> void
{
  /* adjust size of hash table for 2/3 fill rate */
  /* and use a multiple of 2 */

  uint64_t size = 1;
  while (3 * maxelements > 2 * size)
    {
      size <<= 1U;
    }
  mask_ = size - 1;

  table_.resize(size);

  bitmap_ = Bitmap(static_cast<unsigned int>(size));
}


Dbhash::~Dbhash()
{
  clear();
}


auto Dbhash::clear() -> void
{
  /* Release the (potentially large) table and bitmap now rather than holding
     them until the Dbhash is destroyed, and reset to the empty state so the
     index can be reopened. Safe on an instance that was never open()ed. */
  bitmap_ = Bitmap();
  table_.clear();
  table_.shrink_to_fit();
  mask_ = 0;
}


auto Dbhash::search_first(char * seq,
                          uint64_t const seqlen,
                          struct dbhash_search_info_s * info,
                          struct Database const & db) const -> int64_t
{
  auto const hash = hash_cityhash64(seq, seqlen);
  info->hash = hash;
  info->seq = seq;
  info->seqlen = seqlen;
  auto index = hash & mask_;
  auto const * bp = &table_[index];

  while ((bitmap_.is_set(static_cast<unsigned int>(index)))
         and
         ((bp->hash != hash) or
          (seqlen != db.getsequencelen(bp->seqno)) or
          (seqcmp(View<char>{seq, static_cast<std::size_t>(seqlen)}, View<char>{db.getsequence(bp->seqno), static_cast<std::size_t>(seqlen)}) != 0)))
    {
      index = (index + 1) & mask_;
      bp = &table_[index];
    }

  info->index = index;

  if (bitmap_.is_set(static_cast<unsigned int>(index)))
    {
      return static_cast<int64_t>(bp->seqno);
    }
  return -1;
}


auto Dbhash::search_next(struct dbhash_search_info_s * info, struct Database const & db) const -> int64_t
{
  auto const hash = info->hash;
  auto const * seq = info->seq;
  auto const seqlen = info->seqlen;
  auto index = (info->index + 1) & mask_;
  auto const * bp = &table_[index];

  while ((bitmap_.is_set(static_cast<unsigned int>(index)))
         and
         ((bp->hash != hash) or
          (seqlen != db.getsequencelen(bp->seqno)) or
          (seqcmp(View<char>{seq, static_cast<std::size_t>(seqlen)}, View<char>{db.getsequence(bp->seqno), static_cast<std::size_t>(seqlen)}) != 0)))
    {
      index = (index + 1) & mask_;
      bp = &table_[index];
    }

  info->index = index;

  if (bitmap_.is_set(static_cast<unsigned int>(index)))
    {
      return static_cast<int64_t>(bp->seqno);
    }
  return -1;
}


auto Dbhash::add(char * seq, uint64_t const seqlen, uint64_t const seqno, struct Database const & db) -> void
{
  struct dbhash_search_info_s info;

  auto ret = search_first(seq, seqlen, &info, db);
  while (ret >= 0)
    {
      ret = search_next(&info, db);
    }

  bitmap_.set(static_cast<unsigned int>(info.index));
  auto & bucket = table_[info.index];
  bucket.hash = info.hash;
  bucket.seqno = seqno;
}


auto Dbhash::add_all(struct Database const & db, struct Parameters const & parameters) -> void
{
  Progress progress("Hashing database sequences", db.getsequencecount(), parameters);
  std::vector<char> normalized(db.getlongestsequence() + 1);
  for (uint64_t seqno = 0; seqno < db.getsequencecount(); ++seqno)
    {
      auto const * seq = db.getsequence(seqno);
      auto const seqlen = db.getsequencelen(seqno);
      string_normalize(Span<char>{normalized.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq, static_cast<std::size_t>(seqlen)});
      add(normalized.data(), seqlen, seqno, db);
      progress.update(seqno + 1);
    }
}
