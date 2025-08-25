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
#include "bitmap.h"
#include "maps.h"
#include <cstdint>  // int64_t, uint64_t
#include <cstring>  // std::memset
#include <vector>


static struct bitmap_s * dbhash_bitmap;
static uint64_t dbhash_size;
static unsigned int dbhash_shift;
static uint64_t dbhash_mask;
std::vector<struct dbhash_bucket_s> dbhash_table;


auto dbhash_seqcmp(char const * a, char const * b, uint64_t n) -> int
{
  char const * p = a;
  char const * q = b;

  if (n <= 0)
    {
      return 0;
    }

  while ((n-- > 0) and (chrmap_4bit[(int) (*p)] == chrmap_4bit[(int) (*q)]))
    {
      if ((n == 0) or (*p == 0) or (*q == 0))
        {
          break;
        }
      ++p;
      ++q;
    }

  return chrmap_4bit[(int) (*p)] - chrmap_4bit[(int) (*q)];
}


auto dbhash_open(uint64_t const maxelements) -> void
{
  /* adjust size of hash table for 2/3 fill rate */
  /* and use a multiple of 2 */

  dbhash_size = 1;
  dbhash_shift = 0;
  while (3 * maxelements > 2 * dbhash_size)
    {
      dbhash_size <<= 1U;
      ++dbhash_shift;
    }
  dbhash_mask = dbhash_size - 1;

  dbhash_table.resize(dbhash_size);

  dbhash_bitmap = bitmap_init(dbhash_size);
  bitmap_reset_all(dbhash_bitmap);
}


auto dbhash_close() -> void
{
  bitmap_free(dbhash_bitmap);
  dbhash_bitmap = nullptr;
}


auto dbhash_search_first(char * seq,
                         uint64_t const seqlen,
                         struct dbhash_search_info_s * info) -> int64_t
{
  auto const hash = hash_cityhash64(seq, seqlen);
  info->hash = hash;
  info->seq = seq;
  info->seqlen = seqlen;
  auto index = hash & dbhash_mask;
  auto * bp = &dbhash_table[index];

  while ((bitmap_get(dbhash_bitmap, index) != 0U)
         and
         ((bp->hash != hash) or
          (seqlen != db_getsequencelen(bp->seqno)) or
          (dbhash_seqcmp(seq, db_getsequence(bp->seqno), seqlen) != 0)))
    {
      index = (index + 1) & dbhash_mask;
      bp = &dbhash_table[index];
    }

  info->index = index;

  if (bitmap_get(dbhash_bitmap, index) != 0U)
    {
      return bp->seqno;
    }
  return -1;
}


auto dbhash_search_next(struct dbhash_search_info_s * info) -> int64_t
{
  auto const hash = info->hash;
  auto * seq = info->seq;
  auto const seqlen = info->seqlen;
  auto index = (info->index + 1) & dbhash_mask;
  auto * bp = &dbhash_table[index];

  while ((bitmap_get(dbhash_bitmap, index) != 0U)
         and
         ((bp->hash != hash) or
          (seqlen != db_getsequencelen(bp->seqno)) or
          (dbhash_seqcmp(seq, db_getsequence(bp->seqno), seqlen) != 0)))
    {
      index = (index + 1) & dbhash_mask;
      bp = &dbhash_table[index];
    }

  info->index = index;

  if (bitmap_get(dbhash_bitmap, index) != 0U)
    {
      return bp->seqno;
    }
  return -1;
}


auto dbhash_add(char * seq, uint64_t seqlen, uint64_t seqno) -> void
{
  struct dbhash_search_info_s info;

  auto ret = dbhash_search_first(seq, seqlen, &info);
  while (ret >= 0)
    {
      ret = dbhash_search_next(&info);
    }

  bitmap_set(dbhash_bitmap, info.index);
  auto & bucket = dbhash_table[info.index];
  bucket.hash = info.hash;
  bucket.seqno = seqno;
}


auto dbhash_add_all() -> void
{
  progress_init("Hashing database sequences", db_getsequencecount());
  std::vector<char> normalized(db_getlongestsequence() + 1);
  for (uint64_t seqno = 0; seqno < db_getsequencecount(); ++seqno)
    {
      auto const * seq = db_getsequence(seqno);
      auto const seqlen = db_getsequencelen(seqno);
      string_normalize(normalized.data(), seq, seqlen);
      dbhash_add(normalized.data(), seqlen, seqno);
      progress_update(seqno + 1);
    }
  progress_done();
}
