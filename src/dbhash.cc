/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "vsearch5d.h"

static bitmap_t * dbhash_bitmap;
static unsigned long dbhash_size;
static unsigned int dbhash_shift;
static unsigned long dbhash_mask;
static struct dbhash_bucket_s * dbhash_table;

int dbhash_seqcmp(char * a, char * b, unsigned long n)
{
  char * p = a;
  char * q = b;

  if (n <= 0)
    return 0;

  while ((n-- > 0) && (chrmap_4bit[(int)(*p)] == chrmap_4bit[(int)(*q)]))
    {
      if ((n == 0) || (*p == 0) || (*q == 0))
        break;
      p++;
      q++;
    }

  return chrmap_4bit[(int)(*p)] - chrmap_4bit[(int)(*q)];
}

void dbhash_open(unsigned long maxelements)
{
  /* adjust size of hash table for 2/3 fill rate */
  /* and use a multiple of 2 */

  dbhash_size = 1;
  dbhash_shift = 0;
  while (3 * maxelements > 2 * dbhash_size)
    {
      dbhash_size <<= 1;
      dbhash_shift++;
    }
  dbhash_mask = dbhash_size - 1;
  
  dbhash_table = (struct dbhash_bucket_s *)
    xmalloc(sizeof(dbhash_bucket_s) * dbhash_size);
  memset(dbhash_table, 0, sizeof(dbhash_bucket_s) * dbhash_size);

  dbhash_bitmap = bitmap_init(dbhash_size);
  bitmap_reset_all(dbhash_bitmap);
}

void dbhash_close()
{
  bitmap_free(dbhash_bitmap);
  dbhash_bitmap = 0;
  free(dbhash_table);
  dbhash_table = 0;
}

long dbhash_search_first(char * seq,
                        unsigned long seqlen,
                        struct dbhash_search_info_s * info)
{
  
  unsigned long hash = hash_cityhash64(seq, seqlen);
  info->hash = hash;
  info->seq = seq;
  info->seqlen = seqlen;
  unsigned long index = hash & dbhash_mask;
  struct dbhash_bucket_s * bp = dbhash_table + index;

  while (bitmap_get(dbhash_bitmap, index)
         &&
         ((bp->hash != hash) ||
          (seqlen != db_getsequencelen(bp->seqno)) ||
          (dbhash_seqcmp(seq, db_getsequence(bp->seqno), seqlen))))
    {
      index = (index + 1) & dbhash_mask;
      bp = dbhash_table + index;
    }

  info->index = index;
  
  if (bitmap_get(dbhash_bitmap, index))
    return bp->seqno;
  else
    return -1;
}

long dbhash_search_next(struct dbhash_search_info_s * info)
{
  unsigned long hash = info->hash;
  char * seq = info->seq;
  unsigned long seqlen = info->seqlen;
  unsigned long index = (info->index + 1) & dbhash_mask;
  struct dbhash_bucket_s * bp = dbhash_table + index;

  while (bitmap_get(dbhash_bitmap, index)
         &&
         ((bp->hash != hash) ||
          (seqlen != db_getsequencelen(bp->seqno)) ||
          (dbhash_seqcmp(seq, db_getsequence(bp->seqno), seqlen))))
    {
      index = (index + 1) & dbhash_mask;
      bp = dbhash_table + index;
    }

  info->index = index;
  
  if (bitmap_get(dbhash_bitmap, index))
    return bp->seqno;
  else
    return -1;
}

void dbhash_add(char * seq, unsigned long seqlen, unsigned long seqno)
{
  struct dbhash_search_info_s info;
  
  long ret = dbhash_search_first(seq, seqlen, & info);
  while (ret >= 0)
    ret = dbhash_search_next(&info);
  
  bitmap_set(dbhash_bitmap, info.index);
  struct dbhash_bucket_s * bp = dbhash_table + info.index;
  bp->hash = info.hash;
  bp->seqno = seqno;
}

void dbhash_add_one(unsigned long seqno)
{
  char * seq = db_getsequence(seqno);
  unsigned long seqlen = db_getsequencelen(seqno);
  char * normalized = (char*) xmalloc(seqlen+1);
  string_normalize(normalized, seq, seqlen);
  dbhash_add(normalized, seqlen, seqno);
}

void dbhash_add_all()
{
  progress_init("Hashing database sequences", db_getsequencecount());
  char * normalized = (char*) xmalloc(db_getlongestsequence()+1);
  for(unsigned long seqno=0; seqno < db_getsequencecount(); seqno++)
    {
      char * seq = db_getsequence(seqno);
      unsigned long seqlen = db_getsequencelen(seqno);
      string_normalize(normalized, seq, seqlen);
      dbhash_add(normalized, seqlen, seqno);
      progress_update(seqno+1);
    }
  free(normalized);
  progress_done();
}
