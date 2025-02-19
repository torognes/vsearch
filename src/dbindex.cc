/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "maps.h"
#include "unique.h"
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstring>  // std::memset
#include <iterator>  // std::next


unsigned int * kmercount;
uint64_t * kmerhash;
unsigned int * kmerindex;
struct bitmap_s * * kmerbitmap;
unsigned int * dbindex_map;
unsigned int kmerhashsize;
uint64_t kmerindexsize;
unsigned int dbindex_count;
uhandle_s * dbindex_uh;

constexpr unsigned int bitmap_threshold = 8;

static unsigned int bitmap_mincount;


auto dbindex_getbitmap(unsigned int const kmer) -> unsigned char *
{
  auto * a_bitmap_s = *std::next(kmerbitmap, kmer);
  if (a_bitmap_s != nullptr)
    {
      return a_bitmap_s->bitmap;
    }
  return nullptr;
}


auto dbindex_getmatchcount(unsigned int const kmer) -> unsigned int
{
  return *std::next(kmercount, kmer);
}


auto dbindex_getmatchlist(unsigned int const kmer) -> unsigned int *
{
  return std::next(kmerindex, *std::next(kmerhash, kmer));
}


auto dbindex_getmapping(unsigned int const index) -> unsigned int
{
  return *std::next(dbindex_map, index);
}


auto dbindex_getcount() -> unsigned int
{
  return dbindex_count;
}


auto fprint_kmer(std::FILE * output_handle, unsigned int const kmer_length, uint64_t const kmer) -> void
{
  for (auto i = 0U; i < kmer_length; ++i)
    {
      std::fprintf(output_handle, "%c", sym_nt_2bit[(kmer >> (2 * (kmer_length - i - 1))) & 3]);
    }
}

auto dbindex_addsequence(unsigned int seqno, int seqmask) -> void
{
#if 0
  std::printf("Adding seqno %d as index element no %d\n", seqno, dbindex_count);
#endif

  unsigned int uniquecount = 0;
  unsigned int * uniquelist = nullptr;
  unique_count(dbindex_uh, opt_wordlength,
               db_getsequencelen(seqno), db_getsequence(seqno),
               &uniquecount, &uniquelist, seqmask);
  dbindex_map[dbindex_count] = seqno;
  for (unsigned int i = 0; i < uniquecount; i++)
    {
      unsigned int const kmer = uniquelist[i];
      if (kmerbitmap[kmer])
        {
          kmercount[kmer]++;
          bitmap_set(kmerbitmap[kmer], dbindex_count);
        }
      else
        {
          kmerindex[kmerhash[kmer] + (kmercount[kmer]++)] = dbindex_count;
        }
    }
  ++dbindex_count;
}


auto dbindex_addallsequences(int seqmask) -> void
{
  unsigned int const seqcount = db_getsequencecount();
  progress_init("Creating k-mer index", seqcount);
  for (unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      dbindex_addsequence(seqno, seqmask);
      progress_update(seqno);
    }
  progress_done();
}


auto dbindex_prepare(int use_bitmap, int seqmask) -> void
{
  dbindex_uh = unique_init();

  unsigned int const seqcount = db_getsequencecount();
  kmerhashsize = 1U << (2 * opt_wordlength);

  /* allocate memory for kmer count array */
  kmercount = (unsigned int *) xmalloc(kmerhashsize * sizeof(unsigned int));
  std::memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  /* first scan, just count occurences */
  progress_init("Counting k-mers", seqcount);
  for (unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      unsigned int uniquecount = 0;
      unsigned int * uniquelist = nullptr;
      unique_count(dbindex_uh, opt_wordlength,
                   db_getsequencelen(seqno), db_getsequence(seqno),
                   &uniquecount, &uniquelist, seqmask);
      for (unsigned int i = 0; i < uniquecount; i++)
        {
          kmercount[uniquelist[i]]++;
        }
      progress_update(seqno);
    }
  progress_done();

#if 0
  /* dump kmer counts */
  std::FILE * f = fopen_output("kmercounts.txt");
  for (unsigned int kmer=0; kmer < kmerhashsize; kmer++)
    {
      fprint_kmer(f, 8, kmer);
      std::fprintf(f, "\t%d\t%d\n", kmer, kmercount[kmer]);
    }
  std::fclose(f);
#endif

  /* determine minimum kmer count for bitmap usage */
  if (use_bitmap)
    {
      bitmap_mincount = seqcount / bitmap_threshold;
    }
  else
    {
      bitmap_mincount = seqcount + 1;
    }

  /* allocate and zero bitmap pointers */
  kmerbitmap = (struct bitmap_s **) xmalloc(kmerhashsize * sizeof(struct bitmap_s *));
  std::memset(kmerbitmap, 0, kmerhashsize * sizeof(struct bitmap_s *));

  /* hash / bitmap setup */
  /* convert hash counts to position in index */
  kmerhash = (uint64_t *) xmalloc((kmerhashsize + 1) * sizeof(uint64_t));
  uint64_t sum = 0;
  for (unsigned int i = 0; i < kmerhashsize; i++)
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
  kmerindexsize = sum;
  kmerhash[kmerhashsize] = sum;

#if 0
  if (! opt_quiet)
    std::fprintf(stderr, "Unique %ld-mers: %u\n", opt_wordlength, kmerindexsize);
#endif

  /* reset counts */
  std::memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  /* allocate space for actual data */
  kmerindex = (unsigned int *) xmalloc(kmerindexsize * sizeof(unsigned int));

  /* allocate space for mapping from indexno to seqno */
  dbindex_map = (unsigned int *) xmalloc(seqcount * sizeof(unsigned int));

  dbindex_count = 0;

  show_rusage();
}


auto dbindex_free() -> void
{
  xfree(kmerhash);
  xfree(kmerindex);
  xfree(kmercount);
  xfree(dbindex_map);

  for (unsigned int kmer = 0; kmer < kmerhashsize; kmer++)
    {
      if (kmerbitmap[kmer] != nullptr)
        {
          bitmap_free(kmerbitmap[kmer]);
        }
    }
  xfree(kmerbitmap);
  unique_exit(dbindex_uh);
}
