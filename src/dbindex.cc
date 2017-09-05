/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

unsigned int * kmercount;
uint64_t * kmerhash;
unsigned int * kmerindex;
bitmap_t * * kmerbitmap;
unsigned int * dbindex_map;

static unsigned int kmerhashsize;
static uint64_t kmerindexsize;
unsigned int dbindex_count;

static uhandle_s * dbindex_uh;

#define BITMAP_THRESHOLD 8

static unsigned int bitmap_mincount;

void fprint_kmer(FILE * f, unsigned int kk, uint64_t kmer)
{
  uint64_t x = kmer;
  for(unsigned int i=0; i<kk; i++)
    fprintf(f, "%c", sym_nt_2bit[(x >> (2*(kk-i-1))) & 3]);
}

void dbindex_addsequence(unsigned int seqno, int seqmask)
{
#if 0
  printf("Adding seqno %d as index element no %d\n", seqno, dbindex_count);
#endif

  unsigned int uniquecount;
  unsigned int * uniquelist;
  unique_count(dbindex_uh, opt_wordlength,
               db_getsequencelen(seqno), db_getsequence(seqno),
               & uniquecount, & uniquelist, seqmask);
  dbindex_map[dbindex_count] = seqno;
  for(unsigned int i=0; i<uniquecount; i++)
    {
      unsigned int kmer = uniquelist[i];
      if (kmerbitmap[kmer])
        {
          kmercount[kmer]++;
          bitmap_set(kmerbitmap[kmer], dbindex_count);
        }
      else
        kmerindex[kmerhash[kmer]+(kmercount[kmer]++)] = dbindex_count;
    }
  dbindex_count++;
}

void dbindex_addallsequences(int seqmask)
{
  unsigned int seqcount = db_getsequencecount();
  progress_init("Creating index of unique k-mers", seqcount);
  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      dbindex_addsequence(seqno, seqmask);
      progress_update(seqno);
    }
  progress_done();
}

void dbindex_prepare(int use_bitmap, int seqmask)
{
  dbindex_uh = unique_init();

  unsigned int seqcount = db_getsequencecount();
  kmerhashsize = 1 << (2 * opt_wordlength);

  /* allocate memory for kmer count array */
  kmercount = (unsigned int *) xmalloc(kmerhashsize * sizeof(unsigned int));
  memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  /* first scan, just count occurences */
  progress_init("Counting unique k-mers", seqcount);
  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      unsigned int uniquecount;
      unsigned int * uniquelist;
      unique_count(dbindex_uh, opt_wordlength,
                   db_getsequencelen(seqno), db_getsequence(seqno),
                   & uniquecount, & uniquelist, seqmask);
      for(unsigned int i=0; i<uniquecount; i++)
        kmercount[uniquelist[i]]++;
      progress_update(seqno);
    }
  progress_done();

#if 0
  /* dump kmer counts */
  FILE * f = fopen("kmercounts.txt", "w");
  for(unsigned int kmer=0; kmer < kmerhashsize; kmer++)
    {
      fprint_kmer(f, 8, kmer);
      fprintf(f, "\t%d\t%d\n", kmer, kmercount[kmer]);
    }
  fclose(f);
#endif

  /* determine minimum kmer count for bitmap usage */
  if (use_bitmap)
    bitmap_mincount = seqcount / BITMAP_THRESHOLD;
  else
    bitmap_mincount = seqcount + 1;

  /* allocate and zero bitmap pointers */
  kmerbitmap = (bitmap_t **) xmalloc(kmerhashsize * sizeof(bitmap_t *));
  memset(kmerbitmap, 0, kmerhashsize * sizeof(bitmap_t *));

  /* hash / bitmap setup */
  /* convert hash counts to position in index */
  kmerhash = (uint64_t *) xmalloc((kmerhashsize+1) * sizeof(uint64_t));
  uint64_t sum = 0;
  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      kmerhash[i] = sum;
      if (kmercount[i] >= bitmap_mincount)
        {
          kmerbitmap[i] = bitmap_init(seqcount+127); // pad for xmm
          bitmap_reset_all(kmerbitmap[i]);
        }
      else
        sum += kmercount[i];
    }
  kmerindexsize = sum;
  kmerhash[kmerhashsize] = sum;

#if 0
  if (!opt_quiet)
    fprintf(stderr, "Unique %ld-mers: %u\n", opt_wordlength, kmerindexsize);
#endif
  
  /* reset counts */
  memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  /* allocate space for actual data */
  kmerindex = (unsigned int *) xmalloc(kmerindexsize * sizeof(unsigned int));

  /* allocate space for mapping from indexno to seqno */
  dbindex_map = (unsigned int *) xmalloc(seqcount * sizeof(unsigned int));
  
  dbindex_count = 0;
  
  show_rusage();
}

void dbindex_free()
{
  xfree(kmerhash);
  xfree(kmerindex);
  xfree(kmercount);
  xfree(dbindex_map);

  for(unsigned int kmer=0; kmer<kmerhashsize; kmer++)
    if (kmerbitmap[kmer])
      bitmap_free(kmerbitmap[kmer]);
  xfree(kmerbitmap);
  unique_exit(dbindex_uh);
}

void dbindex_udb_write()
{
  int fd_output = 0;

  fd_output = open(opt_output, O_WRONLY | O_CREAT | O_TRUNC);
  if (!fd_output)
    fatal("Unable to open output file for writing");

  db_read(opt_makeudb_usearch, 1);

  if (opt_dbmask == MASK_DUST)
    dust_all();
  else if ((opt_dbmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();

  dbindex_prepare(1, opt_dbmask);
  dbindex_addallsequences(opt_dbmask);

  unsigned int seqcount = db_getsequencecount();
  uint64_t ntcount = db_getnucleotidecount();

  uint64_t header_characters = 0;
  for (unsigned int i=0; i<seqcount; i++)
    header_characters += db_getheaderlen(i) + 1;

  uint64_t kmerhashsize = 1 << (2 * opt_wordlength);

  progress_init("Writing UDB file", kmerhashsize);

  uint64_t buffersize = 4 * MAX(50, seqcount);
  unsigned int * buffer = (unsigned int *) xmalloc(buffersize);
  memset(buffer, 0, buffersize);

  /* Header */
  buffer[0]  = 0x55444246; /* FBDU UDBF */
  buffer[2]  = 32; /* bits */
  buffer[4]  = opt_wordlength; /* default 8 */
  buffer[5]  = 1; /* dbstep */
  buffer[6]  = 100; /* dbaccelpct % */
  buffer[11] = 0; /* slots */
  buffer[13] = (unsigned int) seqcount; /* antall sekvenser */
  buffer[17] = 0x0000746e; /* alphabet: "nt" */
  buffer[49] = 0x55444266; /* fBDU UDBf */
  ssize_t written = write(fd_output, buffer, 50 * 4);
  if (written < 50 * 4)
    fatal("Unable to write to UDB file");

  /* 65536 uint32's with word match counts */
  written = write(fd_output, kmercount, 4 * kmerhashsize);
  if (written < (ssize_t)(4 * kmerhashsize))
    fatal("Unable to write to UDB file");

  /* 3BDU */
  buffer[0] = 0x55444233; /* 3BDU UDB3 */
  written = write(fd_output, buffer, 1 * 4);
  if (written < 1 * 4)
    fatal("Unable to write to UDB file");

  /* lists of sequence no's with matches for all words */
  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      if (kmerbitmap[i])
        {
          memset(buffer, 0, 4 * kmercount[i]);
          unsigned int elements = 0;
          for (unsigned int j = 0; j < seqcount; j++)
            if (bitmap_get(kmerbitmap[i], j))
              buffer[elements++] = j;
          written = write(fd_output, buffer, 4 * elements);
          if (written < 4 * elements)
            fatal("Unable to write to UDB file");
        }
      else
        {
          if (kmercount[i] > 0)
            {
              written = write(fd_output, kmerindex+kmerhash[i], 4*kmercount[i]);
              if (written < 4 * kmercount[i])
                fatal("Unable to write to UDB file");
            }
        }
      progress_update(i);
    }

  /* New header */
  buffer[0] = 0x55444234; /* 4BDU UDB4 */
  /* 0x005e0db3 */
  buffer[1] = 0x005e0db3;
  /* number of sequences, uint32 */
  buffer[2] = (unsigned int) seqcount;
  /* total number of nucleotides, uint64 */
  buffer[3] = (unsigned int)(ntcount & 0xffffffff);
  buffer[4] = (unsigned int)(ntcount >> 32);
  /* total number of header characters, incl zero-terminator, uint64 */
  buffer[5] = (unsigned int)(header_characters & 0xffffffff);
  buffer[6] = (unsigned int)(header_characters >> 32);
  /* 0x005e0db4 */
  buffer[7] = 0x005e0db4;
  written = write(fd_output, buffer, 4 * 8);
  if (written < 4 * 8)
    fatal("Unable to write to UDB file");

  /* indices to headers (uint32) */
  unsigned int sum = 0;
  for (unsigned int i = 0; i < seqcount; i++)
    {
      buffer[i] = sum;
      sum += db_getheaderlen(i) + 1;
    }
  written = write(fd_output, buffer, 4 * seqcount);
  if (written < 4 * seqcount)
    fatal("Unable to write to UDB file");

  /* headers (ascii, zero terminated, not padded) */
  for (unsigned int i = 0; i < seqcount; i++)
    {
      unsigned int len = db_getheaderlen(i);
      written = write(fd_output, db_getheader(i), len+1);
      if (written < len+1)
        fatal("Unable to write to UDB file");
    }

  /* sequence lengths (uint32) */
  for (unsigned int i = 0; i < seqcount; i++)
    buffer[i] = db_getsequencelen(i);
  written = write(fd_output, buffer, 4 * seqcount);
  if (written < 4 * seqcount)
    fatal("Unable to write to UDB file");

  /* sequences (ascii, no term, no pad) */
  for (unsigned int i = 0; i < seqcount; i++)
    {
      unsigned int len = db_getsequencelen(i);
      written = write(fd_output, db_getsequence(i), len);
      if (written < len)
        fatal("Unable to write to UDB file");
    }

  if (close(fd_output) != 0)
    fatal("Unable to close UDB file");

  progress_done();
  dbindex_free();
  db_free();
  xfree(buffer);
}
