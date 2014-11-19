/*
    Copyright (C) 2014 Torbjorn Rognes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "vsearch.h"

unsigned int * kmercount;
unsigned int * kmerhash;
unsigned int * kmerindex;
bitmap_t * * kmerbitmap;
unsigned int * dbindex_map;

unsigned int kmerhashsize;
unsigned int kmerindexsize;
unsigned int dbindex_count;

uhandle_s * dbindex_uh;

#define BITMAP_THRESHOLD 8

unsigned int bitmap_mincount;

void fprint_kmer(FILE * f, unsigned int kk, unsigned long kmer)
{
  unsigned long x = kmer;
  for(unsigned int i=0; i<kk; i++)
    fprintf(f, "%c", sym_nt_2bit[(x >> (2*(kk-i-1))) & 3]);
}

void dbindex_addsequence(unsigned int seqno)
{
#if 0
  printf("Adding seqno %d as index element no %d\n", seqno, dbindex_count);
#endif

  unsigned int uniquecount;
  unsigned int * uniquelist;
  unique_count(dbindex_uh, opt_wordlength,
               db_getsequencelen(seqno), db_getsequence(seqno),
               & uniquecount, & uniquelist);
  dbindex_map[dbindex_count] = seqno;
  for(unsigned int i=0; i<uniquecount; i++)
    {
      unsigned int kmer = uniquelist[i];
      if (kmerbitmap[kmer])
        bitmap_set(kmerbitmap[kmer], dbindex_count);
      else
        kmerindex[kmerhash[kmer]+(kmercount[kmer]++)] = dbindex_count;
    }
  dbindex_count++;
}

void dbindex_addallsequences()
{
  unsigned int seqcount = db_getsequencecount();
  progress_init("Creating index of unique k-mers", seqcount);
  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      dbindex_addsequence(seqno);
      progress_update(seqno);
    }
  progress_done();
}

void dbindex_prepare(int use_bitmap)
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
                   & uniquecount, & uniquelist);
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
  kmerhash = (unsigned int *) xmalloc((kmerhashsize+1) * sizeof(unsigned int));
  unsigned int sum = 0;
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
  free(kmerhash);
  free(kmerindex);
  free(kmercount);
  free(dbindex_map);

  for(unsigned int kmer=0; kmer<kmerhashsize; kmer++)
    if (kmerbitmap[kmer])
      bitmap_free(kmerbitmap[kmer]);
  free(kmerbitmap);
  unique_exit(dbindex_uh);
}
