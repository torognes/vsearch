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

static int k;

unsigned int * kmerhash;
unsigned int * kmerindex;
unsigned int kmerhashsize;
unsigned int kmerindexsize;

inline int dbindex_getkmermatchcount(int kmer)
{
  return kmerhash[kmer+1] - kmerhash[kmer];
}

inline int dbindex_getkmermatch(int kmer, int matchno)
{
  return kmerindex[kmerhash[kmer]+matchno];
}

void fprint_kmer(FILE * f, unsigned int kk, unsigned long kmer)
{
  unsigned long x = kmer;
  for(unsigned int i=0; i<kk; i++)
    fprintf(f, "%c", sym_nt_2bit[(x >> (2*(kk-i-1))) & 3]);
}

void dbindex_build()
{
  k = wordlength;
  unsigned int seqcount = db_getsequencecount();
  kmerhashsize = 1 << (2*k);
  unsigned int * kmercount = (unsigned int *)
    xmalloc(kmerhashsize * sizeof(unsigned int));
  memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));
  struct uhandle_s * uh = unique_init();


  /* first scan, just count occurences */
  
  progress_init("Counting unique k-mers", seqcount);
  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      char * sequence = db_getsequence(seqno);
      long seqlen = db_getsequencelen(seqno);
      unsigned int uniquecount;
      unsigned int * uniquelist;
      unique_count(uh, k, seqlen, sequence, & uniquecount, & uniquelist);
      for(unsigned int i=0; i<uniquecount; i++)
	kmercount[uniquelist[i]]++;
      progress_update(seqno);
    }
  progress_done();


  /* hash setup */
  /* convert hash counts to position in index */

  kmerhash = (unsigned int *) xmalloc((kmerhashsize+1) * sizeof(unsigned int));
  unsigned int sum = 0;
  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      kmerhash[i] = sum;
      sum += kmercount[i];
    }
  kmerindexsize = sum;
  kmerhash[kmerhashsize] = sum;
  free(kmercount);
  fprintf(stderr, "Unique %u-mers: %u\n", k, kmerindexsize);
  show_rusage();


  /* second scan, fill in actual index of the unique kmers */

  kmerindex = (unsigned int *) xmalloc(kmerindexsize * sizeof(unsigned int));
  progress_init("Creating index of unique k-mers", seqcount);
  for(unsigned int seqno = 0; seqno < seqcount ; seqno++)
    {
      char * sequence = db_getsequence(seqno);
      long seqlen = db_getsequencelen(seqno);
      unsigned int uniquecount;
      unsigned int * uniquelist;
      unique_count(uh, k, seqlen, sequence, & uniquecount, & uniquelist);
      for(unsigned int i=0; i<uniquecount; i++)
	kmerindex[kmerhash[uniquelist[i]]++] = seqno;
      progress_update(seqno);
    }
  progress_done();


  /* reset kmerhash pointers (move up) */

  unsigned int temp = 0;
  unsigned int next;
  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      next = kmerhash[i];
      kmerhash[i] = temp;
      temp = next;
    }
  show_rusage();
  unique_exit(uh);
}

void dbindex_free()
{
  free(kmerhash);
  free(kmerindex);
}
