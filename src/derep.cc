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

//#define BITMAP

//#define HASH hash_fnv_1a_64_uc
#define HASH hash_cityhash64

struct bucket
{
  unsigned long seqno_first;
  unsigned long seqno_last;
  unsigned long hash;
  unsigned long size;
};

#ifdef BITMAP
unsigned char * hash_occupied = 0;

void hash_set_occupied(unsigned long hashindex)
{
  hash_occupied[(hashindex) >> 3] |= 1 << (hashindex & 7);
}

int hash_is_occupied(unsigned long hashindex)
{
  return hash_occupied[(hashindex) >> 3] & (1 << (hashindex & 7));
}
#endif

int derep_compare(const void * a, const void * b)
{
  struct bucket * x = (struct bucket *) a;
  struct bucket * y = (struct bucket *) b;

  /* highest abundance first, otherwize keep order */

  if (x->size < y->size)
    return +1;
  else if (x->size > y->size)
    return -1;
  else
    if (x->seqno_first < y->seqno_first)
      return -1;
    else if (x->seqno_first > y->seqno_first)
      return +1;
    else
      return 0;
}

void string_locase(char * locase, char * s, unsigned long len)
{
  /* convert string to lower case, given all chars are in A-Za-z */
  char * p = s;
  char * q = locase;
  for(unsigned long i=0; i<len; i++)
    *q++ = *p++ | 0x20;
}

void derep_fulllength()
{
  FILE * fp_output = 0;
  FILE * fp_uc = 0;

  if (opt_output)
    {
      fp_output = fopen(opt_output, "w");
      if (!fp_output)
	fatal("Unable to open output file for writing");
    }

  if (ucfilename)
    {
      fp_uc = fopen(ucfilename, "w");
      if (!fp_uc)
	fatal("Unable to open output (uc) file for writing");
    }

  db_read(opt_derep_fulllength);

  show_rusage();


  long dbsequencecount = db_getsequencecount();
  
  long hashtablesize = 1;
  int hash_shift = 0;
  while (3 * dbsequencecount > 2 * hashtablesize)
    {
      hashtablesize <<= 1;
      hash_shift++;
    }
  int hash_mask = hashtablesize - 1;

  struct bucket * hashtable =
    (struct bucket *) xmalloc(sizeof(bucket) * hashtablesize);

  memset(hashtable, 0, sizeof(bucket) * hashtablesize);

#ifdef BITMAP
  hash_occupied =
    (unsigned char *) xmalloc(hashtablesize / 8);
  memset(hash_occupied, 0, hashtablesize / 8);
#endif

  long clusters = 0;
  long sumsize = 0;
  unsigned long maxsize = 0;
  double median = 0.0;
  double average = 0.0;

  /* alloc and init table of links to other sequences in cluster */

  long * nextseqtab = (long*) xmalloc(sizeof(long) * dbsequencecount);
  memset(nextseqtab, 0, sizeof(long) * dbsequencecount);

  char * seq_lo = (char*) xmalloc(db_getlongestsequence() + 1);
  char * rc_seq_lo = (char*) xmalloc(db_getlongestsequence() + 1);
  
  progress_init("Dereplicating", dbsequencecount);
  for(long i=0; i<dbsequencecount; i++)
    {
      unsigned long seqlen = db_getsequencelen(i);
      char * seq = db_getsequence(i);

      /* lower case sequence */
      string_locase(seq_lo, seq, seqlen);

      /* reverse complement if necessary */
      if (opt_strand > 1)
	reverse_complement(rc_seq_lo, seq_lo, seqlen);

      /*
	Find free bucket or bucket for identical sequence.
	Make sure sequences are exactly identical
	in case of any hash collision.
	With 64-bit hashes, there is about 50% chance of a
	collision when the number of sequences is about 5e9.
      */

      unsigned long hash = HASH(seq_lo, seqlen);
      unsigned long j = hash & hash_mask;
      struct bucket * bp = hashtable + j;
      
      while (
#ifdef BITMAP
	     (hash_is_occupied(j))
#else
	     (bp->size)
#endif
	     &&
	     ((bp->hash != hash) ||
	      (seqlen != db_getsequencelen(bp->seqno_first)) ||
	      (strncasecmp(seq_lo, db_getsequence(bp->seqno_first), seqlen))))
	{
	  bp++;
	  j++;
	  if (bp >= hashtable + hashtablesize)
	    {
	      bp = hashtable;
	      j = 0;
	    }
	}
	  
      if ((opt_strand > 1) && !bp->size)
	{
	  /* no match on plus strand */
	  /* check minus strand as well */

	  unsigned long rc_hash = HASH(rc_seq_lo, seqlen);
	  struct bucket * rc_bp = hashtable + rc_hash % hashtablesize;
	  unsigned long k = rc_hash & hash_mask;
	  
	  while (
#ifdef BITMAP
		 (hash_is_occupied(j))
#else
		 (rc_bp->size)
#endif
		 &&
		 ((rc_bp->hash != rc_hash) ||
		  (seqlen != db_getsequencelen(rc_bp->seqno_first)) ||
		  (strncasecmp(rc_seq_lo,
			       db_getsequence(rc_bp->seqno_first), 
			       seqlen))))
	    {
	      rc_bp++;
	      k++;
	      if (rc_bp >= hashtable + hashtablesize)
		{
		  rc_bp = hashtable;
		  k++;
		}
	    }

	  if (rc_bp->size)
	    {
	      bp = rc_bp;
	      j = k;
	    }
	}

      long ab = db_getabundance(i); 
      sumsize += ab;

      if (bp->size)
	{
	  /* at least one identical sequence already */
	  bp->size += ab;
	  unsigned long last = bp->seqno_last;
	  nextseqtab[last] = i;
	  bp->seqno_last = i;
	}
      else
	{
	  /* no identical sequences yet */
	  bp->size = ab;
	  bp->hash = hash;
	  bp->seqno_first = i;
	  bp->seqno_last = i;
	  clusters++;
	}

      if (bp->size > maxsize)
	maxsize = bp->size;

#ifdef BITMAP
      hash_set_occupied(j);
#endif

      progress_update(i);
    }
  progress_done();

  free(seq_lo);
  free(rc_seq_lo);
  
  show_rusage();


  progress_init("Sorting", 1);
  qsort(hashtable, hashtablesize, sizeof(bucket), derep_compare);
  progress_done();


  if (clusters > 0)
    {
      if (clusters % 2)
	median = hashtable[(clusters-1)/2].size;
      else
	median = (hashtable[(clusters/2)-1].size +
		  hashtable[clusters/2].size) / 2.0;
    }
  
  average = 1.0 * sumsize / clusters;

  fprintf(stderr,
	  "%ld unique sequences, avg cluster %.1lf, median %.0f, max %ld\n",
	  clusters, average, median, maxsize);

  show_rusage();


  long bigclusters = MIN(clusters, opt_topn);

  long i = 0;
  while ((i<bigclusters) && ((long)(hashtable[i].size) >= opt_minuniquesize))
    i++;
  bigclusters = i;

  if (opt_output)
    {
      progress_init("Writing output file", bigclusters);
      for (long i=0; i<bigclusters; i++)
	{
	  struct bucket * bp = hashtable + i;
	  if (opt_sizeout)
	    db_fprint_fasta_with_size(fp_output, bp->seqno_first, bp->size);
	  else
	    db_fprint_fasta(fp_output, bp->seqno_first);
	  progress_update(i);
	}
      progress_done();
      fclose(fp_output);
      show_rusage();
    }
  
  if (ucfilename)
    {
      progress_init("Writing uc file, first part", clusters);
      for (long i=0; i<clusters; i++)
	{
	  struct bucket * bp = hashtable + i;
	  char * h =  db_getheader(bp->seqno_first);
	  long len = db_getsequencelen(bp->seqno_first);

	  fprintf(fp_uc, "S\t%ld\t%ld\t*\t*\t*\t*\t*\t%s\t*\n",
		  i, len, h);
	  
	  for (unsigned long next = nextseqtab[bp->seqno_first];
	       next;
	       next = nextseqtab[next])
	    fprintf(fp_uc,
		    "H\t%ld\t%ld\t%.1f\t*\t0\t0\t*\t%s\t%s\n",
		    i, len, 100.0, db_getheader(next), h);

	  progress_update(i);
	}
      progress_done();
      show_rusage();
      
      progress_init("Writing uc file, second part", clusters);
      for (long i=0; i<clusters; i++)
	{
	  struct bucket * bp = hashtable + i;
	  fprintf(fp_uc, "C\t%ld\t%ld\t*\t*\t*\t*\t*\t%s\t*\n",
		  i, bp->size, db_getheader(bp->seqno_first));
	  progress_update(i);
	}
      fclose(fp_uc);
      progress_done();
      show_rusage();
    }

  if (bigclusters < clusters)
    fprintf(stderr,
	    "%ld uniques written, %ld clusters < %ld discarded (%.1f%%)\n",
	    bigclusters, clusters - bigclusters, opt_minuniquesize,
	    100.0 * (clusters - bigclusters) / clusters);

  free(nextseqtab);
  free(hashtable);
  db_free();
}
