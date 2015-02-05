/*
    Copyright (C) 2014-2015 Torbjorn Rognes

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

#define HASH hash_cityhash64

struct bucket
{
  unsigned long hash;
  unsigned int seqno_first;
  unsigned int seqno_last;
  unsigned int size;
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

  /* highest abundance first, then by label, otherwise keep order */

  if (x->size < y->size)
    return +1;
  else if (x->size > y->size)
    return -1;
  else
    {
      int r = strcmp(db_getheader(x->seqno_first),
                     db_getheader(y->seqno_first));
      if (r != 0)
        return r;
      else
        {
          if (x->seqno_first < y->seqno_first)
            return -1;
          else if (x->seqno_first > y->seqno_first)
            return +1;
          else
            return 0;
        }
    }
}

void string_normalize(char * normalized, char * s, unsigned int len)
{
  /* convert string to upper case and replace U by T */
  char * p = s;
  char * q = normalized;
  for(unsigned int i=0; i<len; i++)
    *q++ = chrmap_normalize[(int)(*p++)];
}

int seqcmp(char * a, char * b, int n)
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

  if (opt_uc)
    {
      fp_uc = fopen(opt_uc, "w");
      if (!fp_uc)
        fatal("Unable to open output (uc) file for writing");
    }

  db_read(opt_derep_fulllength, 0);

  show_rusage();

  long dbsequencecount = db_getsequencecount();
  
  /* adjust size of hash table for 2/3 fill rate */

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

  unsigned int * nextseqtab = (unsigned int*) xmalloc(sizeof(unsigned int) * dbsequencecount);
  memset(nextseqtab, 0, sizeof(unsigned int) * dbsequencecount);

  char * seq_up = (char*) xmalloc(db_getlongestsequence() + 1);
  char * rc_seq_up = (char*) xmalloc(db_getlongestsequence() + 1);
  
  progress_init("Dereplicating", dbsequencecount);
  for(long i=0; i<dbsequencecount; i++)
    {
      unsigned int seqlen = db_getsequencelen(i);
      char * seq = db_getsequence(i);

      /* normalize sequence: uppercase and replace U by T  */
      string_normalize(seq_up, seq, seqlen);

      /* reverse complement if necessary */
      if (opt_strand > 1)
        reverse_complement(rc_seq_up, seq_up, seqlen);

      /*
        Find free bucket or bucket for identical sequence.
        Make sure sequences are exactly identical
        in case of any hash collision.
        With 64-bit hashes, there is about 50% chance of a
        collision when the number of sequences is about 5e9.
      */

      unsigned long hash = HASH(seq_up, seqlen);
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
              (seqcmp(seq_up, db_getsequence(bp->seqno_first), seqlen))))
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

          unsigned long rc_hash = HASH(rc_seq_up, seqlen);
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
                  (seqcmp(rc_seq_up,
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

      long ab = opt_sizein ? db_getabundance(i) : 1;
      sumsize += ab;

      if (bp->size)
        {
          /* at least one identical sequence already */
          bp->size += ab;
          unsigned int last = bp->seqno_last;
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

  free(seq_up);
  free(rc_seq_up);
  
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
  
  long selected = 0;

  if (opt_output)
    progress_init("Writing output file", clusters);
    
  for (long i=0; i<clusters; i++)
    {
      struct bucket * bp = hashtable + i;
      long size = bp->size;
      if ((size >= opt_minuniquesize) && (size <= opt_maxuniquesize))
        {
          if (opt_output)
            {
              if (opt_sizeout)
                db_fprint_fasta_with_size(fp_output,
                                          bp->seqno_first, bp->size);
              else
                db_fprint_fasta(fp_output,
                                bp->seqno_first);
            }
          selected++;
          if (selected == opt_topn)
            break;
        }
      if (opt_output)
        progress_update(i);
    }

  if (opt_output)
    {
      progress_done();
      fclose(fp_output);
    }

  show_rusage();

  if (opt_uc)
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
          fprintf(fp_uc, "C\t%ld\t%d\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, bp->size, db_getheader(bp->seqno_first));
          progress_update(i);
        }
      fclose(fp_uc);
      progress_done();
      show_rusage();
    }

  if (selected < clusters)
    fprintf(stderr,
            "%ld uniques written, %ld clusters discarded (%.1f%%)\n",
            selected, clusters - selected,
            100.0 * (clusters - selected) / clusters);

  free(nextseqtab);
  free(hashtable);
  db_free();
}
