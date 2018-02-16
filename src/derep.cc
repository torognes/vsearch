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

//#define BITMAP

#define HASH hash_cityhash64

struct bucket
{
  uint64_t hash;
  unsigned int seqno_first;
  unsigned int seqno_last;
  unsigned int size;
  bool deleted;
};

#ifdef BITMAP
unsigned char * hash_occupied = 0;

void hash_set_occupied(uint64_t hashindex)
{
  hash_occupied[(hashindex) >> 3] |= 1 << (hashindex & 7);
}

int hash_is_occupied(uint64_t hashindex)
{
  return hash_occupied[(hashindex) >> 3] & (1 << (hashindex & 7));
}
#endif

int derep_compare(const void * a, const void * b)
{
  struct bucket * x = (struct bucket *) a;
  struct bucket * y = (struct bucket *) b;

  /* highest abundance first, then by label, otherwise keep order */

  if (x->deleted > y->deleted)
    return +1;
  else if (x->deleted < y->deleted)
    return -1;
  else
    {
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
      fp_output = fopen_output(opt_output);
      if (!fp_output)
        fatal("Unable to open output file for writing");
    }

  if (opt_uc)
    {
      fp_uc = fopen_output(opt_uc);
      if (!fp_uc)
        fatal("Unable to open output (uc) file for writing");
    }

  db_read(opt_derep_fulllength, 0);

  show_rusage();

  int64_t dbsequencecount = db_getsequencecount();
  
  /* adjust size of hash table for 2/3 fill rate */

  int64_t hashtablesize = 1;
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

  int64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;
  double median = 0.0;
  double average = 0.0;

  /* alloc and init table of links to other sequences in cluster */

  unsigned int * nextseqtab = (unsigned int*) xmalloc(sizeof(unsigned int) * dbsequencecount);
  const unsigned int terminal = (unsigned int)(-1);
  memset(nextseqtab, -1, sizeof(unsigned int) * dbsequencecount);

  char * match_strand = (char *) xmalloc(dbsequencecount);
  memset(match_strand, 0, dbsequencecount);

  char * seq_up = (char*) xmalloc(db_getlongestsequence() + 1);
  char * rc_seq_up = (char*) xmalloc(db_getlongestsequence() + 1);
  
  progress_init("Dereplicating", dbsequencecount);
  for(int64_t i=0; i<dbsequencecount; i++)
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

      uint64_t hash = HASH(seq_up, seqlen);
      uint64_t j = hash & hash_mask;
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

          uint64_t rc_hash = HASH(rc_seq_up, seqlen);
          struct bucket * rc_bp = hashtable + rc_hash % hashtablesize;
          uint64_t k = rc_hash & hash_mask;
          
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
              match_strand[i] = 1;
            }
        }

      int64_t ab = opt_sizein ? db_getabundance(i) : 1;
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

  xfree(seq_up);
  xfree(rc_seq_up);
  
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

  if (clusters < 1)
    {
      if (!opt_quiet)
        fprintf(stderr,
                "0 unique sequences\n");
      if (opt_log)
        fprintf(fp_log,
                "0 unique sequences\n\n");
    }
  else
    {
      if (!opt_quiet)
        fprintf(stderr,
                "%" PRId64 " unique sequences, avg cluster %.1lf, median %.0f, max %" PRIu64 "\n",
                clusters, average, median, maxsize);
      if (opt_log)
        fprintf(fp_log,
                "%" PRId64 " unique sequences, avg cluster %.1lf, median %.0f, max %" PRIu64 "\n\n",
                clusters, average, median, maxsize);
    }

  show_rusage();
  

  /* count selected */

  int64_t selected = 0;
  for (int64_t i=0; i<clusters; i++)
    {
      struct bucket * bp = hashtable + i;
      int64_t size = bp->size;
      if ((size >= opt_minuniquesize) && (size <= opt_maxuniquesize))
        {
          selected++;
          if (selected == opt_topn)
            break;
        }
    }


  /* write output */

  if (opt_output)
    {
      progress_init("Writing output file", clusters);

      int64_t relabel_count = 0;
      for (int64_t i=0; i<clusters; i++)
        {
          struct bucket * bp = hashtable + i;
          int64_t size = bp->size;
          if ((size >= opt_minuniquesize) && (size <= opt_maxuniquesize))
            {
              relabel_count++;
              fasta_print_general(fp_output,
                                  0,
                                  db_getsequence(bp->seqno_first),
                                  db_getsequencelen(bp->seqno_first),
                                  db_getheader(bp->seqno_first),
                                  db_getheaderlen(bp->seqno_first),
                                  size,
                                  relabel_count,
                                  -1, -1, 0, 0.0);
              if (relabel_count == opt_topn)
                break;
            }
          progress_update(i);
        }

      progress_done();
      fclose(fp_output);
    }

  show_rusage();

  if (opt_uc)
    {
      progress_init("Writing uc file, first part", clusters);
      for (int64_t i=0; i<clusters; i++)
        {
          struct bucket * bp = hashtable + i;
          char * h =  db_getheader(bp->seqno_first);
          int64_t len = db_getsequencelen(bp->seqno_first);

          fprintf(fp_uc, "S\t%" PRId64 "\t%" PRId64 "\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, len, h);
          
          for (unsigned int next = nextseqtab[bp->seqno_first];
               next != terminal;
               next = nextseqtab[next])
            fprintf(fp_uc,
                    "H\t%" PRId64 "\t%" PRId64 "\t%.1f\t%s\t0\t0\t*\t%s\t%s\n",
                    i, len, 100.0,
                    (match_strand[next] ? "-" : "+"),
                    db_getheader(next), h);

          progress_update(i);
        }
      progress_done();
      show_rusage();
      
      progress_init("Writing uc file, second part", clusters);
      for (int64_t i=0; i<clusters; i++)
        {
          struct bucket * bp = hashtable + i;
          fprintf(fp_uc, "C\t%" PRId64 "\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, bp->size, db_getheader(bp->seqno_first));
          progress_update(i);
        }
      fclose(fp_uc);
      progress_done();
      show_rusage();
    }

  if (selected < clusters)
    {
      if (!opt_quiet)
        fprintf(stderr,
                "%" PRId64 " uniques written, %" PRId64 " clusters discarded (%.1f%%)\n",
                selected, clusters - selected,
                100.0 * (clusters - selected) / clusters);

      if (opt_log)
        fprintf(fp_log,
                "%" PRId64 " uniques written, %" PRId64 " clusters discarded (%.1f%%)\n\n",
                selected, clusters - selected,
                100.0 * (clusters - selected) / clusters);
    }
  
  xfree(match_strand);
  xfree(nextseqtab);
  xfree(hashtable);
  db_free();
}


void derep_prefix()
{
  FILE * fp_output = 0;
  FILE * fp_uc = 0;

  if (opt_output)
    {
      fp_output = fopen_output(opt_output);
      if (!fp_output)
        fatal("Unable to open output file for writing");
    }

  if (opt_uc)
    {
      fp_uc = fopen_output(opt_uc);
      if (!fp_uc)
        fatal("Unable to open output (uc) file for writing");
    }

  db_read(opt_derep_prefix, 0);
  
  db_sortbylength_shortest_first();

  show_rusage();

  int64_t dbsequencecount = db_getsequencecount();
  
  /* adjust size of hash table for 2/3 fill rate */

  int64_t hashtablesize = 1;
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

  int64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;
  double median = 0.0;
  double average = 0.0;

  /* alloc and init table of links to other sequences in cluster */

  unsigned int * nextseqtab = (unsigned int*) xmalloc(sizeof(unsigned int) * dbsequencecount);
  const unsigned int terminal = (unsigned int)(-1);
  memset(nextseqtab, -1, sizeof(unsigned int) * dbsequencecount);

  char * seq_up = (char*) xmalloc(db_getlongestsequence() + 1);

  /* make table of hash values of prefixes */

  unsigned int len_longest = db_getlongestsequence();
  unsigned int len_shortest = db_getshortestsequence();
  uint64_t * prefix_hashes = (uint64_t *) 
    xmalloc(sizeof(uint64_t) * (len_longest+1));
  
  progress_init("Dereplicating", dbsequencecount);
  for(int64_t i=0; i<dbsequencecount; i++)
    {
      unsigned int seqlen = db_getsequencelen(i);
      char * seq = db_getsequence(i);

      /* normalize sequence: uppercase and replace U by T  */
      string_normalize(seq_up, seq, seqlen);

      uint64_t ab = opt_sizein ? db_getabundance(i) : 1;
      sumsize += ab;

      /* 
         Look for matching identical or prefix sequences.
         
         Use a hash function that can quickly be applied iteratively on longer
         and longer sequences.
         
         Hash values are generated for all prefixes and saved.

         Should start at exact sequence and then try shorter and shorter
         sequences.
         
         No need to check shorter sequences than the shortest in the database.
         
         Three cases:
         1) Exact match: Update count, point to next
         2) Prefix match: Mark old, insert new, update count, point to next
         3) No match: Insert new entry

      */

      /* compute hashes of all prefixes */

      uint64_t fnv1a_hash = 14695981039346656037ULL;
      prefix_hashes[0] = fnv1a_hash;
      for(unsigned int j = 0; j < seqlen; j++)
        {
          fnv1a_hash ^= seq_up[j];
          fnv1a_hash *= 1099511628211ULL;
          prefix_hashes[j+1] = fnv1a_hash;
        }

      /* first, look for an identical match */

      unsigned int prefix_len = seqlen;

      uint64_t hash = prefix_hashes[prefix_len];
      struct bucket * bp = hashtable + (hash & hash_mask);
      
      while ((bp->size) &&
             ((bp->deleted) ||
              (bp->hash != hash) ||
              (prefix_len != db_getsequencelen(bp->seqno_first)) ||
              (seqcmp(seq_up, db_getsequence(bp->seqno_first), prefix_len))))
        {
          bp++;
          if (bp >= hashtable + hashtablesize)
            bp = hashtable;
        }

      /* at this point, bp points either to (1) a free empty hash bucket, or
         (2) a bucket with an exact match. */

      uint64_t orig_hash = hash;
      struct bucket * orig_bp = bp;

      if (bp->size)
        {
          /* exact match */
          bp->size += ab;
          unsigned int last = bp->seqno_last;
          nextseqtab[last] = i;
          bp->seqno_last = i;
          
          if (bp->size > maxsize)
            maxsize = bp->size;
        }
      else
        {
          /* look for prefix match */
          
          while((! bp->size) && (prefix_len > len_shortest))
            {
              prefix_len--;
              hash = prefix_hashes[prefix_len];
              bp = hashtable + (hash & hash_mask);
              
              while ((bp->size) &&
                     ((bp->deleted) ||
                      (bp->hash != hash) ||
                      (prefix_len != db_getsequencelen(bp->seqno_first)) ||
                      (seqcmp(seq_up,
                              db_getsequence(bp->seqno_first),
                              prefix_len))))
                {
                  bp++;
                  if (bp >= hashtable + hashtablesize)
                    bp = hashtable;
                }
            }
          
          if (bp->size)
            {
              /* prefix match */

              /* get necessary info, then delete prefix from hash */
              unsigned int first = bp->seqno_first;
              unsigned int last = bp->seqno_last;
              unsigned int size = bp->size;
              bp->deleted = true;
              
              /* create new hash entry */
              bp = orig_bp;
              bp->size = size + ab;
              bp->hash = orig_hash;
              bp->seqno_first = i;
              nextseqtab[i] = first;
              bp->seqno_last = last;
              
              if (bp->size > maxsize)
                maxsize = bp->size;
            }
          else
            {
              /* no match */
              orig_bp->size = ab;
              orig_bp->hash = orig_hash;
              orig_bp->seqno_first = i;
              orig_bp->seqno_last = i;
              
              if (ab > maxsize)
                maxsize = ab;
              clusters++;
            }
        }
      
      progress_update(i);
    }
  progress_done();
  
  xfree(prefix_hashes);

  xfree(seq_up);
  
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

  if (clusters < 1)
    {
      if (!opt_quiet)
        fprintf(stderr,
                "0 unique sequences\n");
      if (opt_log)
        fprintf(fp_log,
                "0 unique sequences\n\n");
    }
  else
    {
      if (!opt_quiet)
        fprintf(stderr,
                "%" PRId64 " unique sequences, avg cluster %.1lf, median %.0f, max %" PRIu64 "\n",
                clusters, average, median, maxsize);
      if (opt_log)
        fprintf(fp_log,
                "%" PRId64 " unique sequences, avg cluster %.1lf, median %.0f, max %" PRIu64 "\n\n",
                clusters, average, median, maxsize);
    }

  show_rusage();
  
  /* count selected */

  int64_t selected = 0;
  for (int64_t i=0; i<clusters; i++)
    {
      struct bucket * bp = hashtable + i;
      int64_t size = bp->size;
      if ((size >= opt_minuniquesize) && (size <= opt_maxuniquesize))
        {
          selected++;
          if (selected == opt_topn)
            break;
        }
    }


  /* write output */

  if (opt_output)
    {
      progress_init("Writing output file", clusters);

      int64_t relabel_count = 0;
      for (int64_t i=0; i<clusters; i++)
        {
          struct bucket * bp = hashtable + i;
          int64_t size = bp->size;
          if ((size >= opt_minuniquesize) && (size <= opt_maxuniquesize))
            {
              relabel_count++;
              fasta_print_general(fp_output,
                                  0,
                                  db_getsequence(bp->seqno_first),
                                  db_getsequencelen(bp->seqno_first),
                                  db_getheader(bp->seqno_first),
                                  db_getheaderlen(bp->seqno_first),
                                  size,
                                  relabel_count,
                                  -1, -1, 0, 0.0);
              if (relabel_count == opt_topn)
                break;
            }
          progress_update(i);
        }

      progress_done();
      fclose(fp_output);
    }

  show_rusage();

  if (opt_uc)
    {
      progress_init("Writing uc file, first part", clusters);
      for (int64_t i=0; i<clusters; i++)
        {
          struct bucket * bp = hashtable + i;
          char * h =  db_getheader(bp->seqno_first);
          int64_t len = db_getsequencelen(bp->seqno_first);

          fprintf(fp_uc, "S\t%" PRId64 "\t%" PRId64 "\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, len, h);
          
          for (unsigned int next = nextseqtab[bp->seqno_first];
               next != terminal;
               next = nextseqtab[next])
            fprintf(fp_uc,
                    "H\t%" PRId64 "\t%" PRIu64 "\t%.1f\t+\t0\t0\t*\t%s\t%s\n",
                    i, db_getsequencelen(next), 100.0, db_getheader(next), h);

          progress_update(i);
        }
      progress_done();
      show_rusage();
      
      progress_init("Writing uc file, second part", clusters);
      for (int64_t i=0; i<clusters; i++)
        {
          struct bucket * bp = hashtable + i;
          fprintf(fp_uc, "C\t%" PRId64 "\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, bp->size, db_getheader(bp->seqno_first));
          progress_update(i);
        }
      fclose(fp_uc);
      progress_done();
      show_rusage();
    }

  if (selected < clusters)
    {
      if (!opt_quiet)
        fprintf(stderr,
                "%" PRId64 " uniques written, %" PRId64 " clusters discarded (%.1f%%)\n",
                selected, clusters - selected,
                100.0 * (clusters - selected) / clusters);

      if (opt_log)
        fprintf(fp_log,
                "%" PRId64 " uniques written, %" PRId64 " clusters discarded (%.1f%%)\n\n",
                selected, clusters - selected,
                100.0 * (clusters - selected) / clusters);
    }
  
  xfree(nextseqtab);
  xfree(hashtable);
  db_free();
}
