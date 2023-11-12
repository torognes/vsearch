/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2022, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#define HASH hash_cityhash128

struct sm_bucket
{
  uint128 hash;
  uint64_t size;
};

static struct sm_bucket * hashtable = nullptr;
static uint64_t hashtablesize = 0;

double find_median()
{
  /* find the median size, based on an iterative search starting at e.g. 1 */

  uint64_t cand = 1;    /* candidate for the median */
  uint64_t below = 0;   /* closest value below the candidate */
  uint64_t above = 0;   /* closest value above the candidate */

  uint64_t cand_count;  /* number of clusters with same size as cand */
  uint64_t below_count; /* number of clusters with smaller size than cand */
  uint64_t above_count; /* number of clusters with larger size than cand */

  while (true)
    {
      cand_count = 0;
      below_count = 0;
      above_count = 0;

      for(uint64_t i = 0; i < hashtablesize; i++)
        {
          uint64_t v = hashtable[i].size;
          if (v > 0)
            {
              if (v > cand)
                {
                  if ((above_count == 0) || (v < above))
                    {
                      above = v;
                    }
                  above_count++;
                }
              else if (v < cand)
                {
                  if ((below_count == 0) || (v > below))
                    {
                      below = v;
                    }
                  below_count++;
                }
              else
                {
                  cand_count++;
                }
            }
        }

      if (below_count + cand_count + above_count == 0U) // fix -Wfloat-equal
        return 0;

      if (above_count + cand_count >= below_count)
        // mid >= below_count
        {
          if (above_count <= below_count + cand_count)
            // mid <= below_count + cand_count
            {
              if (above_count == below_count + cand_count)
                // mid == below_count + cand_count
                // same as:
                // (below_count + cand_count + above_count) / 2 == below_count + cand_count
                // which simplifies into:
                // above_count == below_count + cand_count
                {
                  return (cand + above) / 2.0;
                }
              else if (above_count + cand_count == below_count)
                // mid == below_count
                // same as:
                // (below_count + cand_count + above_count) / 2 == below_count
                // which simplifies into:
                // above_count + cand_count == below_count
                {
                  return (below + cand) / 2.0;
                }
              else
                {
                  return cand;
                }
            }
          else
            {
              cand = above;
            }
        }
      else
        {
          cand = below;
        }
    }
}

uint64_t inline hash2bucket(uint128 hash, uint64_t htsize)
{
  return Uint128Low64(hash) % htsize;
}

uint64_t inline next_bucket(uint64_t prev_bucket, uint64_t htsize)
{
  return (prev_bucket + 1) % htsize;
}

void rehash_smallmem()
{
  /* allocate new hash table, 50% larger */
  uint64_t new_hashtablesize = 3 * hashtablesize / 2;
  auto * new_hashtable =
    (struct sm_bucket *) xmalloc(sizeof(struct sm_bucket) * new_hashtablesize);

  /* zero new hash table */
  for(uint64_t j = 0; j < new_hashtablesize; j++)
    {
      new_hashtable[j].hash.first = 0;
      new_hashtable[j].hash.second = 0;
      new_hashtable[j].size = 0;
    }

  /* rehash all from old to new */
  for(uint64_t i = 0; i < hashtablesize; i++)
    {
      struct sm_bucket * old_bp = hashtable + i;
      if (old_bp->size)
        {
          uint64_t k = hash2bucket(old_bp->hash, new_hashtablesize);
          while (new_hashtable[k].size)
            {
              k = next_bucket(k, new_hashtablesize);
            }
          struct sm_bucket * new_bp = new_hashtable + k;
          * new_bp = * old_bp;
        }
    }

  /* free old table */
  xfree(hashtable);

  /* update variables */
  hashtable = new_hashtable;
  hashtablesize = new_hashtablesize;
}

void derep_smallmem(char * input_filename)
{
  /*
    dereplicate full length sequences using a small amount of memory
    output options: --fastaout
  */

  show_rusage();

  fastx_handle h = fastx_open(input_filename);

  if (!h)
    {
      fatal("Unrecognized input file type (not proper FASTA or FASTQ format).");
    }

  if (h->is_pipe)
    {
      fatal("The derep_smallmem command does not support input from a pipe.");
    }

  FILE * fp_fastaout = nullptr;

  if (opt_fastaout)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (!fp_fastaout)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }
  else
    {
      fatal("Ouput file for dereplication must be specified with --fastaout");
    }

  uint64_t filesize = fastx_get_size(h);

  /* allocate initial memory for sequences of length up to 1023 chars */
  int64_t alloc_seqlen = 1024;

  /* allocate initial hashtable with 1024 buckets */

  hashtablesize = 1024;
  hashtable = (struct sm_bucket *) xmalloc(sizeof(struct sm_bucket) * hashtablesize);

  /* zero hash table */
  for(uint64_t j = 0; j < hashtablesize; j++)
    {
      hashtable[j].hash.first = 0;
      hashtable[j].hash.second = 0;
      hashtable[j].size = 0;
    }

  show_rusage();

  char * seq_up = (char*) xmalloc(alloc_seqlen + 1);
  char * rc_seq_up = (char*) xmalloc(alloc_seqlen + 1);

  char * prompt = nullptr;
  if (xsprintf(& prompt, "Dereplicating file %s", input_filename) == -1)
    {
      fatal("Out of memory");
    }

  progress_init(prompt, filesize);

  uint64_t sequencecount = 0;
  uint64_t nucleotidecount = 0;
  int64_t shortest = INT64_MAX;
  int64_t longest = 0;
  uint64_t discarded_short = 0;
  uint64_t discarded_long = 0;
  uint64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;
  double median = 0.0;
  double average = 0.0;

  /* first pass */

  while(fastx_next(h, ! opt_notrunclabels, chrmap_no_change))
    {
      int64_t seqlen = fastx_get_sequence_length(h);

      if (seqlen < opt_minseqlength)
        {
          discarded_short++;
          continue;
        }

      if (seqlen > opt_maxseqlength)
        {
          discarded_long++;
          continue;
        }

      nucleotidecount += seqlen;
      if (seqlen > longest)
        {
          longest = seqlen;
        }
      if (seqlen < shortest)
        {
          shortest = seqlen;
        }

      /* check allocations */

      if (seqlen > alloc_seqlen)
        {
          alloc_seqlen = seqlen;
          seq_up = (char*) xrealloc(seq_up, alloc_seqlen + 1);
          rc_seq_up = (char*) xrealloc(rc_seq_up, alloc_seqlen + 1);

          show_rusage();
        }

      if (100 * (clusters + 1) > 95 * hashtablesize)
        {
          // keep hash table fill rate at max 95% */
          rehash_smallmem();
          show_rusage();
        }

      char * seq = fastx_get_sequence(h);

      /* normalize sequence: uppercase and replace U by T  */
      string_normalize(seq_up, seq, seqlen);

      /* reverse complement if necessary */
      if (opt_strand > 1)
        {
          reverse_complement(rc_seq_up, seq_up, seqlen);
        }

      /*
        Find free bucket or bucket for identical sequence.
        Make sure sequences are exactly identical
        in case of any hash collision.
        With 64-bit hashes, there is about 50% chance of a
        collision when the number of sequences is about 5e9.
      */

      uint128 hash = HASH(seq_up, seqlen);
      uint64_t j =  hash2bucket(hash, hashtablesize);
      struct sm_bucket * bp = hashtable + j;

      while ((bp->size) && (hash != bp->hash))
        {
          j = next_bucket(j, hashtablesize);
          bp = hashtable + j;
        }

      if ((opt_strand > 1) && !bp->size)
        {
          /* no match on plus strand */
          /* check minus strand as well */

          uint128 rc_hash = HASH(rc_seq_up, seqlen);
          uint64_t k =  hash2bucket(rc_hash, hashtablesize);
          struct sm_bucket * rc_bp = hashtable + k;

          while ((rc_bp->size) && (rc_hash != rc_bp->hash))
            {
              k = next_bucket(k, hashtablesize);
              rc_bp = hashtable + k;
            }

          if (rc_bp->size)
            {
              bp = rc_bp;
              j = k;
            }
        }

      int abundance = fastx_get_abundance(h);
      int64_t ab = opt_sizein ? abundance : 1;
      sumsize += ab;

      if (bp->size)
        {
          /* at least one identical sequence already */
          bp->size += ab;
        }
      else
        {
          /* no identical sequences yet */
          bp->size = ab;
          bp->hash = hash;
          clusters++;
        }

      if (bp->size > maxsize)
        {
          maxsize = bp->size;
        }

      sequencecount++;
      progress_update(fastx_get_position(h));
    }
  progress_done();
  xfree(prompt);
  fastx_close(h);

  show_rusage();

  if (!opt_quiet)
    {
      if (sequencecount > 0)
        {
          fprintf(stderr,
                  "%'" PRIu64 " nt in %'" PRIu64 " seqs, min %'" PRIu64
                  ", max %'" PRIu64 ", avg %'.0f\n",
                  nucleotidecount,
                  sequencecount,
                  shortest,
                  longest,
                  nucleotidecount * 1.0 / sequencecount);
        }
      else
        {
          fprintf(stderr,
                  "%'" PRIu64 " nt in %'" PRIu64 " seqs\n",
                  nucleotidecount,
                  sequencecount);
        }
    }

  if (opt_log)
    {
      if (sequencecount > 0)
        {
          fprintf(fp_log,
                  "%'" PRIu64 " nt in %'" PRIu64 " seqs, min %'" PRIu64
                  ", max %'" PRIu64 ", avg %'.0f\n",
                  nucleotidecount,
                  sequencecount,
                  shortest,
                  longest,
                  nucleotidecount * 1.0 / sequencecount);
        }
      else
        {
          fprintf(fp_log,
                  "%'" PRIu64 " nt in %'" PRIu64 " seqs\n",
                  nucleotidecount,
                  sequencecount);
        }
    }

  if (discarded_short)
    {
      fprintf(stderr,
              "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_minseqlength,
              discarded_short,
              (discarded_short == 1 ? "sequence" : "sequences"));

      if (opt_log)
        {
          fprintf(fp_log,
                  "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  opt_minseqlength,
                  discarded_short,
                  (discarded_short == 1 ? "sequence" : "sequences"));
        }
    }

  if (discarded_long)
    {
      fprintf(stderr,
              "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_maxseqlength,
              discarded_long,
              (discarded_long == 1 ? "sequence" : "sequences"));

      if (opt_log)
        {
          fprintf(fp_log,
                  "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  opt_maxseqlength,
                  discarded_long,
                  (discarded_long == 1 ? "sequence" : "sequences"));
        }
    }


  show_rusage();

  average = 1.0 * sumsize / clusters;
  median = find_median();

  if (clusters < 1)
    {
      if (!opt_quiet)
        {
          fprintf(stderr,
                  "0 unique sequences\n");
        }
      if (opt_log)
        {
          fprintf(fp_log,
                  "0 unique sequences\n\n");
        }
    }
  else
    {
      if (!opt_quiet)
        {
          fprintf(stderr,
                  "%" PRId64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n",
                  clusters, average, median, maxsize);
        }
      if (opt_log)
        {
          fprintf(fp_log,
                  "%" PRId64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n\n",
                  clusters, average, median, maxsize);
        }
    }

  show_rusage();

  /* second pass with output */

  fastx_handle h2 = fastx_open(input_filename);
  if (!h2)
    {
      fatal("Cannot open and read from the input file.");
    }

  progress_init("Writing FASTA output file", filesize);

  uint64_t selected = 0;

  while(fastx_next(h2, ! opt_notrunclabels, chrmap_no_change))
    {
      int64_t seqlen = fastx_get_sequence_length(h2);

      if ((seqlen < opt_minseqlength) || (seqlen > opt_maxseqlength))
        {
          continue;
        }

      char * seq = fastx_get_sequence(h2);

      /* normalize sequence: uppercase and replace U by T  */
      string_normalize(seq_up, seq, seqlen);

      /* reverse complement if necessary */
      if (opt_strand > 1)
        {
          reverse_complement(rc_seq_up, seq_up, seqlen);
        }

      uint128 hash = HASH(seq_up, seqlen);
      uint64_t j =  hash2bucket(hash, hashtablesize);
      struct sm_bucket * bp = hashtable + j;

      while ((bp->size) && (hash != bp->hash))
        {
          j = next_bucket(j, hashtablesize);
          bp = hashtable + j;
        }

      if ((opt_strand > 1) && ! bp->size)
        {
          /* no match on plus strand */
          /* check minus strand as well */

          uint128 rc_hash = HASH(rc_seq_up, seqlen);
          uint64_t k =  hash2bucket(rc_hash, hashtablesize);
          struct sm_bucket * rc_bp = hashtable + k;

          while ((rc_bp->size) && (rc_hash != rc_bp->hash))
            {
              k = next_bucket(k, hashtablesize);
              rc_bp = hashtable + k;
            }

          if (rc_bp->size)
            {
              bp = rc_bp;
              j = k;
            }
        }

      int64_t size = bp->size;

      if (size > 0)
        {
          /* print sequence */

          char * header = fastx_get_header(h2);
          int headerlen = fastx_get_header_length(h2);

          if ((size >= opt_minuniquesize) && (size <= opt_maxuniquesize))
            {
              selected++;
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  seq,
                                  seqlen,
                                  header,
                                  headerlen,
                                  size,
                                  selected,
                                  -1.0,
                                  -1, -1, nullptr, 0.0);
            }
          bp->size = -1;
        }

      progress_update(fastx_get_position(h2));
    }
  progress_done();
  fastx_close(h2);
  fclose(fp_fastaout);

  show_rusage();

  if (selected < clusters)
    {
      if (!opt_quiet)
        {
          fprintf(stderr,
                  "%" PRId64 " uniques written, %"
                  PRId64 " clusters discarded (%.1f%%)\n",
                  selected, clusters - selected,
                  100.0 * (clusters - selected) / clusters);
        }

      if (opt_log)
        {
          fprintf(fp_log,
                  "%" PRId64 " uniques written, %"
                  PRId64 " clusters discarded (%.1f%%)\n\n",
                  selected, clusters - selected,
                  100.0 * (clusters - selected) / clusters);
        }
    }

  show_rusage();

  xfree(seq_up);
  xfree(rc_seq_up);
  xfree(hashtable);

  show_rusage();
}
