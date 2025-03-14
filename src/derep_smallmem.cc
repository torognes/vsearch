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
#include "city.h"
#include "maps.h"
// #include "util.h"  // hash_cityhash128, Uint128Low64
#include <algorithm>  // std::min, std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::fclose
#include <cstdlib>  // std::qsort
#include <cstring>  // refactoring: unused?
#include <limits>
#include <string>
#include <vector>


using Hash = decltype(&hash_cityhash128);
static Hash hash_function = hash_cityhash128;


struct sm_bucket
{
  uint128 hash;
  uint64_t size;
};

static struct sm_bucket * hashtable = nullptr;
static uint64_t hashtablesize = 0;


auto find_median() -> double
{
  /* find the median size, based on an iterative search starting at e.g. 1 */

  uint64_t cand = 1;    /* candidate for the median */
  uint64_t below = 0;   /* closest value below the candidate */
  uint64_t above = 0;   /* closest value above the candidate */

  uint64_t cand_count = 0;  /* number of clusters with same size as cand */
  uint64_t below_count = 0; /* number of clusters with smaller size than cand */
  uint64_t above_count = 0; /* number of clusters with larger size than cand */

  while (true)
    {
      cand_count = 0;
      below_count = 0;
      above_count = 0;

      for (uint64_t i = 0; i < hashtablesize; i++)
        {
          auto const v = hashtable[i].size;
          if (v > 0)
            {
              if (v > cand)
                {
                  if ((above_count == 0) or (v < above))
                    {
                      above = v;
                    }
                  ++above_count;
                }
              else if (v < cand)
                {
                  if ((below_count == 0) or (v > below))
                    {
                      below = v;
                    }
                  ++below_count;
                }
              else
                {
                  ++cand_count;
                }
            }
        }

      if (below_count + cand_count + above_count == 0U) { // fix -Wfloat-equal
        return 0;  // unreachable?
      }

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
              if (above_count + cand_count == below_count)
                // mid == below_count
                // same as:
                // (below_count + cand_count + above_count) / 2 == below_count
                // which simplifies into:
                // above_count + cand_count == below_count
                {
                  return (below + cand) / 2.0;  // cannot reach?
                }
              return cand;
            }
          cand = above;
        }
      else
        {
          cand = below;  // cannot reach?
        }
    }
}


inline auto hash2bucket(uint128 hash, uint64_t htsize) -> uint64_t
{
  return Uint128Low64(hash) % htsize;
}


inline auto next_bucket(uint64_t prev_bucket, uint64_t htsize) -> uint64_t
{
  return (prev_bucket + 1) % htsize;
}


auto rehash_smallmem() -> void
{
  /* allocate new hash table, 50% larger */
  auto const new_hashtablesize = 3 * hashtablesize / 2;
  auto * new_hashtable =
    (struct sm_bucket *) xmalloc(sizeof(struct sm_bucket) * new_hashtablesize);

  /* zero new hash table */
  for (uint64_t j = 0; j < new_hashtablesize; j++)
    {
      new_hashtable[j].hash.first = 0;
      new_hashtable[j].hash.second = 0;
      new_hashtable[j].size = 0;
    }

  /* rehash all from old to new */
  for (uint64_t i = 0; i < hashtablesize; i++)
    {
      auto * old_bp = hashtable + i;
      if (old_bp->size != 0U)
        {
          auto k = hash2bucket(old_bp->hash, new_hashtablesize);
          while (new_hashtable[k].size != 0U)
            {
              k = next_bucket(k, new_hashtablesize);
            }
          auto * new_bp = new_hashtable + k;
          * new_bp = * old_bp;
        }
    }

  /* free old table */
  xfree(hashtable);

  /* update variables */
  hashtable = new_hashtable;
  hashtablesize = new_hashtablesize;
}


auto derep_smallmem(struct Parameters const & parameters) -> void
{
  /*
    dereplicate full length sequences using a small amount of memory
    output options: --fastaout
  */

  show_rusage();

  auto * input_filename = parameters.opt_derep_smallmem;
  auto * h = fastx_open(input_filename);

  if (h == nullptr) // refactoring: already checked by fastx_open()?
    {
      fatal("Unrecognized input file type (not proper FASTA or FASTQ format).");
    }

  if (h->is_pipe)
    {
      fatal("The derep_smallmem command does not support input from a pipe.");
    }

  std::FILE * fp_fastaout = nullptr;

  if (parameters.opt_fastaout != nullptr)
    {
      fp_fastaout = fopen_output(parameters.opt_fastaout);
      if (fp_fastaout == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }
  else
    {
      fatal("Output file for dereplication must be specified with --fastaout");
    }

  auto const filesize = fastx_get_size(h);

  /* allocate initial memory for sequences of length up to 1023 chars */
  int64_t alloc_seqlen = 1024;

  /* allocate initial hashtable with 1024 buckets */

  hashtablesize = 1024;
  hashtable = (struct sm_bucket *) xmalloc(sizeof(struct sm_bucket) * hashtablesize);

  /* zero hash table */
  for (uint64_t j = 0; j < hashtablesize; j++)
    {
      hashtable[j].hash.first = 0;
      hashtable[j].hash.second = 0;
      hashtable[j].size = 0;
    }

  show_rusage();

  std::vector<char> seq_up(alloc_seqlen + 1);
  std::vector<char> rc_seq_up(alloc_seqlen + 1);

  std::string prompt = std::string("Dereplicating file ") + input_filename;

  progress_init(prompt.c_str(), filesize);

  uint64_t sequencecount = 0;
  uint64_t nucleotidecount = 0;
  int64_t shortest = std::numeric_limits<long>::max();
  int64_t longest = 0;
  uint64_t discarded_short = 0;
  uint64_t discarded_long = 0;
  uint64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;

  /* first pass */

  while (fastx_next(h, not parameters.opt_notrunclabels, chrmap_no_change))
    {
      int64_t const seqlen = fastx_get_sequence_length(h);

      if (seqlen < parameters.opt_minseqlength)
        {
          ++discarded_short;
          continue;
        }

      if (seqlen > parameters.opt_maxseqlength)
        {
          ++discarded_long;
          continue;
        }

      nucleotidecount += seqlen;
      longest = std::max(seqlen, longest);
      shortest = std::min(seqlen, shortest);

      /* check allocations */

      if (seqlen > alloc_seqlen)
        {
          alloc_seqlen = seqlen;
          seq_up.resize(alloc_seqlen + 1);
          rc_seq_up.resize(alloc_seqlen + 1);

          show_rusage();
        }

      if (100 * (clusters + 1) > 95 * hashtablesize)
        {
          // keep hash table fill rate at max 95% */
          rehash_smallmem();
          show_rusage();
        }

      auto * seq = fastx_get_sequence(h);

      /* normalize sequence: uppercase and replace U by T  */
      string_normalize(seq_up.data(), seq, seqlen);

      /* reverse complement if necessary */
      if (parameters.opt_strand)
        {
          reverse_complement(rc_seq_up.data(), seq_up.data(), seqlen);
        }

      /*
        Find free bucket or bucket for identical sequence.
        Make sure sequences are exactly identical
        in case of any hash collision.
        With 64-bit hashes, there is about 50% chance of a
        collision when the number of sequences is about 5e9.
      */

      auto const hash = hash_function(seq_up.data(), seqlen);
      auto j =  hash2bucket(hash, hashtablesize);
      auto * bp = hashtable + j;

      while ((bp->size != 0U) and (hash != bp->hash))
        {
          j = next_bucket(j, hashtablesize);
          bp = hashtable + j;
        }

      if (parameters.opt_strand and (bp->size == 0U))
        {
          /* no match on plus strand */
          /* check minus strand as well */

          auto const rc_hash = hash_function(rc_seq_up.data(), seqlen);
          auto k =  hash2bucket(rc_hash, hashtablesize);
          auto * rc_bp = hashtable + k;

          while ((rc_bp->size != 0U) and (rc_hash != rc_bp->hash))
            {
              k = next_bucket(k, hashtablesize);
              rc_bp = hashtable + k;
            }

          if (rc_bp->size != 0U)
            {
              bp = rc_bp;
              j = k;  // cppcheck: 'j' is assigned a value that is never used
            }
        }

      int const abundance = fastx_get_abundance(h);
      int64_t const ab = parameters.opt_sizein ? abundance : 1;
      sumsize += ab;

      if (bp->size != 0U)
        {
          /* at least one identical sequence already */
          bp->size += ab;
        }
      else
        {
          /* no identical sequences yet */
          bp->size = ab;
          bp->hash = hash;
          ++clusters;
        }

      maxsize = std::max(bp->size, maxsize);

      ++sequencecount;
      progress_update(fastx_get_position(h));
    }
  progress_done();
  fastx_close(h);

  show_rusage();

  if (not parameters.opt_quiet)
    {
      if (sequencecount > 0)
        {
          fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, min %" PRIu64
                  ", max %" PRIu64 ", avg %.0f\n",
                  nucleotidecount,
                  sequencecount,
                  shortest,
                  longest,
                  nucleotidecount * 1.0 / sequencecount);
        }
      else
        {
          fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  nucleotidecount,
                  sequencecount);
        }
    }

  if (parameters.opt_log != nullptr)
    {
      if (sequencecount > 0)
        {
          fprintf(fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, min %" PRIu64
                  ", max %" PRIu64 ", avg %.0f\n",
                  nucleotidecount,
                  sequencecount,
                  shortest,
                  longest,
                  nucleotidecount * 1.0 / sequencecount);
        }
      else
        {
          fprintf(fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  nucleotidecount,
                  sequencecount);
        }
    }

  if (discarded_short != 0U)
    {
      fprintf(stderr,
              "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              parameters.opt_minseqlength,
              discarded_short,
              (discarded_short == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log != nullptr)
        {
          fprintf(fp_log,
                  "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  parameters.opt_minseqlength,
                  discarded_short,
                  (discarded_short == 1 ? "sequence" : "sequences"));
        }
    }

  if (discarded_long != 0U)
    {
      fprintf(stderr,
              "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              parameters.opt_maxseqlength,
              discarded_long,
              (discarded_long == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log != nullptr)
        {
          fprintf(fp_log,
                  "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  parameters.opt_maxseqlength,
                  discarded_long,
                  (discarded_long == 1 ? "sequence" : "sequences"));
        }
    }


  show_rusage();

  if (clusters < 1)
    {
      if (not parameters.opt_quiet)
        {
          fprintf(stderr,
                  "0 unique sequences\n");
        }
      if (parameters.opt_log != nullptr)
        {
          fprintf(fp_log,
                  "0 unique sequences\n\n");
        }
    }
  else
    {
      auto const average = 1.0 * sumsize / clusters;
      const auto median = find_median();
      if (not parameters.opt_quiet)
        {
          fprintf(stderr,
                  "%" PRId64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n",
                  clusters, average, median, maxsize);
        }
      if (parameters.opt_log != nullptr)
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

  auto * h2 = fastx_open(input_filename);
  if (h2 == nullptr)
    {
      fatal("Cannot open and read from the input file.");
    }

  progress_init("Writing FASTA output file", filesize);

  uint64_t selected = 0;

  while (fastx_next(h2, not parameters.opt_notrunclabels, chrmap_no_change))
    {
      int64_t const seqlen = fastx_get_sequence_length(h2);

      if ((seqlen < parameters.opt_minseqlength) or (seqlen > parameters.opt_maxseqlength))
        {
          continue;
        }

      auto * seq = fastx_get_sequence(h2);

      /* normalize sequence: uppercase and replace U by T  */
      string_normalize(seq_up.data(), seq, seqlen);

      /* reverse complement if necessary */
      if (parameters.opt_strand)
        {
          reverse_complement(rc_seq_up.data(), seq_up.data(), seqlen);
        }

      auto const hash = hash_function(seq_up.data(), seqlen);
      auto j =  hash2bucket(hash, hashtablesize);
      auto * bp = hashtable + j;

      while ((bp->size != 0U) and (hash != bp->hash))
        {
          j = next_bucket(j, hashtablesize);
          bp = hashtable + j;
        }

      if (parameters.opt_strand and (bp->size == 0U))
        {
          /* no match on plus strand */
          /* check minus strand as well */

          auto const rc_hash = hash_function(rc_seq_up.data(), seqlen);
          auto k =  hash2bucket(rc_hash, hashtablesize);
          auto * rc_bp = hashtable + k;

          while ((rc_bp->size != 0U) and (rc_hash != rc_bp->hash))
            {
              k = next_bucket(k, hashtablesize);
              rc_bp = hashtable + k;
            }

          if (rc_bp->size != 0U)
            {
              bp = rc_bp;
              j = k;  // cppcheck: 'j' is assigned a value that is never used
            }
        }

      int64_t const size = bp->size;

      if (size > 0)
        {
          /* print sequence */

          auto * header = fastx_get_header(h2);
          int const headerlen = fastx_get_header_length(h2);

          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++selected;
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
      if (not parameters.opt_quiet)
        {
          fprintf(stderr,
                  "%" PRId64 " uniques written, %"
                  PRId64 " clusters discarded (%.1f%%)\n",
                  selected, clusters - selected,
                  100.0 * (clusters - selected) / clusters);
        }

      if (parameters.opt_log != nullptr)
        {
          fprintf(fp_log,
                  "%" PRId64 " uniques written, %"
                  PRId64 " clusters discarded (%.1f%%)\n\n",
                  selected, clusters - selected,
                  100.0 * (clusters - selected) / clusters);
        }
    }

  show_rusage();

  xfree(hashtable);

  show_rusage();
}
