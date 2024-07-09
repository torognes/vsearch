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
#include "utils/seqcmp.h"
#include <algorithm>  // std::min
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::log10
#include <cstdint> // int64_t, uint64_t
#include <cstdlib>  // std::qsort
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstring>  // std::strlen, std::strcmp, std::memset
#include <limits>
#include <string>
#include <vector>


#define HASH hash_cityhash64

struct bucket
{
  uint64_t hash;
  unsigned int seqno_first;
  unsigned int seqno_last;
  unsigned int size;
  unsigned int count;
  bool deleted;
  char * header;
  char * seq;
  char * qual;
};


auto derep_compare_full(void const * void_lhs, void const * void_rhs) -> int
{
  auto * lhs = (struct bucket *) void_lhs;
  auto * rhs = (struct bucket *) void_rhs;

  /* highest abundance first, then by label, otherwise keep order */

  if (lhs->deleted > rhs->deleted)  // refactoring: deleted is always set to false for derep_fulllength
    {
      return +1;
    }
  else if (lhs->deleted < rhs->deleted)  // refactoring: deleted is always set to false for derep_fulllength
    {
      return -1;
    }
  else
    {
      if (lhs->size < rhs->size)
        {
          return +1;
        }
      else if (lhs->size > rhs->size)
        {
          return -1;
        }
      else
        {
          if (lhs->size == 0)
            {
              return 0;
            }
          auto const result = std::strcmp(lhs->header, rhs->header);
          if (result != 0)
            {
              return result;
            }
          else
            {
              if (lhs->seqno_first < rhs->seqno_first)
                {
                  return -1;
                }
              else if (lhs->seqno_first > rhs->seqno_first)
                {
                  return +1;
                }
              else
                {
                  return 0;  // unreachable
                }
            }
        }
    }
}


auto rehash(struct bucket ** hashtableref, int64_t alloc_clusters) -> void
{
  /*
    double the size of the hash table:

    - allocate the new hash table
    - rehash all entries from the old to the new table
    - free the old table
    - update variables
  */

  struct bucket * old_hashtable = *hashtableref;
  uint64_t const old_hashtablesize = 2 * alloc_clusters;
  uint64_t const new_hashtablesize = 2 * old_hashtablesize;
  uint64_t const new_hash_mask = new_hashtablesize - 1;

  auto * new_hashtable =
    (struct bucket *) xmalloc(sizeof(struct bucket) * new_hashtablesize);
  memset(new_hashtable, 0, sizeof(struct bucket) * new_hashtablesize);

  /* rehash all */
  for (uint64_t i = 0; i < old_hashtablesize; i++)
    {
      struct bucket * old_bp = old_hashtable + i;
      if (old_bp->size)
        {
          uint64_t k = old_bp->hash & new_hash_mask;
          while (new_hashtable[k].size)
            {
              k = (k + 1) & new_hash_mask;
            }
          struct bucket * new_bp = new_hashtable + k;

          * new_bp = * old_bp;
        }
    }

  xfree(old_hashtable);
  * hashtableref = new_hashtable;
}


inline auto convert_quality_to_probability(int const quality_symbol, struct Parameters const & parameters) -> double
{
  int const x = quality_symbol - parameters.opt_fastq_ascii;
  if (x < 2)
    {
      return 0.75;
    }
  else
    {
      return exp10(-x / 10.0);
    }
}


inline auto convert_probability_to_quality(double const p, struct Parameters const & parameters) -> int
{
  // int q = round(-10.0 * log10(p));
  int q = int(-10.0 * log10(p));
  q = MIN(q, parameters.opt_fastq_qmaxout);
  q = MAX(q, parameters.opt_fastq_qminout);
  return parameters.opt_fastq_asciiout + q;
}


auto derep(struct Parameters const & parameters, char * input_filename, bool use_header) -> void
{
  /* dereplicate full length sequences, optionally require identical headers */

  /*
    derep_fulllength output options: --output, --uc (only FASTA, depreciated)
    fastx_uniques output options: --fastaout, --fastqout, --uc, --tabbedout
  */

  show_rusage();

  fastx_handle input_handle = fastx_open(input_filename);

  if (not input_handle)
    {
      fatal("Unrecognized input file type (not proper FASTA or FASTQ format)");  // unreachable? case already handled in fastx_open(), assert(h != nullptr) should always be true
    }

  if (not fastx_is_empty(input_handle))
    {
      if (fastx_is_fastq(input_handle))
        {
          if (not parameters.opt_fastx_uniques) {
            fatal("FASTQ input is only allowed with the fastx_uniques command");
          }
        }
      else
        {
          if (parameters.opt_fastqout) {
            fatal("Cannot write FASTQ output when input file is not in FASTQ "
                  "format");
          }
          if (parameters.opt_tabbedout) {
            fatal("Cannot write tab separated output file when input file is "
                  "not in FASTQ format");
          }
        }
    }

  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_uc = nullptr;
  std::FILE * fp_tabbedout = nullptr;

  if (parameters.opt_fastx_uniques)
    {
      if ((not parameters.opt_uc) and (not parameters.opt_fastaout) and (not parameters.opt_fastqout) and (not parameters.opt_tabbedout)) {
        fatal("Output file for dereplication with fastx_uniques must be "
              "specified with --fastaout, --fastqout, --tabbedout, or --uc");
      }
    } else {
    if ((not parameters.opt_output) and (not parameters.opt_uc)) {
      fatal("Output file for dereplication must be specified with --output "
            "or --uc");
    }
  }

  if (parameters.opt_fastx_uniques)
    {
      if (parameters.opt_fastaout)
        {
          fp_fastaout = fopen_output(parameters.opt_fastaout);
          if (not fp_fastaout)
            {
              fatal("Unable to open FASTA output file for writing");
            }
        }

      if (parameters.opt_fastqout)
        {
          fp_fastqout = fopen_output(parameters.opt_fastqout);
          if (not fp_fastqout)
            {
              fatal("Unable to open FASTQ output file for writing");
            }
        }

      if (parameters.opt_tabbedout)
        {
          fp_tabbedout = fopen_output(parameters.opt_tabbedout);
          if (not fp_tabbedout)
            {
              fatal("Unable to open tab delimited output file for writing");
            }
        }
    }
  else
    {
      if (parameters.opt_output)
        {
          fp_fastaout = fopen_output(parameters.opt_output);
          if (not fp_fastaout)
            {
              fatal("Unable to open FASTA output file for writing");
            }
        }
    }

  if (parameters.opt_uc)
    {
      fp_uc = fopen_output(parameters.opt_uc);
      if (not fp_uc)
        {
          fatal("Unable to open output (uc) file for writing");
        }
    }

  uint64_t const filesize = fastx_get_size(input_handle);


  /* allocate initial memory for 1024 clusters
     with sequences of length 1023 */

  uint64_t alloc_clusters = 1024;
  uint64_t alloc_seqs = 1024;
  int64_t alloc_seqlen = 1023;

  uint64_t hashtablesize = 2 * alloc_clusters;
  uint64_t hash_mask = hashtablesize - 1;
  auto * hashtable =
    (struct bucket *) xmalloc(sizeof(struct bucket) * hashtablesize);
  memset(hashtable, 0, sizeof(struct bucket) * hashtablesize);

  show_rusage();

  constexpr auto terminal = std::numeric_limits<unsigned int>::max();
  std::vector<unsigned int> nextseqtab;
  std::vector<std::string> headertab;
  std::vector<char> match_strand;

  auto const extra_info = parameters.opt_uc or parameters.opt_tabbedout;

  if (extra_info)
    {
      /* If the uc or tabbedout option is in effect,
         we need to keep some extra info.
         Allocate and init memory for this. */

      /* Links to other sequences in cluster */
      nextseqtab.resize(alloc_seqs, terminal);

      /* Pointers to the header strings */
      headertab.resize(alloc_seqs);

      /* Matching strand */
      match_strand.resize(alloc_seqs);
    }

  show_rusage();

  std::vector<char> seq_up(alloc_seqlen + 1);
  std::vector<char> rc_seq_up(alloc_seqlen + 1);
  std::string prompt = std::string("Dereplicating file ") + input_filename;

  progress_init(prompt.c_str(), filesize);

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

  while (fastx_next(input_handle, not parameters.opt_notrunclabels, chrmap_no_change))
    {
      int64_t const seqlen = fastx_get_sequence_length(input_handle);

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
          seq_up.resize(alloc_seqlen + 1);
          rc_seq_up.resize(alloc_seqlen + 1);

          show_rusage();
        }

      if (extra_info and (sequencecount + 1 > alloc_seqs))
        {
          uint64_t const new_alloc_seqs = 2 * alloc_seqs;

          nextseqtab.resize(new_alloc_seqs, terminal);

          headertab.resize(new_alloc_seqs);

          match_strand.resize(new_alloc_seqs);

          alloc_seqs = new_alloc_seqs;

          show_rusage();
        }

      if (clusters + 1 > alloc_clusters)
        {
          uint64_t const new_alloc_clusters = 2 * alloc_clusters;

          rehash(& hashtable, alloc_clusters);

          alloc_clusters = new_alloc_clusters;
          hashtablesize = 2 * alloc_clusters;
          hash_mask = hashtablesize - 1;

          show_rusage();
        }

      char * seq = fastx_get_sequence(input_handle);
      char * header = fastx_get_header(input_handle);
      int64_t const headerlen = fastx_get_header_length(input_handle);
      char * qual = fastx_get_quality(input_handle); // nullptr if FASTA

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

      uint64_t hash_header = 0;
      if (use_header)
        {
          hash_header = HASH(header, headerlen);
        }
      else
        {
          hash_header = 0;
        }

      uint64_t const hash = HASH(seq_up.data(), seqlen) ^ hash_header;
      uint64_t j = hash & hash_mask;
      struct bucket * bp = hashtable + j;

      while ((bp->size) and
             ((hash != bp->hash) or
              (seqcmp(seq_up.data(), bp->seq, seqlen)) or
              (use_header and strcmp(header, bp->header))))
        {
          j = (j + 1) & hash_mask;
          bp = hashtable + j;
        }

      if (parameters.opt_strand and not bp->size)
        {
          /* no match on plus strand */
          /* check minus strand as well */

          uint64_t const rc_hash = HASH(rc_seq_up.data(), seqlen) ^ hash_header;
          uint64_t k = rc_hash & hash_mask;
          struct bucket * rc_bp = hashtable + k;

          while ((rc_bp->size)
                 and
                 ((rc_hash != rc_bp->hash) or
                  (seqcmp(rc_seq_up.data(), rc_bp->seq, seqlen)) or
                  (use_header and strcmp(header, rc_bp->header))))
            {
              k = (k + 1) & hash_mask;
              rc_bp = hashtable + k;
            }

          if (rc_bp->size)
            {
              bp = rc_bp;
              j = k;
              if (extra_info)
                {
                  match_strand[sequencecount] = 1;
                }
            }
        }

      int const abundance = fastx_get_abundance(input_handle);
      int64_t const ab = parameters.opt_sizein ? abundance : 1;
      sumsize += ab;

      if (bp->size)
        {
          /* at least one identical sequence already */
          if (extra_info)
            {
              unsigned int const last = bp->seqno_last;
              nextseqtab[last] = sequencecount;
              bp->seqno_last = sequencecount;
              headertab[sequencecount] = header;
            }

          int64_t const s1 = bp->size;
          int64_t const s2 = ab;
          int64_t const s3 = s1 + s2;

          if (parameters.opt_fastqout)
            {
              /* update quality scores */
              for (int i = 0; i < seqlen; i++)
                {
                  int const q1 = bp->qual[i];
                  int const q2 = qual[i];
                  double const p1 = convert_quality_to_probability(q1, parameters);
                  double const p2 = convert_quality_to_probability(q2, parameters);
                  double p3 = 0.0;

                  /* how to compute the new quality score? */

                  if (parameters.opt_fastq_qout_max)
                    {
                      // fastq_qout_max
                      /* min error prob, highest quality */
                      p3 = std::min(p1, p2);
                    }
                  else
                    {
                      // fastq_qout_avg
                      /* average, as in USEARCH */
                      p3 = (p1 * s1 + p2 * s2) / s3;
                    }

                  // fastq_qout_min
                  /* max error prob, lowest quality */
                  // p3 = MAX(p1, p2);

                  // fastq_qout_first
                  /* keep first */
                  // p3 = p1;

                  // fastq_qout_last
                  /* keep last */
                  // p3 = p2;

                  // fastq_qout_ef
                  /* Compute as multiple independent observations
                     Edgar & Flyvbjerg (2015)
                     But what about s1 and s2? */
                  // p3 = p1 * p2 / 3.0 / (1.0 - p1 - p2 + (4.0 * p1 * p2 / 3.0));

                  /* always worst quality possible, certain error */
                  // p3 = 1.0;

                  // always best quality possible, perfect, no errors */
                  // p3 = 0.0;

                  int const q3 = convert_probability_to_quality(p3, parameters);
                  bp->qual[i] = q3;
                }
            }

          bp->size = s3;
          ++bp->count;
        }
      else
        {
          /* no identical sequences yet */
          bp->size = ab;
          bp->hash = hash;
          bp->seqno_first = sequencecount;
          bp->seqno_last = sequencecount;
          bp->seq = xstrdup(seq);
          bp->header = xstrdup(header);
          bp->count = 1;
          if (qual) {
            bp->qual = xstrdup(qual);
          } else {
            bp->qual = nullptr;
          }
          ++clusters;
        }

      if (bp->size > maxsize)
        {
          maxsize = bp->size;
        }

      ++sequencecount;

      progress_update(fastx_get_position(input_handle));
    }
  progress_done();
  fastx_close(input_handle);

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

  if (parameters.opt_log)
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

  if (discarded_short)
    {
      fprintf(stderr,
              "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              parameters.opt_minseqlength,
              discarded_short,
              (discarded_short == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log)
        {
          fprintf(fp_log,
                  "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  parameters.opt_minseqlength,
                  discarded_short,
                  (discarded_short == 1 ? "sequence" : "sequences"));
        }
    }

  if (discarded_long)
    {
      fprintf(stderr,
              "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              parameters.opt_maxseqlength,
              discarded_long,
              (discarded_long == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log)
        {
          fprintf(fp_log,
                  "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  parameters.opt_maxseqlength,
                  discarded_long,
                  (discarded_long == 1 ? "sequence" : "sequences"));
        }
    }

  show_rusage();

  progress_init("Sorting", 1);
  qsort(hashtable, hashtablesize, sizeof(struct bucket), derep_compare_full);
  progress_done();

  show_rusage();

  if (clusters > 0)
    {
      if (clusters % 2)
        {
          median = hashtable[(clusters - 1) / 2].size;
        }
      else
        {
          median = (hashtable[(clusters / 2) - 1].size +
                    hashtable[clusters / 2].size) / 2.0;
        }
    }

  average = 1.0 * sumsize / clusters;

  if (clusters < 1)
    {
      if (not parameters.opt_quiet)
        {
          fprintf(stderr,
                  "0 unique sequences\n");
        }
      if (parameters.opt_log)
        {
          fprintf(fp_log,
                  "0 unique sequences\n\n");
        }
    }
  else
    {
      if (not parameters.opt_quiet)
        {
          fprintf(stderr,
                  "%" PRId64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n",
                  clusters, average, median, maxsize);
        }
      if (parameters.opt_log)
        {
          fprintf(fp_log,
                  "%" PRId64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n\n",
                  clusters, average, median, maxsize);
        }
    }

  /* count selected */

  uint64_t selected = 0;
  for (uint64_t i = 0; i < clusters; ++i)
    {
      struct bucket * bp = hashtable + i;
      int64_t const size = bp->size;
      if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
        {
          ++selected;
          if (selected == (uint64_t) parameters.opt_topn)
            {
              break;
            }
        }
    }

  show_rusage();

  /* write output */

  if (parameters.opt_output or parameters.opt_fastaout)
    {
      progress_init("Writing FASTA output file", clusters);

      int64_t relabel_count = 0;
      for (uint64_t i = 0; i < clusters; ++i)
        {
          struct bucket * bp = hashtable + i;
          int64_t const size = bp->size;
          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++relabel_count;
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  bp->seq,
                                  strlen(bp->seq),
                                  bp->header,
                                  strlen(bp->header),
                                  size,
                                  relabel_count,
                                  -1.0,
                                  -1, -1, nullptr, 0.0);
              if (relabel_count == parameters.opt_topn)
                {
                  break;
                }
            }
          progress_update(i);
        }

      progress_done();
      fclose(fp_fastaout);
    }

  if (parameters.opt_fastqout)
    {
      progress_init("Writing FASTQ output file", clusters);

      int64_t relabel_count = 0;
      for (uint64_t i = 0; i < clusters; ++i)
        {
          struct bucket * bp = hashtable + i;
          int64_t const size = bp->size;
          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++relabel_count;
              fastq_print_general(fp_fastqout,
                                  bp->seq,
                                  strlen(bp->seq),
                                  bp->header,
                                  strlen(bp->header),
                                  bp->qual,
                                  size,
                                  relabel_count,
                                  -1.0);
              if (relabel_count == parameters.opt_topn)
                {
                  break;
                }
            }
          progress_update(i);
        }

      progress_done();
      fclose(fp_fastqout);
    }

  show_rusage();

  if (parameters.opt_uc)
    {
      progress_init("Writing uc file, first part", clusters);
      for (uint64_t i = 0; i < clusters; ++i)
        {
          struct bucket * bp = hashtable + i;
          char * hh =  bp->header;
          int64_t const len = strlen(bp->seq);

          fprintf(fp_uc, "S\t%" PRId64 "\t%" PRId64 "\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, len, hh);

          for (unsigned int next = nextseqtab[bp->seqno_first];
               next != terminal;
               next = nextseqtab[next])
            {
              fprintf(fp_uc,
                      "H\t%" PRId64 "\t%" PRId64 "\t%.1f\t%s\t0\t0\t*\t%s\t%s\n",
                      i, len, 100.0,
                      (match_strand[next] ? "-" : "+"),
                      headertab[next].c_str(), hh);
            }

          progress_update(i);
        }
      progress_done();

      progress_init("Writing uc file, second part", clusters);
      for (uint64_t i = 0; i < clusters; ++i)
        {
          struct bucket * bp = hashtable + i;
          fprintf(fp_uc, "C\t%" PRId64 "\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, bp->size, bp->header);
          progress_update(i);
        }
      fclose(fp_uc);
      progress_done();
    }

  if (parameters.opt_tabbedout)
    {
      progress_init("Writing tab separated file", clusters);
      for (uint64_t i = 0; i < clusters; ++i)
        {
          struct bucket * bp = hashtable + i;
          char * hh =  bp->header;

          if (parameters.opt_relabel) {
            fprintf(fp_tabbedout,
                    "%s\t%s%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                    hh, parameters.opt_relabel, i + 1, i, (uint64_t) 0, bp->count, hh);
          } else {
            fprintf(fp_tabbedout, "%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                    hh, hh, i, (uint64_t) 0, bp->count, hh);
          }

          uint64_t j = 1;
          for (unsigned int next = nextseqtab[bp->seqno_first];
               next != terminal;
               next = nextseqtab[next])
            {
              if (parameters.opt_relabel) {
                fprintf(fp_tabbedout,
                        "%s\t%s%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                        headertab[next].c_str(), parameters.opt_relabel, i + 1, i, j, bp->count, hh);
              } else {
                fprintf(fp_tabbedout,
                        "%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                        headertab[next].c_str(), hh, i, j, bp->count, hh);
              }
              ++j;
            }

          progress_update(i);
        }
      fclose(fp_tabbedout);
      progress_done();
    }


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

      if (parameters.opt_log)
        {
          fprintf(fp_log,
                  "%" PRId64 " uniques written, %"
                  PRId64 " clusters discarded (%.1f%%)\n\n",
                  selected, clusters - selected,
                  100.0 * (clusters - selected) / clusters);
        }
    }

  show_rusage();

  /* Free all seqs and headers */

  for (uint64_t i = 0; i < clusters; ++i)
    {
      struct bucket * bp = hashtable + i;
      if (bp->size)
        {
          xfree(bp->seq);
          xfree(bp->header);
          if (bp->qual) {
            xfree(bp->qual);
          }
        }
    }

  show_rusage();

  xfree(hashtable);

  show_rusage();
}
