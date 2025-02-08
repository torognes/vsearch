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
#include "utils/seqcmp.h"
#include <algorithm>  // std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdlib>  // std::qsort
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstring>  // std::strcmp, std::memset
#include <limits>
#include <vector>


struct bucket
{
  uint64_t hash = 0;
  unsigned int seqno_first = 0;
  unsigned int seqno_last = 0;
  unsigned int size = 0;
  unsigned int count = 0;
  bool deleted = false;
  char * header = nullptr;
  char * seq = nullptr;
  char * qual = nullptr;
};


auto derep_compare_prefix(const void * a, const void * b) -> int
{
  auto * lhs = (struct bucket *) a;
  auto * rhs = (struct bucket *) b;

  /* deleted(?) first, then by highest abundance, then by label, otherwise keep order */

  if (static_cast<int>(lhs->deleted) > static_cast<int>(rhs->deleted))
    {
      return +1;
    }
  if (static_cast<int>(lhs->deleted) < static_cast<int>(rhs->deleted))
    {
      return -1;
    }

  // both are deleted, compare abundances
  if (lhs->size < rhs->size)
    {
      return +1;
    }
  if (lhs->size > rhs->size)
    {
      return -1;
    }

  // both are deleted, same abundances, compare sequence headers
  auto const result = std::strcmp(db_getheader(lhs->seqno_first),
                                  db_getheader(rhs->seqno_first));
  if (result != 0)
    {
      return result;
    }

  // both are deleted, same abundances, same sequence headers, compare input order
  if (lhs->seqno_first < rhs->seqno_first)
    {
      return -1;
    }
  if (lhs->seqno_first > rhs->seqno_first)
    {
      return +1;
    }
  return 0;
}


auto derep_prefix(struct Parameters const & parameters) -> void
{
  std::FILE * fp_output = nullptr;
  std::FILE * fp_uc = nullptr;

  if (parameters.opt_strand)
    {
      fatal("Option '--strand both' not supported with --derep_prefix");
    }

  if (parameters.opt_output != nullptr)
    {
      fp_output = fopen_output(parameters.opt_output);
      if (fp_output == nullptr)
        {
          fatal("Unable to open output file for writing");
        }
    }

  if (parameters.opt_uc != nullptr)
    {
      fp_uc = fopen_output(parameters.opt_uc);
      if (fp_uc == nullptr)
        {
          fatal("Unable to open output (uc) file for writing");
        }
    }

  db_read(parameters.opt_derep_prefix, 0);

  db_sortbylength_shortest_first();

  show_rusage();

  int64_t const dbsequencecount = db_getsequencecount();

  /* adjust size of hash table for 2/3 fill rate */

  int64_t hashtablesize = 1;
  while (3 * dbsequencecount > 2 * hashtablesize)
    {
      hashtablesize <<= 1U;
    }
  int const hash_mask = hashtablesize - 1;

  std::vector<struct bucket> hashtable(hashtablesize);

  int64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;
  double median = 0.0;
  double average = 0.0;

  /* alloc and init table of links to other sequences in cluster */

  constexpr auto terminal = std::numeric_limits<unsigned int>::max();
  std::vector<unsigned int> nextseqtab(dbsequencecount, terminal);

  std::vector<char> seq_up(db_getlongestsequence() + 1);

  /* make table of hash values of prefixes */

  unsigned int const len_longest = db_getlongestsequence();
  unsigned int const len_shortest = db_getshortestsequence();
  std::vector<uint64_t> prefix_hashes(len_longest + 1);

  progress_init("Dereplicating", dbsequencecount);
  for (int64_t i = 0; i < dbsequencecount; i++)
    {
      unsigned int const seqlen = db_getsequencelen(i);
      char * seq = db_getsequence(i);

      /* normalize sequence: uppercase and replace U by T  */
      string_normalize(seq_up.data(), seq, seqlen);

      uint64_t const ab = parameters.opt_sizein ? db_getabundance(i) : 1;
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
      for (unsigned int j = 0; j < seqlen; j++)
        {
          fnv1a_hash ^= seq_up[j];
          fnv1a_hash *= 1099511628211ULL;
          prefix_hashes[j + 1] = fnv1a_hash;
        }

      /* first, look for an identical match */

      unsigned int prefix_len = seqlen;

      uint64_t hash = prefix_hashes[prefix_len];
      struct bucket * bp = &hashtable[hash & hash_mask];

      while ((bp->size != 0U) and
             ((bp->deleted) or
              (bp->hash != hash) or
              (prefix_len != db_getsequencelen(bp->seqno_first)) or
              (seqcmp(seq_up.data(), db_getsequence(bp->seqno_first), prefix_len) != 0)))
        {
          ++bp;
          if (bp >= &hashtable[hashtablesize])
            {
              bp = hashtable.data();
            }
        }

      /* at this point, bp points either to (1) a free empty hash bucket, or
         (2) a bucket with an exact match. */

      auto const orig_hash = hash;
      struct bucket * orig_bp = bp;

      if (bp->size != 0U)
        {
          /* exact match */
          bp->size += ab;
          auto const last = bp->seqno_last;
          nextseqtab[last] = i;
          bp->seqno_last = i;

          maxsize = std::max<uint64_t>(bp->size, maxsize);
        }
      else
        {
          /* look for prefix match */

          while ((bp->size == 0U) and (prefix_len > len_shortest))
            {
              --prefix_len;
              hash = prefix_hashes[prefix_len];
              bp = &hashtable[hash & hash_mask];

              while ((bp->size != 0U) and
                     ((bp->deleted) or
                      (bp->hash != hash) or
                      (prefix_len != db_getsequencelen(bp->seqno_first)) or
                      (seqcmp(seq_up.data(),
                              db_getsequence(bp->seqno_first),
                              prefix_len) != 0)))
                {
                  ++bp;
                  if (bp >= &hashtable[hashtablesize])
                    {
                      bp = hashtable.data();
                    }
                }
            }

          if (bp->size != 0U)
            {
              /* prefix match */

              /* get necessary info, then delete prefix from hash */
              unsigned int const first = bp->seqno_first;
              unsigned int const last = bp->seqno_last;
              unsigned int const size = bp->size;
              bp->deleted = true;

              /* create new hash entry */
              bp = orig_bp;
              bp->size = size + ab;
              bp->hash = orig_hash;
              bp->seqno_first = i;
              nextseqtab[i] = first;
              bp->seqno_last = last;

              maxsize = std::max<uint64_t>(bp->size, maxsize);
            }
          else
            {
              /* no match */
              orig_bp->size = ab;
              orig_bp->hash = orig_hash;
              orig_bp->seqno_first = i;
              orig_bp->seqno_last = i;

              maxsize = std::max(ab, maxsize);
              ++clusters;
            }
        }

      progress_update(i);
    }
  progress_done();

  show_rusage();

  progress_init("Sorting", 1);
  qsort(hashtable.data(), hashtablesize, sizeof(struct bucket), derep_compare_prefix);
  progress_done();

  if (clusters > 0)
    {
      if ((clusters % 2) != 0)
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
      if (parameters.opt_log != nullptr)
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

  /* count selected */

  int64_t selected = 0;
  for (int64_t i = 0; i < clusters; i++)
    {
      struct bucket * bp = &hashtable[i];
      int64_t const size = bp->size;
      if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
        {
          ++selected;
          if (selected == parameters.opt_topn)
            {
              break;
            }
        }
    }


  /* write output */

  if (parameters.opt_output != nullptr)
    {
      progress_init("Writing output file", clusters);

      int64_t relabel_count = 0;
      for (int64_t i = 0; i < clusters; i++)
        {
          struct bucket * bp = &hashtable[i];
          int64_t const size = bp->size;
          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++relabel_count;
              fasta_print_general(fp_output,
                                  nullptr,
                                  db_getsequence(bp->seqno_first),
                                  db_getsequencelen(bp->seqno_first),
                                  db_getheader(bp->seqno_first),
                                  db_getheaderlen(bp->seqno_first),
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
      fclose(fp_output);
    }

  show_rusage();

  if (parameters.opt_uc != nullptr)
    {
      progress_init("Writing uc file, first part", clusters);
      for (int64_t i = 0; i < clusters; i++)
        {
          struct bucket * bp = &hashtable[i];
          char * h =  db_getheader(bp->seqno_first);
          int64_t const len = db_getsequencelen(bp->seqno_first);

          fprintf(fp_uc, "S\t%" PRId64 "\t%" PRId64 "\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, len, h);

          for (unsigned int next = nextseqtab[bp->seqno_first];
               next != terminal;
               next = nextseqtab[next])
            {
              fprintf(fp_uc,
                      "H\t%" PRId64 "\t%" PRIu64 "\t%.1f\t+\t0\t0\t*\t%s\t%s\n",
                      i, db_getsequencelen(next), 100.0, db_getheader(next), h);
            }

          progress_update(i);
        }
      progress_done();
      show_rusage();

      progress_init("Writing uc file, second part", clusters);
      for (int64_t i = 0; i < clusters; i++)
        {
          struct bucket * bp = &hashtable[i];
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
      if (not parameters.opt_quiet)
        {
          fprintf(stderr,
                  "%" PRId64 " uniques written, %" PRId64
                  " clusters discarded (%.1f%%)\n",
                  selected, clusters - selected,
                  100.0 * (clusters - selected) / clusters);
        }

      if (parameters.opt_log != nullptr)
        {
          fprintf(fp_log,
                  "%" PRId64 " uniques written, %" PRId64
                  " clusters discarded (%.1f%%)\n\n",
                  selected, clusters - selected,
                  100.0 * (clusters - selected) / clusters);
        }
    }

  db_free();
}
