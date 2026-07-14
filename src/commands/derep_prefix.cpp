/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "core/db.hpp"
#include "core/fasta.hpp"
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"
#include "utils/seqcmp.hpp"
#include "utils/span.hpp"
#include "utils/string_normalize.hpp"
#include <algorithm>  // std::max, std::transform
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdlib>  // std::qsort
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstring>  // std::strcmp
#include <iterator>  // std::next
#include <limits>
#include <vector>


// anonymous namespace: 'bucket' is a file-local type; derep.cc defines
// a different struct of the same name, so internal linkage here avoids
// a one-definition-rule violation across translation units
namespace {
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
}


// refactoring: FNV-1A is the hashing function used in std::hash
// (designed for fast hash-table and checksum, not crypto). The
// function below might be redundant with std?

auto compute_hashes_of_all_prefixes(std::vector<uint64_t> & prefix_hashes,
                                    Span<char> const sequence) -> void {
  // Fowler-Noll-Vo (FNV-1A) hash function
  static constexpr auto FNV_offset_basis = uint64_t{14695981039346656037U};
  static constexpr auto FNV_prime = uint64_t{1099511628211U};
  auto FNV1a_hash = FNV_offset_basis;
  prefix_hashes[0] = FNV_offset_basis;
  auto incremental_hash = [&FNV1a_hash](char const nucleotide) -> uint64_t {
    FNV1a_hash ^= static_cast<unsigned char>(nucleotide);
    FNV1a_hash *= FNV_prime;
    return FNV1a_hash;
  };
  std::transform(sequence.cbegin(),
                 sequence.cend(),
                 std::next(prefix_hashes.begin()),
                 incremental_hash);
}


/* The std::qsort comparator below needs the database to compare header strings,
   but a C comparator is a plain function pointer that cannot capture it
   (qsort_r, which passes context, is non-portable). derep_prefix() sets this
   file-scope pointer to its Database immediately before the qsort; the sort is
   single-threaded, so the transient sharing is safe. */
namespace {
  Database const * derep_sort_db = nullptr;
}


auto derep_compare_prefix(const void * a, const void * b) -> int
{
  auto const * lhs = static_cast<struct bucket const *>(a);
  auto const * rhs = static_cast<struct bucket const *>(b);

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
  auto const result = std::strcmp(derep_sort_db->getheader(lhs->seqno_first),
                                  derep_sort_db->getheader(rhs->seqno_first));
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
  if (parameters.opt_strand)
    {
      fatal("Option '--strand both' not supported with --derep_prefix");
    }

  auto output_handle = open_optional_output_file(parameters.opt_output, OutputOption{"--output"});
  std::FILE * const fp_output = output_handle.get();
  auto uc_handle = open_optional_output_file(parameters.opt_uc, OutputOption{"--uc"});
  std::FILE * const fp_uc = uc_handle.get();

  Database db;
  db.read(parameters.opt_derep_prefix, 0, parameters);

  db.sortbylength_shortest_first(parameters);

  // memory-intensive: the entire database is now held in memory

  int64_t const dbsequencecount = static_cast<int64_t>(db.getsequencecount());

  /* adjust size of hash table for 2/3 fill rate */

  int64_t hashtablesize = 1;
  while (3 * dbsequencecount > 2 * hashtablesize)
    {
      hashtablesize <<= 1U;
    }
  uint64_t const hash_mask = static_cast<uint64_t>(hashtablesize - 1);

  std::vector<struct bucket> hashtable(static_cast<std::vector<struct bucket>::size_type>(hashtablesize));

  int64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;
  double median = 0.0;
  double average = 0.0;

  /* alloc and init table of links to other sequences in cluster */

  constexpr auto terminal = std::numeric_limits<unsigned int>::max();
  std::vector<unsigned int> nextseqtab(static_cast<std::vector<unsigned int>::size_type>(dbsequencecount), terminal);

  std::vector<char> seq_up(db.getlongestsequence() + 1);

  /* make table of hash values of prefixes */

  unsigned int const len_longest = static_cast<unsigned int>(db.getlongestsequence());
  unsigned int const len_shortest = static_cast<unsigned int>(db.getshortestsequence());
  std::vector<uint64_t> prefix_hashes(len_longest + 1);

  {
    Progress progress("Dereplicating", static_cast<uint64_t>(dbsequencecount), parameters);
    for (int64_t i = 0; i < dbsequencecount; i++)
      {
        unsigned int const seqlen = static_cast<unsigned int>(db.getsequencelen(static_cast<uint64_t>(i)));
        auto const * seq = db.getsequence(static_cast<uint64_t>(i));

        /* normalize sequence: uppercase and replace U by T  */
        string_normalize(seq_up.data(), seq, seqlen);

        auto const abundance = parameters.opt_sizein ? db.getabundance(static_cast<uint64_t>(i)) : uint64_t{1};
        sumsize += static_cast<int64_t>(abundance);

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

        compute_hashes_of_all_prefixes(prefix_hashes, Span<char>{seq_up.data(), seqlen});

        /* first, look for an identical match */

        auto prefix_len = seqlen;

        uint64_t hash = prefix_hashes[prefix_len];
        auto * bp = &hashtable[hash & hash_mask];

        while ((bp->size != 0U) and
               ((bp->deleted) or
                (bp->hash != hash) or
                (prefix_len != db.getsequencelen(bp->seqno_first)) or
                (seqcmp(seq_up.data(), db.getsequence(bp->seqno_first), prefix_len) != 0)))
          {
            ++bp;
            if (bp > &hashtable.back())
              {
                bp = hashtable.data();
              }
          }

        /* at this point, bp points either to (1) a free empty hash bucket, or
           (2) a bucket with an exact match. */

        auto const orig_hash = hash;
        auto * orig_bp = bp;

        if (bp->size != 0U)
          {
            /* exact match */
            bp->size += static_cast<unsigned int>(abundance);
            auto const last = bp->seqno_last;
            nextseqtab[last] = static_cast<unsigned int>(i);
            bp->seqno_last = static_cast<unsigned int>(i);

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
                        (prefix_len != db.getsequencelen(bp->seqno_first)) or
                        (seqcmp(seq_up.data(),
                                db.getsequence(bp->seqno_first),
                                prefix_len) != 0)))
                  {
                    ++bp;
                    if (bp > &hashtable.back())
                      {
                        bp = hashtable.data();
                      }
                  }
              }

            if (bp->size != 0U)
              {
                /* prefix match */

                /* get necessary info, then delete prefix from hash */
                auto const first = bp->seqno_first;
                auto const last = bp->seqno_last;
                auto const size = bp->size;
                bp->deleted = true;

                /* create new hash entry */
                bp = orig_bp;
                bp->size = static_cast<unsigned int>(size + abundance);
                bp->hash = orig_hash;
                bp->seqno_first = static_cast<unsigned int>(i);
                nextseqtab[static_cast<std::vector<unsigned int>::size_type>(i)] = first;
                bp->seqno_last = last;

                maxsize = std::max<uint64_t>(bp->size, maxsize);
              }
            else
              {
                /* no match */
                orig_bp->size = static_cast<unsigned int>(abundance);
                orig_bp->hash = orig_hash;
                orig_bp->seqno_first = static_cast<unsigned int>(i);
                orig_bp->seqno_last = static_cast<unsigned int>(i);

                maxsize = std::max(abundance, maxsize);
                ++clusters;
              }
          }

        progress.update(static_cast<uint64_t>(i));
      }
  }

  {
    Progress const progress("Sorting", 1, parameters);
    derep_sort_db = &db;
    std::qsort(hashtable.data(), static_cast<size_t>(hashtablesize), sizeof(struct bucket), derep_compare_prefix);
  }

  if (clusters > 0)
    {
      if ((clusters % 2) != 0)
        {
          median = hashtable[static_cast<std::vector<struct bucket>::size_type>((clusters - 1) / 2)].size;
        }
      else
        {
          median = (hashtable[static_cast<std::vector<struct bucket>::size_type>((clusters / 2) - 1)].size +
                    hashtable[static_cast<std::vector<struct bucket>::size_type>(clusters / 2)].size) / 2.0;
        }
    }

  average = 1.0 * static_cast<double>(sumsize) / static_cast<double>(clusters);

  if (clusters < 1)
    {
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr,
                  "0 unique sequences\n");
        }
      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "0 unique sequences\n\n");
        }
    }
  else
    {
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr,
                  "%" PRId64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n",
                  clusters, average, median, maxsize);
        }
      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "%" PRId64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n\n",
                  clusters, average, median, maxsize);
        }
    }

  /* count selected */

  int64_t selected = 0;
  for (int64_t i = 0; i < clusters; i++)
    {
      int64_t const size = hashtable[static_cast<std::vector<struct bucket>::size_type>(i)].size;
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

      int64_t relabel_count = 0;
      {
        Progress progress("Writing output file", static_cast<uint64_t>(clusters), parameters);
        for (int64_t i = 0; i < clusters; i++)
          {
            auto const & bp = hashtable[static_cast<std::vector<struct bucket>::size_type>(i)];
            int64_t const size = bp.size;
            if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
              {
                ++relabel_count;
                fasta_print_general(fp_output,
                                    nullptr,
                                    db.getsequence(bp.seqno_first),
                                    static_cast<int>(db.getsequencelen(bp.seqno_first)),
                                    db.getheader(bp.seqno_first),
                                    static_cast<int>(db.getheaderlen(bp.seqno_first)),
                                    static_cast<uint64_t>(size),
                                    relabel_count,
                                    -1.0,
                                    -1, -1, nullptr, 0.0,
                                    0,
                                    parameters);
                if (relabel_count == parameters.opt_topn)
                  {
                    break;
                  }
              }
            progress.update(static_cast<uint64_t>(i));
          }
      }

      output_handle.reset();
    }

  if (parameters.opt_uc != nullptr)
    {
      {
        Progress progress("Writing uc file, first part", static_cast<uint64_t>(clusters), parameters);
        for (int64_t i = 0; i < clusters; i++)
          {
            auto const & bp = hashtable[static_cast<std::vector<struct bucket>::size_type>(i)];
            auto const * h =  db.getheader(bp.seqno_first);
            int64_t const len = static_cast<int64_t>(db.getsequencelen(bp.seqno_first));

            std::fprintf(fp_uc, "S\t%" PRId64 "\t%" PRId64 "\t*\t*\t*\t*\t*\t%s\t*\n",
                    i, len, h);

            for (auto next = nextseqtab[bp.seqno_first];
                 next != terminal;
                 next = nextseqtab[next])
              {
                std::fprintf(fp_uc,
                        "H\t%" PRId64 "\t%" PRIu64 "\t%.1f\t+\t0\t0\t*\t%s\t%s\n",
                        i, db.getsequencelen(next), 100.0, db.getheader(next), h);
              }

            progress.update(static_cast<uint64_t>(i));
          }
      }

      {
        Progress progress("Writing uc file, second part", static_cast<uint64_t>(clusters), parameters);
        for (int64_t i = 0; i < clusters; i++)
          {
            auto const & bp = hashtable[static_cast<std::vector<struct bucket>::size_type>(i)];
            std::fprintf(fp_uc, "C\t%" PRId64 "\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
                    i, bp.size, db.getheader(bp.seqno_first));
            progress.update(static_cast<uint64_t>(i));
          }
        uc_handle.reset();
      }
    }

  if (selected < clusters)
    {
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr,
                  "%" PRId64 " uniques written, %" PRId64
                  " clusters discarded (%.1f%%)\n",
                  selected, clusters - selected,
                  100.0 * static_cast<double>(clusters - selected) / static_cast<double>(clusters));
        }

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "%" PRId64 " uniques written, %" PRId64
                  " clusters discarded (%.1f%%)\n\n",
                  selected, clusters - selected,
                  100.0 * static_cast<double>(clusters - selected) / static_cast<double>(clusters));
        }
    }

  db.clear();
}
