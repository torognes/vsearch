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

#include "utils/view.hpp"
#include "vsearch.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/db.hpp"
#include "core/derep_stats.hpp"  // Derep_stats, report_*
#include "core/fasta.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/hash_table_size.hpp"  // table_size_two_thirds
#include "utils/median.hpp"  // median_of_descending
#include "utils/open_file.hpp"
#include "utils/seqcmp.hpp"
#include "utils/span.hpp"
#include "utils/string_normalize.hpp"
#include <algorithm>  // std::max, std::sort, std::transform
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
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
    /* total abundance of the cluster, as wide as the ;size= annotation it comes
       from (header_get_size returns int64_t); see the same field in
       core/derep.cpp for what a 32-bit field did to a large annotation */
    uint64_t size = 0;
    unsigned int count = 0;
    bool deleted = false;
  };


  /* the empty-bucket sentinel, distinct from 'deleted': a slot that was
     superseded by a longer prefix is still claimed, and the probe has to step
     over it rather than stop. */
  auto is_occupied(struct bucket const & entry) noexcept -> bool
  {
    return entry.size != 0U;
  }
}


// refactoring: FNV-1A is the hashing function used in std::hash
// (designed for fast hash-table and checksum, not crypto). The
// function below might be redundant with std?

namespace {
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
}  // anonymous namespace


auto derep_prefix(struct Parameters const & parameters) -> void
{
  if (parameters.opt_strand)
    {
      fatal("Option '--strand both' not supported with --derep_prefix");
    }

  /* same output requirement as derep_fulllength and derep_id (the manual
     lists --output as mandatory); without it a run with no output option
     was a silent no-op that exited 0 */
  if ((parameters.opt_output == nullptr) and (parameters.opt_uc == nullptr)) {
    fatal("Output file for dereplication must be specified with --output "
          "or --uc");
  }

  auto output_handle = open_optional_output_file(parameters.opt_output, OutputOption{"--output"});
  auto uc_handle = open_optional_output_file(parameters.opt_uc, OutputOption{"--uc"});
  std::FILE * const fp_uc = uc_handle.get();

  Database db;
  db.read(parameters.opt_derep_prefix, 0, parameters);

  db.sortbylength_shortest_first(parameters);

  // memory-intensive: the entire database is now held in memory

  auto const dbsequencecount = static_cast<int64_t>(db.getsequencecount());

  /* adjust size of hash table for 2/3 fill rate */

  auto const hashtablesize =
    vsearch::table_size_two_thirds(static_cast<uint64_t>(dbsequencecount));
  auto const hash_mask = hashtablesize - 1;

  std::vector<struct bucket> hashtable(static_cast<std::vector<struct bucket>::size_type>(hashtablesize));

  /* only clusters, sumsize and maxsize apply here: this engine dereplicates a
     database already read into memory, so it has no input-scanning phase to
     report on and no length filtering of its own. The other six members keep
     their defaults and are never read. */
  Derep_stats stats;
  double median = 0.0;
  double average = 0.0;

  /* alloc and init table of links to other sequences in cluster */

  constexpr auto terminal = std::numeric_limits<unsigned int>::max();
  std::vector<unsigned int> nextseqtab(static_cast<std::vector<unsigned int>::size_type>(dbsequencecount), terminal);

  std::vector<char> seq_up(db.getlongestsequence());

  /* make table of hash values of prefixes */

  auto const len_longest = static_cast<unsigned int>(db.getlongestsequence());
  auto const len_shortest = static_cast<unsigned int>(db.getshortestsequence());
  std::vector<uint64_t> prefix_hashes(len_longest + 1);

  {
    Progress progress("Dereplicating", static_cast<uint64_t>(dbsequencecount), parameters);
    for (int64_t i = 0; i < dbsequencecount; i++)
      {
        auto const sequence = db.sequence_view(static_cast<uint64_t>(i));
        auto const seqlen = static_cast<unsigned int>(sequence.size());

        /* normalize sequence: uppercase and replace U by T  */
        normalize_into(seq_up, sequence);

        auto const abundance = parameters.opt_sizein ? db.getabundance(static_cast<uint64_t>(i)) : uint64_t{1};
        stats.sumsize += static_cast<int64_t>(abundance);

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

        compute_hashes_of_all_prefixes(prefix_hashes, make_span(seq_up).first(seqlen));

        /* first, look for an identical match */

        auto prefix_len = seqlen;

        uint64_t hash = prefix_hashes[prefix_len];
        auto * bp = &hashtable[hash & hash_mask];

        while (is_occupied(*bp) and
               (bp->deleted or
                (bp->hash != hash) or
                (prefix_len != db.getsequencelen(bp->seqno_first)) or
                (seqcmp(make_view(seq_up).first(static_cast<std::size_t>(prefix_len)),
                        db.sequence_view(bp->seqno_first).first(static_cast<std::size_t>(prefix_len))) != 0)))
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

        if (is_occupied(*bp))
          {
            /* exact match */
            bp->size += abundance;
            auto const last = bp->seqno_last;
            nextseqtab[last] = static_cast<unsigned int>(i);
            bp->seqno_last = static_cast<unsigned int>(i);

            stats.maxsize = std::max<uint64_t>(bp->size, stats.maxsize);
          }
        else
          {
            /* look for prefix match */

            while ((not is_occupied(*bp)) and (prefix_len > len_shortest))
              {
                --prefix_len;
                hash = prefix_hashes[prefix_len];
                bp = &hashtable[hash & hash_mask];

                while (is_occupied(*bp) and
                       (bp->deleted or
                        (bp->hash != hash) or
                        (prefix_len != db.getsequencelen(bp->seqno_first)) or
                        (seqcmp(make_view(seq_up).first(static_cast<std::size_t>(prefix_len)),
                                db.sequence_view(bp->seqno_first).first(static_cast<std::size_t>(prefix_len))) != 0)))
                  {
                    ++bp;
                    if (bp > &hashtable.back())
                      {
                        bp = hashtable.data();
                      }
                  }
              }

            if (is_occupied(*bp))
              {
                /* prefix match */

                /* get necessary info, then delete prefix from hash */
                auto const first = bp->seqno_first;
                auto const last = bp->seqno_last;
                auto const size = bp->size;
                bp->deleted = true;

                /* create new hash entry */
                bp = orig_bp;
                bp->size = size + abundance;
                bp->hash = orig_hash;
                bp->seqno_first = static_cast<unsigned int>(i);
                nextseqtab[static_cast<std::vector<unsigned int>::size_type>(i)] = first;
                bp->seqno_last = last;

                stats.maxsize = std::max<uint64_t>(bp->size, stats.maxsize);
              }
            else
              {
                /* no match */
                orig_bp->size = abundance;
                orig_bp->hash = orig_hash;
                orig_bp->seqno_first = static_cast<unsigned int>(i);
                orig_bp->seqno_last = static_cast<unsigned int>(i);

                stats.maxsize = std::max(abundance, stats.maxsize);
                ++stats.clusters;
              }
          }

        progress.update(static_cast<uint64_t>(i));
      }
  }

  {
    Progress const progress("Sorting", 1, parameters);

    /* deleted(?) first, then by highest abundance, then by label, otherwise keep order */
    auto const compare_prefix = [&db](struct bucket const & lhs, struct bucket const & rhs) -> bool
    {
      if (static_cast<int>(lhs.deleted) != static_cast<int>(rhs.deleted))
        {
          return static_cast<int>(lhs.deleted) < static_cast<int>(rhs.deleted);
        }

      // both are deleted, compare abundances
      if (lhs.size != rhs.size)
        {
          return lhs.size > rhs.size;
        }

      // both are deleted, same abundances, compare sequence headers
      auto const result = db.header_view(lhs.seqno_first)
                            .compare(db.header_view(rhs.seqno_first));
      if (result != 0)
        {
          return result < 0;
        }

      // both are deleted, same abundances, same sequence headers, compare input order
      return lhs.seqno_first < rhs.seqno_first;
    };

    /* skip when the database is empty: the lone empty bucket has no valid
       header for the comparator to dereference */
    if (dbsequencecount > 0)
      {
        std::sort(hashtable.begin(), hashtable.end(), compare_prefix);
      }
  }

  /* the live clusters are the leading 'clusters' entries: compare_prefix puts
     the not-deleted buckets first and orders them by decreasing size, and
     'stats.clusters' counts exactly those (it is incremented only when no prefix
     match was found, never on the path that marks a bucket deleted) */
  median = median_of_descending(make_view(hashtable).first(stats.clusters),
                                [](struct bucket const & entry) { return entry.size; });

  average = 1.0 * static_cast<double>(stats.sumsize) / static_cast<double>(stats.clusters);

  report_unique_summary(stats, average, median, parameters);

  /* count selected */

  uint64_t selected = 0;
  for (uint64_t i = 0; i < stats.clusters; i++)
    {
      auto const size = static_cast<int64_t>(hashtable[static_cast<std::vector<struct bucket>::size_type>(i)].size);
      if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
        {
          ++selected;
          if (selected == static_cast<uint64_t>(parameters.opt_topn))
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
        Progress progress("Writing output file", stats.clusters, parameters);
        for (uint64_t i = 0; i < stats.clusters; i++)
          {
            auto const & bp = hashtable[static_cast<std::vector<struct bucket>::size_type>(i)];
            auto const size = static_cast<int64_t>(bp.size);
            if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
              {
                ++relabel_count;
                fasta_print_general(output_handle.get(),
                                    nullptr,
                                    db.record(bp.seqno_first),
                                    OutputAnnotations{static_cast<uint64_t>(size), relabel_count},
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
        Progress progress("Writing uc file, first part", stats.clusters, parameters);
        for (uint64_t i = 0; i < stats.clusters; i++)
          {
            auto const & bp = hashtable[static_cast<std::vector<struct bucket>::size_type>(i)];
            auto const header = db.header_view(bp.seqno_first);
            auto const len = static_cast<int64_t>(db.getsequencelen(bp.seqno_first));

            fprint(fp_uc, "S\t");
            fprint_integer(fp_uc, i);
            fprint(fp_uc, '\t');
            fprint_integer(fp_uc, len);
            fprint(fp_uc, "\t*\t*\t*\t*\t*\t");
            fprint(fp_uc, header);
            fprint(fp_uc, "\t*\n");

            for (auto next = nextseqtab[bp.seqno_first];
                 next != terminal;
                 next = nextseqtab[next])
              {
                fprint(fp_uc, "H\t");
                fprint_integer(fp_uc, i);
                fprint(fp_uc, '\t');
                fprint_integer(fp_uc, db.getsequencelen(next));
                fprint(fp_uc, '\t');
                std::fprintf(fp_uc, "%.1f", 100.0);
                fprint(fp_uc, "\t+\t0\t0\t*\t");
                fprint(fp_uc, db.header_view(next));
                fprint(fp_uc, '\t');
                fprint(fp_uc, header);
                fprint(fp_uc, '\n');
              }

            progress.update(static_cast<uint64_t>(i));
          }
      }

      {
        Progress progress("Writing uc file, second part", stats.clusters, parameters);
        for (uint64_t i = 0; i < stats.clusters; i++)
          {
            auto const & bp = hashtable[static_cast<std::vector<struct bucket>::size_type>(i)];
            fprint(fp_uc, "C\t");
            fprint_integer(fp_uc, i);
            fprint(fp_uc, '\t');
            fprint_integer(fp_uc, bp.size);
            fprint(fp_uc, "\t*\t*\t*\t*\t*\t");
            fprint(fp_uc, db.header_view(bp.seqno_first));
            fprint(fp_uc, "\t*\n");
            progress.update(static_cast<uint64_t>(i));
          }
        uc_handle.reset();
      }
    }

  report_selected(selected, stats, parameters);

  db.clear();
}
