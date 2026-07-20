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

#include "utils/span.hpp"
#include "utils/view.hpp"
#include "vsearch.hpp"
#include "core/fasta.hpp"
#include "core/fastx.hpp"
#include "utils/progress.hpp"
#include "vendored/city.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
// #include "util.h"  // hash_cityhash128
#include "utils/cityhash.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/string_normalize.hpp"
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
static constexpr Hash hash_function = hash_cityhash128;


struct sm_bucket
{
  uint128 hash;
  uint64_t size;
};


auto find_median(std::vector<struct sm_bucket> const & hashtable) -> double
{
  /* find the median size, based on an iterative search starting at e.g. 1 */

  uint64_t cand = 1;    /* candidate for the median */
  uint64_t below = 0;   /* closest value below the candidate */
  uint64_t above = 0;   /* closest value above the candidate */

  while (true)
    {
      uint64_t cand_count = 0;  /* number of clusters with same size as cand */
      uint64_t below_count = 0; /* number of clusters with smaller size than cand */
      uint64_t above_count = 0; /* number of clusters with larger size than cand */

      for (uint64_t i = 0; i < hashtable.size(); i++)
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
                  return static_cast<double>(cand + above) / 2.0;
                }
              if (above_count + cand_count == below_count)
                // mid == below_count
                // same as:
                // (below_count + cand_count + above_count) / 2 == below_count
                // which simplifies into:
                // above_count + cand_count == below_count
                {
                  return static_cast<double>(below + cand) / 2.0;  // cannot reach?
                }
              return static_cast<double>(cand);
            }
          cand = above;
        }
      else
        {
          cand = below;  // cannot reach?
        }
    }
}


inline auto hash2bucket(uint128 hash, uint64_t const htsize) -> uint64_t
{
  // extract hash's first uint64_t, cast to the size of the hash table
  return hash.first % htsize;
}


inline auto next_bucket(uint64_t const prev_bucket, uint64_t const htsize) -> uint64_t
{
  return (prev_bucket + 1) % htsize;
}


auto rehash_smallmem(std::vector<struct sm_bucket> & hashtable) -> void
{
  /* allocate new hash table, 50% larger */
  auto const new_hashtablesize = 3 * hashtable.size() / 2;
  std::vector<struct sm_bucket> new_hashtable(new_hashtablesize);

  /* rehash all from old to new */
  for (auto const & old_bucket : hashtable)
    {
      if (old_bucket.size != 0U)
        {
          auto k = hash2bucket(old_bucket.hash, new_hashtablesize);
          while (new_hashtable[k].size != 0U)
            {
              k = next_bucket(k, new_hashtablesize);
            }
          new_hashtable[k] = old_bucket;
        }
    }

  hashtable.swap(new_hashtable);
}


auto derep_smallmem(struct Parameters const & parameters) -> void
{
  /*
    dereplicate full length sequences using a small amount of memory
    output options: --fastaout
  */

  auto * input_filename = parameters.opt_derep_smallmem;
  auto h = fastx_open(input_filename, parameters);

  if (h->is_pipe_input())
    {
      fatal("The derep_smallmem command does not support input from a pipe.");
    }

  auto const output_handle = open_mandatory_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  std::FILE * const fp_fastaout = output_handle.get();

  auto const filesize = h->get_size();

  /* allocate initial memory for sequences of length up to 1023 chars */
  int64_t alloc_seqlen = 1024;

  /* allocate initial hashtable with 1024 buckets */

  std::vector<struct sm_bucket> hashtable(1024);

  // memory-intensive: the hash table has been allocated

  std::vector<char> seq_up(static_cast<size_t>(alloc_seqlen) + 1);
  std::vector<char> rc_seq_up(static_cast<size_t>(alloc_seqlen) + 1);

  std::string const prompt = std::string("Dereplicating file ") + input_filename;


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

  {
    Progress progress(prompt.c_str(), filesize, parameters);
    while (h->next(not parameters.opt_notrunclabels, chrmap_no_change()))
      {
        int64_t const seqlen = static_cast<int64_t>(h->get_sequence_length());

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

        nucleotidecount += static_cast<uint64_t>(seqlen);
        longest = std::max(seqlen, longest);
        shortest = std::min(seqlen, shortest);

        /* check allocations */

        if (seqlen > alloc_seqlen)
          {
            alloc_seqlen = seqlen;
            seq_up.resize(static_cast<size_t>(alloc_seqlen) + 1);
            rc_seq_up.resize(static_cast<size_t>(alloc_seqlen) + 1);

            // memory-intensive: sequence buffers grown to fit the longest sequence
          }

        if (100 * (clusters + 1) > 95 * hashtable.size())
          {
            // keep hash table fill rate at max 95% */
            rehash_smallmem(hashtable);
            // memory-intensive: the hash table has been resized (rehash)
          }

        auto const * seq = h->get_sequence();

        /* normalize sequence: uppercase and replace U by T  */
        string_normalize(Span<char>{seq_up.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq, static_cast<std::size_t>(seqlen)});

        /* reverse complement if necessary */
        if (parameters.opt_strand)
          {
            reverse_complement(Span<char>{rc_seq_up.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq_up.data(), static_cast<std::size_t>(seqlen)});
          }

        /*
          Find a free bucket, or the bucket holding this sequence. Sequences
          are matched by their 128-bit CityHash alone — there is no byte-wise
          comparison here (unlike derep_fulllength and derep_prefix), a
          deliberate memory tradeoff. A 128-bit hash collision would merge two
          distinct sequences, but the probability only approaches 50% near
          2^64 (~1.8e19) sequences.
        */

        auto const hash = hash_function(seq_up.data(), static_cast<uint64_t>(seqlen));
        auto j =  hash2bucket(hash, hashtable.size());
        auto * bp = &hashtable[j];

        while ((bp->size != 0U) and (hash != bp->hash))
          {
            j = next_bucket(j, hashtable.size());
            bp = &hashtable[j];
          }

        if (parameters.opt_strand and (bp->size == 0U))
          {
            /* no match on plus strand */
            /* check minus strand as well */

            auto const rc_hash = hash_function(rc_seq_up.data(), static_cast<uint64_t>(seqlen));
            auto k =  hash2bucket(rc_hash, hashtable.size());
            auto * rc_bp = &hashtable[k];

            while ((rc_bp->size != 0U) and (rc_hash != rc_bp->hash))
              {
                k = next_bucket(k, hashtable.size());
                rc_bp = &hashtable[k];
              }

            if (rc_bp->size != 0U)
              {
                bp = rc_bp;
                j = k;  // cppcheck: 'j' is assigned a value that is never used
              }
          }

        int64_t const abundance = h->get_abundance();
        int64_t const ab = parameters.opt_sizein ? abundance : 1;
        sumsize += ab;

        if (bp->size != 0U)
          {
            /* at least one identical sequence already */
            bp->size += static_cast<uint64_t>(ab);
          }
        else
          {
            /* no identical sequences yet */
            bp->size = static_cast<uint64_t>(ab);
            bp->hash = hash;
            ++clusters;
          }

        maxsize = std::max(bp->size, maxsize);

        ++sequencecount;
        progress.update(h->get_position());
      }
  }
  h->report_stripped_warning(parameters);

  if (not parameters.opt_quiet)
    {
      if (sequencecount > 0)
        {
          std::fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, min %" PRId64
                  ", max %" PRId64 ", avg %.0f\n",
                  nucleotidecount,
                  sequencecount,
                  shortest,
                  longest,
                  static_cast<double>(nucleotidecount) / static_cast<double>(sequencecount));
        }
      else
        {
          std::fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  nucleotidecount,
                  sequencecount);
        }
    }

  if (parameters.opt_log != nullptr)
    {
      if (sequencecount > 0)
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, min %" PRId64
                  ", max %" PRId64 ", avg %.0f\n",
                  nucleotidecount,
                  sequencecount,
                  shortest,
                  longest,
                  static_cast<double>(nucleotidecount) / static_cast<double>(sequencecount));
        }
      else
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  nucleotidecount,
                  sequencecount);
        }
    }

  if (discarded_short != 0U)
    {
      std::fprintf(stderr,
              "minseqlength %" PRId64 ": %" PRIu64 " %s discarded.\n",
              parameters.opt_minseqlength,
              discarded_short,
              (discarded_short == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "minseqlength %" PRId64 ": %" PRIu64 " %s discarded.\n\n",
                  parameters.opt_minseqlength,
                  discarded_short,
                  (discarded_short == 1 ? "sequence" : "sequences"));
        }
    }

  if (discarded_long != 0U)
    {
      std::fprintf(stderr,
              "maxseqlength %" PRId64 ": %" PRIu64 " %s discarded.\n",
              parameters.opt_maxseqlength,
              discarded_long,
              (discarded_long == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "maxseqlength %" PRId64 ": %" PRIu64 " %s discarded.\n\n",
                  parameters.opt_maxseqlength,
                  discarded_long,
                  (discarded_long == 1 ? "sequence" : "sequences"));
        }
    }


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
      auto const average = static_cast<double>(sumsize) / static_cast<double>(clusters);
      const auto median = find_median(hashtable);
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr,
                  "%" PRIu64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n",
                  clusters, average, median, maxsize);
        }
      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n\n",
                  clusters, average, median, maxsize);
        }
    }

  /* second pass with output */

  auto h2 = fastx_open(input_filename, parameters);


  uint64_t selected = 0;

  {
    Progress progress("Writing FASTA output file", filesize, parameters);
    while (h2->next(not parameters.opt_notrunclabels, chrmap_no_change()))
      {
        int64_t const seqlen = static_cast<int64_t>(h2->get_sequence_length());

        if ((seqlen < parameters.opt_minseqlength) or (seqlen > parameters.opt_maxseqlength))
          {
            continue;
          }

        auto const * seq = h2->get_sequence();

        /* normalize sequence: uppercase and replace U by T  */
        string_normalize(Span<char>{seq_up.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq, static_cast<std::size_t>(seqlen)});

        /* reverse complement if necessary */
        if (parameters.opt_strand)
          {
            reverse_complement(Span<char>{rc_seq_up.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq_up.data(), static_cast<std::size_t>(seqlen)});
          }

        auto const hash = hash_function(seq_up.data(), static_cast<uint64_t>(seqlen));
        auto j =  hash2bucket(hash, hashtable.size());
        auto * bp = &hashtable[j];

        while ((bp->size != 0U) and (hash != bp->hash))
          {
            j = next_bucket(j, hashtable.size());
            bp = &hashtable[j];
          }

        if (parameters.opt_strand and (bp->size == 0U))
          {
            /* no match on plus strand */
            /* check minus strand as well */

            auto const rc_hash = hash_function(rc_seq_up.data(), static_cast<uint64_t>(seqlen));
            auto k =  hash2bucket(rc_hash, hashtable.size());
            auto * rc_bp = &hashtable[k];

            while ((rc_bp->size != 0U) and (rc_hash != rc_bp->hash))
              {
                k = next_bucket(k, hashtable.size());
                rc_bp = &hashtable[k];
              }

            if (rc_bp->size != 0U)
              {
                bp = rc_bp;
                j = k;  // cppcheck: 'j' is assigned a value that is never used
              }
          }

        int64_t const size = static_cast<int64_t>(bp->size);

        if (size > 0)
          {
            /* print sequence */

            auto const * header = h2->get_header();
            int const headerlen = static_cast<int>(h2->get_header_length());

            if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
              {
                ++selected;
                fasta_print_general(fp_fastaout,
                                    nullptr,
                                    seq,
                                    static_cast<int>(seqlen),
                                    header,
                                    headerlen,
                                    static_cast<uint64_t>(size),
                                    static_cast<int64_t>(selected),
                                    -1.0,
                                    -1, -1, nullptr, 0.0,
                                    0,
                                    parameters);
              }
            bp->size = static_cast<uint64_t>(-1);
          }

        progress.update(h2->get_position());
      }
  }
  h2->report_stripped_warning(parameters);

  if (selected < clusters)
    {
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr,
                  "%" PRIu64 " uniques written, %"
                  PRIu64 " clusters discarded (%.1f%%)\n",
                  selected, clusters - selected,
                  100.0 * static_cast<double>(clusters - selected) / static_cast<double>(clusters));
        }

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64 " uniques written, %"
                  PRIu64 " clusters discarded (%.1f%%)\n\n",
                  selected, clusters - selected,
                  100.0 * static_cast<double>(clusters - selected) / static_cast<double>(clusters));
        }
    }
}
