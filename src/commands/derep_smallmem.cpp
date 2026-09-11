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
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/derep_stats.hpp"  // Derep_stats, report_*
#include "core/fasta.hpp"
#include "core/fastx.hpp"
#include "utils/progress.hpp"
#include "vendored/city.h"
#include "utils/base_mapping.hpp"
#include "utils/fatal.hpp"
#include "utils/grow_to_fit.hpp"  // vsearch::grow_to_fit
#include "utils/open_file.hpp"
#include "utils/print_view.hpp"  // fprint
// #include "util.h"  // hash_cityhash128
#include "utils/cityhash.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/string_normalize.hpp"
#include <algorithm>  // std::min, std::max
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf
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


namespace {

/* the empty-bucket sentinel. This table stores no sequence at all (a deliberate
   memory tradeoff, see the probe below), so 'size' is the only field that can
   say whether a slot is claimed. */
auto is_occupied(struct sm_bucket const & entry) noexcept -> bool
{
  return entry.size != 0U;
}


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

      for (auto const & bucket : hashtable)
        {
          auto const v = bucket.size;
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


inline auto hash2bucket(uint128 const hash, uint64_t const htsize) -> uint64_t
{
  // extract hash's first uint64_t, cast to the size of the hash table
  return hash.first % htsize;
}


inline auto next_bucket(uint64_t const prev_bucket, uint64_t const htsize) -> uint64_t
{
  return (prev_bucket + 1) % htsize;
}


/*
  The bucket holding this hash, or the free bucket where it belongs: linear
  probing from hash2bucket(), stepping with next_bucket().

  Sequences are matched by their 128-bit CityHash alone -- there is no
  byte-wise comparison here (unlike derep_fulllength and derep_prefix), a
  deliberate memory tradeoff. A 128-bit hash collision would merge two
  distinct sequences, but the probability only approaches 50% near
  2^64 (~1.8e19) sequences.
*/
inline auto find_bucket(std::vector<struct sm_bucket> & hashtable,
                        uint128 const hash) -> struct sm_bucket &
{
  auto const htsize = hashtable.size();
  auto index = hash2bucket(hash, htsize);
  /* the one raw pointer left in the probe, walking in step with the index so
     that the bucket's address does not have to be recomputed on return; it
     used to be four, one per copy of this loop */
  auto * bucket = &hashtable[index];
  while (is_occupied(*bucket) and (hash != bucket->hash))
    {
      index = next_bucket(index, htsize);
      bucket = &hashtable[index];
    }
  return *bucket;
}


/* the table is grown before it passes this fill rate */
constexpr auto max_fill_rate_percent = uint64_t{95};
constexpr auto percent_scale = uint64_t{100};


/* Whether adding one more cluster would take the table past its fill-rate cap.
   Written as a cross-multiplication so that it stays integer arithmetic. */
inline auto is_too_full(uint64_t const clusters, std::size_t const table_size) -> bool
{
  return percent_scale * (clusters + 1) > max_fill_rate_percent * table_size;
}


auto rehash_smallmem(std::vector<struct sm_bucket> & hashtable) -> void
{
  /* allocate new hash table, 50% larger */
  auto const new_hashtablesize = 3 * hashtable.size() / 2;
  std::vector<struct sm_bucket> new_hashtable(new_hashtablesize);

  /* rehash all from old to new */
  for (auto const & old_bucket : hashtable)
    {
      if (is_occupied(old_bucket))
        {
          auto & new_bucket = find_bucket(new_hashtable, old_bucket.hash);
          /* the table holds one bucket per distinct hash, so find_bucket can
             only have stopped on a free one */
          assert(not is_occupied(new_bucket));
          new_bucket = old_bucket;
        }
    }

  hashtable.swap(new_hashtable);
}
}  // anonymous namespace


auto derep_smallmem(struct Parameters const & parameters) -> void
{
  /*
    dereplicate full length sequences using a small amount of memory
    output options: --fastaout
  */

  auto * input_filename = parameters.input_filename;
  auto h = fastx_open(input_filename, parameters);

  if (h->is_pipe_input())
    {
      fatal("The derep_smallmem command does not support input from a pipe.");
    }

  auto const output_handle = open_mandatory_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});

  auto const filesize = h->get_size();

  /* allocate initial memory for sequences of length up to 1023 chars */
  int64_t const initial_seqlen = 1024;

  /* allocate initial hashtable with 1024 buckets */

  static constexpr auto initial_bucket_count = std::size_t{1024};
  std::vector<struct sm_bucket> hashtable(initial_bucket_count);

  // memory-intensive: the hash table has been allocated

  std::vector<char> seq_up(static_cast<size_t>(initial_seqlen) + 1);
  std::vector<char> rc_seq_up(static_cast<size_t>(initial_seqlen) + 1);

  std::string const prompt = std::string("Dereplicating file ") + input_filename;


  Derep_stats stats;

  /* first pass */

  {
    Progress progress(prompt, filesize, parameters);
    while (h->next(not parameters.opt_notrunclabels, Mapping::none))
      {
        auto const sequence = h->sequence_view();
        auto const seqlen = static_cast<int64_t>(sequence.size());

        if (seqlen < parameters.opt_minseqlength)
          {
            ++stats.discarded_short;
            continue;
          }

        if (seqlen > parameters.opt_maxseqlength)
          {
            ++stats.discarded_long;
            continue;
          }

        stats.nucleotidecount += static_cast<uint64_t>(seqlen);
        stats.longest = std::max(seqlen, stats.longest);
        stats.shortest = std::min(seqlen, stats.shortest);

        /* check allocations */

        // memory-intensive: sequence buffers grown to fit the longest sequence
        vsearch::grow_to_fit(seq_up, static_cast<size_t>(seqlen));
        if (parameters.opt_strand)
          {
            vsearch::grow_to_fit(rc_seq_up, static_cast<size_t>(seqlen));
          }

        if (is_too_full(stats.clusters, hashtable.size()))
          {
            // keep hash table fill rate at max 95% */
            rehash_smallmem(hashtable);
            // memory-intensive: the hash table has been resized (rehash)
          }

        /* normalize sequence: uppercase and replace U by T  */
        auto const seq_up_v = normalize_into(seq_up, sequence);

        auto const hash = hash_function(seq_up_v);
        auto * bucket = &find_bucket(hashtable, hash);

        if (parameters.opt_strand and not is_occupied(*bucket))
          {
            /* no match on plus strand */
            /* check minus strand as well */

            /* the reverse complement is only ever read here, so it is only
               computed here: a record that matched a cluster on the plus
               strand does not need one, and on dereplication input most
               records do */
            reverse_complement(make_span(rc_seq_up).first(static_cast<std::size_t>(seqlen)), seq_up_v);
            auto const rc_hash = hash_function(make_view(rc_seq_up).first(static_cast<std::size_t>(seqlen)));
            auto & rc_bucket = find_bucket(hashtable, rc_hash);

            if (is_occupied(rc_bucket))
              {
                bucket = &rc_bucket;
              }
          }

        int64_t const abundance = h->get_abundance();
        int64_t const ab = parameters.opt_sizein ? abundance : 1;
        stats.sumsize += ab;

        if (is_occupied(*bucket))
          {
            /* at least one identical sequence already */
            bucket->size += static_cast<uint64_t>(ab);
          }
        else
          {
            /* no identical sequences yet */
            bucket->size = static_cast<uint64_t>(ab);
            bucket->hash = hash;
            ++stats.clusters;
          }

        stats.maxsize = std::max(bucket->size, stats.maxsize);

        ++stats.sequencecount;
        progress.update(h->get_position());
      }
  }
  h->report_stripped_warning(parameters);

  report_input_stats(stats, parameters);

  report_length_filtered(parameters, "minseqlength", parameters.opt_minseqlength, stats.discarded_short);
  report_length_filtered(parameters, "maxseqlength", parameters.opt_maxseqlength, stats.discarded_long);


  {
    /* both are unused when there is no cluster to report, and the median
       costs a full pass over the hash table, so neither is computed then */
    auto average = 0.0;
    auto median = 0.0;
    if (stats.clusters >= 1)
      {
        average = static_cast<double>(stats.sumsize) / static_cast<double>(stats.clusters);
        median = find_median(hashtable);
      }
    report_unique_summary(stats, average, median, parameters);
  }

  /* second pass with output */

  auto h2 = fastx_open(input_filename, parameters);


  uint64_t selected = 0;

  {
    Progress progress("Writing FASTA output file", filesize, parameters);
    while (h2->next(not parameters.opt_notrunclabels, Mapping::none))
      {
        auto const sequence = h2->sequence_view();
        auto const seqlen = static_cast<int64_t>(sequence.size());

        if ((seqlen < parameters.opt_minseqlength) or (seqlen > parameters.opt_maxseqlength))
          {
            continue;
          }

        /* normalize sequence: uppercase and replace U by T  */
        auto const seq_up_v = normalize_into(seq_up, sequence);

        auto const hash = hash_function(seq_up_v);
        auto * bucket = &find_bucket(hashtable, hash);

        if (parameters.opt_strand and not is_occupied(*bucket))
          {
            /* no match on plus strand */
            /* check minus strand as well */

            /* the reverse complement is only ever read here, so it is only
               computed here: a record that matched a cluster on the plus
               strand does not need one, and on dereplication input most
               records do */
            reverse_complement(make_span(rc_seq_up).first(static_cast<std::size_t>(seqlen)), seq_up_v);
            auto const rc_hash = hash_function(make_view(rc_seq_up).first(static_cast<std::size_t>(seqlen)));
            auto & rc_bucket = find_bucket(hashtable, rc_hash);

            if (is_occupied(rc_bucket))
              {
                bucket = &rc_bucket;
              }
          }

        auto const size = static_cast<int64_t>(bucket->size);

        if (size > 0)
          {
            /* print sequence */

            if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
              {
                ++selected;
                fasta_print_general(output_handle.get(),
                                    sequence,
                                    h2->header_view(),
                                    OutputAnnotations{static_cast<uint64_t>(size),
                                                      static_cast<int64_t>(selected)},
                                    parameters);
              }
            bucket->size = static_cast<uint64_t>(-1);
          }

        progress.update(h2->get_position());
      }
  }
  h2->report_stripped_warning(parameters);

  report_selected(selected, stats, parameters);
}
