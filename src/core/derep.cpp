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
#include "utils/progress.hpp"
#include "core/derep.hpp"
#include "core/derep_internal.hpp"
#include "core/fasta.hpp"  // fasta_print_general
#include "core/fastq.hpp"  // fastq_print_general
#include "core/fastx.hpp"  // fastx_open, fastx_next, fastx_get_*
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/seqcmp.hpp"
#include "utils/cityhash.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/string_normalize.hpp"
#include <algorithm>  // std::count_if, std::min, std::sort
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::log10, std::pow
#include <cstdint> // int64_t, uint64_t
#include <cstdlib>  // std::ldiv
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstring>  // std::strcmp
#include <limits>
#include <memory>  // std::unique_ptr
#include <string>
#include <utility>  // std::move
#include <vector>


// refactoring:
// replace with std::unordered_map (default hashing)
// if performance are bad, see Victor_Ciura's Cpp Talk "So You Think You Can Hash"
// then make a CityHash hasher object and use it with std::unordered_map
using Hash = decltype(&hash_cityhash64);
static constexpr Hash hash_function = hash_cityhash64;


// anonymous namespace: 'bucket' is a file-local type; derep_prefix.cc
// defines a different struct of the same name, so internal linkage
// here avoids a one-definition-rule violation across translation units
namespace {
  struct bucket
  {
    uint64_t hash = 0;
    unsigned int seqno_first = 0;
    unsigned int seqno_last = 0;
    unsigned int size = 0;
    unsigned int count = 0;
    unsigned int seqlen = 0;  /* sequence length (used by API to avoid strlen) */
    bool deleted = false;
    std::string header;
    std::string seq;
    std::string qual;  /* empty when FASTA (no quality) */
  };
}


static auto rehash(std::vector<struct bucket> & hashtable) -> void
{
  // new double-size hash table
  uint64_t const new_hashtable_size = 2 * hashtable.size();
  uint64_t const new_hash_mask = new_hashtable_size - 1;
  std::vector<struct bucket> new_hashtable(new_hashtable_size);

  // rehash all entries from the old to the new table (move, not copy, so the
  // std::string members are not deep-copied on every rehash)
  for (auto & old_bucket : hashtable) {
    if (old_bucket.size != 0U) {
      auto new_index = old_bucket.hash & new_hash_mask;
      while (new_hashtable[new_index].size != 0U) {
        new_index = (new_index + 1) & new_hash_mask;
      }
      new_hashtable[new_index] = std::move(old_bucket);
    }
  }
  hashtable.swap(new_hashtable);
}


// anonymous namespace: limit visibility and usage to this translation unit
namespace {


  constexpr auto terminal = std::numeric_limits<unsigned int>::max();


  auto count_selected(std::vector<struct bucket> const & hashtable,
                      struct Parameters const & parameters) -> uint64_t {
    auto size_in_range = [&](struct bucket const & bucket) -> bool {
      auto const size = bucket.size;
      return ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize));
    };
    auto const selected = std::count_if(hashtable.begin(), hashtable.end(),
                                        size_in_range);
    return std::min(static_cast<uint64_t>(selected),
                    static_cast<uint64_t>(parameters.opt_topn));
  }


  // refactoring: same as find_median_abundance()
  auto find_median_size(std::vector<struct bucket> const & hashtable, uint64_t const num_used) -> double {
    static constexpr auto half = 0.5;
    if (num_used == 0) {
      return 0.0;
    }
    auto const midpoint = std::ldiv(static_cast<long>(num_used), 2L);
    auto const is_odd = ((num_used % 2) != 0U);
    if (is_odd) {
      // index is zero-based, so if size == 3, midpoint == 1
      auto const index = static_cast<std::size_t>(midpoint.quot);
      return hashtable[index].size;
    }
    // pair number of elements:
    // index is zero-based, so if size == 4, midpoint == 2, lhs index == 1
    auto const lhs_index = static_cast<std::size_t>(midpoint.quot - 1);
    auto const rhs_index = static_cast<std::size_t>(midpoint.quot);
    auto const lhs_size = hashtable[lhs_index].size;
    auto const rhs_size = hashtable[rhs_index].size;
    // sorted by decreasing abundance: lhs size > rhs size
    // limit risk of integer additon overflow:
    // a >= b ; (a + b) / 2 == b + (a - b) / 2
    return rhs_size + ((lhs_size - rhs_size) * half);
  }


  // refactorig: duplicate of q2p()?
  inline auto convert_quality_symbol_to_probability(int const quality_symbol, struct Parameters const & parameters) -> double
  {
    static constexpr auto minimal_quality_value = 2;
    static constexpr auto maximal_probability = 0.75;
    auto const quality_value = quality_symbol - static_cast<int>(parameters.opt_fastq_ascii);
    if (quality_value < minimal_quality_value)
      {
        return maximal_probability;
      }
    static constexpr auto base = 10.0;
    return std::pow(base, -quality_value / base);
  }


  inline auto convert_probability_to_quality_symbol(double const probability, struct Parameters const & parameters) -> int
  {
    static constexpr auto base = 10.0;
    auto quality_value = static_cast<int64_t>(std::trunc(-base * std::log10(probability)));
    quality_value = std::min(quality_value, parameters.opt_fastq_qmaxout);
    quality_value = std::max(quality_value, parameters.opt_fastq_qminout);
    return static_cast<int>(quality_value + parameters.opt_fastq_asciiout);
  }

}  // end of anonymous namespace


// strict-weak-ordering predicate for std::sort: true iff lhs sorts before rhs.
// A faithful translation of the former derep_compare_full C-comparator
// (positive -> lhs after -> false; negative -> lhs before -> true; zero ->
// equivalent -> fall through). static because it takes the anonymous-namespace
// bucket type by reference, like rehash().
static auto derep_bucket_before(struct bucket const & lhs, struct bucket const & rhs) -> bool
{
  /* highest abundance first, then by label, otherwise keep order */

  if (lhs.deleted and not rhs.deleted)  // refactoring: deleted is always set to false for derep_fulllength
    {
      return false;
    }
  if (not lhs.deleted and rhs.deleted)  // refactoring: deleted is always set to false for derep_fulllength
    {
      return true;
    }
  // same status
  if (lhs.size < rhs.size)
    {
      return false;
    }
  if (lhs.size > rhs.size)
    {
      return true;
    }
  // same abundance
  if (lhs.size == 0)
    {
      return false;
    }
  auto const result = std::strcmp(lhs.header.c_str(), rhs.header.c_str());
  if (result != 0)
    {
      return result < 0;
    }
  // same header (label)
  if (lhs.seqno_first < rhs.seqno_first)
    {
      return true;
    }
  if (lhs.seqno_first > rhs.seqno_first)
    {
      return false;
    }
  // same ordinal value (impossible)
  return false;  // unreachable
}


// output writers: one per format, each iterating the sorted clusters. Placed
// in an anonymous namespace (internal linkage) like the other file-local
// helpers, and called through output_results() below.
namespace {

  auto write_fasta_output(std::FILE * fp_fastaout,
                          std::vector<struct bucket> const & hashtable,
                          uint64_t const clusters,
                          struct Parameters const & parameters) -> void
  {
    {
      Progress progress("Writing FASTA output file", clusters, parameters);
      int64_t relabel_count = 0;
      for (uint64_t i = 0; i < clusters; ++i)
        {
          auto const & cluster = hashtable[i];
          int64_t const size = cluster.size;
          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++relabel_count;
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  cluster.seq.c_str(),
                                  static_cast<int>(cluster.seq.size()),
                                  cluster.header.c_str(),
                                  static_cast<int>(cluster.header.size()),
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
          progress.update(i);
        }
    }
  }


  auto write_fastq_output(std::FILE * fp_fastqout,
                          std::vector<struct bucket> const & hashtable,
                          uint64_t const clusters,
                          struct Parameters const & parameters) -> void
  {
    {
      Progress progress("Writing FASTQ output file", clusters, parameters);
      int64_t relabel_count = 0;
      for (uint64_t i = 0; i < clusters; ++i)
        {
          auto const & cluster = hashtable[i];
          int64_t const size = cluster.size;
          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++relabel_count;
              fastq_print_general(fp_fastqout,
                                  cluster.seq.c_str(),
                                  static_cast<int>(cluster.seq.size()),
                                  cluster.header.c_str(),
                                  static_cast<int>(cluster.header.size()),
                                  cluster.qual.c_str(),
                                  static_cast<uint64_t>(size),
                                  relabel_count,
                                  -1.0,
                                  parameters);
              if (relabel_count == parameters.opt_topn)
                {
                  break;
                }
            }
          progress.update(i);
        }
    }
  }


  auto write_uc_output(std::FILE * fp_uc,
                       std::vector<struct bucket> const & hashtable,
                       uint64_t const clusters,
                       std::vector<unsigned int> const & nextseqtab,
                       std::vector<std::string> const & headertab,
                       std::vector<char> const & match_strand,
                       struct Parameters const & parameters) -> void
  {
    {
      Progress progress("Writing uc file, first part", clusters, parameters);
      for (uint64_t i = 0; i < clusters; ++i)
        {
          auto const & cluster = hashtable[i];
          int64_t const len = static_cast<int64_t>(cluster.seq.size());

          std::fprintf(fp_uc, "S\t%" PRIu64 "\t%" PRId64 "\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, len, cluster.header.c_str());

          for (auto next = nextseqtab[cluster.seqno_first];
               next != terminal;
               next = nextseqtab[next])
            {
              std::fprintf(fp_uc,
                      "H\t%" PRIu64 "\t%" PRId64 "\t%.1f\t%s\t0\t0\t*\t%s\t%s\n",
                      i, len, 100.0,
                      ((match_strand[next] != 0) ? "-" : "+"),
                      headertab[next].c_str(), cluster.header.c_str());
            }

          progress.update(i);
        }
    }

    {
      Progress progress("Writing uc file, second part", clusters, parameters);
      for (uint64_t i = 0; i < clusters; ++i)
        {
          auto const & cluster = hashtable[i];
          std::fprintf(fp_uc, "C\t%" PRIu64 "\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
                  i, cluster.size, cluster.header.c_str());
          progress.update(i);
        }
    }
  }


  auto write_tabbedout_output(std::FILE * fp_tabbedout,
                              std::vector<struct bucket> const & hashtable,
                              uint64_t const clusters,
                              std::vector<unsigned int> const & nextseqtab,
                              std::vector<std::string> const & headertab,
                              struct Parameters const & parameters) -> void
  {
    {
      Progress progress("Writing tab separated file", clusters, parameters);
      for (uint64_t i = 0; i < clusters; ++i)
        {
          auto const & cluster = hashtable[i];

          if (parameters.opt_relabel != nullptr) {
            std::fprintf(fp_tabbedout,
                    "%s\t%s%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                    cluster.header.c_str(), parameters.opt_relabel, i + 1, i, static_cast<uint64_t>(0), cluster.count, cluster.header.c_str());
          } else {
            std::fprintf(fp_tabbedout, "%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                    cluster.header.c_str(), cluster.header.c_str(), i, static_cast<uint64_t>(0), cluster.count, cluster.header.c_str());
          }

          uint64_t j = 1;
          for (auto next = nextseqtab[cluster.seqno_first];
               next != terminal;
               next = nextseqtab[next])
            {
              if (parameters.opt_relabel != nullptr) {
                std::fprintf(fp_tabbedout,
                        "%s\t%s%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                        headertab[next].c_str(), parameters.opt_relabel, i + 1, i, j, cluster.count, cluster.header.c_str());
              } else {
                std::fprintf(fp_tabbedout,
                        "%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t%u\t%s\n",
                        headertab[next].c_str(), cluster.header.c_str(), i, j, cluster.count, cluster.header.c_str());
              }
              ++j;
            }

          progress.update(i);
        }
    }
  }


  // dispatcher: write every requested output format, closing each handle
  // immediately after its writer (RAII handles are owned by derep()).
  auto output_results(struct Parameters const & parameters,
                      std::vector<struct bucket> const & hashtable,
                      uint64_t const clusters,
                      std::vector<unsigned int> const & nextseqtab,
                      std::vector<std::string> const & headertab,
                      std::vector<char> const & match_strand,
                      OutputFileHandle & fastaout_handle,
                      OutputFileHandle & fastqout_handle,
                      OutputFileHandle & uc_handle,
                      OutputFileHandle & tabbedout_handle) -> void
  {
    if ((parameters.opt_output != nullptr) or (parameters.opt_fastaout != nullptr))
      {
        write_fasta_output(fastaout_handle.get(), hashtable, clusters, parameters);
        fastaout_handle.reset();
      }

    if (parameters.opt_fastqout != nullptr)
      {
        write_fastq_output(fastqout_handle.get(), hashtable, clusters, parameters);
        fastqout_handle.reset();
      }

    if (parameters.opt_uc != nullptr)
      {
        write_uc_output(uc_handle.get(), hashtable, clusters, nextseqtab, headertab, match_strand, parameters);
        uc_handle.reset();
      }

    if (parameters.opt_tabbedout != nullptr)
      {
        write_tabbedout_output(tabbedout_handle.get(), hashtable, clusters, nextseqtab, headertab, parameters);
        tabbedout_handle.reset();
      }
  }

}  // end of anonymous namespace (output writers)


namespace {
  struct Derep_stats
  {
    uint64_t sequencecount = 0;
    uint64_t nucleotidecount = 0;
    int64_t shortest = std::numeric_limits<int64_t>::max();
    int64_t longest = 0;
    uint64_t discarded_short = 0;
    uint64_t discarded_long = 0;
    uint64_t clusters = 0;
    int64_t sumsize = 0;
    uint64_t maxsize = 0;
  };
}  // end of anonymous namespace


// streams every input record and builds the cluster hash table: normalize,
// find-or-create the cluster, accumulate abundance (and merge FASTQ quality).
// Grows the table and side tables as needed. Returns the run statistics; the
// filled hashtable and side tables are left in the caller-owned arguments.
// static because it takes the anonymous-namespace bucket type, like rehash().
// Not noexcept: it reads input, allocates, and may call fatal().
static auto dereplicating(std::unique_ptr<fastx_s> const & input_handle,
                          struct Parameters const & parameters,
                          Derep_mode const mode,
                          char const * input_filename,
                          std::vector<struct bucket> & hashtable,
                          std::vector<unsigned int> & nextseqtab,
                          std::vector<std::string> & headertab,
                          std::vector<char> & match_strand) -> Derep_stats
{
  /* derep_id is the only command that also requires identical headers to
     collapse two sequences into one */
  bool const use_header = (mode == Derep_mode::id);

  auto const filesize = input_handle->get_size();

  /* allocate initial memory for 1024 clusters
     with sequences of length 1023 */

  uint64_t alloc_clusters = 1024;
  uint64_t alloc_seqs = 1024;
  int64_t alloc_seqlen = 1023;

  uint64_t hashtablesize = 2 * alloc_clusters;
  uint64_t hash_mask = hashtablesize - 1;
  hashtable.resize(hashtablesize);

  // memory-intensive: the hash table has been allocated

  auto const extra_info = (parameters.opt_uc != nullptr) or (parameters.opt_tabbedout != nullptr);

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

  // memory-intensive: per-sequence buffers have been allocated

  std::vector<char> seq_up(static_cast<std::size_t>(alloc_seqlen) + 1);
  std::vector<char> rc_seq_up(static_cast<std::size_t>(alloc_seqlen) + 1);
  std::string const prompt = std::string("Dereplicating file ") + input_filename;


  uint64_t sequencecount = 0;
  uint64_t nucleotidecount = 0;
  auto shortest = std::numeric_limits<int64_t>::max();
  int64_t longest = 0;
  uint64_t discarded_short = 0;
  uint64_t discarded_long = 0;
  uint64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;

  {
    Progress progress(prompt.c_str(), filesize, parameters);
    while (input_handle->next(not parameters.opt_notrunclabels, chrmap_no_change()))
      {
        int64_t const seqlen = static_cast<int64_t>(input_handle->get_sequence_length());

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
            seq_up.resize(static_cast<std::size_t>(alloc_seqlen) + 1);
            rc_seq_up.resize(static_cast<std::size_t>(alloc_seqlen) + 1);

            // memory-intensive: sequence buffers grown to fit the longest sequence
          }

        if (extra_info and (sequencecount + 1 > alloc_seqs))
          {
            uint64_t const new_alloc_seqs = 2 * alloc_seqs;

            nextseqtab.resize(new_alloc_seqs, terminal);

            headertab.resize(new_alloc_seqs);

            match_strand.resize(new_alloc_seqs);

            alloc_seqs = new_alloc_seqs;

            // memory-intensive: per-sequence buffers have been grown
          }

        if (clusters + 1 > alloc_clusters)
          {
            uint64_t const new_alloc_clusters = 2 * alloc_clusters;

            rehash(hashtable);

            alloc_clusters = new_alloc_clusters;
            hashtablesize = 2 * alloc_clusters;
            hash_mask = hashtablesize - 1;

            // memory-intensive: the hash table has been resized (rehash)
          }

        auto const * seq = input_handle->get_sequence();
        auto const * header = input_handle->get_header();
        auto const headerlen = input_handle->get_header_length();
        auto const * qual = input_handle->get_quality(); // nullptr if FASTA

        /* normalize sequence: uppercase and replace U by T  */
        string_normalize(Span<char>{seq_up.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq, static_cast<std::size_t>(seqlen)});

        /* reverse complement if necessary */
        if (parameters.opt_strand)
          {
            reverse_complement(Span<char>{rc_seq_up.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq_up.data(), static_cast<std::size_t>(seqlen)});
          }

        /*
          Find free bucket or bucket for identical sequence.
          Make sure sequences are exactly identical
          in case of any hash collision.
          With 64-bit hashes, there is about 50% chance of a
          collision when the number of sequences is about 5e9.
        */

        auto const hash_header = use_header ? hash_function(header, headerlen) : uint64_t{0};

        auto const hash = hash_function(seq_up.data(), static_cast<uint64_t>(seqlen)) ^ hash_header;
        auto j = hash & hash_mask;
        auto * bp = &hashtable[j];  // refactoring: rename to "cluster"

        while ((bp->size != 0U) and
               ((hash != bp->hash) or
                (seqcmp(View<char>{seq_up.data(), static_cast<std::size_t>(seqlen)}, View<char>{bp->seq.data(), static_cast<std::size_t>(seqlen)}) != 0) or
                (use_header and (std::strcmp(header, bp->header.c_str()) != 0))))
          {
            j = (j + 1) & hash_mask;
            bp = &hashtable[j];
          }

        if (parameters.opt_strand and (bp->size == 0U))
          {
            /* no match on plus strand */
            /* check minus strand as well */

            auto const rc_hash = hash_function(rc_seq_up.data(), static_cast<uint64_t>(seqlen)) ^ hash_header;
            auto k = rc_hash & hash_mask;
            auto * rc_bp = &hashtable[k];

            while ((rc_bp->size != 0U)
                   and
                   ((rc_hash != rc_bp->hash) or
                    (seqcmp(View<char>{rc_seq_up.data(), static_cast<std::size_t>(seqlen)}, View<char>{rc_bp->seq.data(), static_cast<std::size_t>(seqlen)}) != 0) or
                    (use_header and (std::strcmp(header, rc_bp->header.c_str()) != 0))))
              {
                k = (k + 1) & hash_mask;
                rc_bp = &hashtable[k];
              }

            if (rc_bp->size != 0U)
              {
                bp = rc_bp;
                j = k;
                if (extra_info)
                  {
                    match_strand[sequencecount] = 1;
                  }
              }
          }

        auto const abundance = parameters.opt_sizein ? input_handle->get_abundance() : int64_t{1};
        sumsize += abundance;

        if (bp->size != 0U)
          {
            /* at least one identical sequence already */
            if (extra_info)
              {
                unsigned int const last = bp->seqno_last;
                nextseqtab[last] = static_cast<unsigned int>(sequencecount);
                bp->seqno_last = static_cast<unsigned int>(sequencecount);
                headertab[sequencecount] = header;
              }

            int64_t const s1 = bp->size;
            int64_t const s2 = abundance;
            int64_t const s3 = s1 + s2;

            if (parameters.opt_fastqout != nullptr)
              {
                /* update quality scores */
                for (int i = 0; i < seqlen; i++)
                  {
                    int const q1 = bp->qual[static_cast<std::size_t>(i)];
                    int const q2 = qual[i];
                    auto const p1 = convert_quality_symbol_to_probability(q1, parameters);
                    auto const p2 = convert_quality_symbol_to_probability(q2, parameters);
                    auto p3 = 0.0;

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
                        p3 = ((p1 * static_cast<double>(s1)) + (p2 * static_cast<double>(s2))) / static_cast<double>(s3);
                      }

                    // fastq_qout_min
                    /* max error prob, lowest quality */
                    // p3 = std::max(p1, p2);

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

                    int const q3 = convert_probability_to_quality_symbol(p3, parameters);
                    bp->qual[static_cast<std::size_t>(i)] = static_cast<char>(q3);
                  }
              }

            bp->size = static_cast<unsigned int>(s3);
            ++bp->count;
          }
        else
          {
            /* no identical sequences yet */
            bp->size = static_cast<unsigned int>(abundance);
            bp->hash = hash;
            bp->seqno_first = static_cast<unsigned int>(sequencecount);
            bp->seqno_last = static_cast<unsigned int>(sequencecount);
            bp->seq = seq;
            bp->header = header;
            bp->count = 1;
            if (qual != nullptr) {
              bp->qual = qual;
            } else {
              bp->qual.clear();
            }
            ++clusters;
          }

        maxsize = std::max<uint64_t>(bp->size, maxsize);

        ++sequencecount;

        progress.update(input_handle->get_position());
      }
  }

  Derep_stats stats;
  stats.sequencecount = sequencecount;
  stats.nucleotidecount = nucleotidecount;
  stats.shortest = shortest;
  stats.longest = longest;
  stats.discarded_short = discarded_short;
  stats.discarded_long = discarded_long;
  stats.clusters = clusters;
  stats.sumsize = sumsize;
  stats.maxsize = maxsize;
  return stats;
}


// statistics / summary reporting: each helper folds the former
// stderr-then-log duplicate blocks. The message text is written verbatim to
// stderr and, when --log is in effect, to the log; the log copy is followed
// by the extra blank line the original code emitted (the "\n\n" endings).
namespace {

  auto report_input_stats(Derep_stats const & stats,
                          struct Parameters const & parameters) -> void
  {
    auto emit = [&](std::FILE * fp) -> void {
      if (stats.sequencecount > 0)
        {
          std::fprintf(fp,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, min %" PRId64
                  ", max %" PRId64 ", avg %.0f\n",
                  stats.nucleotidecount,
                  stats.sequencecount,
                  stats.shortest,
                  stats.longest,
                  static_cast<double>(stats.nucleotidecount) * 1.0 / static_cast<double>(stats.sequencecount));
        }
      else
        {
          std::fprintf(fp,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  stats.nucleotidecount,
                  stats.sequencecount);
        }
    };
    if (not parameters.opt_quiet)
      {
        emit(stderr);
      }
    if (parameters.opt_log != nullptr)
      {
        emit(parameters.fp_log);
      }
  }


  auto report_length_filtered(struct Parameters const & parameters,
                              char const * option_name,
                              int64_t const length_limit,
                              uint64_t const discarded) -> void
  {
    if (discarded == 0U)
      {
        return;
      }
    auto emit = [&](std::FILE * fp) -> void {
      std::fprintf(fp,
              "%s %" PRId64 ": %" PRIu64 " %s discarded.\n",
              option_name,
              length_limit,
              discarded,
              (discarded == 1 ? "sequence" : "sequences"));
    };
    emit(stderr);
    if (parameters.opt_log != nullptr)
      {
        emit(parameters.fp_log);
        std::fputc('\n', parameters.fp_log);
      }
  }


  auto report_unique_summary(Derep_stats const & stats,
                             double const average,
                             double const median,
                             struct Parameters const & parameters) -> void
  {
    auto emit = [&](std::FILE * fp) -> void {
      if (stats.clusters < 1)
        {
          std::fprintf(fp,
                  "0 unique sequences\n");
        }
      else
        {
          std::fprintf(fp,
                  "%" PRIu64
                  " unique sequences, avg cluster %.1lf, median %.0f, max %"
                  PRIu64 "\n",
                  stats.clusters, average, median, stats.maxsize);
        }
    };
    if (not parameters.opt_quiet)
      {
        emit(stderr);
      }
    if (parameters.opt_log != nullptr)
      {
        emit(parameters.fp_log);
        std::fputc('\n', parameters.fp_log);
      }
  }


  auto report_selected(uint64_t const selected,
                       Derep_stats const & stats,
                       struct Parameters const & parameters) -> void
  {
    if (selected >= stats.clusters)
      {
        return;
      }
    auto emit = [&](std::FILE * fp) -> void {
      std::fprintf(fp,
              "%" PRIu64 " uniques written, %"
              PRIu64 " clusters discarded (%.1f%%)\n",
              selected, stats.clusters - selected,
              100.0 * static_cast<double>(stats.clusters - selected) / static_cast<double>(stats.clusters));
    };
    if (not parameters.opt_quiet)
      {
        emit(stderr);
      }
    if (parameters.opt_log != nullptr)
      {
        emit(parameters.fp_log);
        std::fputc('\n', parameters.fp_log);
      }
  }

}  // end of anonymous namespace (statistics reporting)


// used by --derep_fulllength, --derep_id, and --fastx_uniques
auto derep(struct Parameters const & parameters, char const * input_filename, Derep_mode const mode) -> void
{
  /* dereplicate full length sequences, optionally require identical headers */

  /*
    derep_fulllength output options: --output, --uc (only FASTA, depreciated)
    fastx_uniques output options: --fastaout, --fastqout, --uc, --tabbedout
  */

  auto input_handle = fastx_open(input_filename, parameters);

  if (not input_handle->is_empty_input())
    {
      if (input_handle->is_fastq_input())
        {
          if (mode != Derep_mode::uniques) {
            fatal("FASTQ input is only allowed with the fastx_uniques command");
          }
        }
      else
        {
          if (parameters.opt_fastqout != nullptr) {
            fatal("Cannot write FASTQ output when input file is not in FASTQ "
                  "format");
          }
          if (parameters.opt_tabbedout != nullptr) {
            fatal("Cannot write tab separated output file when input file is "
                  "not in FASTQ format");
          }
        }
    }

  OutputFileHandle fastaout_handle;
  OutputFileHandle fastqout_handle;
  OutputFileHandle uc_handle;
  OutputFileHandle tabbedout_handle;

  if (mode == Derep_mode::uniques)
    {
      if ((parameters.opt_uc == nullptr) and (parameters.opt_fastaout == nullptr) and (parameters.opt_fastqout == nullptr) and (parameters.opt_tabbedout == nullptr)) {
        fatal("Output file for dereplication with fastx_uniques must be "
              "specified with --fastaout, --fastqout, --tabbedout, or --uc");
      }
    } else {
    if ((parameters.opt_output == nullptr) and (parameters.opt_uc == nullptr)) {
      fatal("Output file for dereplication must be specified with --output "
            "or --uc");
    }
  }

  if (mode == Derep_mode::uniques)
    {
      fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
      fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
      tabbedout_handle = open_optional_output_file(parameters.opt_tabbedout, OutputOption{"--tabbedout"});
    }
  else
    {
      fastaout_handle = open_optional_output_file(parameters.opt_output, OutputOption{"--output"});
    }

  uc_handle = open_optional_output_file(parameters.opt_uc, OutputOption{"--uc"});


  std::vector<struct bucket> hashtable;
  std::vector<unsigned int> nextseqtab;
  std::vector<std::string> headertab;
  std::vector<char> match_strand;

  auto const stats = dereplicating(input_handle, parameters, mode, input_filename,
                                   hashtable, nextseqtab, headertab, match_strand);
  input_handle->report_stripped_warning(parameters);

  report_input_stats(stats, parameters);
  report_length_filtered(parameters, "minseqlength", parameters.opt_minseqlength, stats.discarded_short);
  report_length_filtered(parameters, "maxseqlength", parameters.opt_maxseqlength, stats.discarded_long);

  {
    Progress const progress("Sorting", 1, parameters);
    std::sort(hashtable.begin(), hashtable.end(), derep_bucket_before);
  }

  auto const median = find_median_size(hashtable, stats.clusters);
  auto const average = 1.0 * static_cast<double>(stats.sumsize) / static_cast<double>(stats.clusters);
  report_unique_summary(stats, average, median, parameters);

  /* count selected */

  auto const selected = count_selected(hashtable, parameters);

  /* write output */

  output_results(parameters, hashtable, stats.clusters,
                 nextseqtab, headertab, match_strand,
                 fastaout_handle, fastqout_handle, uc_handle, tabbedout_handle);

  report_selected(selected, stats, parameters);

  /* the buckets own their seq/header/qual as std::string; the hashtable's
     destruction releases them (RAII) */
}


/* === Library API implementation === */

struct derep_session_s {
  std::vector<struct bucket> hashtable;
  uint64_t hashtablesize = 0;
  uint64_t hash_mask = 0;
  uint64_t alloc_clusters = 0;
  uint64_t clusters = 0;
  unsigned int next_seqno = 0;  /* insertion counter for deterministic sort order */
  std::vector<char> seq_up;
  bool finalized = false;
};


auto derep_session_alloc() -> struct derep_session_s *
{
  return new derep_session_s {};
}


auto derep_session_free(struct derep_session_s * ds) -> void
{
  if (ds != nullptr)
    {
      derep_session_cleanup(ds);
      delete ds;
    }
}


auto derep_session_init(struct derep_session_s * ds) -> void
{
  /* Release any state from a previous session on the same handle first, so an
     init -> add* -> get_results -> init reuse starts from a clean slate. The
     buckets own their seq/header/qual as std::string (RAII), so cleanup just
     clears the hashtable and resets the counters. derep_session_cleanup() is
     idempotent and a no-op on a freshly allocated (empty) session (L2e). */
  derep_session_cleanup(ds);

  ds->alloc_clusters = 1024;
  ds->hashtablesize = 2 * ds->alloc_clusters;
  ds->hash_mask = ds->hashtablesize - 1;
  ds->hashtable.resize(ds->hashtablesize);
  ds->seq_up.resize(1024);
  ds->clusters = 0;
  ds->next_seqno = 0;
  ds->finalized = false;
}


auto derep_add_sequence(struct derep_session_s * ds,
                        const char * header,
                        const char * sequence,
                        int const seqlen,
                        int64_t const abundance) -> void
{
  if (seqlen <= 0)
    {
      return;
    }

  if (ds->finalized)
    {
      ds->finalized = false;  /* allow re-sort on next get_results */
    }

  /* Grow seq_up buffer if needed */
  if (seqlen + 1 > static_cast<int>(ds->seq_up.size()))
    {
      ds->seq_up.resize(static_cast<std::size_t>(seqlen) + 1);
    }

  /* Normalize: uppercase, U→T */
  string_normalize(Span<char>{ds->seq_up.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{sequence, static_cast<std::size_t>(seqlen)});

  /* Rehash if needed */
  if (ds->clusters + 1 > ds->alloc_clusters)
    {
      ds->alloc_clusters *= 2;
      rehash(ds->hashtable);
      ds->hashtablesize = ds->hashtable.size();
      ds->hash_mask = ds->hashtablesize - 1;
    }

  /* Hash and probe */
  auto const hash = hash_cityhash64(ds->seq_up.data(), static_cast<uint64_t>(seqlen));
  auto j = hash & ds->hash_mask;
  auto * bp = &ds->hashtable[j];

  while ((bp->size != 0U) and
         ((hash != bp->hash) or
          (seqcmp(View<char>{ds->seq_up.data(), static_cast<std::size_t>(seqlen)}, View<char>{bp->seq.data(), static_cast<std::size_t>(seqlen)}) != 0)))
    {
      j = (j + 1) & ds->hash_mask;
      bp = &ds->hashtable[j];
    }

  if (bp->size != 0U)
    {
      /* Existing unique sequence — merge */
      bp->size += static_cast<unsigned int>(abundance);
      ++bp->count;
    }
  else
    {
      /* New unique sequence */
      bp->size = static_cast<unsigned int>(abundance);
      bp->hash = hash;
      bp->seq = ds->seq_up.data();
      bp->header = header;
      bp->count = 1;
      bp->seqlen = static_cast<unsigned int>(seqlen);
      bp->qual.clear();
      bp->seqno_first = ds->next_seqno;
      bp->seqno_last = ds->next_seqno;
      ++ds->clusters;
    }
  ++ds->next_seqno;
}


auto derep_get_results(struct derep_session_s * ds,
                       struct derep_result_s * results,
                       int const max_results,
                       int * result_count) -> void
{
  /* Guard against a null output array: without this, a populated session with
     max_results > 0 would write results[0] through a null pointer (L2e). */
  if (results == nullptr)
    {
      if (result_count != nullptr)
        {
          *result_count = 0;
        }
      return;
    }

  if (!ds->finalized)
    {
      /* Sort the hashtable — same comparator as CLI */
      std::sort(ds->hashtable.begin(), ds->hashtable.end(), derep_bucket_before);
      ds->finalized = true;
    }

  int count = 0;
  for (uint64_t i = 0; i < ds->hashtablesize and count < max_results; ++i)
    {
      auto const & b = ds->hashtable[i];
      if (b.size == 0U)
        {
          break;  /* sorted: all empty buckets are at the end */
        }
      results[count].header = b.header.c_str();
      results[count].sequence = b.seq.c_str();
      results[count].abundance = b.size;
      results[count].seqlen = b.seqlen;
      results[count].count = static_cast<int>(b.count);
      ++count;
    }
  *result_count = count;
}


auto derep_session_cleanup(struct derep_session_s * ds) -> void
{
  /* the buckets own their seq/header/qual as std::string; clearing the
     hashtable destroys them (RAII) */
  ds->hashtable.clear();
  ds->clusters = 0;
  ds->finalized = false;
}
