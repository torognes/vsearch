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
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/derep.hpp"
#include "core/derep_internal.hpp"
#include "core/derep_stats.hpp"
#include "utils/quality_table.hpp"  // vsearch::QualityTable
#include "core/fasta.hpp"  // fasta_print_general
#include "core/fastq.hpp"  // fastq_print_general
#include "core/fastx.hpp"  // fastx_open, fastx_next, fastx_get_*
#include "core/quality_range.hpp"  // vsearch::check_quality_score
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/median.hpp"
#include "utils/open_file.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/cityhash.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/string_normalize.hpp"
#include <algorithm>  // std::count_if, std::equal, std::max, std::min, std::minmax_element, std::sort, std::transform
#include <array>  // std::array
#include <cassert>  // assert
#include <cmath>  // std::log10, std::pow
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <limits>
#include <memory>  // std::unique_ptr
#include <string>
#include <utility>  // std::move
#include <vector>


// refactoring: deliberately not std::unordered_map. This table is sorted in
// place once the run is over (derep_bucket_before, below), so by the end it is
// not a map at all; and the six hand-rolled open-addressing tables in the tree
// differ on seven independent axes, which is why they share sizing arithmetic
// (utils/hash_table_size.hpp) and a named occupancy predicate rather than a
// container. Reasoning and inventory: DONE_20260825_flatmap_helpers.md.
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
    /* total abundance of the cluster, as wide as the ;size= annotation it comes
       from (header_get_size returns int64_t). A 32-bit field truncated any
       annotation above 4294967295 silently, and an abundance that was an exact
       multiple of 2^32 truncated to zero -- which is this table's empty-bucket
       sentinel, so the record was dropped from the output altogether. */
    uint64_t size = 0;
    unsigned int count = 0;
    unsigned int seqlen = 0;  /* sequence length (used by API to avoid strlen) */
    std::string header;
    std::string seq;
    std::string qual;  /* empty when FASTA (no quality) */
  };


  /* the empty-bucket sentinel. Naming it spares every probe from having to
     remember which field carries that meaning -- the abundance bug recorded on
     'size' above lived exactly in that gap, and the tables in this family each
     pick a different field for the job ('size' here, 'count' in
     core/unique.cpp, 'pos' in core/kmerhash.cpp). */
  auto is_occupied(struct bucket const & entry) noexcept -> bool
  {
    return entry.size != 0U;
  }


  /*
    A linear probe stops at the first bucket that is free, or that already holds
    this exact record; it steps over every other one. Records are matched by
    their hash first, then byte by byte, because two distinct sequences can land
    in the same bucket: with 64-bit hashes, there is about 50% chance of a
    collision when the number of sequences is about 5e9.

    The length test is what makes the byte comparison safe. Without it, a
    collision between sequences of different lengths compared the shorter one
    against a window running past the end of the longer -- and reported them
    identical, merging a record into the wrong cluster.

    The comparison folds through chrmap_4bit rather than comparing raw bytes:
    the CLI stores each representative exactly as it was read (see the output,
    which preserves the case of the first occurrence), while the incoming
    sequence has been normalized, so this compares normalized against raw and
    must be blind to case and to U versus T.
  */
  auto holds_another_record(struct bucket const & candidate,
                            uint64_t const hash,
                            View<char> const seq,
                            bool const use_header = false,
                            View<char> const header = View<char>{}) -> bool
  {
    if (not is_occupied(candidate)) { return false; }  // free bucket
    if (candidate.hash != hash) { return true; }
    if (candidate.seq.size() != seq.size()) { return true; }
    auto const * const map_4bit_table = chrmap_4bit();
    auto const same_nucleotide = [map_4bit_table](char const lhs, char const rhs) -> bool {
      // the table has 256 entries, but sequence data is ASCII: this is the
      // range check map_4bit() made through to_uchar() before the table was
      // hoisted out of the loop (compiled out in release, where NDEBUG is set)
      assert(static_cast<unsigned char>(lhs) < 128U);
      assert(static_cast<unsigned char>(rhs) < 128U);
      return map_4bit_table[static_cast<unsigned char>(lhs)] ==
             map_4bit_table[static_cast<unsigned char>(rhs)];
    };
    if (not std::equal(seq.cbegin(), seq.cend(), candidate.seq.cbegin(),
                       same_nucleotide)) { return true; }
    return use_header and (header != make_view(candidate.header));
  }
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
    if (is_occupied(old_bucket)) {
      auto new_index = old_bucket.hash & new_hash_mask;
      while (is_occupied(new_hashtable[new_index])) {
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
      /* compared as int64_t, the type of the two bounds */
      auto const size = static_cast<int64_t>(bucket.size);
      return ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize));
    };
    auto const selected = std::count_if(hashtable.begin(), hashtable.end(),
                                        size_in_range);
    return std::min(static_cast<uint64_t>(selected),
                    static_cast<uint64_t>(parameters.opt_topn));
  }


  /* Reject a quality string that steps outside [--fastq_qmin, --fastq_qmax].

     fastx_uniques is the only derep mode that accepts FASTQ input, and it
     used to accept the two bounds without ever applying them, while still
     decoding every symbol to merge it. A Q50 symbol went in unremarked and
     came out silently clamped to Q41 by the qmaxout rule, and a phred+64
     file read at the default offset 33 dereplicated to fabricated qualities
     with exit 0 -- where filter, eestats, fastq_stats and fastq_convert all
     stop and name the option to raise. The wording is theirs.

     The hot loop only asks whether anything was out of range, and answers
     without branching: a legal symbol lands in the window
     [ascii + qmin, ascii + qmax], and unsigned wrap-around folds the two
     bound comparisons into one subtraction and one comparison. It is a
     reduction rather than an early-exit search so that the compiler may
     vectorize it; this walks every quality byte of every record. Measured on
     100k reads of 283 bp (--fastqout, 15 runs, one hyperfine session): 216.5
     ms unchanged, 222.9 ms here, 230.9 ms for the std::minmax_element form
     this replaced -- 3.0% against 6.7%.

     Naming the offending value is the cold path, and pays for its second
     pass only on a record that is about to be rejected anyway. The score is
     a strictly increasing function of the symbol, so the two extremes are
     the only candidates. The parser has already rejected every byte outside
     33-126 (quality_policy in core/fastq.cpp), so the plain char comparison
     there cannot go negative.

     The message is the shared one (vsearch::check_quality_score below); what
     is local here is only how the two extremes are found.

     It keeps its own pass rather than reading the parser's per-record range
     (fastx.hpp's track_quality_range(), which commands/fastq_stats.cpp does
     use). Measured on 100k records with --fastqout, -O2
     -falign-functions=64, pinned core, 10 runs: 193.9 ms with the loop
     below, 197.4 ms reading the parser's range -- 1.8% slower, not faster.
     The loop below is branchless and over a contiguous buffer, so it
     vectorizes; the parser's two comparisons sit inside a loop that already
     does a table lookup, a branch and a store per byte, and cannot. The
     expensive std::minmax_element runs only once the window test has already
     failed, which is the last record of the run. fastq_stats had no such
     guard -- it called std::minmax_element on every record -- which is why
     the same move made it 3.4% faster. */
  auto check_quality_range(View<char> const quality_symbols,
                           struct Parameters const & parameters,
                           vsearch::QualityLocation const location) -> void
  {
    if (quality_symbols.empty()) { return; }

    auto const ascii_offset = static_cast<int>(parameters.opt_fastq_ascii);
    /* cli.cc keeps ascii + qmin at 33 or more and ascii + qmax at 126 or
       less, so neither the lowest legal symbol nor the window width can
       leave the unsigned char domain */
    assert(ascii_offset + parameters.opt_fastq_qmin >= 33);
    assert(ascii_offset + parameters.opt_fastq_qmax <= 126);
    auto const lowest_legal =
      static_cast<unsigned char>(ascii_offset + parameters.opt_fastq_qmin);
    auto const window_width =
      static_cast<unsigned char>(parameters.opt_fastq_qmax - parameters.opt_fastq_qmin);

    auto any_outside = false;
    for (auto const symbol : quality_symbols)
      {
        auto const distance_from_lowest =
          static_cast<unsigned char>(static_cast<unsigned char>(symbol) - lowest_legal);
        any_outside |= (distance_from_lowest > window_width);
      }
    if (not any_outside) { return; }

    auto const extremes = std::minmax_element(quality_symbols.begin(),
                                              quality_symbols.end());
    vsearch::check_quality_score(*std::get<0>(extremes) - ascii_offset, parameters, location);
    vsearch::check_quality_score(*std::get<1>(extremes) - ascii_offset, parameters, location);
  }


  // refactoring: duplicate of q2p()?
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
    /* encoded with the INPUT offset: stored quality strings keep the input
       encoding whether or not they were touched by a merge, and the output
       offset is applied exactly once, by write_fastq_output(). Encoding
       merged symbols with the output offset here corrupted clusters of
       abundance three or more when the two offsets differed: the next merge
       decoded the stored symbol with the input offset again. */
    return static_cast<int>(quality_value + parameters.opt_fastq_ascii);
  }


  /* The output symbol standing for the error probability of one input
     symbol. write_fastq_output() sends every stored symbol through this
     table exactly once -- stored quality strings always keep the input
     encoding -- so it is both the --fastq_asciiout re-encoding step and
     the qminout/qmaxout clamp, for merged and unmerged records alike.

     Going through 10^-(q/10) and back does not survive the round trip.
     -10 * log10(10^-0.2) evaluates to 1.999999999999999778, so trunc() sent
     quality 2 back as quality 1; and at qualities 3, 5, 6, 8, 10 and 17 the
     answer moved with the abundances, because the weighted mean of two equal
     probabilities is not bit-exactly that probability and trunc()'s
     threshold sits exactly where the table's values land. Two identical Q10
     reads merged to Q10 at abundances 1+1 but to Q9 at 1+2.

     There is no logarithm to get wrong here. The forward table caps every
     probability at 0.75 (the chance of guessing one of four bases wrong), so
     qualities below 2 all share that entry and are only representable as the
     1 that 0.75 floors to; every other quality represents itself. Verified
     against an exact reference over both ASCII offsets, three qminout and
     four qmaxout settings, all 94 legal symbols: 2256 combinations, no
     disagreement. */
  auto merged_symbol_table(struct Parameters const & parameters) -> std::array<char, vsearch::quality_symbol_count>
  {
    static constexpr auto uninformative_quality = int64_t{1};
    std::array<char, vsearch::quality_symbol_count> table {{}};
    for (std::size_t symbol = 0; symbol < table.size(); ++symbol)
      {
        auto quality_value = static_cast<int64_t>(symbol) - parameters.opt_fastq_ascii;
        quality_value = std::max(quality_value, uninformative_quality);
        quality_value = std::min(quality_value, parameters.opt_fastq_qmaxout);
        quality_value = std::max(quality_value, parameters.opt_fastq_qminout);
        table[symbol] = static_cast<char>(quality_value + parameters.opt_fastq_asciiout);
      }
    return table;
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

  if (lhs.size < rhs.size)
    {
      return false;
    }
  if (lhs.size > rhs.size)
    {
      return true;
    }
  // same abundance
  if (not is_occupied(lhs))
    {
      return false;
    }
  auto const result = make_view(lhs.header)
                        .compare(make_view(rhs.header));
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
          auto const size = static_cast<int64_t>(cluster.size);
          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++relabel_count;
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  make_view(cluster.seq),
                                  make_view(cluster.header),
                                  OutputAnnotations{static_cast<uint64_t>(size), relabel_count},
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
      /* stored quality strings keep the input encoding, merged or not; this
         table is the one place the qminout/qmaxout clamp and the
         --fastq_asciiout offset are applied */
      auto const merged_symbol = merged_symbol_table(parameters);
      std::string encoded_quality;
      int64_t relabel_count = 0;
      for (uint64_t i = 0; i < clusters; ++i)
        {
          auto const & cluster = hashtable[i];
          auto const size = static_cast<int64_t>(cluster.size);
          if ((size >= parameters.opt_minuniquesize) and (size <= parameters.opt_maxuniquesize))
            {
              ++relabel_count;
              encoded_quality.resize(cluster.qual.size());
              std::transform(cluster.qual.cbegin(), cluster.qual.cend(),
                             encoded_quality.begin(),
                             [&merged_symbol](char const symbol) -> char
                {
                  /* the FASTQ readers reject DEL and every byte above it */
                  assert(static_cast<unsigned char>(symbol) < vsearch::quality_symbol_count);
                  return merged_symbol[static_cast<unsigned char>(symbol)];
                });
              fastq_print_general(fp_fastqout,
                                  make_view(cluster.seq),
                                  make_view(cluster.header),
                                  make_view(encoded_quality),
                                  OutputAnnotations{static_cast<uint64_t>(size), relabel_count},
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
          auto const len = static_cast<int64_t>(cluster.seq.size());

          fprint(fp_uc, "S\t");
          fprint_integer(fp_uc, i);
          fprint(fp_uc, '\t');
          fprint_integer(fp_uc, len);
          fprint(fp_uc, "\t*\t*\t*\t*\t*\t");
          fprint(fp_uc, make_view(cluster.header));
          fprint(fp_uc, "\t*\n");

          for (auto next = nextseqtab[cluster.seqno_first];
               next != terminal;
               next = nextseqtab[next])
            {
              fprint(fp_uc, "H\t");
              fprint_integer(fp_uc, i);
              fprint(fp_uc, '\t');
              fprint_integer(fp_uc, len);
              fprint(fp_uc, '\t');
              std::fprintf(fp_uc, "%.1f", 100.0);
              fprint(fp_uc, '\t');
              std::fputs(((match_strand[next] != 0) ? "-" : "+"), fp_uc);
              fprint(fp_uc, "\t0\t0\t*\t");
              fprint(fp_uc, make_view(headertab[next]));
              fprint(fp_uc, '\t');
              fprint(fp_uc, make_view(cluster.header));
              fprint(fp_uc, '\n');
            }

          progress.update(i);
        }
    }

    {
      Progress progress("Writing uc file, second part", clusters, parameters);
      for (uint64_t i = 0; i < clusters; ++i)
        {
          auto const & cluster = hashtable[i];
          fprint(fp_uc, "C\t");
          fprint_integer(fp_uc, i);
          fprint(fp_uc, '\t');
          fprint_integer(fp_uc, cluster.size);
          fprint(fp_uc, "\t*\t*\t*\t*\t*\t");
          fprint(fp_uc, make_view(cluster.header));
          fprint(fp_uc, "\t*\n");
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
            fprint(fp_tabbedout, make_view(cluster.header));
            fprint(fp_tabbedout, '\t');
            std::fputs(parameters.opt_relabel, fp_tabbedout);
            fprint_integer(fp_tabbedout, i + 1);
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, i);
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, static_cast<uint64_t>(0));
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, cluster.count);
            fprint(fp_tabbedout, '\t');
            fprint(fp_tabbedout, make_view(cluster.header));
            fprint(fp_tabbedout, '\n');
          } else {
            fprint(fp_tabbedout, make_view(cluster.header));
            fprint(fp_tabbedout, '\t');
            fprint(fp_tabbedout, make_view(cluster.header));
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, i);
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, static_cast<uint64_t>(0));
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, cluster.count);
            fprint(fp_tabbedout, '\t');
            fprint(fp_tabbedout, make_view(cluster.header));
            fprint(fp_tabbedout, '\n');
          }

          uint64_t j = 1;
          for (auto next = nextseqtab[cluster.seqno_first];
               next != terminal;
               next = nextseqtab[next])
            {
              if (parameters.opt_relabel != nullptr) {
                fprint(fp_tabbedout, make_view(headertab[next]));
                fprint(fp_tabbedout, '\t');
                std::fputs(parameters.opt_relabel, fp_tabbedout);
                fprint_integer(fp_tabbedout, i + 1);
                fprint(fp_tabbedout, '\t');
                fprint_integer(fp_tabbedout, i);
                fprint(fp_tabbedout, '\t');
                fprint_integer(fp_tabbedout, j);
                fprint(fp_tabbedout, '\t');
                fprint_integer(fp_tabbedout, cluster.count);
                fprint(fp_tabbedout, '\t');
                fprint(fp_tabbedout, make_view(cluster.header));
                fprint(fp_tabbedout, '\n');
              } else {
                fprint(fp_tabbedout, make_view(headertab[next]));
                fprint(fp_tabbedout, '\t');
                fprint(fp_tabbedout, make_view(cluster.header));
                fprint(fp_tabbedout, '\t');
                fprint_integer(fp_tabbedout, i);
                fprint(fp_tabbedout, '\t');
                fprint_integer(fp_tabbedout, j);
                fprint(fp_tabbedout, '\t');
                fprint_integer(fp_tabbedout, cluster.count);
                fprint(fp_tabbedout, '\t');
                fprint(fp_tabbedout, make_view(cluster.header));
                fprint(fp_tabbedout, '\n');
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

  uint64_t hashtablesize = 2 * alloc_clusters;
  uint64_t hash_mask = hashtablesize - 1;
  hashtable.resize(hashtablesize);

  // memory-intensive: the hash table has been allocated

  auto const extra_info = (parameters.opt_uc != nullptr) or (parameters.opt_tabbedout != nullptr);

  /* one 10^-(q/10) per quality symbol, not per base; the cap at 0.75 is the
     "quality below 2 carries no information" rule this loop applied per call */
  vsearch::QualityTable const quality_table(static_cast<int>(parameters.opt_fastq_ascii),
                                            vsearch::ProbabilityCap::random_guess);

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

  /* deliberately declared outside the record loop: the buffers grow to
     the longest record seen and are reused, so after the first few
     records no allocation happens at all (cppcheck's variableScope
     hint would reallocate per record) */
  std::vector<char> seq_up(1024);
  std::vector<char> rc_seq_up(1024);
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
    Progress progress(prompt, filesize, parameters);
    while (input_handle->next(not parameters.opt_notrunclabels, chrmap_no_change()))
      {
        auto const sequence = input_handle->sequence_view();
        auto const seqlen = static_cast<int64_t>(sequence.size());

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

        if (seq_up.size() < static_cast<std::size_t>(seqlen) + 1)
          {
            seq_up.resize(static_cast<std::size_t>(seqlen) + 1);
            rc_seq_up.resize(static_cast<std::size_t>(seqlen) + 1);

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

        auto const header_v = input_handle->header_view();
        auto const quality = input_handle->quality_view(); // empty if FASTA

        /* self-guarded: an empty quality view returns immediately */
        check_quality_range(quality, parameters, input_handle->quality_location());

        /* normalize sequence: uppercase and replace U by T  */
        auto const seq_up_v = normalize_into(seq_up, sequence);

        /* reverse complement if necessary */
        if (parameters.opt_strand)
          {
            reverse_complement(make_span(rc_seq_up).first(static_cast<std::size_t>(seqlen) + 1), seq_up_v);
          }

        /* Find free bucket or bucket for identical sequence (see
           holds_another_record) */

        auto const hash_header = use_header ? hash_function(header_v) : uint64_t{0};

        auto const hash = hash_function(seq_up_v) ^ hash_header;
        auto j = hash & hash_mask;
        auto * bp = &hashtable[j];

        while (holds_another_record(*bp, hash, seq_up_v, use_header, header_v))
          {
            j = (j + 1) & hash_mask;
            bp = &hashtable[j];
          }

        auto matched_minus_strand = false;

        if (parameters.opt_strand and not is_occupied(*bp))
          {
            /* no match on plus strand */
            /* check minus strand as well */

            auto const rc_seq_up_v = make_view(rc_seq_up).first(static_cast<std::size_t>(seqlen));
            auto const rc_hash = hash_function(rc_seq_up_v) ^ hash_header;
            auto k = rc_hash & hash_mask;
            auto * rc_bp = &hashtable[k];

            while (holds_another_record(*rc_bp, rc_hash, rc_seq_up_v, use_header, header_v))
              {
                k = (k + 1) & hash_mask;
                rc_bp = &hashtable[k];
              }

            if (is_occupied(*rc_bp))
              {
                bp = rc_bp;
                j = k;
                matched_minus_strand = true;
                if (extra_info)
                  {
                    match_strand[sequencecount] = 1;
                  }
              }
          }

        auto const abundance = parameters.opt_sizein ? input_handle->get_abundance() : int64_t{1};
        sumsize += abundance;

        if (is_occupied(*bp))
          {
            /* at least one identical sequence already */
            if (extra_info)
              {
                unsigned int const last = bp->seqno_last;
                nextseqtab[last] = static_cast<unsigned int>(sequencecount);
                bp->seqno_last = static_cast<unsigned int>(sequencecount);
                headertab[sequencecount].assign(header_v.data(), header_v.size());
              }

            auto const s1 = static_cast<int64_t>(bp->size);
            int64_t const s2 = abundance;
            int64_t const s3 = s1 + s2;

            if (parameters.opt_fastqout != nullptr)
              {
                /* update quality scores; a member grouped through its
                   reverse complement aligns position i of the stored (plus
                   strand) sequence with position seqlen - 1 - i of its own
                   quality string, so read that string reversed (quality
                   symbols carry no complement, only position) */
                for (int i = 0; i < seqlen; i++)
                  {
                    auto const member_position = matched_minus_strand ? (seqlen - 1 - i) : i;
                    auto const symbol1 = bp->qual[static_cast<std::size_t>(i)];
                    auto const symbol2 = quality[static_cast<std::size_t>(member_position)];
                    /* the FASTQ readers reject DEL and every byte above it,
                       so both symbols are printable and compare as expected */
                    assert(symbol1 > 0);
                    assert(symbol2 > 0);

                    /* how to compute the new quality score? */

                    if (parameters.opt_fastq_qout_max)
                      {
                        // fastq_qout_max
                        /* min error prob, highest quality */
                        /* the merged probability is always one of the two
                           inputs, so the merged symbol is the higher-quality
                           input symbol and there is no round trip to make */
                        bp->qual[static_cast<std::size_t>(i)] =
                          std::max(symbol1, symbol2);
                        continue;
                      }

                    if (symbol1 == symbol2)
                      {
                        /* both members agree on this base, so the merged
                           quality is that quality, whatever the abundances */
                        bp->qual[static_cast<std::size_t>(i)] = symbol1;
                        continue;
                      }

                    auto const p1 = quality_table[symbol1];
                    auto const p2 = quality_table[symbol2];

                    // fastq_qout_avg
                    /* average, as in USEARCH */
                    auto const p3 = ((p1 * static_cast<double>(s1)) + (p2 * static_cast<double>(s2))) / static_cast<double>(s3);

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

            bp->size = static_cast<uint64_t>(s3);
            ++bp->count;
          }
        else
          {
            /* no identical sequences yet */
            bp->size = static_cast<uint64_t>(abundance);
            bp->hash = hash;
            bp->seqno_first = static_cast<unsigned int>(sequencecount);
            bp->seqno_last = static_cast<unsigned int>(sequencecount);
            bp->seq.assign(sequence.data(), sequence.size());
            bp->header.assign(header_v.data(), header_v.size());
            bp->count = 1;
            if (quality.empty()) {
              bp->qual.clear();
            } else {
              bp->qual.assign(quality.data(), quality.size());
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




namespace {

// The input file is the argument of the very option that selects the mode, so
// the two are one choice: derep() used to take both, leaving three call sites
// free to pass a file name and a mode that disagree, and nothing to catch it.
// Same derivation derep_smallmem.cpp:225 already makes for its own option.
auto derep_input_filename(struct Parameters const & parameters,
                          Derep_mode const mode) -> char const *
{
  switch (mode)
    {
    case Derep_mode::fulllength:
      return parameters.opt_derep_fulllength;
    case Derep_mode::id:
      return parameters.opt_derep_id;
    case Derep_mode::uniques:
      return parameters.opt_fastx_uniques;
    }
  assert(false);  // unreachable: -Wswitch checks that every Derep_mode is listed
  return nullptr;
}

}  // end of anonymous namespace


// used by --derep_fulllength, --derep_id, and --fastx_uniques
auto derep(struct Parameters const & parameters, Derep_mode const mode) -> void
{
  auto * const input_filename = derep_input_filename(parameters, mode);

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
          /* --tabbedout used to be rejected here too, but none of its six
             columns depends on quality data, so fasta input is accepted */
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

  auto const median = median_of_descending(
      make_view(hashtable).first(static_cast<std::size_t>(stats.clusters)),
      [](struct bucket const & entry) { return entry.size; });
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
  auto const seq_up_v = normalize_into(ds->seq_up, View<char>{sequence, static_cast<std::size_t>(seqlen)});

  /* Rehash if needed */
  if (ds->clusters + 1 > ds->alloc_clusters)
    {
      ds->alloc_clusters *= 2;
      rehash(ds->hashtable);
      ds->hashtablesize = ds->hashtable.size();
      ds->hash_mask = ds->hashtablesize - 1;
    }

  /* Hash and probe */
  auto const hash = hash_cityhash64(seq_up_v);
  auto j = hash & ds->hash_mask;
  auto * bp = &ds->hashtable[j];

  while (holds_another_record(*bp, hash, seq_up_v))
    {
      j = (j + 1) & ds->hash_mask;
      bp = &ds->hashtable[j];
    }

  if (is_occupied(*bp))
    {
      /* Existing unique sequence — merge */
      bp->size += static_cast<uint64_t>(abundance);
      ++bp->count;
    }
  else
    {
      /* New unique sequence */
      bp->size = static_cast<uint64_t>(abundance);
      bp->hash = hash;
      /* seq_up_v, not seq_up.data(): the normalized sequence's length is
         already known here, so the assignment does not have to find the
         terminator again for every unique record */
      bp->seq.assign(seq_up_v.cbegin(), seq_up_v.cend());
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
      if (not is_occupied(b))
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
