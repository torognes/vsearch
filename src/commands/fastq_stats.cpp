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

#include "vsearch.hpp"
#include "core/fastq.hpp"
#include "core/quality_range.hpp"  // vsearch::check_quality_score
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/view.hpp"
#include <array>
#include <algorithm>  // std::max, std::min, std::find_if, std::transform, std::minmax_element, std::for_each
#include <cassert>
#include <cmath>  // std::pow
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <functional>  // std::plus
#include <iterator>  // std::distance, std::next
#include <limits>
#include <numeric>  // std::partial_sum, std::inner_product, std::iota, std::accumulate
#include <vector>


constexpr auto n_eight_bit_values = std::size_t{256};
using Length_vs_Quality_counts = std::vector<std::vector<uint64_t>>;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct Distributions {
    double avgq = 0.0;
    double avgp = 0.0;
    double avgee = 0.0;
    double rate = 0.0;
  };

  struct Stats {
    uint64_t len_min;
    uint64_t len_max;
    double n_symbols;
    uint64_t seq_count;
    double n_sequences;
    std::vector<uint64_t> length_dist;
    std::vector<uint64_t> quality_dist;
    std::vector<struct Distributions> distributions;
  };


  auto q2p(double const quality_score) -> double {
    static constexpr auto base = 10.0;
    return std::pow(base, -quality_score / base);
  }



  /* The range comes from the parser, which recorded it while copying the
     quality line, so this no longer walks the string a second time.
     symbol_to_score is std::iota from -opt_fastq_ascii and therefore
     increases with the symbol, so the lowest symbol carries the lowest score
     and the highest the highest -- the same two values std::minmax_element
     used to return here. It stays signed throughout: --fastq_qmin may be
     negative (cli.cc only requires fastq_ascii + fastq_qmin >= 33, so down
     to -31 at offset 64), and so may the decoded score; comparing through
     unsigned int once turned every negative bound into a value near 2^32,
     which failed every quality value in every file.

     This command used to phrase its own message ("... out of range (0-41)",
     plus a four-line hint pointing at --fastq_chars). It now uses the shared
     one, so that a user who moves between commands on the same bad file
     reads the same sentence -- see DONE_20260825_quality_range.md phase 5. */
  auto check_minmax_scores(QualitySymbolRange const range,
                           std::vector<int64_t> const & symbol_to_score,
                           struct Parameters const & parameters,
                           vsearch::QualityLocation const location) -> void {
    if (not range.seen()) { return; }
    vsearch::check_quality_score(symbol_to_score[range.lowest], parameters, location);
    vsearch::check_quality_score(symbol_to_score[range.highest], parameters, location);
  }


  auto const is_observed = [](uint64_t const count) -> bool { return count != 0UL; };


  auto find_smallest(std::vector<uint64_t> const & observables) -> unsigned long {
    auto const first_hit =
      std::find_if(observables.begin(), observables.end(),
                   is_observed);
    if (first_hit == observables.end()) {
      return 0UL;
    }
    return static_cast<unsigned long>(
                                      std::distance(observables.begin(), first_hit));
  }


  auto find_largest(std::vector<uint64_t> const & observables) -> unsigned long {
    auto const last_hit =
      std::find_if(observables.rbegin(), observables.rend(),
                   is_observed);
    if (last_hit == observables.rend()) {
      return 0UL;
    }
    return static_cast<unsigned long>(
                                      std::distance(last_hit, observables.rend()) - 1);
  }


  auto compute_cumulative_sum(std::vector<uint64_t> const & read_length_table)
    -> std::vector<uint64_t> {
    std::vector<uint64_t> cumulative_sum_of_lengths(read_length_table.size());
    std::partial_sum(read_length_table.cbegin(), read_length_table.cend(),
                     cumulative_sum_of_lengths.begin());
    return cumulative_sum_of_lengths;
  }


  auto compute_number_of_symbols(std::vector<uint64_t> const & n_reads_per_length)
    -> double {
    // total number of nucleotides = sum(read_length * n_reads_with_that_length)
    std::vector<uint64_t> read_lengths(n_reads_per_length.size());
    std::iota(read_lengths.begin(), read_lengths.end(), 0UL);
    return std::inner_product(read_lengths.begin(), read_lengths.end(),
                              n_reads_per_length.begin(), double{0});
  }


  auto compute_n_symbols_per_length(Length_vs_Quality_counts const & qual_length_table) -> std::vector<uint64_t> {
    // sum_counts is the sum of observed valid symbols for each length
    // (invalid symbols are guaranteed to be set to zero)
    std::vector<uint64_t> sum_counts(qual_length_table.size());
    std::transform(
                   qual_length_table.begin(), qual_length_table.end(), sum_counts.begin(),
                   [](std::vector<uint64_t> const & quality_symbols) -> std::uint64_t {
                     return std::accumulate(quality_symbols.begin(), quality_symbols.end(),
                                            std::uint64_t{0});
                   });
    return sum_counts;
  }


  auto precompute_quality_scores(struct Parameters const & parameters) -> std::vector<int64_t> {
    // quality score = quality symbol - opt_fastq_ascii
    //
    // opt_fastq_ascii is a fix value and quality_symbol increases
    // linearly, so a vector of quality_scores can be generated by
    // std::iota
    //
    // the walk starts at symbol zero, not at the offset: symbols below
    // opt_fastq_ascii decode to a negative score, which is what a Solexa
    // quality string carries. Leaving those entries at their zero-initialized
    // value made every such symbol read as quality zero, so they escaped the
    // --fastq_qmin check and were tallied as Q=0 / Pe=1.0 in the report.
    std::vector<int64_t> quality_scores(n_eight_bit_values);
    std::iota(quality_scores.begin(), quality_scores.end(), -parameters.opt_fastq_ascii);
    return quality_scores;
  }


  auto compute_sum_quality_scores_per_length(Length_vs_Quality_counts const & qual_length_table,
                                             struct Parameters const & parameters) -> std::vector<int64_t> {
    // sum_quality_scores is the sum of observed scores for each length
    std::vector<int64_t> sum_quality_scores(qual_length_table.size());
    auto const quality_scores = precompute_quality_scores(parameters);
    std::transform(
                   qual_length_table.begin(), qual_length_table.end(),
                   sum_quality_scores.begin(),
                   [&quality_scores](std::vector<uint64_t> const & quality_symbols)
                   -> std::int64_t {
                     /* explicit product: the default one would multiply an
                        unsigned count by a signed score and settle on the
                        unsigned type, turning a negative score into a huge
                        positive contribution */
                     return std::inner_product(quality_symbols.begin(),
                                               quality_symbols.end(),
                                               quality_scores.begin(),
                                               std::int64_t{0},
                                               std::plus<std::int64_t>{},
                                               [](std::uint64_t const count,
                                                  std::int64_t const score) -> std::int64_t {
                                                 return static_cast<std::int64_t>(count) * score;
                                               });
                   });
    return sum_quality_scores;
  }


  auto precompute_probability_values(struct Parameters const & parameters) -> std::vector<double> {
    // probability value = 10 ^ - (quality score / 10)
    // quality score = quality symbol - opt_fastq_ascii
    //
    // opt_fastq_ascii is a fix value and quality_symbol increases
    // linearly, so a vector of quality_scores can be generated by
    // std::iota and then transformed into probability values
    auto const quality_scores = precompute_quality_scores(parameters);
    std::vector<double> probability_values(n_eight_bit_values);
    std::transform(quality_scores.cbegin(), quality_scores.cend(),
                   probability_values.begin(),
                   [](int64_t const quality_score) -> double {
                     return q2p(static_cast<double>(quality_score));
                   });
    return probability_values;
  }


  auto compute_sum_error_probabilities_per_length(Length_vs_Quality_counts const & qual_length_table, struct Parameters const & parameters) -> std::vector<double> {
    // sum_error_probabilities is the sum of observed probabilities for each length
    std::vector<double> sum_error_probabilities(qual_length_table.size());
    auto const probability_values = precompute_probability_values(parameters);
    std::transform(
                   qual_length_table.begin(), qual_length_table.end(),
                   sum_error_probabilities.begin(),
                   [&probability_values](std::vector<uint64_t> const & quality_symbols)
                   -> double {
                     return std::inner_product(quality_symbols.begin(),
                                               quality_symbols.end(),
                                               probability_values.begin(),
                                               double{0});
                   });
    return sum_error_probabilities;
  }


  auto compute_distribution_of_quality_symbols(Length_vs_Quality_counts const & length_vs_quality) -> std::vector<uint64_t> {
    // for each quality symbol: sum symbol observations for each position
    std::vector<uint64_t> distribution(n_eight_bit_values);
    std::for_each(length_vs_quality.begin(), length_vs_quality.end(),
                  [& distribution](std::vector<uint64_t> const & observations) -> void {
                    std::transform(observations.begin(),
                                   observations.end(),
                                   distribution.begin(),
                                   distribution.begin(),
                                   std::plus<uint64_t>{});
                  });
    return distribution;
  }


  auto compute_distributions(
                             unsigned int const len_max,
                             Length_vs_Quality_counts const & qual_length_table,
                             std::vector<double> const & sumee_length_table,
                             struct Parameters const & parameters) -> std::vector<struct Distributions> {
    std::vector<struct Distributions> distributions(len_max + 1);

    auto const sum_counts = compute_n_symbols_per_length(qual_length_table);
    auto const sum_quality_scores = compute_sum_quality_scores_per_length(qual_length_table, parameters);
    auto const sum_error_probabilities = compute_sum_error_probabilities_per_length(qual_length_table, parameters);

    auto position = std::size_t{0};
    for (auto & distribution: distributions) {
      auto const n_symbols = static_cast<double>(sum_counts[position]);
      auto const length = static_cast<double>(position + 1);
      auto const sum_quality_score = static_cast<double>(sum_quality_scores[position]);
      distribution.avgq = sum_quality_score / n_symbols;
      distribution.avgp = sum_error_probabilities[position] / n_symbols;
      distribution.avgee = sumee_length_table[position] / n_symbols;
      distribution.rate = distribution.avgee / length;
      ++position;
    }
    return distributions;
  }


  auto find_first_complete_EE_filtering(struct Stats const & stats,
                                        std::vector<std::array<uint64_t, 4>> const & ee_length_table) -> uint64_t {
    // - find the first position where no reads remain when filtering
    //   with the threshold EE > 1.0
    // - by construction, more stringent filtering thresholds also have no reads,
    // - downstream positions also contain no reads
    auto read_count_is_null = [](std::array<uint64_t, 4> const & read_counts) -> bool {
      return read_counts.front() == 0;
    };
    auto const iterator = std::find_if(ee_length_table.cbegin(), ee_length_table.cend(), read_count_is_null);
    if (iterator == ee_length_table.end()) {
      return stats.len_max;
    }
    return static_cast<uint64_t>(std::distance(ee_length_table.begin(), iterator));
  }


  // section 1
  auto report_read_length_distribution(std::FILE * log_handle,
                                       struct Stats const & stats,
                                       std::vector<uint64_t> const & read_length_table) -> void {
    assert(log_handle != nullptr);
    fprint(log_handle, '\n');
    fprint(log_handle, "Read length distribution");
    fprint(log_handle, '\n');
    fprint(log_handle, "      L           N      Pct   AccPct");
    fprint(log_handle, '\n');
    fprint(log_handle, "-------  ----------  -------  -------");
    fprint(log_handle, '\n');
    for (auto length = stats.len_max; length >= stats.len_min; --length)
      {
        if (read_length_table[length] != 0) {
          auto const previous_count = (length != 0) ? static_cast<double>(stats.length_dist[length - 1]) : 0;
          /* the "%2s" this replaces padded a field whose two branches are both
             already two characters wide, so the padding never applied */
          std::fputs((length == stats.len_max ? ">=" : "  "), log_handle);
          fprint_integer(log_handle, length, 5);
          fprint(log_handle, "  ");
          fprint_integer(log_handle, read_length_table[length], 10);
          fprint(log_handle, "   ");
          std::fprintf(log_handle, "%5.1lf",
                       static_cast<double>(read_length_table[length]) * 100.0 / stats.n_sequences);
          fprint(log_handle, "%   ");
          std::fprintf(log_handle, "%5.1lf",
                       100.0 * (stats.n_sequences - previous_count) / stats.n_sequences);
          fprint(log_handle, "%\n");
        }
        if (length == 0UL) { break; }
      }
  }


  // section 2
  auto report_q_score_distribution(
                                   std::FILE * log_handle,
                                   struct Stats const & stats,
                                   std::vector<double> const & symbol_to_probability,
                                   std::vector<int64_t> const & symbol_to_score) -> void {
    assert(log_handle != nullptr);
    auto const qmin = static_cast<int>(find_smallest(stats.quality_dist));
    auto const qmax = static_cast<int>(find_largest(stats.quality_dist));

    fprint(log_handle, '\n');
    fprint(log_handle, "Q score distribution");
    fprint(log_handle, '\n');
    fprint(log_handle, "ASCII    Q       Pe           N      Pct   AccPct");
    fprint(log_handle, '\n');
    fprint(log_handle, "-----  ---  -------  ----------  -------  -------");
    fprint(log_handle, '\n');
    uint64_t qual_accum = 0;
    for (auto quality_symbol = qmax ; quality_symbol >= qmin ; --quality_symbol)
      {
        auto const symbol_index = static_cast<std::size_t>(quality_symbol);
        if (stats.quality_dist[symbol_index] == 0) { continue; }

        qual_accum += stats.quality_dist[symbol_index];
        fprint(log_handle, "    ");
        fprint(log_handle, static_cast<char>(quality_symbol));
        fprint(log_handle, "  ");
        fprint_integer(log_handle, symbol_to_score[symbol_index], 3);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%7.5lf", symbol_to_probability[symbol_index]);
        fprint(log_handle, "  ");
        fprint_integer(log_handle, stats.quality_dist[symbol_index], 10);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%6.1lf", 100.0 * static_cast<double>(stats.quality_dist[symbol_index]) / stats.n_symbols);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%6.1lf", 100.0 * static_cast<double>(qual_accum) / stats.n_symbols);
        fprint(log_handle, "%\n");
      }
  }


  // section 3
  auto report_length_vs_quality_distribution(std::FILE * log_handle,
                                             struct Stats const & stats) -> void {
    assert(log_handle != nullptr);
    fprint(log_handle, '\n');
    fprint(log_handle, "    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct");
    fprint(log_handle, '\n');
    fprint(log_handle, "-----  -------  ----  -------  --------  -----  ---------  --------");
    fprint(log_handle, '\n');

    for (auto length = uint64_t{2}; length <= stats.len_max; ++length)
      {
        auto const previous_count = static_cast<double>(stats.length_dist[length - 1]);
        auto const & distribution = stats.distributions[length - 1];
        auto const PctRecs = 100.0 * (stats.n_sequences - previous_count) / stats.n_sequences;
        auto const AvgQ = distribution.avgq;
        auto const AvgP = distribution.avgp;
        auto const AvgEE = distribution.avgee;
        auto const Rate = distribution.rate;

        fprint_integer(log_handle, length, 5);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%6.1lf", PctRecs);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%4.1lf", AvgQ);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%7.5lf", q2p(AvgQ));
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%8.6lf", AvgP);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%5.2lf", AvgEE);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%9.6lf", Rate);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%7.3lf", 100.0 * Rate);
        fprint(log_handle, "%\n");
      }
  }


  // section 4
  auto report_expected_error_and_length_filtering(std::FILE * log_handle,
                                                  struct Stats const & stats,
                                                  std::vector<std::array<uint64_t, 4>> const & ee_length_table) -> void {
    assert(log_handle != nullptr);
    fprint(log_handle, '\n');
    fprint(log_handle, "    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000");
    fprint(log_handle, '\n');
    fprint(log_handle, "-----  -------  -------  -------  -------  -------  -------  -------  -------");
    fprint(log_handle, '\n');

    std::vector<double> read_percentage(ee_length_table[0].size());
    auto const max_length = find_first_complete_EE_filtering(stats, ee_length_table);
    for (auto length = max_length; length >= 1UL; --length)
      {
        auto const & read_count = ee_length_table[length - 1];
        std::transform(
            read_count.begin(), read_count.end(), read_percentage.begin(),
            [&stats](unsigned long const count) -> double {
              return 100.0 * static_cast<double>(count) / stats.n_sequences;
            });

        fprint_integer(log_handle, length, 5);
        fprint(log_handle, "  ");
        fprint_integer(log_handle, read_count[0], 7);
        fprint(log_handle, "  ");
        fprint_integer(log_handle, read_count[1], 7);
        fprint(log_handle, "  ");
        fprint_integer(log_handle, read_count[2], 7);
        fprint(log_handle, "  ");
        fprint_integer(log_handle, read_count[3], 7);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%6.2lf", read_percentage[0]);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%6.2lf", read_percentage[1]);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%6.2lf", read_percentage[2]);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%6.2lf", read_percentage[3]);
        fprint(log_handle, "%\n");
      }
  }


  // section 5
  auto report_minimum_quality_and_length_filtering(std::FILE * log_handle,
                                                   struct Stats const & stats,
                                                   std::vector<std::array<uint64_t, 4>> const & q_length_table) -> void {
    assert(log_handle != nullptr);
    fprint(log_handle, '\n');
    fprint(log_handle, "Truncate at first Q");
    fprint(log_handle, '\n');
    fprint(log_handle, "  Len     Q=5    Q=10    Q=15    Q=20");
    fprint(log_handle, '\n');
    fprint(log_handle, "-----  ------  ------  ------  ------");
    fprint(log_handle, '\n');
    auto const mid_length = std::max(uint64_t{1}, stats.len_max / 2);
    std::vector<double> read_percentage(q_length_table[0].size());
    for (auto length = stats.len_max; length >= mid_length; --length)
      {
        auto const & read_count = q_length_table[length - 1];
        std::transform(
            read_count.begin(), read_count.end(), read_percentage.begin(),
            [&stats](unsigned long const count) -> double {
              return 100.0 * static_cast<double>(count) / stats.n_sequences;
            });

        fprint_integer(log_handle, length, 5);
        fprint(log_handle, "  ");
        std::fprintf(log_handle, "%5.1lf", read_percentage[0]);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%5.1lf", read_percentage[1]);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%5.1lf", read_percentage[2]);
        fprint(log_handle, "%  ");
        std::fprintf(log_handle, "%5.1lf", read_percentage[3]);
        fprint(log_handle, "%\n");
      }
  }


  // closing section
  auto report_sequence_stats(std::FILE * log_handle, struct Stats const & stats) -> void {
    assert(log_handle != nullptr);
    static constexpr auto a_million = double{1000000};
    auto const n_sequences = static_cast<double>(stats.seq_count);
    fprint(log_handle, '\n');
    fprint_integer(log_handle, stats.seq_count, 10);
    fprint(log_handle, "  Recs (");
    std::fprintf(log_handle, "%.1lf", n_sequences / a_million);
    fprint(log_handle, "M), 0 too long\n");
    if (stats.seq_count != 0)
      {
        std::fprintf(log_handle, "%10.1lf", 1.0 * stats.n_symbols / n_sequences);
        fprint(log_handle, "  Avg length\n");
      }
    std::fprintf(log_handle, "%9.1lf", stats.n_symbols / a_million);
    fprint(log_handle, "M  Bases\n");
  }

}  // end of anonymous namespace


auto fastq_stats(struct Parameters const & parameters) -> void
{
  constexpr auto initial_memory_allocation = std::size_t{512};
  constexpr std::array<int64_t, 4> quality_thresholds = {5, 10, 15, 20};
  constexpr std::array<double, 4> ee_thresholds = { 1.0, 0.5, 0.25, 0.1 };

  /* the whole report is written to the log file; the manual lists --log
     as mandatory, and without it the command silently discarded every
     computed statistic and exited 0 */
  if (parameters.fp_log == nullptr)
    {
      fatal("Output file for fastq_stats must be specified with --log");
    }

  auto input_handle = fastq_open(parameters.opt_fastq_stats, parameters);
  /* the parser records the quality-symbol range of each record as it copies
     the quality line, which is what check_minmax_scores() reads below */
  input_handle->track_quality_range();

  auto const filesize = input_handle->get_size();


  auto const symbol_to_score = precompute_quality_scores(parameters);
  auto const symbol_to_probability = precompute_probability_values(parameters);
  std::vector<uint64_t> read_length_table(initial_memory_allocation);
  Length_vs_Quality_counts qual_length_table(initial_memory_allocation, std::vector<uint64_t>(n_eight_bit_values));
  /* ee_length_table and q_length_table hold, per position and threshold, the
     count of reads whose running minimum quality (resp. running expected
     error) still passes the threshold at that position. During the scan they
     record only each read's LAST passing position; a read passes at every
     position up to that one, so a suffix sum after the scan turns them into
     the per-position counts the report expects. */
  std::vector<std::array<uint64_t, 4>> ee_length_table(initial_memory_allocation);
  std::vector<std::array<uint64_t, 4>> q_length_table(initial_memory_allocation);
  std::vector<double> sumee_length_table(initial_memory_allocation);

  // refactoring: separate parse_fastq() and compute_distributions()
  // note: fastq parsing represents 99% of total wallclock time
  {
    Progress progress("Reading FASTQ file", filesize, parameters);
    while (input_handle->next(false, chrmap_upcase()))
      {

        /* update length statistics */

        auto const length = input_handle->get_sequence_length();

        if (length + 1 > read_length_table.size())
          {
            read_length_table.resize(length + 1);
            qual_length_table.resize(length + 1, std::vector<uint64_t>(n_eight_bit_values));
            ee_length_table.resize(length + 1);
            q_length_table.resize(length + 1);
            sumee_length_table.resize(length + 1);
          }

        ++read_length_table[length];  // can NOT be derived from qual_length_table


        /* update quality statistics */

        auto const * quality_symbols = input_handle->get_quality();
        auto expected_error = 0.0;
        auto qmin = std::numeric_limits<int64_t>::max();  // lowest Q value observed so far in this read

        check_minmax_scores(input_handle->quality_symbol_range(), symbol_to_score,
                            parameters, input_handle->quality_location());

        /* Within a read, qmin only decreases and expected_error only
           increases, so each threshold flips exactly once, from passing to
           failing: "qmin > 5, 10, 15 or 20" and "EE <= 1.0, 0.5, 0.25 or
           0.1" both hold on a prefix of the read. Instead of testing the
           four thresholds at every base, track how many still pass
           (quality_thresholds ascending, so the passing ones are its first
           n_passing_q; ee_thresholds descending, likewise) and record only
           the position where each one flips; the suffix sum after the scan
           spreads that record over the whole prefix. */
        auto n_passing_q = quality_thresholds.size();
        auto n_passing_ee = ee_thresholds.size();

        for (auto i = 0UL; i < length; ++i)
          {
            auto const quality_symbol = static_cast<unsigned char>(quality_symbols[i]);
            auto const quality_score = symbol_to_score[quality_symbol];

            ++qual_length_table[i][quality_symbol];

            if (quality_score < qmin)
              {
                qmin = quality_score;
                while ((n_passing_q != 0) and (qmin <= quality_thresholds[n_passing_q - 1]))
                  {
                    --n_passing_q;
                    if (i != 0) {
                      ++q_length_table[i - 1][n_passing_q];  // last passing position
                    }
                  }
              }

            expected_error += symbol_to_probability[quality_symbol];

            sumee_length_table[i] += expected_error;  // can NOT be derived from qual_length_table

            while ((n_passing_ee != 0) and (expected_error > ee_thresholds[n_passing_ee - 1]))
              {
                --n_passing_ee;
                if (i != 0) {
                  ++ee_length_table[i - 1][n_passing_ee];  // last passing position
                }
              }
          }

        // thresholds still passing at the end of the read pass everywhere
        if (length != 0)
          {
            for (auto threshold = 0UL; threshold < n_passing_q; ++threshold) {
              ++q_length_table[length - 1][threshold];
            }
            for (auto threshold = 0UL; threshold < n_passing_ee; ++threshold) {
              ++ee_length_table[length - 1][threshold];
            }
          }

        progress.update(input_handle->get_position());
      }
  }
  input_handle->report_stripped_warning(parameters);


  // note: operations below represent 1% of total wallclock time

  /* a read whose last passing position is p passes at every position up to
     p: suffix sums turn the last-passing-position records into the
     per-position read counts the report tables hold (see the scan above) */
  auto const spread_over_prefixes = [](std::vector<std::array<uint64_t, 4>> & table) -> void {
    std::array<uint64_t, 4> reads_passing {{}};
    for (auto position = table.size(); position != 0; --position)
      {
        auto & row = table[position - 1];
        std::transform(row.begin(), row.end(),
                       reads_passing.begin(), reads_passing.begin(),
                       std::plus<uint64_t>{});
        row = reads_passing;
      }
  };
  spread_over_prefixes(q_length_table);
  spread_over_prefixes(ee_length_table);

  /* compute various distributions */

  auto const n_reads =
    std::accumulate(read_length_table.begin(), read_length_table.end(), std::uint64_t{0});
  auto const stats = Stats{
    find_smallest(read_length_table), find_largest(read_length_table),
    compute_number_of_symbols(read_length_table),
    n_reads,
    static_cast<double>(n_reads),
    compute_cumulative_sum(read_length_table),
    compute_distribution_of_quality_symbols(qual_length_table),
    compute_distributions(static_cast<unsigned int>(find_largest(read_length_table)), qual_length_table, sumee_length_table, parameters),
  };


  /* print report */

  if (parameters.fp_log != nullptr)
    {
      report_read_length_distribution(parameters.fp_log, stats, read_length_table);  // section 1
      report_q_score_distribution(parameters.fp_log, stats, symbol_to_probability, symbol_to_score);  // section 2
      report_length_vs_quality_distribution(parameters.fp_log, stats);  // section 3
      report_expected_error_and_length_filtering(parameters.fp_log, stats, ee_length_table);  // section 4
      report_minimum_quality_and_length_filtering(parameters.fp_log, stats, q_length_table);  // section 5
      report_sequence_stats(parameters.fp_log, stats);  // closing section
    }

  if (not parameters.opt_quiet)
    {
      fprint(stderr, "Read ");
      fprint_integer(stderr, stats.seq_count);
      fprint(stderr, " sequences.\n");
    }
}
