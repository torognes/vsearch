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
#include "utils/maps.hpp"
#include <array>
#include <algorithm>  // std::max, std::find_if, std::transform, std::minmax_element
#include <cassert>
#include <cinttypes>  // macros PRIu64 (for uint64_t) and PRId64 (for int64_t)
#include <cmath>  // std::pow
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <functional>  // std::plus
#include <iterator>  // std::distance, std::next
#include <limits>
#include <numeric>  // std::partial_sum
#include <string>
#include <vector>


constexpr auto initial_memory_allocation = std::size_t{512};
constexpr std::array<uint64_t, 4> quality_thresholds = {5, 10, 15, 20};
constexpr std::array<double, 4> ee_thresholds = { 1.0, 0.5, 0.25, 0.1 };
constexpr auto n_eight_bit_values = std::size_t{256};

struct Span {
  char * start;
  std::size_t n_elements;
};

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


auto q2p(double quality_score) -> double {
  static constexpr auto base = 10.0;
  return std::pow(base, -quality_score / base);
}


auto q2p(uint64_t quality_score) -> double {
  static constexpr auto base = 10.0;
  auto const quality_score_double = static_cast<double>(quality_score);
  return std::pow(base, -quality_score_double / base);
}


auto check_quality_score(struct Parameters const & parameters, unsigned int const quality_score) -> void {
  auto const is_in_accepted_range =
    (quality_score >= static_cast<unsigned int>(parameters.opt_fastq_qmin)) and
    (quality_score <= static_cast<unsigned int>(parameters.opt_fastq_qmax));

  if (is_in_accepted_range) {
    return;
  }

  std::string const message =
    std::string("FASTQ quality value (") + std::to_string(quality_score) +
    ") out of range (" + std::to_string(parameters.opt_fastq_qmin) + "-" +
    std::to_string(parameters.opt_fastq_qmax) + ").\n" +
    "Please adjust the FASTQ quality base character or range with the\n" +
    "--fastq_ascii, --fastq_qmin or --fastq_qmax options. For a complete\n" +
    "diagnosis with suggested values, please run vsearch --fastq_chars file.";
  fatal(message.c_str());
}


auto check_minmax_scores(struct Span const a_span,
                         std::vector<uint64_t> const & symbol_to_score,
                         struct Parameters const & parameters) -> void {
  if (a_span.n_elements == 0) { return; }
  assert(a_span.n_elements <= std::numeric_limits<int64_t>::max());
  auto * const end = std::next(a_span.start, static_cast<int64_t>(a_span.n_elements));
  auto const minmax_scores =
    std::minmax_element(a_span.start, end);
  auto const qmin = symbol_to_score[*std::get<0>(minmax_scores)];
  auto const qmax = symbol_to_score[*std::get<1>(minmax_scores)];
  check_quality_score(parameters, qmin);
  check_quality_score(parameters, qmax);
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


auto compute_n_symbols_per_length(std::vector<std::array<uint64_t, n_eight_bit_values>> const & qual_length_table) -> std::vector<uint64_t> {
  // sum_counts is the sum of observed valid symbols for each length
  // (invalid symbols are guaranteed to be set to zero)
  std::vector<uint64_t> sum_counts(qual_length_table.size());
  std::transform(
      qual_length_table.begin(), qual_length_table.end(), sum_counts.begin(),
      [](std::array<uint64_t, n_eight_bit_values> const & quality_symbols) -> std::uint64_t {
        return std::accumulate(quality_symbols.begin(), quality_symbols.end(),
                               std::uint64_t{0});
      });
  return sum_counts;
}


auto precompute_quality_scores(struct Parameters const & parameters) -> std::vector<uint64_t> {
  // quality score = quality symbol - opt_fastq_ascii
  //
  // opt_fastq_ascii is a fix value and quality_symbol increases
  // linearly, so a vector of quality_scores can be generated by
  // std::iota
  std::vector<uint64_t> quality_scores(n_eight_bit_values);
  auto starting_position = std::next(quality_scores.begin(), parameters.opt_fastq_ascii);
  std::iota(starting_position, quality_scores.end(), 0UL);
  return quality_scores;
}


auto compute_sum_quality_scores_per_length(std::vector<std::array<uint64_t, n_eight_bit_values>> const & qual_length_table, struct Parameters const & parameters) -> std::vector<uint64_t> {
  // sum_quality_scores is the sum of observed scores for each length
  std::vector<uint64_t> sum_quality_scores(qual_length_table.size());
  auto const quality_scores = precompute_quality_scores(parameters);
  std::transform(
      qual_length_table.begin(), qual_length_table.end(),
      sum_quality_scores.begin(),
      [&quality_scores](std::array<uint64_t, n_eight_bit_values> const & quality_symbols)
          -> std::uint64_t {
        return std::inner_product(quality_symbols.begin(),
                                  quality_symbols.end(),
                                  quality_scores.begin(),
                                  std::uint64_t{0});
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
                 [](uint64_t const quality_score) -> double {
                   return q2p(quality_score);
                 });
  return probability_values;
}


auto compute_sum_error_probabilities_per_length(std::vector<std::array<uint64_t, n_eight_bit_values>> const & qual_length_table, struct Parameters const & parameters) -> std::vector<double> {
  // sum_error_probabilities is the sum of observed probabilities for each length
  std::vector<double> sum_error_probabilities(qual_length_table.size());
  auto const probability_values = precompute_probability_values(parameters);
  std::transform(
      qual_length_table.begin(), qual_length_table.end(),
      sum_error_probabilities.begin(),
      [&probability_values](std::array<uint64_t, n_eight_bit_values> const & quality_symbols)
          -> double {
        return std::inner_product(quality_symbols.begin(),
                                  quality_symbols.end(),
                                  probability_values.begin(),
                                  double{0});
      });
  return sum_error_probabilities;
}


auto compute_distribution_of_quality_symbols(std::vector<std::array<uint64_t, n_eight_bit_values>> const & length_vs_quality) -> std::vector<uint64_t> {
  // for each quality symbol: sum symbol observations for each position
  std::vector<uint64_t> distribution(n_eight_bit_values);
  std::for_each(length_vs_quality.begin(), length_vs_quality.end(),
                [& distribution](std::array<uint64_t, n_eight_bit_values> const & observations) {
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
    std::vector<std::array<uint64_t, n_eight_bit_values>> const & qual_length_table,
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
    distribution.rate = distributions[position].avgee / length;
    ++position;
  }
  return distributions;
}


// section 1
auto report_read_length_distribution(std::FILE * log_handle,
                                     struct Stats const & stats,
                                     std::vector<uint64_t> const & read_length_table) -> void {
  assert(log_handle != nullptr);
  std::fprintf(log_handle, "\n%s\n%s\n%s\n",
               "Read length distribution",
               "      L           N      Pct   AccPct",
               "-------  ----------  -------  -------");
  for (auto length = stats.len_max; length >= stats.len_min; --length)
    {
      if (read_length_table[length] != 0) {
        auto const previous_count = (length != 0) ? static_cast<double>(stats.length_dist[length - 1]) : 0;
        std::fprintf(log_handle, "%2s%5" PRIu64 "  %10" PRIu64 "   %5.1lf%%   %5.1lf%%\n",
                     (length == stats.len_max ? ">=" : "  "),
                     length,
                     read_length_table[length],
                     static_cast<double>(read_length_table[length]) * 100.0 / stats.n_sequences,
                     100.0 * (stats.n_sequences - previous_count) / stats.n_sequences);
      }
      if (length == 0UL) { break; }
    }
}


// section 2
auto report_q_score_distribution(
    std::FILE * log_handle,
    struct Stats const & stats,
    std::vector<double> const & symbol_to_probability,
    std::vector<uint64_t> const & symbol_to_score) -> void {
  assert(log_handle != nullptr);
  auto const qmin = static_cast<int>(find_smallest(stats.quality_dist));
  auto const qmax = static_cast<int>(find_largest(stats.quality_dist));

  std::fprintf(log_handle, "\n%s\n%s\n%s\n",
               "Q score distribution",
               "ASCII    Q       Pe           N      Pct   AccPct",
               "-----  ---  -------  ----------  -------  -------");
  uint64_t qual_accum = 0;
  for (auto quality_symbol = qmax ; quality_symbol >= qmin ; --quality_symbol)
    {
      if (stats.quality_dist[quality_symbol] == 0) { continue; }

      qual_accum += stats.quality_dist[quality_symbol];
      std::fprintf(log_handle,
                   "    %c  %3" PRIu64 "  %7.5lf  %10" PRIu64 "  %6.1lf%%  %6.1lf%%\n",
                   quality_symbol,
                   symbol_to_score[quality_symbol],
                   symbol_to_probability[quality_symbol],
                   stats.quality_dist[quality_symbol],
                   100.0 * static_cast<double>(stats.quality_dist[quality_symbol]) / stats.n_symbols,
                   100.0 * static_cast<double>(qual_accum) / stats.n_symbols);
    }
}


// section 3
auto report_length_vs_quality_distribution(std::FILE * log_handle,
                                           struct Stats const & stats) -> void {
  assert(log_handle != nullptr);
  std::fprintf(log_handle, "\n%s\n%s\n",
               "    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct",
               "-----  -------  ----  -------  --------  -----  ---------  --------");

  for (auto length = 2UL; length <= stats.len_max; ++length)
    {
      auto const previous_count = static_cast<double>(stats.length_dist[length - 1]);
      auto const & distribution = stats.distributions[length - 1];
      auto const PctRecs = 100.0 * (stats.n_sequences - previous_count) / stats.n_sequences;
      auto const AvgQ = distribution.avgq;
      auto const AvgP = distribution.avgp;
      auto const AvgEE = distribution.avgee;
      auto const Rate = distribution.rate;

      std::fprintf(log_handle,
                   "%5" PRIu64 "  %6.1lf%%  %4.1lf  %7.5lf  %8.6lf  %5.2lf  %9.6lf  %7.3lf%%\n",
                   length,
                   PctRecs,
                   AvgQ,
                   q2p(AvgQ),
                   AvgP,
                   AvgEE,
                   Rate,
                   100.0 * Rate);
    }
}


// section 4
auto report_expected_error_and_length_filtering(std::FILE * log_handle,
                                                struct Stats const & stats,
                                                std::vector<std::array<uint64_t, 4>> const & ee_length_table) -> void {
  assert(log_handle != nullptr);
  std::fprintf(log_handle, "\n%s\n%s\n",
               "    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000",
               "-----  -------  -------  -------  -------  -------  -------  -------  -------");

  std::vector<double> read_percentage;
  read_percentage.reserve(ee_length_table[0].size());
  for (auto length = stats.len_max; length >= 1UL; --length)
    {
      auto const & read_count = ee_length_table[length - 1];
      for (auto const count : read_count) {
        read_percentage.emplace_back(100.0 * static_cast<double>(count) / stats.n_sequences);
      }

      if (read_count[0] != 0)
        {
          std::fprintf(log_handle,
                       "%5" PRIu64 "  %7" PRIu64 "  %7" PRIu64 "  %7" PRIu64 "  %7" PRIu64 "  "
                       "%6.2lf%%  %6.2lf%%  %6.2lf%%  %6.2lf%%\n",
                       length,
                       read_count[0], read_count[1],
                       read_count[2], read_count[3],
                       read_percentage[0], read_percentage[1],
                       read_percentage[2], read_percentage[3]);
        }
      read_percentage.clear();
    }
}


// section 5
auto report_minimum_quality_and_length_filtering(std::FILE * log_handle,
                                                 struct Stats const & stats,
                                                 std::vector<std::array<uint64_t, 4>> const & q_length_table) -> void {
  assert(log_handle != nullptr);
  std::fprintf(log_handle, "\n%s\n%s\n%s\n",
               "Truncate at first Q",
               "  Len     Q=5    Q=10    Q=15    Q=20",
               "-----  ------  ------  ------  ------");
  auto const mid_length = std::max(1UL, stats.len_max / 2);
  std::vector<double> read_percentage;
  read_percentage.reserve(q_length_table[0].size());
  for (auto length = stats.len_max; length >= mid_length; --length)
    {
      for (auto const count : q_length_table[length - 1]) {
        read_percentage.emplace_back(100.0 * static_cast<double>(count) / stats.n_sequences);
      }

      std::fprintf(log_handle, "%5" PRIu64 "  %5.1lf%%  %5.1lf%%  %5.1lf%%  %5.1lf%%\n",
                   length, read_percentage[0], read_percentage[1],
                   read_percentage[2], read_percentage[3]);
      read_percentage.clear();
    }
}


// closing section
auto report_sequence_stats(std::FILE * log_handle, struct Stats const & stats) -> void {
  assert(log_handle != nullptr);
  static constexpr auto a_million = double{1000000};
  auto const n_sequences = static_cast<double>(stats.seq_count);
  std::fprintf(log_handle, "\n%10" PRIu64 "  Recs (%.1lfM), 0 too long\n",
               stats.seq_count, n_sequences / a_million);
  if (stats.seq_count != 0)
    {
      std::fprintf(log_handle, "%10.1lf  Avg length\n", 1.0 * stats.n_symbols / n_sequences);
    }
  std::fprintf(log_handle, "%9.1lfM  Bases\n", stats.n_symbols / a_million);
}


auto fastq_stats(struct Parameters const & parameters) -> void
{
  auto * input_handle = fastq_open(parameters.opt_fastq_stats);

  auto const filesize = fastq_get_size(input_handle);

  progress_init("Reading FASTQ file", filesize);

  auto const symbol_to_score = precompute_quality_scores(parameters);
  auto const symbol_to_probability = precompute_probability_values(parameters);
  std::vector<uint64_t> read_length_table(initial_memory_allocation);
  std::vector<std::array<uint64_t, n_eight_bit_values>> qual_length_table(initial_memory_allocation);
  std::vector<std::array<uint64_t, 4>> ee_length_table(initial_memory_allocation);
  std::vector<std::array<uint64_t, 4>> q_length_table(initial_memory_allocation);
  std::vector<double> sumee_length_table(initial_memory_allocation);

  // note: fastq parsing represents 99% of total wallclock time
  while (fastq_next(input_handle, false, chrmap_upcase_vector.data()))
    {

      /* update length statistics */

      auto const length = fastq_get_sequence_length(input_handle);

      if (length + 1 > read_length_table.size())
        {
          read_length_table.resize(length + 1);
          qual_length_table.resize(length + 1);
          ee_length_table.resize(length + 1);
          q_length_table.resize(length + 1);
          sumee_length_table.resize(length + 1);
        }

      ++read_length_table[length];  // can NOT be derived from qual_length_table


      /* update quality statistics */

      auto * quality_symbols = fastq_get_quality(input_handle);
      auto expected_error = 0.0;
      auto qmin = std::numeric_limits<uint64_t>::max();  // lowest Q value observed so far in this read

      check_minmax_scores(Span{quality_symbols, length}, symbol_to_score, parameters);

      // search for first Q <= 20, increment q_length_table[i][3] for i = 0 to
      // position of first Q
      // search for the 4 positions at once?
      for (auto i = 0UL; i < length; ++i)
        {
          auto const quality_symbol = static_cast<unsigned char>(quality_symbols[i]);
          auto const quality_score = symbol_to_score[quality_symbol];

          ++qual_length_table[i][quality_symbol];

          qmin = std::min(quality_score, qmin);

          // increment quality observations if the current Q > 5, 10, 15, or 20
          std::transform(quality_thresholds.begin(), quality_thresholds.end(),
                         q_length_table[i].begin(), q_length_table[i].begin(),
                         [qmin](uint64_t const threshold, uint64_t current_value) -> uint64_t {
                           return current_value + (qmin > threshold ? 1 : 0);
                         });

          expected_error += symbol_to_probability[quality_symbol];

          sumee_length_table[i] += expected_error;  // can NOT be derived from qual_length_table

          // increment EE observations if the current EE <= 1.0, 0.5, 0.25, or 0.1
          std::transform(ee_thresholds.begin(), ee_thresholds.end(),
                         ee_length_table[i].begin(), ee_length_table[i].begin(),
                         [expected_error](double const threshold, uint64_t current_value) -> uint64_t {
                           return current_value + (expected_error <= threshold ? 1 : 0);
                         });
        }

      progress_update(fastq_get_position(input_handle));
    }
  progress_done();
  fastq_close(input_handle);


  // note: operations below represent 1% of total wallclock time

  /* compute various distributions */

  auto const stats = Stats{
    find_smallest(read_length_table), find_largest(read_length_table),
        compute_number_of_symbols(read_length_table),
        std::accumulate(read_length_table.begin(), read_length_table.end(),
                        std::uint64_t{0}),
        static_cast<double>(std::accumulate(read_length_table.begin(),
                                            read_length_table.end(),
                                            std::uint64_t{0})),
        compute_cumulative_sum(read_length_table),
        compute_distribution_of_quality_symbols(qual_length_table),
        compute_distributions(find_largest(read_length_table), qual_length_table, sumee_length_table, parameters)
  };


  /* print report */

  if (fp_log != nullptr)
    {
      report_read_length_distribution(fp_log, stats, read_length_table);  // section 1
      report_q_score_distribution(fp_log, stats, symbol_to_probability, symbol_to_score);  // section 2
      report_length_vs_quality_distribution(fp_log, stats);  // section 3
      report_expected_error_and_length_filtering(fp_log, stats, ee_length_table);  // section 4
      report_minimum_quality_and_length_filtering(fp_log, stats, q_length_table);  // section 5
      report_sequence_stats(fp_log, stats);  // closing section
    }

  if (not parameters.opt_quiet)
    {
      std::fprintf(stderr, "Read %" PRIu64 " sequences.\n", stats.seq_count);
    }
}
