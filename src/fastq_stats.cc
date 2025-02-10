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
#include <algorithm>  // std::max, std::find_if, std::transform
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <iterator>  // std::distance
#include <limits>
#include <numeric>  // std::partial_sum
#include <string>
#include <vector>


constexpr auto initial_memory_allocation = std::size_t{512};
constexpr std::array<int, 4> quality_thresholds = {5, 10, 15, 20};
constexpr std::array<double, 4> ee_thresholds = { 1.0, 0.5, 0.25, 0.1 };
constexpr auto n_eight_bit_values = std::size_t{256};

struct Distributions {
  double avgq = 0.0;
  double avgp = 0.0;
  double avgee = 0.0;
  double rate = 0.0;
};


auto q2p(double quality_value) -> double {
  static constexpr auto base = 10.0;
  return std::pow(base, -quality_value / base);
}

// refactoring: overload, or template, or reduce types for quality values?
// auto q2p(long int quality_value) -> double {
//   static constexpr auto base = 10.0;
//   return std::pow(base, -quality_value / base);
// }

// auto q2p(int quality_value) -> double {
//   static constexpr auto base = 10.0;
//   return std::pow(base, -quality_value / base);
// }

// auto q2p(uint64_t quality_value) -> double {
//   static constexpr auto base = 10.0;
//   return std::pow(base, -quality_value / base);
// }


auto check_quality_score(struct Parameters const & parameters, int const quality_score) -> void {
  auto const is_in_accepted_range =
      (quality_score >= parameters.opt_fastq_qmin) and
      (quality_score <= parameters.opt_fastq_qmax);

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
  -> uint64_t {
  // total number of nucleotides = sum(read_length * n_reads_with_that_length)
  std::vector<uint64_t> read_lengths(n_reads_per_length.size());
  std::iota(read_lengths.begin(), read_lengths.end(), 0UL);
  return std::inner_product(read_lengths.begin(), read_lengths.end(),
                            n_reads_per_length.begin(), std::uint64_t{0});
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


auto compute_distributions(unsigned int const len_max, std::vector<std::array<uint64_t, n_eight_bit_values>> const & qual_length_table, std::vector<double> const & sumee_length_table, struct Parameters const & parameters) -> std::vector<struct Distributions> {
  std::vector<struct Distributions> distributions(len_max + 1);

  auto const sum_counts = compute_n_symbols_per_length(qual_length_table);
  auto const sum_quality_scores = compute_sum_quality_scores_per_length(qual_length_table, parameters);
  auto const sum_error_probabilities = compute_sum_error_probabilities_per_length(qual_length_table, parameters);

  for (auto i = 0UL; i <= len_max; i++)
    {
      distributions[i].avgq = 1.0 * sum_quality_scores[i] / sum_counts[i];
      distributions[i].avgp = sum_error_probabilities[i] / sum_counts[i];
      distributions[i].avgee = sumee_length_table[i] / sum_counts[i];
      distributions[i].rate = distributions[i].avgee / (i + 1);
  }
  return distributions;
}


auto fastq_stats(struct Parameters const & parameters) -> void
{
  static constexpr auto a_million = double{1000000};
  auto * input_handle = fastq_open(parameters.opt_fastq_stats);

  auto const filesize = fastq_get_size(input_handle);

  progress_init("Reading FASTQ file", filesize);

  std::vector<uint64_t> read_length_table(initial_memory_allocation);
  std::vector<std::array<uint64_t, n_eight_bit_values>> qual_length_table(initial_memory_allocation);
  std::vector<std::array<uint64_t, 4>> ee_length_table(initial_memory_allocation);
  std::vector<std::array<uint64_t, 4>> q_length_table(initial_memory_allocation);
  std::vector<double> sumee_length_table(initial_memory_allocation);
  std::vector<uint64_t> quality_chars(n_eight_bit_values);

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

      ++read_length_table[length];  // refactoring: could be derived from qual_length_table with std::adjacent_difference(v.rbegin(), v.rend(), read_length_table.rbegin());


      /* update quality statistics */

      auto * quality_symbols = fastq_get_quality(input_handle);
      auto expected_error = 0.0;
      auto qmin_this = std::numeric_limits<int>::max();  // lowest Q value observed in this read
      for (auto i = 0UL; i < length; i++)
        {
          auto const quality_symbol = static_cast<int>(quality_symbols[i]);
          int const quality_score = quality_symbol - parameters.opt_fastq_ascii;
          check_quality_score(parameters, quality_score);  // refactoring: perform after parsing the read

          ++quality_chars[quality_symbol];  // refactoring: could be derived from qual_length_table

          ++qual_length_table[i][quality_symbol];

          expected_error += q2p(quality_score);  // refactoring: replace with lookup

          sumee_length_table[i] += expected_error;

          // increment EE observations if the current EE <= 1.0, 0.5, 0.25, or 0.1
          std::transform(ee_thresholds.begin(), ee_thresholds.end(),
                         ee_length_table[i].begin(), ee_length_table[i].begin(),
                         [expected_error](double const threshold, uint64_t current_value) -> uint64_t {
                           return current_value + (expected_error <= threshold ? 1 : 0);
                         });

          qmin_this = std::min(quality_score, qmin_this);

          // increment quality observations if the current Q > 5, 10, 15, or 20
          std::transform(quality_thresholds.begin(), quality_thresholds.end(),
                         q_length_table[i].begin(), q_length_table[i].begin(),
                         [qmin_this](int const threshold, uint64_t current_value) -> uint64_t {
                           return current_value + (qmin_this > threshold ? 1 : 0);
                         });
        }

      progress_update(fastq_get_position(input_handle));
    }
  progress_done();
  fastq_close(input_handle);

  // note: operations below represent 1% of total wallclock time

  /* compute various distributions */

  auto const n_symbols = compute_number_of_symbols(read_length_table);
  auto const seq_count = std::accumulate(read_length_table.begin(), read_length_table.end(), std::uint64_t{0});
  auto const len_min = find_smallest(read_length_table);
  auto const len_max = find_largest(read_length_table);
  auto const qmin = static_cast<int>(find_smallest(quality_chars));
  auto const qmax = static_cast<int>(find_largest(quality_chars));
  auto const length_dist = compute_cumulative_sum(read_length_table);
  auto const distributions = compute_distributions(len_max, qual_length_table, sumee_length_table, parameters);


  /* print report */


  if (fp_log != nullptr)
    {
      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "Read length distribution\n");
      std::fprintf(fp_log, "      L           N      Pct   AccPct\n");
      std::fprintf(fp_log, "-------  ----------  -------  -------\n");
      // refactoring: std::for_each(read_length_table.rbegin(), read_length_table.rend(), [](uint64_t const read_count) { ... });
      for (auto i = len_max; i >= len_min; i--)
        {
          if (read_length_table[i] > 0)
            {
              std::fprintf(fp_log, "%2s%5" PRId64 "  %10" PRIu64 "   %5.1lf%%   %5.1lf%%\n",
                      (i == len_max ? ">=" : "  "),
                      i,
                      read_length_table[i],
                      read_length_table[i] * 100.0 / seq_count,
                      100.0 * (seq_count - (i > 0 ? length_dist[i - 1] : 0)) / seq_count);
            }
          if (i == 0UL) { break; }
        }

      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "Q score distribution\n");
      std::fprintf(fp_log, "ASCII    Q       Pe           N      Pct   AccPct\n");
      std::fprintf(fp_log, "-----  ---  -------  ----------  -------  -------\n");

      int64_t qual_accum = 0;
      for (auto quality_symbol = qmax ; quality_symbol >= qmin ; quality_symbol--)
        {
          if (quality_chars[quality_symbol] > 0)
            {
              qual_accum += quality_chars[quality_symbol];
              std::fprintf(fp_log,
                      "    %c  %3" PRId64 "  %7.5lf  %10" PRIu64 "  %6.1lf%%  %6.1lf%%\n",
                      quality_symbol,
                      quality_symbol - parameters.opt_fastq_ascii,
                      q2p(quality_symbol - parameters.opt_fastq_ascii),
                      quality_chars[quality_symbol],
                      100.0 * quality_chars[quality_symbol] / n_symbols,
                      100.0 * qual_accum / n_symbols);
            }
        }

      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct\n");
      std::fprintf(fp_log, "-----  -------  ----  -------  --------  -----  ---------  --------\n");

      for (auto i = 2UL; i <= len_max; i++)
        {
          auto const & distribution = distributions[i - 1];
          auto const PctRecs = 100.0 * (seq_count - length_dist[i - 1]) / seq_count;
          auto const AvgQ = distribution.avgq;
          auto const AvgP = distribution.avgp;
          auto const AvgEE = distribution.avgee;
          auto const Rate = distribution.rate;

          std::fprintf(fp_log,
                  "%5" PRId64 "  %6.1lf%%  %4.1lf  %7.5lf  %8.6lf  %5.2lf  %9.6lf  %7.3lf%%\n",
                  i,
                  PctRecs,
                  AvgQ,
                  q2p(AvgQ),
                  AvgP,
                  AvgEE,
                  Rate,
                  100.0 * Rate);
        }

      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000\n");
      std::fprintf(fp_log, "-----  -------  -------  -------  -------  -------  -------  -------  -------\n");

      for (auto i = len_max; i >= 1UL; i--)
        {
          std::array<int64_t, 4> read_count {{}};
          std::array<double, 4> read_percentage {{}};

          for (auto j = 0; j < 4; j++)
            {
              read_count[j] = ee_length_table[i - 1][j];
              read_percentage[j] = 100.0 * read_count[j] / seq_count;
            }

          if (read_count[0] > 0)
            {
              std::fprintf(fp_log,
                      "%5" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  "
                      "%6.2lf%%  %6.2lf%%  %6.2lf%%  %6.2lf%%\n",
                      i,
                      read_count[0], read_count[1],
                      read_count[2], read_count[3],
                      read_percentage[0], read_percentage[1],
                      read_percentage[2], read_percentage[3]);
            }
        }


      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "Truncate at first Q\n");
      std::fprintf(fp_log, "  Len     Q=5    Q=10    Q=15    Q=20\n");
      std::fprintf(fp_log, "-----  ------  ------  ------  ------\n");

      auto const mid_length = std::max(1UL, len_max / 2);
      for (auto i = len_max; i >= mid_length; i--)
        {
          std::array<double, 4> read_percentage {{}};

          for (auto j = 0; j < 4; j++)
            {
              read_percentage[j] = 100.0 * q_length_table[i - 1][j] / seq_count;
            }

          std::fprintf(fp_log, "%5" PRId64 "  %5.1lf%%  %5.1lf%%  %5.1lf%%  %5.1lf%%\n",
                  i,
                  read_percentage[0], read_percentage[1],
                  read_percentage[2], read_percentage[3]);
        }

      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "%10" PRIu64 "  Recs (%.1lfM), 0 too long\n",
              seq_count, seq_count / a_million);
      if (seq_count > 0)
        {
          std::fprintf(fp_log, "%10.1lf  Avg length\n", 1.0 * n_symbols / seq_count);
        }
      std::fprintf(fp_log, "%9.1lfM  Bases\n", n_symbols / a_million);
    }

  if (not parameters.opt_quiet)
    {
      std::fprintf(stderr, "Read %" PRIu64 " sequences.\n", seq_count);
    }
}
