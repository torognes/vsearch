/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <algorithm>  // std::max, std::find_if
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <limits>
#include <numeric>  // std::partial_sum
#include <string>
#include <vector>


constexpr auto initial_memory_allocation = std::size_t{512};
constexpr std::array<int, 4> quality_thresholds = {5, 10, 15, 20};
constexpr std::array<double, 4> ee_limits = { 1.0, 0.5, 0.25, 0.1 };


auto q2p(double quality_value) -> double {
  static constexpr auto base = 10.0;
  return std::pow(base, -quality_value / base);
}


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


auto is_observed = [](uint64_t const count) { return count != 0UL; };


auto find_smallest(std::vector<uint64_t> const & observables) -> unsigned long {
  assert(observables.size() != 0U);
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
  assert(observables.size() != 0U);
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


auto fastq_stats(struct Parameters const & parameters) -> void
{
  static constexpr auto n_eight_bit_values = std::size_t{256};
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

      ++read_length_table[length];


      /* update quality statistics */

      auto * quality_symbols = fastq_get_quality(input_handle);
      auto expected_error = 0.0;
      auto qmin_this = std::numeric_limits<int>::max();
      for (auto i = 0UL; i < length; i++)
        {
          auto const quality_symbol = static_cast<int>(quality_symbols[i]);

          int const quality_score = quality_symbol - parameters.opt_fastq_ascii;
          check_quality_score(parameters, quality_score);

          ++quality_chars[quality_symbol];

          ++qual_length_table[i][quality_symbol];

          expected_error += q2p(quality_score);

          sumee_length_table[i] += expected_error;

          for (auto j = 0; j < 4; j++)
            {
              if (expected_error <= ee_limits[j])
                {
                  ++ee_length_table[i][j];
                }
              else
                {
                  break;
                }
            }

          qmin_this = std::min(quality_score, qmin_this);

          for (auto j = 0; j < 4; j++)
            {
              if (qmin_this > quality_thresholds[j])  // Q > 5, 10, 15, 20
                {
                  ++q_length_table[i][j];
                }
              else
                {
                  break;
                }
            }
        }

      progress_update(fastq_get_position(input_handle));
    }
  progress_done();
  fastq_close(input_handle);

  /* compute various distributions */

  auto const symbols = compute_number_of_symbols(read_length_table);
  auto const seq_count = std::accumulate(read_length_table.begin(), read_length_table.end(), std::uint64_t{0});
  auto const len_min = find_smallest(read_length_table);
  auto const len_max = find_largest(read_length_table);
  auto const qmin = static_cast<int>(find_smallest(quality_chars));
  auto const qmax = static_cast<int>(find_largest(quality_chars));
  auto const length_dist = compute_cumulative_sum(read_length_table);
  std::vector<double> rate_dist(len_max + 1);
  std::vector<double> avgq_dist(len_max + 1);
  std::vector<double> avgee_dist(len_max + 1);
  std::vector<double> avgp_dist(len_max + 1);
  // std::vector<struct Stats> distributions(len_max + 1);  // refactoring above vectors

  for (auto i = 0UL; i <= len_max; i++)
    {
      int64_t sum_counts = 0;
      int64_t sum_quality_scores = 0;
      auto sum_error_probabilities = 0.0;
      for (auto quality_symbol = qmin; quality_symbol <= qmax; quality_symbol++)
        {
          int const quality_score = quality_symbol - parameters.opt_fastq_ascii;
          int const count = qual_length_table[i][quality_symbol];
          sum_counts += count;
          sum_quality_scores += count * quality_score;
          sum_error_probabilities += count * q2p(quality_score);
        }
      avgq_dist[i] = 1.0 * sum_quality_scores / sum_counts;
      avgp_dist[i] = sum_error_probabilities / sum_counts;
      avgee_dist[i] = sumee_length_table[i] / sum_counts;
      rate_dist[i] = avgee_dist[i] / (i + 1);
    }

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
                      100.0 * quality_chars[quality_symbol] / symbols,
                      100.0 * qual_accum / symbols);
            }
        }

      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct\n");
      std::fprintf(fp_log, "-----  -------  ----  -------  --------  -----  ---------  --------\n");

      for (auto i = 2UL; i <= len_max; i++)
        {
          auto const PctRecs = 100.0 * (seq_count - length_dist[i - 1]) / seq_count;
          auto const AvgQ = avgq_dist[i - 1];
          auto const AvgP = avgp_dist[i - 1];
          auto const AvgEE = avgee_dist[i - 1];
          auto const Rate = rate_dist[i - 1];

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
          std::fprintf(fp_log, "%10.1lf  Avg length\n", 1.0 * symbols / seq_count);
        }
      std::fprintf(fp_log, "%9.1lfM  Bases\n", symbols / a_million);
    }

  if (not parameters.opt_quiet)
    {
      std::fprintf(stderr, "Read %" PRIu64 " sequences.\n", seq_count);
    }
}
