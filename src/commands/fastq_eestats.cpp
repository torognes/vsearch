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

#include "commands/fastq_eestats.hpp"
#include "core/eestats.hpp"
#include "utils/quality_table.hpp"  // vsearch::QualityTable
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "vsearch.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include <algorithm>  // std::max, std::min, std::sort
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <limits>
#include <utility>  // std::pair, std::move
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* Sparse histogram of expected-error bins for one read position.

     Its predecessor, one std::map<int64_t, uint64_t> per position, already
     kept memory proportional to the bins actually observed (E12), but paid a
     tree walk through pointer-chased nodes plus a node allocation per base:
     gprof attributes 96% of the command's runtime to that loop on a real
     MiSeq file. Same sparse behaviour, contiguous storage: add() appends to
     a small buffer, and once the buffer has grown as large as the merged
     histogram it is sorted and folded in with one linear pass, so a base
     costs an append plus O(log buffer) comparisons in cache-resident
     memory, and a bin costs 16 bytes instead of a tree node. */
  class BinHistogram {
  public:
    auto add(int64_t const bin) -> void {
      pending_.push_back(bin);
      if (pending_.size() >= std::max(min_buffer_length, bins_.size())) {
        merge_pending();
      }
    }

    /* (bin, count) pairs, sorted by bin, each count non-zero; the reader's
       ascending walk over e_int is thus the same as with the former
       std::map. Not const: remaining pending values are folded in first. */
    auto bins() -> std::vector<std::pair<int64_t, uint64_t>> const & {
      merge_pending();
      return bins_;
    }

  private:
    auto merge_pending() -> void {
      if (pending_.empty()) {
        return;
      }
      std::sort(pending_.begin(), pending_.end());
      std::vector<std::pair<int64_t, uint64_t>> merged;
      merged.reserve(bins_.size() + pending_.size());
      auto old_bin = bins_.cbegin();
      auto new_value = pending_.cbegin();
      while (new_value != pending_.cend()) {
        for (; (old_bin != bins_.cend()) and (old_bin->first < *new_value); ++old_bin) {
          merged.push_back(*old_bin);
        }
        auto const value = *new_value;
        auto count = uint64_t{0};
        for (; (new_value != pending_.cend()) and (*new_value == value); ++new_value) {
          ++count;
        }
        if ((old_bin != bins_.cend()) and (old_bin->first == value)) {
          count += old_bin->second;
          ++old_bin;
        }
        merged.emplace_back(value, count);
      }
      merged.insert(merged.end(), old_bin, bins_.cend());
      bins_ = std::move(merged);
      pending_.clear();  // keeps its capacity, so no reallocation next round
    }

    static constexpr auto min_buffer_length = std::size_t{1024};
    std::vector<std::pair<int64_t, uint64_t>> bins_;
    std::vector<int64_t> pending_;
  };

}  // end of anonymous namespace


auto fastq_eestats(struct Parameters const & parameters) -> void
{
  if (parameters.opt_output == nullptr) {
    fatal("Output file for fastq_eestats must be specified with --output");
  }

  auto h = fastq_open(parameters.opt_fastq_eestats, parameters);

  uint64_t const filesize = h->get_size();

  auto const output_handle = open_optional_output_file(parameters.opt_output, OutputOption{"--output"});
  std::FILE * const fp_output = output_handle.get();


  uint64_t seq_count = 0;

  int64_t len_alloc = 10;

  const int resolution = 1000;
  /* rows of qual_length_table are indexed by the raw quality value (0..qmax),
     so size them by qmax rather than (qmax - qmin): subtracting qmin here
     while indexing by the unshifted value overflowed the row for qmin >= 2 */
  int const max_quality = static_cast<int>(parameters.opt_fastq_qmax + 1);

  std::vector<uint64_t> read_length_table(static_cast<size_t>(len_alloc));
  std::vector<uint64_t> qual_length_table(static_cast<size_t>(len_alloc * (max_quality + 1)));
  /* Sparse per-position expected-error histogram: ee_length_table[pos] holds
     one (e_int bin, count) pair per observed bin. This replaces a dense
     triangular table of size ~resolution*len^2/2 (almost entirely zeros) that
     OOM'd on long reads; memory scales with the cells actually observed, and
     BinHistogram::bins() keeps the reader's ascending walk over e_int
     unchanged (E12). */
  std::vector<BinHistogram> ee_length_table(static_cast<size_t>(len_alloc));
  std::vector<double> sum_ee_length_table(static_cast<size_t>(len_alloc));
  std::vector<double> sum_pe_length_table(static_cast<size_t>(len_alloc));

  int64_t len_min = std::numeric_limits<long>::max();

  /* one 10^-(q/10) per quality symbol, not per base; the cap at 1.0 is the
     std::max(qual, 0) the decoded quality value still goes through below */
  vsearch::QualityTable const quality_table(static_cast<int>(parameters.opt_fastq_ascii),
                                            vsearch::ProbabilityCap::certain_error);

  /* per-symbol decoded scores and range verdicts; a rejected symbol goes back
     through fastq_get_qual_eestats for the unchanged fatal message */
  vsearch::QualityScoreTable const score_table(parameters);

  int64_t len_max = 0;

  {
    Progress progress("Reading FASTQ file", filesize, parameters);
    while (h->next(false, chrmap_upcase()))
      {
        ++seq_count;

        auto const len = static_cast<int64_t>(h->get_sequence_length());
        char const * q = h->get_quality();

        /* update length statistics */

        int64_t const new_alloc = len + 1;

        if (new_alloc > len_alloc)
          {
            read_length_table.resize(static_cast<size_t>(new_alloc));
            qual_length_table.resize(static_cast<size_t>(new_alloc * (max_quality + 1)));
            ee_length_table.resize(static_cast<size_t>(new_alloc));
            sum_ee_length_table.resize(static_cast<size_t>(new_alloc));
            sum_pe_length_table.resize(static_cast<size_t>(new_alloc));

            len_alloc = new_alloc;
          }

        len_min = std::min(len, len_min);
        len_max = std::max(len, len_max);

        /* update quality statistics */

        double ee = 0.0;

        for (int64_t i = 0; i < len; i++)
          {
            ++read_length_table[static_cast<size_t>(i)];

            /* quality score */

            if (not score_table.accepts(q[i]))
              {
                static_cast<void>(fastq_get_qual_eestats(q[i], parameters));  // fatal
              }
            auto const qual = score_table.score(q[i]);
            ++qual_length_table[static_cast<size_t>(((max_quality + 1) * i) + qual)];


            /* probability of error (Pe) */

            auto const probability_of_error = quality_table[q[i]];
            sum_pe_length_table[static_cast<size_t>(i)] += probability_of_error;


            /* expected number of errors */

            ee += probability_of_error;

            auto const e_int = std::min<int64_t>(resolution * (i + 1), static_cast<int>(resolution * ee));
            ee_length_table[static_cast<size_t>(i)].add(e_int);

            sum_ee_length_table[static_cast<size_t>(i)] += ee;
          }
        progress.update(h->get_position());
      }
  }

  fprint(fp_output, "Pos\tRecs\tPctRecs\t"
                    "Min_Q\tLow_Q\tMed_Q\tMean_Q\tHi_Q\tMax_Q\t"
                    "Min_Pe\tLow_Pe\tMed_Pe\tMean_Pe\tHi_Pe\tMax_Pe\t"
                    "Min_EE\tLow_EE\tMed_EE\tMean_EE\tHi_EE\tMax_EE\n");

  for (int64_t i = 0; i < len_max; i++)
    {
      auto const reads = static_cast<int64_t>(read_length_table[static_cast<size_t>(i)]);
      double const pctrecs = 100.0 * static_cast<double>(reads) / static_cast<double>(seq_count);


      /* q */

      double min_q = -1.0;
      double low_q = -1.0;
      double med_q = -1.0;
      double hi_q  = -1.0;
      double max_q = -1.0;

      double qsum = 0;
      double n = 0;
      for (int q = 0; q <= max_quality; q++)
        {
          auto const x = static_cast<double>(qual_length_table[static_cast<size_t>(((max_quality + 1) * i) + q)]);

          if (x > 0)
            {
              qsum += q * x;
              n += x;

              if (min_q < 0)
                {
                  min_q = q;
                }

              if ((low_q < 0) && (n >= 0.25 * static_cast<double>(reads)))
                {
                  low_q = q;
                }

              if ((med_q < 0) && (n >= 0.50 * static_cast<double>(reads)))
                {
                  med_q = q;
                }

              if ((hi_q < 0)  && (n >= 0.75 * static_cast<double>(reads)))
                {
                  hi_q = q;
                }

              max_q = q;
            }
        }

      double const mean_q = 1.0 * qsum / static_cast<double>(reads);


      /* pe */

      double min_pe = -1.0;
      double low_pe = -1.0;
      double med_pe = -1.0;
      double hi_pe  = -1.0;
      double max_pe = -1.0;

      double pesum = 0;
      n = 0;
      for (int q = max_quality; q >= 0; q--)
        {
          auto const x = static_cast<double>(qual_length_table[static_cast<size_t>(((max_quality + 1) * i) + q)]);

          if (x > 0)
            {
              double const pe = q2p(q);
              pesum += pe * x;
              n += x;

              if (min_pe < 0)
                {
                  min_pe = pe;
                }

              if ((low_pe < 0) && (n >= 0.25 * static_cast<double>(reads)))
                {
                  low_pe = pe;
                }

              if ((med_pe < 0) && (n >= 0.50 * static_cast<double>(reads)))
                {
                  med_pe = pe;
                }

              if ((hi_pe < 0) && (n >= 0.75 * static_cast<double>(reads)))
                {
                  hi_pe = pe;
                }

              max_pe = pe;
            }
        }

      double const mean_pe = 1.0 * pesum / static_cast<double>(reads);


      /* expected errors */

      double min_ee = -1.0;
      double low_ee = -1.0;
      double med_ee = -1.0;
      double hi_ee  = -1.0;
      double max_ee = -1.0;

      /* Walk the observed e_int bins for this position in ascending order
         (bins() returns them sorted); every stored bin has a non-zero count,
         so this is exactly the non-zero subset the former dense loop acted on. */
      n = 0;
      for (auto const & bin : ee_length_table[static_cast<size_t>(i)].bins())
        {
          int64_t const e = bin.first;
          n += static_cast<double>(bin.second);

          if (min_ee < 0)
            {
              min_ee = static_cast<double>(e);
            }

          if ((low_ee < 0) && (n >= 0.25 * static_cast<double>(reads)))
            {
              low_ee = static_cast<double>(e);
            }

          if ((med_ee < 0) && (n >= 0.50 * static_cast<double>(reads)))
            {
              med_ee = static_cast<double>(e);
            }

          if ((hi_ee < 0)  && (n >= 0.75 * static_cast<double>(reads)))
            {
              hi_ee = static_cast<double>(e);
            }

          max_ee = static_cast<double>(e);
        }

      double const mean_ee = sum_ee_length_table[static_cast<size_t>(i)] / static_cast<double>(reads);

      min_ee  = (min_ee  + 0.5) / resolution;
      low_ee  = (low_ee  + 0.5) / resolution;
      med_ee  = (med_ee  + 0.5) / resolution;
      hi_ee   = (hi_ee   + 0.5) / resolution;
      max_ee  = (max_ee  + 0.5) / resolution;

      fprint_integer(fp_output, i + 1);
      fprint(fp_output, '\t');
      fprint_integer(fp_output, reads);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.1lf", pctrecs);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.1lf", min_q);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.1lf", low_q);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.1lf", med_q);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.1lf", mean_q);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.1lf", hi_q);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.1lf", max_q);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lg", min_pe);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lg", low_pe);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lg", med_pe);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lg", mean_pe);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lg", hi_pe);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lg", max_pe);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lf", min_ee);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lf", low_ee);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lf", med_ee);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lf", mean_ee);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lf", hi_ee);
      fprint(fp_output, '\t');
      std::fprintf(fp_output, "%.2lf", max_ee);
      fprint(fp_output, '\n');
    }

  h->report_stripped_warning(parameters);
}
