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
#include "maps.h"
#include <algorithm>  // std::max, std::min
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::memset
#include <limits>
#include <vector>


inline auto fastq_get_qual_eestats(char q) -> int
{
  int const qual = q - opt_fastq_ascii;

  if (qual < opt_fastq_qmin)
    {
      fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
              PRId64 ")\n",
              qual, opt_fastq_qmin);
      if (fp_log)
        {
          fprintf(stderr,
                  "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
                  PRId64 ")\n",
                  qual, opt_fastq_qmin);
        }
      exit(EXIT_FAILURE);
    }
  else if (qual > opt_fastq_qmax)
    {
      fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
              PRId64 ")\n",
              qual, opt_fastq_qmax);
      fprintf(stderr,
              "By default, quality values range from 0 to 41.\n"
              "To allow higher quality values, "
              "please use the option --fastq_qmax %d\n", qual);
      if (fp_log)
        {
          fprintf(fp_log,
                  "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
                  PRId64 ")\n",
                  qual, opt_fastq_qmax);
          fprintf(fp_log,
                  "By default, quality values range from 0 to 41.\n"
                  "To allow higher quality values, "
                  "please use the option --fastq_qmax %d\n", qual);
        }
      exit(EXIT_FAILURE);
    }
  return qual;
}


auto q2p(int quality_value) -> double
{
  static constexpr auto base = 10.0;
  return std::pow(base, -quality_value / base);
}

auto ee_start(int pos, int resolution) -> int64_t
{
  return pos * (resolution * (pos + 1) + 2) / 2;
}


auto fastq_eestats() -> void
{
  if (not opt_output) {
    fatal("Output file for fastq_eestats must be specified with --output");
  }

  fastx_handle h = fastq_open(opt_fastq_eestats);

  uint64_t const filesize = fastq_get_size(h);

  std::FILE * fp_output = nullptr;

  if (opt_output)
    {
      fp_output = fopen_output(opt_output);
      if (not fp_output)
        {
          fatal("Unable to open output file for writing");
        }
    }

  progress_init("Reading FASTQ file", filesize);

  uint64_t seq_count = 0;

  int64_t len_alloc = 10;

  const int resolution = 1000;
  int const max_quality = opt_fastq_qmax - opt_fastq_qmin + 1;

  int64_t ee_size = ee_start(len_alloc, resolution);

  std::vector<uint64_t> read_length_table(len_alloc);
  std::vector<uint64_t> qual_length_table(len_alloc * (max_quality + 1));
  std::vector<uint64_t> ee_length_table(ee_size);
  std::vector<double> sum_ee_length_table(len_alloc);
  std::vector<double> sum_pe_length_table(len_alloc);

  int64_t len_min = std::numeric_limits<long>::max();
  int64_t len_max = 0;

  while (fastq_next(h, false, chrmap_upcase))
    {
      ++seq_count;

      int64_t const len = fastq_get_sequence_length(h);
      char * q = fastq_get_quality(h);

      /* update length statistics */

      int64_t const new_alloc = len + 1;

      if (new_alloc > len_alloc)
        {
          int64_t const new_ee_size = ee_start(new_alloc, resolution);

          read_length_table.resize(new_alloc);
          qual_length_table.resize(new_alloc * (max_quality + 1));
          ee_length_table.resize(new_ee_size);
          sum_ee_length_table.resize(new_alloc);
          sum_pe_length_table.resize(new_alloc);

          len_alloc = new_alloc;
          ee_size = new_ee_size;
        }

      len_min = std::min(len, len_min);
      len_max = std::max(len, len_max);

      /* update quality statistics */

      double ee = 0.0;

      for (int64_t i = 0; i < len; i++)
        {
          ++read_length_table[i];

          /* quality score */

          auto const qual = std::max(fastq_get_qual_eestats(q[i]), 0);
          ++qual_length_table[((max_quality + 1) * i) + qual];


          /* probability of error (Pe) */

          auto const probability_of_error = q2p(qual);
          sum_pe_length_table[i] += probability_of_error;


          /* expected number of errors */

          ee += probability_of_error;

          auto const e_int = std::min<int64_t>(resolution * (i + 1), (int) (resolution * ee));
          ++ee_length_table[ee_start(i, resolution) + e_int];

          sum_ee_length_table[i] += ee;
        }
      progress_update(fastq_get_position(h));
    }
  progress_done();

  fprintf(fp_output,
          "Pos\tRecs\tPctRecs\t"
          "Min_Q\tLow_Q\tMed_Q\tMean_Q\tHi_Q\tMax_Q\t"
          "Min_Pe\tLow_Pe\tMed_Pe\tMean_Pe\tHi_Pe\tMax_Pe\t"
          "Min_EE\tLow_EE\tMed_EE\tMean_EE\tHi_EE\tMax_EE\n");

  for (int64_t i = 0; i < len_max; i++)
    {
      int64_t const reads = read_length_table[i];
      double const pctrecs = 100.0 * reads / seq_count;


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
          double const x = qual_length_table[((max_quality + 1) * i) + q];

          if (x > 0)
            {
              qsum += q * x;
              n += x;

              if (min_q < 0)
                {
                  min_q = q;
                }

              if ((low_q < 0) && (n >= 0.25 * reads))
                {
                  low_q = q;
                }

              if ((med_q < 0) && (n >= 0.50 * reads))
                {
                  med_q = q;
                }

              if ((hi_q < 0)  && (n >= 0.75 * reads))
                {
                  hi_q = q;
                }

              max_q = q;
            }
        }

      double const mean_q = 1.0 * qsum / reads;


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
          double const x = qual_length_table[((max_quality + 1) * i) + q];

          if (x > 0)
            {
              double const pe = q2p(q);
              pesum += pe * x;
              n += x;

              if (min_pe < 0)
                {
                  min_pe = pe;
                }

              if ((low_pe < 0) && (n >= 0.25 * reads))
                {
                  low_pe = pe;
                }

              if ((med_pe < 0) && (n >= 0.50 * reads))
                {
                  med_pe = pe;
                }

              if ((hi_pe < 0) && (n >= 0.75 * reads))
                {
                  hi_pe = pe;
                }

              max_pe = pe;
            }
        }

      double const mean_pe = 1.0 * pesum / reads;


      /* expected errors */

      double min_ee = -1.0;
      double low_ee = -1.0;
      double med_ee = -1.0;
      double hi_ee  = -1.0;
      double max_ee = -1.0;

      int64_t const ee_offset = ee_start(i, resolution);
      int64_t const max_errors = resolution * (i + 1);

      n = 0;
      for (int64_t e = 0; e <= max_errors; e++)
        {
          int64_t const x = ee_length_table[ee_offset + e];

          if (x > 0)
            {
              n += x;

              if (min_ee < 0)
                {
                  min_ee = e;
                }

              if ((low_ee < 0) && (n >= 0.25 * reads))
                {
                  low_ee = e;
                }

              if ((med_ee < 0) && (n >= 0.50 * reads))
                {
                  med_ee = e;
                }

              if ((hi_ee < 0)  && (n >= 0.75 * reads))
                {
                  hi_ee = e;
                }

              max_ee = e;
            }
        }

      double const mean_ee = sum_ee_length_table[i] / reads;

      min_ee  = (min_ee  + 0.5) / resolution;
      low_ee  = (low_ee  + 0.5) / resolution;
      med_ee  = (med_ee  + 0.5) / resolution;
      hi_ee   = (hi_ee   + 0.5) / resolution;
      max_ee  = (max_ee  + 0.5) / resolution;

      fprintf(fp_output,
              "%" PRId64 "\t%" PRId64 "\t%.1lf"
              "\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf"
              "\t%.2lg\t%.2lg\t%.2lg\t%.2lg\t%.2lg\t%.2lg"
              "\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
              i + 1, reads, pctrecs,
              min_q,  low_q,  med_q,  mean_q,  hi_q,  max_q,
              min_pe, low_pe, med_pe, mean_pe, hi_pe, max_pe,
              min_ee, low_ee, med_ee, mean_ee, hi_ee, max_ee);
    }

  fclose(fp_output);

  fastq_close(h);
}


auto fastq_eestats2() -> void
{
  if (! opt_output) {
    fatal("Output file for fastq_eestats2 must be specified with --output");
  }

  fastx_handle h = fastq_open(opt_fastq_eestats2);

  uint64_t const filesize = fastq_get_size(h);

  std::FILE * fp_output = nullptr;

  if (opt_output)
    {
      fp_output = fopen_output(opt_output);
      if (! fp_output)
        {
          fatal("Unable to open output file for writing");
        }
    }

  progress_init("Reading FASTQ file", filesize);

  uint64_t seq_count = 0;
  uint64_t symbols = 0;
  uint64_t longest = 0;

  int len_steps = 0;

  std::vector<uint64_t> count_table;

  while (fastq_next(h, false, chrmap_upcase))
    {
      ++seq_count;

      uint64_t const len = fastq_get_sequence_length(h);
      char * q = fastq_get_quality(h);

      /* update length statistics */

      if (len > longest)
        {
          longest = len;
          // opt_length_cutoffs_longest is an int between 1 and INT_MAX
          int const high = MIN(longest, (uint64_t) (opt_length_cutoffs_longest));
          int const new_len_steps = 1 + MAX(0, ((high - opt_length_cutoffs_shortest)
                                          / opt_length_cutoffs_increment));

          if (new_len_steps > len_steps)
            {
              count_table.resize(new_len_steps * opt_ee_cutoffs_count);
              len_steps = new_len_steps;
            }
        }

      /* update quality statistics */

      symbols += len;

      double ee = 0.0;

      for (uint64_t i = 0; i < len; i++)
        {
          /* quality score */

          auto const qual = std::max(fastq_get_qual_eestats(q[i]), 0);

          auto const pe = q2p(qual);

          ee += pe;

          for (int x = 0; x < len_steps; x++)
            {
              uint64_t const len_cutoff = opt_length_cutoffs_shortest + (x * opt_length_cutoffs_increment);
              if (i + 1 == len_cutoff)
                {
                  for (int y = 0; y < opt_ee_cutoffs_count; y++)
                    {
                      if (ee <= opt_ee_cutoffs_values[y])
                        {
                          ++count_table[(x * opt_ee_cutoffs_count) + y];
                        }
                    }
                }
            }
        }

      progress_update(fastq_get_position(h));
    }
  progress_done();

  fprintf(fp_output,
          "%" PRIu64 " reads",
          seq_count);

  if (seq_count > 0)
    {
      fprintf(fp_output,
              ", max len %" PRIu64 ", avg %.1f",
              longest, 1.0 * symbols / seq_count);
    }
  fprintf(fp_output, "\n\n");

  fprintf(fp_output, "Length");
  for (int y = 0; y < opt_ee_cutoffs_count; y++)
    {
      fprintf(fp_output, "         MaxEE %.2f", opt_ee_cutoffs_values[y]);
    }
  fprintf(fp_output, "\n");
  fprintf(fp_output, "------");
  for (int y = 0; y < opt_ee_cutoffs_count; y++)
    {
      fprintf(fp_output, "   ----------------");
    }
  fprintf(fp_output, "\n");

  for (int x = 0; x < len_steps; x++)
    {
      int const len_cutoff = opt_length_cutoffs_shortest + (x * opt_length_cutoffs_increment);

      if (len_cutoff > opt_length_cutoffs_longest)
        {
          break;
        }

      fprintf(fp_output, "%6d", len_cutoff);

      for (int y = 0; y < opt_ee_cutoffs_count; y++)
        {
          fprintf(fp_output,
                  "   %8" PRIu64 "(%5.1f%%)",
                  count_table[(x * opt_ee_cutoffs_count) + y],
                  100.0 * count_table[(x * opt_ee_cutoffs_count) + y] / seq_count);
        }
      fprintf(fp_output, "\n");
    }

  if (fp_log)
    {
      fprintf(fp_log,
              "%" PRIu64 " reads, max len %" PRIu64 ", avg %.1f\n\n",
              seq_count, longest, 1.0 * symbols / seq_count);

      fprintf(fp_log, "Length");
      for (int y = 0; y < opt_ee_cutoffs_count; y++)
        {
          fprintf(fp_log, "         MaxEE %.2f", opt_ee_cutoffs_values[y]);
        }
      fprintf(fp_log, "\n");
      fprintf(fp_log, "------");
      for (int y = 0; y < opt_ee_cutoffs_count; y++)
        {
          fprintf(fp_log, "   ----------------");
        }
      fprintf(fp_log, "\n");

      for (int x = 0; x < len_steps; x++)
        {
          int const len_cutoff = opt_length_cutoffs_shortest + (x * opt_length_cutoffs_increment);

          if (len_cutoff > opt_length_cutoffs_longest)
            {
              break;
            }

          fprintf(fp_log, "%6d", len_cutoff);

          for (int y = 0; y < opt_ee_cutoffs_count; y++)
            {
              fprintf(fp_log,
                      "   %8" PRIu64 "(%5.1f%%)",
                      count_table[(x * opt_ee_cutoffs_count) + y],
                      100.0 * count_table[(x * opt_ee_cutoffs_count) + y] / seq_count);
            }
          fprintf(fp_log, "\n");
        }
    }

  fclose(fp_output);

  fastq_close(h);
}
