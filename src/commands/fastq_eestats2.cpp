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

#include "commands/fastq_eestats2.hpp"
#include "core/eestats.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "vsearch.hpp"
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include <algorithm>  // std::max, std::min
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <vector>


auto fastq_eestats2(struct Parameters const & parameters) -> void
{
  if (parameters.opt_output == nullptr) {
    fatal("Output file for fastq_eestats2 must be specified with --output");
  }

  /* read the expected-error cutoffs directly from the configured Parameters
     vector (E1). */
  auto const & ee_cutoffs = parameters.opt_ee_cutoffs;
  auto const ee_cutoffs_count = static_cast<int>(ee_cutoffs.size());

  fastx_handle h = fastq_open(parameters.opt_fastq_eestats2, parameters);

  uint64_t const filesize = fastq_get_size(h);

  auto const output_handle = open_optional_output_file(parameters.opt_output, OutputOption{"--output"});
  std::FILE * const fp_output = output_handle.get();


  uint64_t seq_count = 0;
  uint64_t symbols = 0;
  uint64_t longest = 0;

  int len_steps = 0;

  std::vector<uint64_t> count_table;

  {
    Progress progress("Reading FASTQ file", filesize, parameters);
    while (fastq_next(h, false, chrmap_upcase()))
      {
        ++seq_count;

        auto const len = fastq_get_sequence_length(h);
        auto const * q = fastq_get_quality(h);

        /* update length statistics */

        if (len > longest)
          {
            longest = len;
            // parameters.opt_length_cutoffs_longest is an int between 1 and INT_MAX
            int const high = static_cast<int>(std::min(longest, static_cast<uint64_t>(parameters.opt_length_cutoffs_longest)));
            auto const new_len_steps = 1 + std::max(0, ((high - parameters.opt_length_cutoffs_shortest)
                                            / parameters.opt_length_cutoffs_increment));

            if (new_len_steps > len_steps)
              {
                count_table.resize(static_cast<size_t>(new_len_steps) * static_cast<size_t>(ee_cutoffs_count));
                len_steps = new_len_steps;
              }
          }

        /* update quality statistics */

        symbols += len;

        auto ee = 0.0;

        for (uint64_t i = 0; i < len; i++)
          {
            /* quality score */

            auto const qual = std::max(fastq_get_qual_eestats(q[i], parameters), 0);

            auto const pe = q2p(qual);

            ee += pe;

            for (int x = 0; x < len_steps; x++)
              {
                uint64_t const len_cutoff = static_cast<uint64_t>(parameters.opt_length_cutoffs_shortest + (x * parameters.opt_length_cutoffs_increment));
                if (i + 1 == len_cutoff)
                  {
                    for (int y = 0; y < ee_cutoffs_count; y++)
                      {
                        if (ee <= ee_cutoffs[static_cast<size_t>(y)])
                          {
                            ++count_table[((static_cast<size_t>(x) * static_cast<size_t>(ee_cutoffs_count)) + static_cast<size_t>(y))];
                          }
                      }
                  }
              }
          }

        progress.update(fastq_get_position(h));
      }
  }

  std::fprintf(fp_output,
          "%" PRIu64 " reads",
          seq_count);

  if (seq_count > 0)
    {
      std::fprintf(fp_output,
              ", max len %" PRIu64 ", avg %.1f",
              longest, 1.0 * static_cast<double>(symbols) / static_cast<double>(seq_count));
    }
  std::fprintf(fp_output, "\n\n");

  std::fprintf(fp_output, "Length");
  for (int y = 0; y < ee_cutoffs_count; y++)
    {
      std::fprintf(fp_output, "         MaxEE %.2f", ee_cutoffs[static_cast<size_t>(y)]);
    }
  std::fprintf(fp_output, "\n");
  std::fprintf(fp_output, "------");
  for (int y = 0; y < ee_cutoffs_count; y++)
    {
      std::fprintf(fp_output, "   ----------------");
    }
  std::fprintf(fp_output, "\n");

  for (int x = 0; x < len_steps; x++)
    {
      int const len_cutoff = parameters.opt_length_cutoffs_shortest + (x * parameters.opt_length_cutoffs_increment);

      if (len_cutoff > parameters.opt_length_cutoffs_longest)
        {
          break;
        }

      std::fprintf(fp_output, "%6d", len_cutoff);

      for (int y = 0; y < ee_cutoffs_count; y++)
        {
          std::fprintf(fp_output,
                  "   %8" PRIu64 "(%5.1f%%)",
                  count_table[((static_cast<size_t>(x) * static_cast<size_t>(ee_cutoffs_count)) + static_cast<size_t>(y))],
                  100.0 * static_cast<double>(count_table[((static_cast<size_t>(x) * static_cast<size_t>(ee_cutoffs_count)) + static_cast<size_t>(y))]) / static_cast<double>(seq_count));
        }
      std::fprintf(fp_output, "\n");
    }

  if (parameters.fp_log != nullptr)
    {
      std::fprintf(parameters.fp_log,
              "%" PRIu64 " reads, max len %" PRIu64 ", avg %.1f\n\n",
              seq_count, longest, 1.0 * static_cast<double>(symbols) / static_cast<double>(seq_count));

      std::fprintf(parameters.fp_log, "Length");
      for (int y = 0; y < ee_cutoffs_count; y++)
        {
          std::fprintf(parameters.fp_log, "         MaxEE %.2f", ee_cutoffs[static_cast<size_t>(y)]);
        }
      std::fprintf(parameters.fp_log, "\n");
      std::fprintf(parameters.fp_log, "------");
      for (int y = 0; y < ee_cutoffs_count; y++)
        {
          std::fprintf(parameters.fp_log, "   ----------------");
        }
      std::fprintf(parameters.fp_log, "\n");

      for (int x = 0; x < len_steps; x++)
        {
          int const len_cutoff = parameters.opt_length_cutoffs_shortest + (x * parameters.opt_length_cutoffs_increment);

          if (len_cutoff > parameters.opt_length_cutoffs_longest)
            {
              break;
            }

          std::fprintf(parameters.fp_log, "%6d", len_cutoff);

          for (int y = 0; y < ee_cutoffs_count; y++)
            {
              std::fprintf(parameters.fp_log,
                      "   %8" PRIu64 "(%5.1f%%)",
                      count_table[((static_cast<size_t>(x) * static_cast<size_t>(ee_cutoffs_count)) + static_cast<size_t>(y))],
                      100.0 * static_cast<double>(count_table[((static_cast<size_t>(x) * static_cast<size_t>(ee_cutoffs_count)) + static_cast<size_t>(y))]) / static_cast<double>(seq_count));
            }
          std::fprintf(parameters.fp_log, "\n");
        }
    }

  fastq_close(h, parameters);
}
