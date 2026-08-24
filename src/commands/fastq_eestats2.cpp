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
#include "utils/quality_table.hpp"  // vsearch::QualityTable
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "vsearch.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include <algorithm>  // std::max, std::min
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <initializer_list>
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* what the reading loop accumulated over the whole input */
  struct ReadTotals {
    uint64_t seq_count;
    uint64_t symbols;
    uint64_t longest;
  };


  /* the report goes verbatim to --output and, when set, to --log; when
     there are no reads at all, "max len" and "avg" are omitted from the
     "N reads" line (avg would divide by zero) */
  auto report_eestats2(std::FILE * output_stream,
                       struct Parameters const & parameters,
                       struct ReadTotals const & totals,
                       std::vector<uint64_t> const & count_table,
                       int const len_steps) -> void
  {
    auto const & ee_cutoffs = parameters.opt_ee_cutoffs;
    auto const ee_cutoffs_count = static_cast<int>(ee_cutoffs.size());
    auto const seq_count = totals.seq_count;

    fprint_integer(output_stream, seq_count);
    fprint(output_stream, " reads");

    if (seq_count > 0)
      {
        fprint(output_stream, ", max len ");
        fprint_integer(output_stream, totals.longest);
        fprint(output_stream, ", avg ");
        std::fprintf(output_stream, "%.1f", 1.0 * static_cast<double>(totals.symbols) / static_cast<double>(seq_count));
      }
    fprint(output_stream, "\n\n");

    fprint(output_stream, "Length");
    for (int y = 0; y < ee_cutoffs_count; y++)
      {
        fprint(output_stream, "         MaxEE ");
        std::fprintf(output_stream, "%.2f", ee_cutoffs[static_cast<size_t>(y)]);
      }
    fprint(output_stream, '\n');
    fprint(output_stream, "------");
    for (int y = 0; y < ee_cutoffs_count; y++)
      {
        fprint(output_stream, "   ----------------");
      }
    fprint(output_stream, '\n');

    for (int x = 0; x < len_steps; x++)
      {
        int const len_cutoff = parameters.opt_length_cutoffs_shortest + (x * parameters.opt_length_cutoffs_increment);

        if (len_cutoff > parameters.opt_length_cutoffs_longest)
          {
            break;
          }

        fprint_integer(output_stream, len_cutoff, 6);

        for (int y = 0; y < ee_cutoffs_count; y++)
          {
            auto const count = count_table[((static_cast<size_t>(x) * static_cast<size_t>(ee_cutoffs_count)) + static_cast<size_t>(y))];
            fprint(output_stream, "   ");
            fprint_integer(output_stream, count, 8);
            fprint(output_stream, '(');
            std::fprintf(output_stream, "%5.1f", 100.0 * static_cast<double>(count) / static_cast<double>(seq_count));
            fprint(output_stream, "%)");
          }
        fprint(output_stream, '\n');
      }
  }

}  // end of anonymous namespace


auto fastq_eestats2(struct Parameters const & parameters) -> void
{
  if (parameters.opt_output == nullptr) {
    fatal("Output file for fastq_eestats2 must be specified with --output");
  }

  /* read the expected-error cutoffs directly from the configured Parameters
     vector (E1). */
  auto const & ee_cutoffs = parameters.opt_ee_cutoffs;
  auto const ee_cutoffs_count = static_cast<int>(ee_cutoffs.size());

  auto h = fastq_open(parameters.opt_fastq_eestats2, parameters);

  uint64_t const filesize = h->get_size();

  auto const output_handle = open_optional_output_file(parameters.opt_output, OutputOption{"--output"});
  std::FILE * const fp_output = output_handle.get();


  uint64_t seq_count = 0;
  uint64_t symbols = 0;
  uint64_t longest = 0;

  int len_steps = 0;

  std::vector<uint64_t> count_table;

  /* one 10^-(q/10) per quality symbol, not per base; the cap at 1.0 is the
     std::max(qual, 0) the decoded quality value still goes through below */
  vsearch::QualityTable const quality_table(static_cast<int>(parameters.opt_fastq_ascii),
                                            vsearch::ProbabilityCap::certain_error);

  /* per-symbol range verdicts; a rejected symbol goes back through
     fastq_get_qual_eestats for the unchanged fatal message */
  vsearch::QualityScoreTable const score_table(parameters);


  {
    Progress progress("Reading FASTQ file", filesize, parameters);
    while (h->next(false, chrmap_upcase()))
      {
        ++seq_count;

        auto const len = h->get_sequence_length();
        auto const * q = h->get_quality();

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

        /* the length cutoffs (shortest + x * increment, with increment >= 1,
           enforced at option parsing) strictly increase with x, and so do the
           positions, so one pass through the cutoffs is enough per read: x
           always points at the next cutoff this read can still reach */
        int x = 0;

        for (uint64_t i = 0; i < len; i++)
          {
            /* quality score */

            /* range-check only: this command never uses the decoded value,
               it indexes the table by the symbol. The check still has to run
               per base, so that an out-of-range quality is reported at the
               same position and with the same message as before. */
            if (not score_table.accepts(q[i]))
              {
                static_cast<void>(fastq_get_qual_eestats(q[i], parameters));  // fatal
              }

            auto const pe = quality_table[q[i]];

            ee += pe;

            if (x < len_steps)
              {
                auto const len_cutoff = static_cast<uint64_t>(parameters.opt_length_cutoffs_shortest + (x * parameters.opt_length_cutoffs_increment));
                if (i + 1 == len_cutoff)
                  {
                    for (int y = 0; y < ee_cutoffs_count; y++)
                      {
                        if (ee <= ee_cutoffs[static_cast<size_t>(y)])
                          {
                            ++count_table[((static_cast<size_t>(x) * static_cast<size_t>(ee_cutoffs_count)) + static_cast<size_t>(y))];
                          }
                      }
                    ++x;
                  }
              }
          }

        progress.update(h->get_position());
      }
  }

  ReadTotals const totals = {seq_count, symbols, longest};

  // same report to each requested destination (--log is optional)
  for (auto * output_stream : {fp_output, parameters.fp_log})
    {
      if (output_stream != nullptr)
        {
          report_eestats2(output_stream, parameters, totals, count_table, len_steps);
        }
    }

  h->report_stripped_warning(parameters);
}
