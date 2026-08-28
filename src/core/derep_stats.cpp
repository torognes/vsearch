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

#include "core/derep_stats.hpp"
#include "core/discarded_message.hpp"  // vsearch::print_discarded
#include "vsearch.hpp"  // struct Parameters
#include "utils/print_view.hpp"  // fprint, fprint_integer
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fputs


/* One writer per report, taking the destination as its first argument, so
   that the report_* entry points below are left with nothing but the choice
   of destinations: stderr unless --quiet, and the --log file when one is
   open. Same shape as report_orient() (commands/orient.cpp), write_report()
   (commands/sff_convert.cpp) and stats_message()
   (commands/fastx_syncpairs.cpp); it replaces a lambda per entry point.
   See TBD_20260824_report_destinations.md. */
namespace {

  auto print_input_stats(std::FILE * output_stream,
                         Derep_stats const & stats) -> void
  {
    if (stats.sequencecount > 0)
      {
        fprint_integer(output_stream, stats.nucleotidecount);
        fprint(output_stream, " nt in ");
        fprint_integer(output_stream, stats.sequencecount);
        fprint(output_stream, " seqs, min ");
        fprint_integer(output_stream, stats.shortest);
        fprint(output_stream, ", max ");
        fprint_integer(output_stream, stats.longest);
        fprint(output_stream, ", avg ");
        std::fprintf(output_stream, "%.0f", static_cast<double>(stats.nucleotidecount) * 1.0 / static_cast<double>(stats.sequencecount));
        fprint(output_stream, '\n');
      }
    else
      {
        fprint_integer(output_stream, stats.nucleotidecount);
        fprint(output_stream, " nt in ");
        fprint_integer(output_stream, stats.sequencecount);
        fprint(output_stream, " seqs\n");
      }
  }




  auto print_unique_summary(std::FILE * output_stream,
                            Derep_stats const & stats,
                            double const average,
                            double const median) -> void
  {
    if (stats.clusters < 1)
      {
        fprint(output_stream, "0 unique sequences\n");
      }
    else
      {
        fprint_integer(output_stream, stats.clusters);
        fprint(output_stream, " unique sequences, avg cluster ");
        std::fprintf(output_stream, "%.1lf", average);
        fprint(output_stream, ", median ");
        std::fprintf(output_stream, "%.0f", median);
        fprint(output_stream, ", max ");
        fprint_integer(output_stream, stats.maxsize);
        fprint(output_stream, '\n');
      }
  }


  auto print_selected(std::FILE * output_stream,
                      uint64_t const selected,
                      Derep_stats const & stats) -> void
  {
    fprint_integer(output_stream, selected);
    fprint(output_stream, " uniques written, ");
    fprint_integer(output_stream, stats.clusters - selected);
    fprint(output_stream, " clusters discarded (");
    std::fprintf(output_stream, "%.1f", 100.0 * static_cast<double>(stats.clusters - selected) / static_cast<double>(stats.clusters));
    fprint(output_stream, "%)\n");
  }

}  // anonymous namespace


auto report_input_stats(Derep_stats const & stats,
                        struct Parameters const & parameters) -> void
{
  if (not parameters.opt_quiet)
    {
      print_input_stats(stderr, stats);
    }
  if (parameters.fp_log != nullptr)
    {
      print_input_stats(parameters.fp_log, stats);
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
  /* no --quiet gate: this is a warning, and --quiet suppresses "messages to
     stdout and stderr, except for warnings and error messages" (see
     man/commands/fragments/option_quiet.md, and vsearch::warn() which
     encodes the same contract) */
  vsearch::print_discarded(stderr, option_name, length_limit, discarded);
  if (parameters.fp_log != nullptr)
    {
      vsearch::print_discarded(parameters.fp_log, option_name, length_limit, discarded);
      fprint(parameters.fp_log, '\n');
    }
}


auto report_unique_summary(Derep_stats const & stats,
                           double const average,
                           double const median,
                           struct Parameters const & parameters) -> void
{
  if (not parameters.opt_quiet)
    {
      print_unique_summary(stderr, stats, average, median);
    }
  if (parameters.fp_log != nullptr)
    {
      print_unique_summary(parameters.fp_log, stats, average, median);
      fprint(parameters.fp_log, '\n');
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
  if (not parameters.opt_quiet)
    {
      print_selected(stderr, selected, stats);
    }
  if (parameters.fp_log != nullptr)
    {
      print_selected(parameters.fp_log, selected, stats);
      fprint(parameters.fp_log, '\n');
    }
}
