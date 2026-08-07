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
#include "vsearch.hpp"  // struct Parameters
#include "utils/print_view.hpp"  // fprint, fprint_integer
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fputs


auto report_input_stats(Derep_stats const & stats,
                        struct Parameters const & parameters) -> void
{
  auto emit = [&](std::FILE * fp) -> void {
    if (stats.sequencecount > 0)
      {
        fprint_integer(fp, stats.nucleotidecount);
        fprint(fp, " nt in ");
        fprint_integer(fp, stats.sequencecount);
        fprint(fp, " seqs, min ");
        fprint_integer(fp, stats.shortest);
        fprint(fp, ", max ");
        fprint_integer(fp, stats.longest);
        fprint(fp, ", avg ");
        std::fprintf(fp, "%.0f", static_cast<double>(stats.nucleotidecount) * 1.0 / static_cast<double>(stats.sequencecount));
        fprint(fp, '\n');
      }
    else
      {
        fprint_integer(fp, stats.nucleotidecount);
        fprint(fp, " nt in ");
        fprint_integer(fp, stats.sequencecount);
        fprint(fp, " seqs\n");
      }
  };
  if (not parameters.opt_quiet)
    {
      emit(stderr);
    }
  if (parameters.opt_log != nullptr)
    {
      emit(parameters.fp_log);
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
  auto emit = [&](std::FILE * fp) -> void {
    std::fputs(option_name, fp);
    fprint(fp, ' ');
    fprint_integer(fp, length_limit);
    fprint(fp, ": ");
    fprint_integer(fp, discarded);
    fprint(fp, ' ');
    std::fputs((discarded == 1 ? "sequence" : "sequences"), fp);
    fprint(fp, " discarded.\n");
  };
  emit(stderr);
  if (parameters.opt_log != nullptr)
    {
      emit(parameters.fp_log);
      fprint(parameters.fp_log, '\n');
    }
}


auto report_unique_summary(Derep_stats const & stats,
                           double const average,
                           double const median,
                           struct Parameters const & parameters) -> void
{
  auto emit = [&](std::FILE * fp) -> void {
    if (stats.clusters < 1)
      {
        fprint(fp, "0 unique sequences\n");
      }
    else
      {
        fprint_integer(fp, stats.clusters);
        fprint(fp, " unique sequences, avg cluster ");
        std::fprintf(fp, "%.1lf", average);
        fprint(fp, ", median ");
        std::fprintf(fp, "%.0f", median);
        fprint(fp, ", max ");
        fprint_integer(fp, stats.maxsize);
        fprint(fp, '\n');
      }
  };
  if (not parameters.opt_quiet)
    {
      emit(stderr);
    }
  if (parameters.opt_log != nullptr)
    {
      emit(parameters.fp_log);
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
  auto emit = [&](std::FILE * fp) -> void {
    fprint_integer(fp, selected);
    fprint(fp, " uniques written, ");
    fprint_integer(fp, stats.clusters - selected);
    fprint(fp, " clusters discarded (");
    std::fprintf(fp, "%.1f", 100.0 * static_cast<double>(stats.clusters - selected) / static_cast<double>(stats.clusters));
    fprint(fp, "%)\n");
  };
  if (not parameters.opt_quiet)
    {
      emit(stderr);
    }
  if (parameters.opt_log != nullptr)
    {
      emit(parameters.fp_log);
      fprint(parameters.fp_log, '\n');
    }
}
