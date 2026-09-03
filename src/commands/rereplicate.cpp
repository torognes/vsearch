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
#include <memory>  // std::unique_ptr
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/warn.hpp"  // vsearch::warn
#include <cstdio>  // std::FILE
#include <cstdint>  // int64_t


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* The two counts travel together and are both int64_t, so they are named in
     a struct rather than passed as two adjacent swappable arguments. No
     default member initializers: they would make this a non-aggregate before
     C++14, and the call site brace-initializes it. */
  struct RereplicationCounts
  {
    int64_t reads;
    int64_t amplicons;
  };


  auto print_rereplicated(std::FILE * output_stream,
                          RereplicationCounts const counts) -> void
  {
    fprint(output_stream, "Rereplicated ");
    fprint_integer(output_stream, counts.reads);
    fprint(output_stream, " reads from ");
    fprint_integer(output_stream, counts.amplicons);
    fprint(output_stream, " amplicons\n");
  }

}  // end of anonymous namespace


auto rereplicate(struct Parameters const & parameters) -> void
{
  auto const output_handle = open_mandatory_output_file(parameters.opt_output, OutputOption{"--output"});
  auto input_handle = fasta_open(parameters.input_filename, parameters);
  auto const filesize = static_cast<int64_t>(input_handle->get_size());

  int64_t n_amplicons = 0;
  auto missing_abundance = false;
  int64_t n_reads = 0;
  auto const truncateatspace = not parameters.opt_notrunclabels;
  {
    Progress progress("Rereplicating", static_cast<uint64_t>(filesize), parameters);
    while (input_handle->next(truncateatspace, chrmap_no_change()))
      {
        ++n_amplicons;
        auto abundance = input_handle->get_abundance_and_presence();
        if (abundance == 0)
          {
            missing_abundance = true;
            abundance = 1;
          }

        for (int64_t i = 0; i < abundance; ++i)
          {
            ++n_reads;
            fasta_print_general(output_handle.get(),
                                nullptr,
                                input_handle->record(),
                                OutputAnnotations{1, n_reads},
                                parameters);
          }

        progress.update(input_handle->get_position());
      }
  }

  /* Outside the --quiet block below, and emitted once for both destinations:
     warnings survive --quiet by documented contract, and the shared reporter
     writes stderr and the log itself. */
  if (missing_abundance)
    {
      vsearch::warn("Missing abundance information for some input sequences, assumed 1");
    }

  auto const totals = RereplicationCounts{n_reads, n_amplicons};

  if (not parameters.opt_quiet)
    {
      print_rereplicated(stderr, totals);
    }

  if (parameters.fp_log != nullptr)
    {
      print_rereplicated(parameters.fp_log, totals);
    }

  input_handle->report_stripped_warning(parameters);
}
