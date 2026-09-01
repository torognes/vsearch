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

#include "commands/fastx_revcomp.hpp"
#include "utils/span.hpp"
#include "utils/view.hpp"
#include "vsearch.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/reverse_complement.hpp"
#include <algorithm>  // std::reverse_copy
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fclose
#include <vector>


constexpr auto initial_memory_allocation = 512;


auto fastx_revcomp(struct Parameters const & parameters) -> void
{
  /* deliberately declared outside the record loop: the buffers grow to
     the longest record seen and are reused, so after the first few
     records no allocation happens at all (cppcheck's variableScope
     hint would reallocate per record) */
  std::vector<char> seq_buffer(initial_memory_allocation);
  std::vector<char> qual_buffer(initial_memory_allocation);

  if ((parameters.opt_fastaout == nullptr) && (parameters.opt_fastqout == nullptr)) {
    fatal("No output files specified");
  }

  auto input_handle = fastx_open(parameters.opt_fastx_revcomp, parameters);

  // if (input_handle == nullptr)
  //   {
  //     fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
  //   }

  if ((parameters.opt_fastqout != nullptr) && not input_handle->is_fastq_input())
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  auto const filesize = input_handle->get_size();

  auto fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  auto fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});

  {
    int64_t count = 0;  // the ordinal fed to --relabel; int would wrap at 2^31 records
    Progress progress(input_handle->is_fastq_format() ? "Reading FASTQ file" : "Reading FASTA file", filesize, parameters);
    while (input_handle->next(false, chrmap_no_change()))
      {
        ++count;

        /* header */

        auto const header = input_handle->header_view();
        auto const abundance = input_handle->get_abundance();


        /* sequence */

        auto const sequence = input_handle->sequence_view();
        auto const length = sequence.size();

        if (seq_buffer.size() < length + 1)
          {
            seq_buffer.resize(length + 1);
          }
        if (qual_buffer.size() < length + 1)
          {
            qual_buffer.resize(length + 1);
          }

        reverse_complement(make_span(seq_buffer).first(length + 1), sequence);


        /* quality values */

        /* reverse quality values (the view is empty for fasta input) */
        auto const quality = input_handle->quality_view();
        std::reverse_copy(quality.cbegin(), quality.cend(), qual_buffer.begin());
        qual_buffer[quality.size()] = '\0';

        if (parameters.opt_fastaout != nullptr)
          {
            fasta_print_general(fastaout_handle.get(),
                                nullptr,
                                make_view(seq_buffer).first(length),
                                header,
                                OutputAnnotations{static_cast<uint64_t>(abundance), count},
                                parameters);
          }

        if (parameters.opt_fastqout != nullptr)
          {
            fastq_print_general(fastqout_handle.get(),
                                make_view(seq_buffer).first(length),
                                header,
                                make_view(qual_buffer).first(length),
                                OutputAnnotations{static_cast<uint64_t>(abundance), count},
                                parameters);
          }

        progress.update(input_handle->get_position());
      }
  }

  /* close in declaration order at a defined point (destructors would run in
     reverse, changing the flush order when both streams share stdout);
     reset() is a no-op on an empty handle, so unopened outputs need no
     guard. */
  fastaout_handle.reset();
  fastqout_handle.reset();

  input_handle->report_stripped_warning(parameters);
}
