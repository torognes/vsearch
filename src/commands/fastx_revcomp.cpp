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
#include "vsearch.hpp"
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/reverse_complement.hpp"
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fclose, std::size_t
#include <vector>


constexpr auto initial_memory_allocation = 512;


auto fastx_revcomp(struct Parameters const & parameters) -> void
{
  uint64_t buffer_alloc = initial_memory_allocation;
  std::vector<char> seq_buffer(buffer_alloc);
  std::vector<char> qual_buffer(buffer_alloc);

  if ((parameters.opt_fastaout == nullptr) && (parameters.opt_fastqout == nullptr)) {
    fatal("No output files specified");
  }

  auto * input_handle = fastx_open(parameters.opt_fastx_revcomp, parameters);

  // if (input_handle == nullptr)
  //   {
  //     fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
  //   }

  if ((parameters.opt_fastqout != nullptr) && ! (input_handle->is_fastq || input_handle->is_empty))
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  auto const filesize = fastx_get_size(input_handle);

  auto fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  auto fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
  std::FILE * const fp_fastaout = fastaout_handle.get();
  std::FILE * const fp_fastqout = fastqout_handle.get();

  auto count = 0;

  {
    Progress progress(input_handle->is_fastq ? "Reading FASTQ file" : "Reading FASTA file", filesize, parameters);
    while (fastx_next(input_handle, false, chrmap_no_change()))
      {
        ++count;

        /* header */

        auto const hlen = fastx_get_header_length(input_handle);
        auto const * header = fastx_get_header(input_handle);
        auto const abundance = fastx_get_abundance(input_handle);


        /* sequence */

        auto const length = fastx_get_sequence_length(input_handle);

        if (length + 1 > buffer_alloc)
          {
            buffer_alloc = length + 1;
            seq_buffer.resize(buffer_alloc);
            qual_buffer.resize(buffer_alloc);
          }

        auto const * p = fastx_get_sequence(input_handle);
        reverse_complement(seq_buffer.data(), p, static_cast<int64_t>(length));


        /* quality values */

        auto const * q = fastx_get_quality(input_handle);

        if (fastx_is_fastq(input_handle))
          {
            /* reverse quality values */
            for (uint64_t i = 0; i < length; i++)
              {
                qual_buffer[i] = q[length - 1 - i];
              }
            qual_buffer[length] = 0;
          }

        if (parameters.opt_fastaout != nullptr)
          {
            fasta_print_general(fp_fastaout,
                                nullptr,
                                seq_buffer.data(),
                                static_cast<int>(length),
                                header,
                                static_cast<int>(hlen),
                                static_cast<uint64_t>(abundance),
                                count,
                                -1.0,
                                -1, -1, nullptr, 0.0,
                                0,
                                parameters);
          }

        if (parameters.opt_fastqout != nullptr)
          {
            fastq_print_general(fp_fastqout,
                                seq_buffer.data(),
                                static_cast<int>(length),
                                header,
                                static_cast<int>(hlen),
                                qual_buffer.data(),
                                static_cast<uint64_t>(abundance),
                                count,
                                -1.0,
                                parameters);
          }

        progress.update(fastx_get_position(input_handle));
      }
  }

  if (parameters.opt_fastaout != nullptr)
    {
      fastaout_handle.reset();
    }

  if (parameters.opt_fastqout != nullptr)
    {
      fastqout_handle.reset();
    }

  fastx_close(input_handle, parameters);
}
