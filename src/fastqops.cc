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
#include <algorithm>  // std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t
#include <vector>


constexpr auto initial_memory_allocation = 512;


auto fastx_revcomp() -> void
{
  uint64_t buffer_alloc = initial_memory_allocation;
  std::vector<char> seq_buffer(buffer_alloc);
  std::vector<char> qual_buffer(buffer_alloc);

  if ((opt_fastaout == nullptr) && (opt_fastqout == nullptr)) {
    fatal("No output files specified");
  }

  auto * input_handle = fastx_open(opt_fastx_revcomp);

  if (input_handle == nullptr)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  if ((opt_fastqout != nullptr) && ! (input_handle->is_fastq || input_handle->is_empty))
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  auto const filesize = fastx_get_size(input_handle);

  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastqout = nullptr;

  if (opt_fastaout != nullptr)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (fp_fastaout == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout != nullptr)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (fp_fastqout == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (input_handle->is_fastq)
    {
      progress_init("Reading FASTQ file", filesize);
    }
  else
    {
      progress_init("Reading FASTA file", filesize);
    }

  auto count = 0;

  while (fastx_next(input_handle, false, chrmap_no_change))
    {
      ++count;

      /* header */

      auto const hlen = fastx_get_header_length(input_handle);
      auto * header = fastx_get_header(input_handle);
      auto const abundance = fastx_get_abundance(input_handle);


      /* sequence */

      auto const length = fastx_get_sequence_length(input_handle);

      if (length + 1 > buffer_alloc)
        {
          buffer_alloc = length + 1;
          seq_buffer.resize(buffer_alloc);
          qual_buffer.resize(buffer_alloc);
        }

      auto * p = fastx_get_sequence(input_handle);
      reverse_complement(seq_buffer.data(), p, length);


      /* quality values */

      auto * q = fastx_get_quality(input_handle);

      if (fastx_is_fastq(input_handle))
        {
          /* reverse quality values */
          for (uint64_t i = 0; i < length; i++)
            {
              qual_buffer[i] = q[length - 1 - i];
            }
          qual_buffer[length] = 0;
        }

      if (opt_fastaout != nullptr)
        {
          fasta_print_general(fp_fastaout,
                              nullptr,
                              seq_buffer.data(),
                              length,
                              header,
                              hlen,
                              abundance,
                              count,
                              -1.0,
                              -1, -1, nullptr, 0.0);
        }

      if (opt_fastqout != nullptr)
        {
          fastq_print_general(fp_fastqout,
                              seq_buffer.data(),
                              length,
                              header,
                              hlen,
                              qual_buffer.data(),
                              abundance,
                              count,
                              -1.0);
        }

      progress_update(fastx_get_position(input_handle));
    }
  progress_done();

  if (opt_fastaout != nullptr)
    {
      fclose(fp_fastaout);
    }

  if (opt_fastqout != nullptr)
    {
      fclose(fp_fastqout);
    }

  fastx_close(input_handle);
}


auto fastq_convert() -> void
{
  if (opt_fastqout == nullptr) {
    fatal("No output file specified with --fastqout");
  }

  auto * input_handle = fastq_open(opt_fastq_convert);

  if (input_handle == nullptr)
    {
      fatal("Unable to open FASTQ file");
    }

  auto const filesize = fastq_get_size(input_handle);

  std::FILE * fp_fastqout = nullptr;

  fp_fastqout = fopen_output(opt_fastqout);
  if (fp_fastqout == nullptr)
    {
      fatal("Unable to open FASTQ output file for writing");
    }

  progress_init("Reading FASTQ file", filesize);

  auto n_entries = 1;
  static constexpr auto default_expected_error = -1.0;  // refactoring: print no ee value?
  while (fastq_next(input_handle, false, chrmap_no_change))
    {
      /* header */

      auto * header = fastq_get_header(input_handle);
      auto const abundance = fastq_get_abundance(input_handle);

      /* sequence */

      auto const length = fastq_get_sequence_length(input_handle);
      auto * sequence = fastq_get_sequence(input_handle);

      /* convert quality values */

      auto * quality = fastq_get_quality(input_handle);
      for (uint64_t i = 0; i < length; i++)
        {
          int q = quality[i] - opt_fastq_ascii;
          if (q < opt_fastq_qmin)
            {
              fprintf(stderr,
                      "\nFASTQ quality score (%d) below minimum (%" PRId64
                      ") in entry no %" PRIu64
                      " starting on line %" PRIu64 "\n",
                      q,
                      opt_fastq_qmin,
                      fastq_get_seqno(input_handle) + 1,
                      fastq_get_lineno(input_handle));
              fatal("FASTQ quality score too low");
            }
          if (q > opt_fastq_qmax)
            {
              fprintf(stderr,
                      "\nFASTQ quality score (%d) above maximum (%" PRId64
                      ") in entry no %" PRIu64
                      " starting on line %" PRIu64 "\n",
                      q,
                      opt_fastq_qmax,
                      fastq_get_seqno(input_handle) + 1,
                      fastq_get_lineno(input_handle));
              fatal("FASTQ quality score too high");
            }
          q = std::max<int64_t>(q, opt_fastq_qminout);
          q = std::min<int64_t>(q, opt_fastq_qmaxout);
          q += opt_fastq_asciiout;
          q = std::max(q, 33);
          q = std::min(q, 126);
          quality[i] = q;
        }
      quality[length] = 0;

      int const hlen = fastq_get_header_length(input_handle);
      fastq_print_general(fp_fastqout,
                          sequence,
                          length,
                          header,
                          hlen,
                          quality,
                          abundance,
                          n_entries,
                          default_expected_error);  // refactoring: prefer function overload?

      ++n_entries;
      progress_update(fastq_get_position(input_handle));
    }

  progress_done();

  fclose(fp_fastqout);
  fastq_close(input_handle);
}
