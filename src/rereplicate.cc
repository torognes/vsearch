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
#include "utils/maps.hpp"
#include <cinttypes> // macros PRIu64 and PRId64
#include <cstdio>  // std::FILE, std::fprintf
#include <cstdint>  // int64_t


auto rereplicate(struct Parameters & parameters) -> void
{
  if (parameters.opt_output == nullptr) {
    fatal("FASTA output file for rereplicate must be specified with --output");
  }

  auto * fp_output = fopen_output(parameters.opt_output);
  if (fp_output == nullptr) {
    fatal("Unable to open FASTA output file for writing");
  }

  opt_xsize = true;
  parameters.opt_xsize = true;
  fastx_handle file_handle = fasta_open(parameters.opt_rereplicate);
  auto const filesize = static_cast<int64_t>(fasta_get_size(file_handle));

  progress_init("Rereplicating", filesize);

  int64_t n_amplicons = 0;
  int64_t missing = 0;
  int64_t n_reads = 0;
  auto const truncateatspace = not parameters.opt_notrunclabels;
  while (fasta_next(file_handle, truncateatspace, chrmap_no_change_vector.data()))
    {
      ++n_amplicons;
      int64_t abundance = fasta_get_abundance_and_presence(file_handle);
      if (abundance == 0)
        {
          ++missing;
          abundance = 1;
        }

      for(int64_t i = 0; i < abundance; ++i)
        {
          ++n_reads;
          if (parameters.opt_output != nullptr)
            {
              fasta_print_general(fp_output,
                                  nullptr,
                                  fasta_get_sequence(file_handle),
                                  static_cast<int>(fasta_get_sequence_length(file_handle)),
                                  fasta_get_header(file_handle),
                                  static_cast<int>(fasta_get_header_length(file_handle)),
                                  1,
                                  static_cast<int>(n_reads),
                                  -1.0,
                                  -1, -1, nullptr, 0.0);
            }
        }

      progress_update(fasta_get_position(file_handle));
    }
  progress_done();

  if (not parameters.opt_quiet)
    {
      if (missing != 0)
        {
          std::fprintf(stderr, "WARNING: Missing abundance information for some input sequences, assumed 1\n");
        }
      std::fprintf(stderr, "Rereplicated %" PRId64 " reads from %" PRId64 " amplicons\n", n_reads, n_amplicons);
    }

  if (parameters.opt_log != nullptr)
    {
      if (missing != 0)
        {
          std::fprintf(stderr, "WARNING: Missing abundance information for some input sequences, assumed 1\n");
        }
      std::fprintf(fp_log, "Rereplicated %" PRId64 " reads from %" PRId64 " amplicons\n", n_reads, n_amplicons);
    }

  fasta_close(file_handle);
  if (fp_output != nullptr) {
    static_cast<void>(std::fclose(fp_output));
  }
}
