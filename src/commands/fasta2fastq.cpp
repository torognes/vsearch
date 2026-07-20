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
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/progress.hpp"
#include <cassert>
#include <cstdint>
#include <cstdio>  // std::FILE, std::size_t, std::fclose
#include <vector>


auto fasta2fastq(struct Parameters const & parameters) -> void
{
  auto const max_ascii_value = static_cast<char>(parameters.opt_fastq_asciiout + parameters.opt_fastq_qmaxout);

  auto fp_input = fasta_open(parameters.opt_fasta2fastq, parameters);
  assert(fp_input != nullptr);  // check performed in fasta_open(fastx_open())

  auto const output_handle = open_mandatory_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
  assert(parameters.opt_fastqout != nullptr);  // check performed above

  static constexpr auto initial_length = 1024U;
  std::vector<char> quality(initial_length, max_ascii_value);

  Progress progress("Converting FASTA file to FASTQ",
                    fp_input->get_size(),
                    parameters);

  auto counter = 0;
  while (fp_input->next(false, chrmap_no_change()))
    {
      /* get sequence length and allocate more mem if necessary */

      auto const length = fp_input->get_sequence_length();

      if (quality.size() < length + 1)
        {
          quality.resize(length + 1, max_ascii_value);
        }

      // note: adding '\0' and the end of the quality string is not necessary,
      // fastq_print_general() uses 'length' for both sequence and quality

      ++counter;

      /* write to fastq file */
      fastq_print_general(output_handle.get(),
                          fp_input->get_sequence(),
                          static_cast<int>(length),
                          fp_input->get_header(),
                          static_cast<int>(fp_input->get_header_length()),
                          quality.data(),
                          static_cast<uint64_t>(fp_input->get_abundance()),
                          counter,
                          -1.0,
                          parameters);

      progress.update(fp_input->get_position());
    }

  fp_input->report_stripped_warning(parameters);
}
