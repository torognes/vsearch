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
#include <algorithm>  // std::transform
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <string>


struct input_file {
  char * name = nullptr;
  fastx_handle handle = nullptr;
};

struct input_files {
  input_file forward;
  input_file reverse;
};

struct output_file {
  char * name = nullptr;
  std::FILE * handle = nullptr;
};

struct output_files {
  output_file fasta;
  output_file fastq;
};


auto check_parameters(struct Parameters const & parameters) -> void {
  if (parameters.opt_reverse == nullptr) {
    fatal("No reverse reads file specified with --reverse");
  }

  if ((parameters.opt_fastqout == nullptr) and (parameters.opt_fastaout == nullptr)) {
    fatal("No output files specified");
  }

  if (parameters.opt_join_padgap.length() != parameters.opt_join_padgapq.length()) {
    fatal("Strings given by --join_padgap and --join_padgapq differ in length");
  }
}


auto open_input_files(struct Parameters const & parameters) -> struct input_files {
  struct input_files infiles;
  infiles.forward.name = parameters.opt_fastq_join;
  infiles.reverse.name = parameters.opt_reverse;
  if (infiles.forward.name != nullptr) {
    infiles.forward.handle = fastq_open(infiles.forward.name);
  }
  if (infiles.reverse.name != nullptr) {
    infiles.reverse.handle = fastq_open(infiles.reverse.name);
  }
  return infiles;
}


auto open_output_files(struct Parameters const & parameters) -> struct output_files {
  struct output_files outfiles;
  outfiles.fasta.name = parameters.opt_fastaout;
  outfiles.fastq.name = parameters.opt_fastqout;
  if (outfiles.fasta.name != nullptr) {
    outfiles.fasta.handle = fopen_output(outfiles.fasta.name);
  }
  if (outfiles.fastq.name != nullptr) {
    outfiles.fastq.handle = fopen_output(outfiles.fastq.name);
  }
  return outfiles;
}


auto check_output_files(struct output_files const & outfiles) -> void {
  if (outfiles.fasta.name != nullptr) {
    if (outfiles.fasta.handle == nullptr) {
      fatal("Unable to open file for writing (%s)", outfiles.fasta.name);
    }
  }
  if (outfiles.fastq.name != nullptr) {
    if (outfiles.fastq.handle == nullptr) {
      fatal("Unable to open file for writing (%s)", outfiles.fastq.name);
    }
  }
}


auto close_output_files(struct output_files const & outfiles) -> void {
  for (auto * fp_outputfile : {outfiles.fasta.handle, outfiles.fastq.handle}) {
    if (fp_outputfile != nullptr) {
      static_cast<void>(std::fclose(fp_outputfile));
    }
  }
}


auto close_input_files(struct input_files const & infiles) -> void {
  for (auto * fp_inputfile : {infiles.forward.handle, infiles.reverse.handle}) {
    if (fp_inputfile != nullptr) {
        fastq_close(fp_inputfile);
    }
  }
}


auto stats_message(std::FILE * output_stream,
                   uint64_t const total) -> void {
  static_cast<void>(std::fprintf(output_stream,
                                 "%" PRIu64 " pairs joined\n",
                                 total));
}


auto output_stats_message(struct Parameters const & parameters,
                          uint64_t const total,
                          char const * log_filename) -> void {
  if (log_filename == nullptr) {
    return;
  }
  stats_message(parameters.fp_log, total);
}


auto output_stats_message(struct Parameters const & parameters,
                          uint64_t const total) -> void {
  if (parameters.opt_quiet) {
    return;
  }
  stats_message(stderr, total);
}


auto fastq_join(struct Parameters const & parameters) -> void
{

  /* check parameters */

  check_parameters(parameters);

  /* open and check input and output files */

  auto const infiles = open_input_files(parameters);
  // check_input_files(infiles)? already done by the function fastq_open()
  auto const outfiles = open_output_files(parameters);
  check_output_files(outfiles);

  /* main */

  auto const filesize = fastq_get_size(infiles.forward.handle);
  progress_init("Joining reads", filesize);

  /* do it */

  constexpr auto bufferlength = 1024U;
  auto const padlen = parameters.opt_join_padgap.length();
  uint64_t total = 0;
  std::string final_sequence;
  final_sequence.reserve(bufferlength + padlen + bufferlength);
  std::string final_quality;
  final_quality.reserve(final_sequence.capacity());
  std::string reverse_sequence;
  reverse_sequence.reserve(bufferlength);
  std::string reverse_quality;
  reverse_quality.reserve(bufferlength);

  while (fastq_next(infiles.forward.handle, false, chrmap_no_change_array.data()))
    {
      if (not fastq_next(infiles.reverse.handle, false, chrmap_no_change_array.data()))
        {
          fatal("More forward reads than reverse reads");
        }

      final_sequence.clear();
      final_quality.clear();
      reverse_sequence.clear();
      reverse_quality.clear();

      auto const fwd_seq_length = fastq_get_sequence_length(infiles.forward.handle);
      auto const rev_seq_length = fastq_get_sequence_length(infiles.reverse.handle);
      auto const needed = fwd_seq_length + padlen + rev_seq_length;

      /* allocate enough memory */
      if (rev_seq_length > reverse_sequence.capacity()) {
        reverse_sequence.reserve(rev_seq_length);
      }
      if (rev_seq_length > reverse_quality.capacity()) {
        reverse_quality.reserve(rev_seq_length);
      }
      if (needed > final_sequence.capacity()) {
        final_sequence.reserve(needed);
      }

      /* reverse read: reverse-complement sequence */

      reverse_sequence.assign(fastq_get_sequence(infiles.reverse.handle), fwd_seq_length);
      std::reverse(reverse_sequence.begin(), reverse_sequence.end());
      std::transform(reverse_sequence.begin(),
                     reverse_sequence.end(),
                     reverse_sequence.begin(),
                     [](char const & lhs) -> char {
                       auto const unsigned_lhs = static_cast<unsigned char>(lhs);
                       auto const complement_lhs = chrmap_complement_vector[unsigned_lhs];
                       return static_cast<char>(complement_lhs);
                     });

      /* reverse read: reverse quality */

      reverse_quality.assign(fastq_get_quality(infiles.reverse.handle), fwd_seq_length);
      std::reverse(reverse_quality.begin(), reverse_quality.end());

      /* join them */

      final_sequence = std::string{fastq_get_sequence(infiles.forward.handle), fwd_seq_length} + parameters.opt_join_padgap + reverse_sequence;
      final_quality = std::string{fastq_get_quality(infiles.forward.handle), fwd_seq_length} + parameters.opt_join_padgapq + reverse_quality;

      /* write output */

      if (parameters.opt_fastqout != nullptr)
        {
          fastq_print_general(outfiles.fastq.handle,
                              const_cast<char *>(final_sequence.c_str()),
                              static_cast<int>(needed),
                              fastq_get_header(infiles.forward.handle),
                              static_cast<int>(fastq_get_header_length(infiles.forward.handle)),
                              const_cast<char *>(final_quality.c_str()),
                              0,
                              static_cast<int>(total + 1),
                              -1.0);
        }

      if (parameters.opt_fastaout != nullptr)
        {
          fasta_print_general(outfiles.fasta.handle,
                              nullptr,
                              const_cast<char *>(final_sequence.c_str()),
                              static_cast<int>(needed),
                              fastq_get_header(infiles.forward.handle),
                              static_cast<int>(fastq_get_header_length(infiles.forward.handle)),
                              0,
                              static_cast<int>(total + 1),
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0);
        }

      ++total;
      progress_update(fastq_get_position(infiles.forward.handle));
    }

  progress_done();

  if (fastq_next(infiles.reverse.handle, false, chrmap_no_change_array.data()))
    {
      fatal("More reverse reads than forward reads");
    }

  output_stats_message(parameters, total);
  output_stats_message(parameters, total, parameters.opt_log);

  /* clean up */

  close_output_files(outfiles);
  close_input_files(infiles);
}
