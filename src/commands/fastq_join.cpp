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

#include "commands/fastq_join.hpp"
#include "vsearch.hpp"
#include <memory>  // std::unique_ptr
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"  // fasta_print_general, fasta_get_abundance
#include "core/fastq.hpp"  // fastq_open, fastq_get_sequence, fastq_get_quality
#include "core/fastx.hpp"  // fastx_handle
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/grow_to_fit.hpp"  // vsearch::grow_to_fit
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/span.hpp"  // make_span
#include "utils/view.hpp"
#include <cassert>
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE
#include <iterator>  // std::reverse_iterator
#include <string>
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct input_file {
    char * name = nullptr;
    std::unique_ptr<fastx_s> handle;
  };

  struct input_files {
    input_file forward;
    input_file reverse;
  };

  struct output_file {
    char * name = nullptr;
    OutputFileHandle handle;
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
      infiles.forward.handle = fastq_open(infiles.forward.name, parameters);
    }
    if (infiles.reverse.name != nullptr) {
      infiles.reverse.handle = fastq_open(infiles.reverse.name, parameters);
    }
    return infiles;
  }


  auto open_output_files(struct Parameters const & parameters) -> struct output_files {
    struct output_files outfiles;
    outfiles.fasta.name = parameters.opt_fastaout;
    outfiles.fastq.name = parameters.opt_fastqout;
    outfiles.fasta.handle = open_optional_output_file(outfiles.fasta.name, OutputOption{"--fastaout"});
    outfiles.fastq.handle = open_optional_output_file(outfiles.fastq.name, OutputOption{"--fastqout"});
    return outfiles;
  }


  auto close_output_files(struct output_files & outfiles) -> void {
    /* reset in this fixed order (scope-exit destruction runs in reverse) so
       that outputs sharing stdout flush as they did before RAII */
    outfiles.fasta.handle.reset();
    outfiles.fastq.handle.reset();
  }


  auto close_input_files(struct input_files const & infiles, struct Parameters const & parameters) -> void {
    /* emit the stripped-character warning for each open input; the handles
       themselves are released when 'infiles' (which owns them) is destroyed. */
    if (infiles.forward.handle != nullptr) {
      infiles.forward.handle->report_stripped_warning(parameters);
    }
    if (infiles.reverse.handle != nullptr) {
      infiles.reverse.handle->report_stripped_warning(parameters);
    }
  }


  auto stats_message(std::FILE * output_stream,
                     uint64_t const total) -> void {
    fprint_integer(output_stream, total);
    fprint(output_stream, " pairs joined\n");
  }
}  // end of anonymous namespace


auto fastq_join(struct Parameters const & parameters) -> void
{

  /* check parameters */

  check_parameters(parameters);

  /* open and check input and output files */

  auto const infiles = open_input_files(parameters);
  // check_input_files(infiles)? already done by the function fastq_open()
  auto outfiles = open_output_files(parameters);

  /* main */

  auto const filesize = infiles.forward.handle->get_size();

  /* do it */

  constexpr auto bufferlength = 1024U;
  auto const padlen = parameters.opt_join_padgap.length();
  uint64_t total = 0;
  std::string final_sequence;
  final_sequence.reserve(bufferlength + padlen + bufferlength);
  std::string final_quality;
  final_quality.reserve(final_sequence.capacity());
  /* scratch for the reverse-complemented reverse read; the shared helper
     terminates it with a '\0', hence the extra byte. Deliberately hoisted out
     of the read loop so its allocation is reused across records. */
  std::vector<char> rc_buffer(bufferlength + 1);

  {
    Progress progress("Joining reads", filesize, parameters);
    while (infiles.forward.handle->next(false, chrmap_no_change()))
      {
        if (not infiles.reverse.handle->next(false, chrmap_no_change()))
          {
            fatal("More forward reads than reverse reads");
          }

        auto const fwd_sequence = infiles.forward.handle->sequence_view();
        auto const rev_sequence = infiles.reverse.handle->sequence_view();
        auto const fwd_seq_length = fwd_sequence.size();
        auto const rev_seq_length = rev_sequence.size();
        auto const needed = fwd_seq_length + padlen + rev_seq_length;
        /* parsed from the forward header once; both writers below reuse it */
        auto const abundance = static_cast<uint64_t>(infiles.forward.handle->get_abundance());

        /* reverse read: reverse-complement sequence (the shared one-pass
           helper; the former in-place std::reverse + std::transform called
           the cross-TU map_complement() once per base) */

        vsearch::grow_to_fit(rc_buffer, rev_seq_length);
        reverse_complement(make_span(rc_buffer), rev_sequence);

        /* join them: forward read, pad gap, reverse-complemented reverse
           read -- appended in place, where building the concatenation from
           std::string temporaries re-allocated on every record. The reverse
           quality is appended back to front, which is all its reversal is. */

        final_sequence.assign(fwd_sequence.cbegin(), fwd_sequence.cend());
        final_sequence.append(parameters.opt_join_padgap);
        final_sequence.append(rc_buffer.data(), rev_seq_length);

        auto const fwd_quality = infiles.forward.handle->quality_view();
        auto const rev_quality = infiles.reverse.handle->quality_view();
        final_quality.assign(fwd_quality.cbegin(), fwd_quality.cend());
        final_quality.append(parameters.opt_join_padgapq);
        final_quality.append(std::reverse_iterator<char const *>{rev_quality.cend()},
                             std::reverse_iterator<char const *>{rev_quality.cbegin()});

        assert(final_sequence.size() == needed);
        assert(final_quality.size() == needed);

        /* write output */

        if (parameters.opt_fastqout != nullptr)
          {
            fastq_print_general(outfiles.fastq.handle.get(),
                                make_view(final_sequence).first(needed),
                                infiles.forward.handle->header_view(),
                                make_view(final_quality).first(needed),
                                OutputAnnotations{abundance,
                                                  static_cast<int64_t>(total + 1)},
                                parameters);
          }

        if (parameters.opt_fastaout != nullptr)
          {
            fasta_print_general(outfiles.fasta.handle.get(),
                                nullptr,
                                make_view(final_sequence).first(needed),
                                infiles.forward.handle->header_view(),
                                OutputAnnotations{abundance,
                                                  static_cast<int64_t>(total + 1)},
                                parameters);
          }

        ++total;
        progress.update(infiles.forward.handle->get_position());
      }
  }

  if (infiles.reverse.handle->next(false, chrmap_no_change()))
    {
      fatal("More reverse reads than forward reads");
    }

  if (not parameters.opt_quiet) {
    stats_message(stderr, total);
  }
  if (parameters.fp_log != nullptr) {
    stats_message(parameters.fp_log, total);
  }

  /* clean up */

  close_output_files(outfiles);
  close_input_files(infiles, parameters);
}
