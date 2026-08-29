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

#include "commands/fastq_convert.hpp"
#include "vsearch.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fastq.hpp"
#include "core/fastx.hpp"  // fastx_s, byte_range
#include "core/quality_range.hpp"  // vsearch::check_quality_score
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/view.hpp"
#include <algorithm>  // std::max, std::min, std::transform
#include <array>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::size_t
#include <vector>


/* the printable ASCII range a quality character must end up in, whatever
   --fastq_asciiout and the clamps above it asked for */
constexpr auto printable_low = 33;
constexpr auto printable_high = 126;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  constexpr auto n_symbols = std::size_t{byte_range};  // byte_range: core/fastx.hpp

  /* Precompute the whole per-character conversion once: input quality symbol
     -> output quality symbol. A '\0' marks a symbol whose score falls outside
     --fastq_qmin/--fastq_qmax (it cannot collide with a converted symbol,
     which is confined to the printable range above). */
  auto make_quality_mapping(struct Parameters const & parameters)
    -> std::array<char, n_symbols>
  {
    std::array<char, n_symbols> mapping {{}};
    for (auto symbol = std::size_t{0}; symbol < n_symbols; ++symbol)
      {
        /* the same arithmetic as the former per-character conversion: the
           symbol is a char (signed on x86-64 and Windows), so a byte beyond
           127 yields a negative score and is rejected, exactly as before */
        auto const as_char = static_cast<char>(static_cast<unsigned char>(symbol));
        auto const score = static_cast<int64_t>(as_char) - parameters.opt_fastq_ascii;
        if ((score < parameters.opt_fastq_qmin) or (score > parameters.opt_fastq_qmax))
          {
            continue;  // stays '\0': rejected below, with the score recomputed
          }
        auto rebased = std::max(score, parameters.opt_fastq_qminout);
        rebased = std::min(rebased, parameters.opt_fastq_qmaxout);
        rebased += parameters.opt_fastq_asciiout;
        rebased = std::max<int64_t>(rebased, printable_low);
        rebased = std::min<int64_t>(rebased, printable_high);
        mapping[symbol] = static_cast<char>(rebased);
      }
    return mapping;
  }


}  // end of anonymous namespace


auto fastq_convert(struct Parameters const & parameters) -> void
{
  if (parameters.opt_fastqout == nullptr) {
    fatal("No output file specified with --fastqout");
  }

  auto input_handle = fastq_open(parameters.opt_fastq_convert, parameters);

  auto const filesize = input_handle->get_size();

  auto fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
  std::FILE * const fp_fastqout = fastqout_handle.get();


  static constexpr auto default_expected_error = -1.0;  // refactoring: print no ee value?
  auto const quality_mapping = make_quality_mapping(parameters);
  {
    std::vector<char> normalized_quality;
    auto n_entries = 1;
    Progress progress("Reading FASTQ file", filesize, parameters);
    while (input_handle->next(false, chrmap_no_change()))
      {
        /* header */

        auto const * header = input_handle->get_header();
        auto const abundance = input_handle->get_abundance();

        /* sequence */

        auto const length = input_handle->get_sequence_length();
        auto const * sequence = input_handle->get_sequence();

        /* convert quality values */

        /* Rebase each score onto the requested output offset: subtract the
           input offset, refuse anything outside --fastq_qmin/--fastq_qmax,
           clamp to --fastq_qminout/--fastq_qmaxout, add the output offset, and
           clamp to the printable range. All of it is precomputed per symbol in
           quality_mapping (see make_quality_mapping above), so the per-character
           work is a single table load. */
        normalized_quality.resize(length + 1);
        auto const quality = input_handle->quality_view();
        std::transform(quality.begin(), quality.end(), normalized_quality.begin(),
                       [&](char const quality_char) -> char
          {
            auto const mapped = quality_mapping[static_cast<unsigned char>(quality_char)];
            if (mapped == '\0')
              {
                /* This command used to print its own two-part message
                   ("FASTQ quality score (X) below minimum (Y) in entry no N
                   starting on line L", then a shorter fatal()). It was the
                   only site that named the record and the line; the shared
                   message now carries those for every command that can
                   supply them -- see DONE_20260825_quality_range.md phase 5. */
                auto const score = static_cast<int>(quality_char - parameters.opt_fastq_ascii);
                vsearch::check_quality_score(score, parameters,
                                             input_handle->quality_location());
              }
            return mapped;
          });

        int const hlen = static_cast<int>(input_handle->get_header_length());
        OutputAnnotations annotations {static_cast<uint64_t>(abundance), n_entries};
        annotations.expected_error = default_expected_error;
        fastq_print_general(fp_fastqout,
                            View<char>{sequence, static_cast<std::size_t>(length)},
                            View<char>{header, static_cast<std::size_t>(hlen)},
                            make_view(normalized_quality).first(static_cast<std::size_t>(length)),
                            annotations,
                            parameters);

        ++n_entries;
        normalized_quality.clear();
        progress.update(input_handle->get_position());
      }
  }


  fastqout_handle.reset();
  input_handle->report_stripped_warning(parameters);
}
