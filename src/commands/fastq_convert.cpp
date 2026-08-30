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
#include <cmath>  // std::log10, std::lround, std::pow
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

  /* A Solexa score is not a Phred score: Solexa is -10 log10(p / (1 - p)),
     Phred is -10 log10(p). The two share the offset 64 and nothing else, so
     --fastq_solexa converts the score before any output clamp is applied
     (see TBD_20260825_solexa_conversion.md). The formula is

         Q_phred = 10 log10(10^(Q_solexa / 10) + 1)

     taken from Biopython's phred_quality_from_solexa(); an independent
     implementation reproduces all five of that function's doctest values.

     The mapping is the identity from Solexa 10 upward, so all of the
     substance is in the fifteen symbols below it -- which are exactly the
     ones a quality filter acts on. It is lossy on purpose: six Solexa pairs
     collapse onto one Phred score each ({-5,-4} -> 1, {-3,-2} -> 2,
     {-1,0} -> 3, {1,2} -> 4, {3,4} -> 5, {9,10} -> 10). That is inherent to
     the two scales, not to this implementation, and it is why there is no
     reverse conversion.

     There are no exact .5 ties anywhere in the legal range (-5..62), so
     round-half-even, round-half-up and std::lround all agree on every input;
     the choice of rounding rule cannot move an output symbol here. */
  auto solexa_to_phred(int64_t const solexa_score) noexcept -> int64_t
  {
    static constexpr auto phred_scale_factor = 10.0;  // the 10 in -10 log10(p)
    auto const odds =
      std::pow(phred_scale_factor,
               static_cast<double>(solexa_score) / phred_scale_factor);
    /* std::lround returns long: identical to int64_t where long is 64 bits,
       a widening conversion where it is 32 (mingw), never a narrowing one */
    return std::lround(phred_scale_factor * std::log10(odds + 1.0));
  }


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
        /* after the --fastq_qmin/--fastq_qmax test, which is stated on the
           Solexa score (that is what --fastq_qmin -5 means), and before the
           output clamps, or the six collapsing pairs land in the wrong
           place */
        auto const converted = parameters.opt_fastq_solexa
          ? solexa_to_phred(score) : score;
        auto rebased = std::max(converted, parameters.opt_fastq_qminout);
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
