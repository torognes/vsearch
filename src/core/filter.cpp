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

#include "core/filter.hpp"
#include "vsearch.hpp"
#include <memory>  // std::unique_ptr
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/view.hpp"  // View<char>
#include <algorithm>  // std::count_if, std::min, std::max
#include <cstddef>  // std::size_t
#include <cmath>  // std::pow, std::signbit
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE
#include <limits>
#include <string>  // std::string, std::to_string


namespace {
inline auto fastq_get_qual(char const quality_symbol, struct Parameters const & parameters) -> int
{
  int const quality_score = quality_symbol - static_cast<int>(parameters.opt_fastq_ascii);

  // route through fatal() (which writes stderr and the --log file, and in a
  // library session throws instead of std::exit()ing); the printed text is
  // unchanged.
  if (quality_score < parameters.opt_fastq_qmin)
    {
      fatal("FASTQ quality value (" + std::to_string(quality_score) + ") below qmin ("
            + std::to_string(parameters.opt_fastq_qmin) + ")");
    }
  else if (quality_score > parameters.opt_fastq_qmax)
    {
      fatal("FASTQ quality value (" + std::to_string(quality_score) + ") above qmax ("
            + std::to_string(parameters.opt_fastq_qmax) + ")\n"
            "By default, quality values range from 0 to 41.\n"
            "To allow higher quality values, "
            "please use the option --fastq_qmax " + std::to_string(quality_score));
    }
  return quality_score;
}
}  // anonymous namespace


struct analysis_res
{
  bool discarded = false;
  bool truncated = false;
  /* the kept part of the record, and the matching part of its quality string
     (empty for a FASTA record). Subspans of the reader's own buffers, so they
     stay valid until the next read; holding them removes the start/length pair
     every consumer had to re-slice the handle with. */
  View<char> sequence;
  View<char> quality;
  double ee = -1.0;
};


namespace {
auto analyse(fastx_handle input_handle, struct Parameters const & parameters) -> struct analysis_res
{
  auto const fastq_trunclen = static_cast<int>(parameters.opt_fastq_trunclen);
  auto const fastq_trunclen_keep = static_cast<int>(parameters.opt_fastq_trunclen_keep);
  struct analysis_res res;
  int start = 0;
  int length = static_cast<int>(input_handle->get_sequence_length());
  auto const old_length = length;

  /* strip left (5') end */
  if (parameters.opt_fastq_stripleft < length)
    {
      start += static_cast<int>(parameters.opt_fastq_stripleft);
      length -= static_cast<int>(parameters.opt_fastq_stripleft);
    }
  else
    {
      start = length;
      length = 0;
    }

  /* strip right (3') end */
  if (parameters.opt_fastq_stripright < length)
    {
      length -= static_cast<int>(parameters.opt_fastq_stripright);
    }
  else
    {
      length = 0;
    }

  /* truncate trailing (3') part */
  if (parameters.opt_fastq_trunclen >= 0)
    {
      length = std::min(length, fastq_trunclen);
    }

  /* truncate trailing (3') part, but keep if short */
  if (parameters.opt_fastq_trunclen_keep >= 0)
    {
      length = std::min(length, fastq_trunclen_keep);
    }

  if (input_handle->is_fastq_format())
    {
      /* truncate by quality and expected errors (ee) */
      res.ee = 0.0;
      static constexpr auto base = 10.0;
      auto const * quality_symbols = input_handle->get_quality() + start;
      for (auto i = 0; i < length; ++i)
        {
          auto const quality_score = fastq_get_qual(quality_symbols[i], parameters);
          auto const expected_error = std::pow(base, -quality_score / base);
          res.ee += expected_error;

          if ((quality_score <= parameters.opt_fastq_truncqual) or
              (res.ee > parameters.opt_fastq_truncee) or
              (res.ee > parameters.opt_fastq_truncee_rate * (i + 1)))
            {
              res.ee -= expected_error;
              length = i;
              break;
            }

          if (quality_score < parameters.opt_fastq_minqual)
            {
              res.discarded = true;
            }
        }

      /* filter by expected errors (ee) */
      if (res.ee > parameters.opt_fastq_maxee)
        {
          res.discarded = true;
        }
      if ((length > 0) and ((res.ee / length) > parameters.opt_fastq_maxee_rate))
        {
          res.discarded = true;
        }
    }

  /* filter by length */
  if ((parameters.opt_fastq_trunclen >= 0) and (length < parameters.opt_fastq_trunclen))
    {
      res.discarded = true;
    }
  if (length < parameters.opt_fastq_minlen)
    {
      res.discarded = true;
    }
  if (length > parameters.opt_fastq_maxlen)
    {
      res.discarded = true;
    }

  /* filter by n's, over the kept part of the sequence: the reader has had a
     sequence_view() since the streaming reader became a class, so the window
     is a subspan rather than a raw get_sequence() + start */
  auto const kept = input_handle->sequence_view()
    .subspan(static_cast<std::size_t>(start), static_cast<std::size_t>(length));
  auto const ncount = std::count_if(kept.begin(), kept.end(),
                                    [](char const nucleotide) -> bool {
                                      return (nucleotide == 'N') or (nucleotide == 'n');
                                    });
  if (ncount > parameters.opt_fastq_maxns)
    {
      res.discarded = true;
    }

  /* filter by abundance */
  auto const abundance = input_handle->get_abundance();
  if (abundance < parameters.opt_minsize)
    {
      res.discarded = true;
    }
  if (abundance > parameters.opt_maxsize)
    {
      res.discarded = true;
    }

  res.truncated = length < old_length;

  /* publish the kept window as views into the reader's buffers */
  auto const window_start = static_cast<std::size_t>(start);
  auto const window_length = static_cast<std::size_t>(length);
  res.sequence = input_handle->sequence_view().subspan(window_start, window_length);
  if (input_handle->is_fastq_format())
    {
      res.quality = input_handle->quality_view().subspan(window_start, window_length);
    }

  return res;
}
}  // anonymous namespace


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // defined further below, in the anonymous namespace
  auto check_parameters(struct Parameters const & parameters) -> void;

}  // end of anonymous namespace


auto filter(bool const fastq_only, char const * filename, struct Parameters const & parameters) -> void
{
  check_parameters(parameters);

  static constexpr auto dbl_max_local = std::numeric_limits<double>::max();  // refactoring: redundant?
  static constexpr auto long_min = std::numeric_limits<long>::min();

  if ((parameters.opt_fastqout == nullptr) and (parameters.opt_fastaout == nullptr) and
      (parameters.opt_fastqout_discarded == nullptr) and (parameters.opt_fastaout_discarded == nullptr) and
      (parameters.opt_fastqout_rev == nullptr) and (parameters.opt_fastaout_rev == nullptr) and
      (parameters.opt_fastqout_discarded_rev == nullptr) and (parameters.opt_fastaout_discarded_rev == nullptr))
    {
      fatal("No output files specified");
    }

  auto forward_handle = fastx_open(filename, parameters);
  std::unique_ptr<fastx_s> reverse_handle;

  if (not forward_handle->is_fastq_input())
    {
      if (fastq_only)
        {
          fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");
      }
      else if (parameters.opt_eeout or (parameters.opt_fastq_ascii != 33) or parameters.opt_fastq_eeout or
               (parameters.opt_fastq_maxee < dbl_max_local) or
               (parameters.opt_fastq_maxee_rate < dbl_max_local) or (parameters.opt_fastqout != nullptr) or
               (parameters.opt_fastq_qmax < 41) or (parameters.opt_fastq_qmin > 0) or
               (parameters.opt_fastq_truncee < dbl_max_local) or
               (parameters.opt_fastq_truncee_rate < dbl_max_local) or
               (parameters.opt_fastq_truncqual < long_min) or
               (parameters.opt_fastq_minqual > 0) or
               (parameters.opt_fastqout_discarded != nullptr) or
               (parameters.opt_fastqout_discarded_rev != nullptr) or
               (parameters.opt_fastqout_rev != nullptr))
        {
          fatal("The following options are not accepted with the fastx_filter command when the input is a FASTA file, because quality scores are not available: eeout, fastq_ascii, fastq_eeout, fastq_maxee, fastq_maxee_rate, fastq_minqual, fastq_out, fastq_qmax, fastq_qmin, fastq_truncee, fastq_truncee_rate, fastq_truncqual,  fastqout_discarded, fastqout_discarded_rev, fastqout_rev");
        }
    }

  auto const filesize = forward_handle->get_size();

  if (parameters.opt_reverse != nullptr)
    {
      reverse_handle = fastx_open(parameters.opt_reverse, parameters);

      if (forward_handle->is_fastq_format() != reverse_handle->is_fastq_format())
        {
          fatal("The forward and reverse input sequence must in the same format, either FASTA or FASTQ");
        }
    }

  OutputFileHandle fp_fastaout;
  OutputFileHandle fp_fastqout;
  OutputFileHandle fp_fastaout_discarded;
  OutputFileHandle fp_fastqout_discarded;

  OutputFileHandle fp_fastaout_rev;
  OutputFileHandle fp_fastqout_rev;
  OutputFileHandle fp_fastaout_discarded_rev;
  OutputFileHandle fp_fastqout_discarded_rev;

  fp_fastaout = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  fp_fastqout = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
  fp_fastaout_discarded = open_optional_output_file(parameters.opt_fastaout_discarded, OutputOption{"--fastaout_discarded"});
  fp_fastqout_discarded = open_optional_output_file(parameters.opt_fastqout_discarded, OutputOption{"--fastqout_discarded"});

  if (reverse_handle != nullptr)
    {
      fp_fastaout_rev = open_optional_output_file(parameters.opt_fastaout_rev, OutputOption{"--fastaout_rev"});
      fp_fastqout_rev = open_optional_output_file(parameters.opt_fastqout_rev, OutputOption{"--fastqout_rev"});
      fp_fastaout_discarded_rev = open_optional_output_file(parameters.opt_fastaout_discarded_rev, OutputOption{"--fastaout_discarded_rev"});
      fp_fastqout_discarded_rev = open_optional_output_file(parameters.opt_fastqout_discarded_rev, OutputOption{"--fastqout_discarded_rev"});
    }

  int64_t kept = 0;
  int64_t discarded = 0;
  int64_t truncated = 0;

  {
    Progress progress("Reading input file", filesize, parameters);
    while (forward_handle->next(false, chrmap_no_change()))
      {
        if ((reverse_handle != nullptr) and not reverse_handle->next(false, chrmap_no_change()))
          {
            fatal("More forward reads than reverse reads");
          }

        struct analysis_res res1;
        res1.ee = 0.0;
        struct analysis_res res2;

        res1 = analyse(forward_handle.get(), parameters);
        if (reverse_handle != nullptr)
          {
            res2 = analyse(reverse_handle.get(), parameters);
          }

        if (res1.discarded or res2.discarded)
          {
            /* discard the sequence(s) */

            ++discarded;

            OutputAnnotations forward_annotations {
              static_cast<uint64_t>(forward_handle->get_abundance()), discarded};
            forward_annotations.expected_error = res1.ee;

            if (parameters.opt_fastaout_discarded != nullptr)
              {
                fasta_print_general(fp_fastaout_discarded.get(),
                                    nullptr,
                                    res1.sequence,
                                    forward_handle->header_view(),
                                    forward_annotations,
                                    parameters);
              }

            if (parameters.opt_fastqout_discarded != nullptr)
              {
                fastq_print_general(fp_fastqout_discarded.get(),
                                    res1.sequence,
                                    forward_handle->header_view(),
                                    res1.quality,
                                    forward_annotations,
                                    parameters);
              }

            if (reverse_handle != nullptr)
              {
                OutputAnnotations reverse_annotations {
                  static_cast<uint64_t>(reverse_handle->get_abundance()), discarded};
                reverse_annotations.expected_error = res2.ee;

                if (parameters.opt_fastaout_discarded_rev != nullptr)
                  {
                    fasta_print_general(fp_fastaout_discarded_rev.get(),
                                        nullptr,
                                        res2.sequence,
                                        reverse_handle->header_view(),
                                        reverse_annotations,
                                        parameters);
                  }

                if (parameters.opt_fastqout_discarded_rev != nullptr)
                  {
                    fastq_print_general(fp_fastqout_discarded_rev.get(),
                                        res2.sequence,
                                        reverse_handle->header_view(),
                                        res2.quality,
                                        reverse_annotations,
                                        parameters);
                  }
              }
          }
        else
          {
            /* keep the sequence(s) */

            ++kept;

            if (res1.truncated or res2.truncated)
              {
                ++truncated;
              }

            OutputAnnotations forward_annotations {
              static_cast<uint64_t>(forward_handle->get_abundance()), kept};
            forward_annotations.expected_error = res1.ee;

            if (parameters.opt_fastaout != nullptr)
              {
                fasta_print_general(fp_fastaout.get(),
                                    nullptr,
                                    res1.sequence,
                                    forward_handle->header_view(),
                                    forward_annotations,
                                    parameters);
              }

            if (parameters.opt_fastqout != nullptr)
              {
                fastq_print_general(fp_fastqout.get(),
                                    res1.sequence,
                                    forward_handle->header_view(),
                                    res1.quality,
                                    forward_annotations,
                                    parameters);
              }

            if (reverse_handle != nullptr)
              {
                OutputAnnotations reverse_annotations {
                  static_cast<uint64_t>(reverse_handle->get_abundance()), kept};
                reverse_annotations.expected_error = res2.ee;

                if (parameters.opt_fastaout_rev != nullptr)
                  {
                    fasta_print_general(fp_fastaout_rev.get(),
                                        nullptr,
                                        res2.sequence,
                                        reverse_handle->header_view(),
                                        reverse_annotations,
                                        parameters);
                  }

                if (parameters.opt_fastqout_rev != nullptr)
                  {
                    fastq_print_general(fp_fastqout_rev.get(),
                                        res2.sequence,
                                        reverse_handle->header_view(),
                                        res2.quality,
                                        reverse_annotations,
                                        parameters);
                  }
              }
          }

        progress.update(forward_handle->get_position());
      }
  }

  if ((reverse_handle != nullptr) and reverse_handle->next(false, chrmap_no_change()))
    {
      fatal("More reverse reads than forward reads");
    }

  if (not parameters.opt_quiet)
    {
      fprint_integer(stderr, kept);
      fprint(stderr, " sequences kept (of which ");
      fprint_integer(stderr, truncated);
      fprint(stderr, " truncated), ");
      fprint_integer(stderr, discarded);
      fprint(stderr, " sequences discarded.\n");
    }

  if (parameters.fp_log != nullptr)
    {
      fprint_integer(parameters.fp_log, kept);
      fprint(parameters.fp_log, " sequences kept (of which ");
      fprint_integer(parameters.fp_log, truncated);
      fprint(parameters.fp_log, " truncated), ");
      fprint_integer(parameters.fp_log, discarded);
      fprint(parameters.fp_log, " sequences discarded.\n");
    }

  if (reverse_handle != nullptr)
    {
      if (parameters.opt_fastaout_rev != nullptr)
        {
          fp_fastaout_rev.reset();
        }

      if (parameters.opt_fastqout_rev != nullptr)
        {
          fp_fastqout_rev.reset();
        }

      if (parameters.opt_fastaout_discarded_rev != nullptr)
        {
          fp_fastaout_discarded_rev.reset();
        }

      if (parameters.opt_fastqout_discarded_rev != nullptr)
        {
          fp_fastqout_discarded_rev.reset();
        }

      reverse_handle->report_stripped_warning(parameters);
    }

  if (parameters.opt_fastaout != nullptr)
    {
      fp_fastaout.reset();
    }

  if (parameters.opt_fastqout != nullptr)
    {
      fp_fastqout.reset();
    }

  if (parameters.opt_fastaout_discarded != nullptr)
    {
      fp_fastaout_discarded.reset();
    }

  if (parameters.opt_fastqout_discarded != nullptr)
    {
      fp_fastqout_discarded.reset();
    }

  forward_handle->report_stripped_warning(parameters);
}


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto check_parameters(struct Parameters const & parameters) -> void {
    static constexpr auto long_min = std::numeric_limits<long>::min();

    auto const is_negative = std::signbit(parameters.opt_fastq_truncee_rate);
    if (is_negative) {
      fatal("--fastq_truncee_rate cannot be negative");
    }

    if (parameters.opt_fastq_minqual < 0) {
      fatal("--fastq_minqual cannot be negative");
    }

    /* Reject out-of-domain values, matching the "positive integer" /
       "real" contract stated in the vsearch manual page. */

    if (parameters.opt_fastq_maxee <= 0.0) {
      /* expected error is the sum of per-base error probabilities;
         probabilities are strictly positive (min quality score is 93,
         giving ~10^-9.3 > 0), so EE > 0. A threshold of 0.0 or below
         would reject every read -- almost certainly a user mistake. */
      fatal("Argument to --fastq_maxee must be positive");
    }

    if (std::signbit(parameters.opt_fastq_maxee_rate)) {
      fatal("Argument to --fastq_maxee_rate cannot be negative");
    }

    if (std::signbit(parameters.opt_fastq_truncee)) {
      fatal("Argument to --fastq_truncee cannot be negative");
    }

    if (parameters.opt_fastq_maxlen < 1) {
      fatal("Argument to --fastq_maxlen must be a positive integer");
    }

    if (parameters.opt_fastq_maxns < 0) {
      fatal("Argument to --fastq_maxns must be a non-negative integer");
    }

    if (parameters.opt_fastq_minlen < 1) {
      fatal("Argument to --fastq_minlen must be a positive integer");
    }

    /* Default sentinel is -1 ("not set"); preserve it but reject any
       other non-positive value. */
    if ((parameters.opt_fastq_trunclen != -1) and (parameters.opt_fastq_trunclen < 1)) {
      fatal("Argument to --fastq_trunclen must be a positive integer");
    }

    if ((parameters.opt_fastq_trunclen_keep != -1) and (parameters.opt_fastq_trunclen_keep < 1)) {
      fatal("Argument to --fastq_trunclen_keep must be a positive integer");
    }

    /* Quality score range: 0..93 (Phred scores encoded in fastq).
       The default value is std::numeric_limits<long>::min(), meaning
       "no truncation"; skip the range check in that case so the
       default is preserved. */
    if ((parameters.opt_fastq_truncqual != long_min) and
        ((parameters.opt_fastq_truncqual < 0) or (parameters.opt_fastq_truncqual > 93))) {
      fatal("Argument to --fastq_truncqual must be in range 0..93");
    }

    if (parameters.opt_fastq_stripleft < 0) {
      fatal("Argument to --fastq_stripleft must be a non-negative integer");
    }

    if (parameters.opt_fastq_stripright < 0) {
      fatal("Argument to --fastq_stripright must be a non-negative integer");
    }
  }

}  // end of anonymous namespace
