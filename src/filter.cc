/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <algorithm>  // std::min, std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <limits>


inline auto fastq_get_qual(char const quality_symbol) -> int
{
  int const quality_score = quality_symbol - opt_fastq_ascii;

  if (quality_score < opt_fastq_qmin)
    {
      std::fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
              PRId64 ")\n",
              quality_score, opt_fastq_qmin);
      if (fp_log != nullptr)
        {
          std::fprintf(stderr,
                  "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
                  PRId64 ")\n",
                  quality_score, opt_fastq_qmin);
        }
      std::exit(EXIT_FAILURE);
    }
  else if (quality_score > opt_fastq_qmax)
    {
      std::fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
              PRId64 ")\n",
              quality_score, opt_fastq_qmax);
      std::fprintf(stderr,
              "By default, quality values range from 0 to 41.\n"
              "To allow higher quality values, "
              "please use the option --fastq_qmax %d\n", quality_score);
      if (fp_log != nullptr)
        {
          std::fprintf(fp_log,
                  "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
                  PRId64 ")\n",
                  quality_score, opt_fastq_qmax);
          std::fprintf(fp_log,
                  "By default, quality values range from 0 to 41.\n"
                  "To allow higher quality values, "
                  "please use the option --fastq_qmax %d\n", quality_score);
        }
      std::exit(EXIT_FAILURE);
    }
  return quality_score;
}


struct analysis_res
{
  bool discarded = false;
  bool truncated = false;
  int start = 0;
  int length = 0;
  double ee = -1.0;
};


auto analyse(fastx_handle input_handle) -> struct analysis_res
{
  auto const fastq_trunclen = static_cast<int>(opt_fastq_trunclen);
  auto const fastq_trunclen_keep = static_cast<int>(opt_fastq_trunclen_keep);
  struct analysis_res res;
  res.length = static_cast<int>(fastx_get_sequence_length(input_handle));
  auto const old_length = res.length;

  /* strip left (5') end */
  if (opt_fastq_stripleft < res.length)
    {
      res.start += opt_fastq_stripleft;
      res.length -= opt_fastq_stripleft;
    }
  else
    {
      res.start = res.length;
      res.length = 0;
    }

  /* strip right (3') end */
  if (opt_fastq_stripright < res.length)
    {
      res.length -= opt_fastq_stripright;
    }
  else
    {
      res.length = 0;
    }

  /* truncate trailing (3') part */
  if (opt_fastq_trunclen >= 0)
    {
      res.length = std::min(res.length, fastq_trunclen);
    }

  /* truncate trailing (3') part, but keep if short */
  if (opt_fastq_trunclen_keep >= 0)
    {
      res.length = std::min(res.length, fastq_trunclen_keep);
    }

  if (input_handle->is_fastq)
    {
      /* truncate by quality and expected errors (ee) */
      res.ee = 0.0;
      static constexpr auto base = 10.0;
      auto * quality_symbols = fastx_get_quality(input_handle) + res.start;
      for (auto i = 0; i < res.length; ++i)
        {
          auto const quality_score = fastq_get_qual(quality_symbols[i]);
          auto const expected_error = std::pow(base, -quality_score / base);
          res.ee += expected_error;

          if ((quality_score <= opt_fastq_truncqual) or
              (res.ee > opt_fastq_truncee))
            {
              res.ee -= expected_error;
              res.length = i;
              break;
            }
        }

      /* filter by expected errors (ee) */
      if (res.ee > opt_fastq_maxee)
        {
          res.discarded = true;
        }
      if ((res.length > 0) and ((res.ee / res.length) > opt_fastq_maxee_rate))
        {
          res.discarded = true;
        }
    }

  /* filter by length */
  if ((opt_fastq_trunclen >= 0) and (res.length < opt_fastq_trunclen))
    {
      res.discarded = true;
    }
  if (res.length < opt_fastq_minlen)
    {
      res.discarded = true;
    }
  if (res.length > opt_fastq_maxlen)
    {
      res.discarded = true;
    }

  /* filter by n's */  // refactoring: std::count_if()
  int64_t ncount = 0;
  auto * nucleotides = fastx_get_sequence(input_handle) + res.start;
  for (auto i = 0; i < res.length; ++i)
    {
      auto const nucleotide = nucleotides[i];
      if ((nucleotide == 'N') or (nucleotide == 'n'))
        {
          ++ncount;
        }
    }
  if (ncount > opt_fastq_maxns)
    {
      res.discarded = true;
    }

  /* filter by abundance */
  auto const abundance = fastx_get_abundance(input_handle);
  if (abundance < opt_minsize)
    {
      res.discarded = true;
    }
  if (abundance > opt_maxsize)
    {
      res.discarded = true;
    }

  res.truncated = res.length < old_length;

  return res;
}


auto filter(bool const fastq_only, char * filename) -> void
{
  static constexpr auto dbl_max = std::numeric_limits<double>::max();
  static constexpr auto long_min = std::numeric_limits<long>::min();

  if ((opt_fastqout == nullptr) and (opt_fastaout == nullptr) and
      (opt_fastqout_discarded == nullptr) and (opt_fastaout_discarded == nullptr) and
      (opt_fastqout_rev == nullptr) and (opt_fastaout_rev == nullptr) and
      (opt_fastqout_discarded_rev == nullptr) and (opt_fastaout_discarded_rev == nullptr))
    {
      fatal("No output files specified");
    }

  auto * forward_handle = fastx_open(filename);
  fastx_handle reverse_handle = nullptr;  // refactoring: direct initialization

  if (forward_handle == nullptr)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  if (not (forward_handle->is_fastq or forward_handle->is_empty))
    {
      if (fastq_only)
        {
          fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");
        }
      else if (opt_eeout or
               (opt_fastq_ascii != 33) or
               opt_fastq_eeout or
               (opt_fastq_maxee < dbl_max) or
               (opt_fastq_maxee_rate < dbl_max) or
               (opt_fastqout != nullptr) or
               (opt_fastq_qmax < 41) or
               (opt_fastq_qmin > 0) or
               (opt_fastq_truncee < dbl_max) or
               (opt_fastq_truncqual < long_min) or
               (opt_fastqout_discarded != nullptr) or
               (opt_fastqout_discarded_rev != nullptr) or
               (opt_fastqout_rev != nullptr))
        {
          fatal("The following options are not accepted with the fastx_filter command when the input is a FASTA file, because quality scores are not available: eeout, fastq_ascii, fastq_eeout, fastq_maxee, fastq_maxee_rate, fastq_out, fastq_qmax, fastq_qmin, fastq_truncee, fastq_truncqual,  fastqout_discarded, fastqout_discarded_rev, fastqout_rev");
        }
    }

  auto const filesize = fastx_get_size(forward_handle);

  if (opt_reverse != nullptr)
    {
      reverse_handle = fastx_open(opt_reverse);

      if (reverse_handle == nullptr)
        {
          fatal("Unrecognized file type (not proper FASTA or FASTQ format) for reverse reads");
        }

      if (forward_handle->is_fastq != reverse_handle->is_fastq)
        {
          fatal("The forward and reverse input sequence must in the same format, either FASTA or FASTQ");
        }

      if (not (reverse_handle->is_fastq or reverse_handle->is_empty))
        {
          if (fastq_only)
            {
              fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");
            }
          else if (opt_eeout or
                   (opt_fastq_ascii != 33) or
                   opt_fastq_eeout or
                   (opt_fastq_maxee < dbl_max) or
                   (opt_fastq_maxee_rate < dbl_max) or
                   (opt_fastqout != nullptr) or
                   (opt_fastq_qmax < 41) or
                   (opt_fastq_qmin > 0) or
                   (opt_fastq_truncee < dbl_max) or
                   (opt_fastq_truncqual < long_min) or
                   (opt_fastqout_discarded != nullptr) or
                   (opt_fastqout_discarded_rev != nullptr) or
                   (opt_fastqout_rev != nullptr))
            {
              fatal("The following options are not accepted with the fastx_filter command when the input is a FASTA file, because quality scores are not available: eeout, fastq_ascii, fastq_eeout, fastq_maxee, fastq_maxee_rate, fastq_out, fastq_qmax, fastq_qmin, fastq_truncee, fastq_truncqual,  fastqout_discarded, fastqout_discarded_rev, fastqout_rev");
            }
        }
    }

  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_fastaout_discarded = nullptr;
  std::FILE * fp_fastqout_discarded = nullptr;

  std::FILE * fp_fastaout_rev = nullptr;
  std::FILE * fp_fastqout_rev = nullptr;
  std::FILE * fp_fastaout_discarded_rev = nullptr;
  std::FILE * fp_fastqout_discarded_rev = nullptr;

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

  if (opt_fastaout_discarded != nullptr)
    {
      fp_fastaout_discarded = fopen_output(opt_fastaout_discarded);
      if (fp_fastaout_discarded == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout_discarded != nullptr)
    {
      fp_fastqout_discarded = fopen_output(opt_fastqout_discarded);
      if (fp_fastqout_discarded == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (reverse_handle != nullptr)
    {
      if (opt_fastaout_rev != nullptr)
        {
          fp_fastaout_rev = fopen_output(opt_fastaout_rev);
          if (fp_fastaout_rev == nullptr)
            {
              fatal("Unable to open FASTA output file for writing");
            }
        }

      if (opt_fastqout_rev != nullptr)
        {
          fp_fastqout_rev = fopen_output(opt_fastqout_rev);
          if (fp_fastqout_rev == nullptr)
            {
              fatal("Unable to open FASTQ output file for writing");
            }
        }

      if (opt_fastaout_discarded_rev != nullptr)
        {
          fp_fastaout_discarded_rev = fopen_output(opt_fastaout_discarded_rev);
          if (fp_fastaout_discarded_rev == nullptr)
            {
              fatal("Unable to open FASTA output file for writing");
            }
        }

      if (opt_fastqout_discarded_rev != nullptr)
        {
          fp_fastqout_discarded_rev = fopen_output(opt_fastqout_discarded_rev);
          if (fp_fastqout_discarded_rev == nullptr)
            {
              fatal("Unable to open FASTQ output file for writing");
            }
        }
    }

  progress_init("Reading input file", filesize);

  int64_t kept = 0;
  int64_t discarded = 0;
  int64_t truncated = 0;

  while (fastx_next(forward_handle, false, chrmap_no_change))
    {
      if ((reverse_handle != nullptr) and not fastx_next(reverse_handle, false, chrmap_no_change))
        {
          fatal("More forward reads than reverse reads");
        }

      struct analysis_res res1;
      res1.ee = 0.0;
      struct analysis_res res2;

      res1 = analyse(forward_handle);
      if (reverse_handle != nullptr)
        {
          res2 = analyse(reverse_handle);
        }

      if (res1.discarded or res2.discarded)
        {
          /* discard the sequence(s) */

          ++discarded;

          if (opt_fastaout_discarded != nullptr)
            {
              fasta_print_general(fp_fastaout_discarded,
                                  nullptr,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_abundance(forward_handle),
                                  discarded,
                                  res1.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout_discarded != nullptr)
            {
              fastq_print_general(fp_fastqout_discarded,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_quality(forward_handle) + res1.start,
                                  fastx_get_abundance(forward_handle),
                                  discarded,
                                  res1.ee);
            }

          if (reverse_handle != nullptr)
            {
              if (opt_fastaout_discarded_rev != nullptr)
                {
                  fasta_print_general(fp_fastaout_discarded_rev,
                                      nullptr,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_abundance(reverse_handle),
                                      discarded,
                                      res2.ee,
                                      -1,
                                      -1,
                                      nullptr,
                                      0.0);
                }

              if (opt_fastqout_discarded_rev != nullptr)
                {
                  fastq_print_general(fp_fastqout_discarded_rev,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_quality(reverse_handle) + res2.start,
                                      fastx_get_abundance(reverse_handle),
                                      discarded,
                                      res2.ee);
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

          if (opt_fastaout != nullptr)
            {
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_abundance(forward_handle),
                                  kept,
                                  res1.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout != nullptr)
            {
              fastq_print_general(fp_fastqout,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_quality(forward_handle) + res1.start,
                                  fastx_get_abundance(forward_handle),
                                  kept,
                                  res1.ee);
            }

          if (reverse_handle != nullptr)
            {
              if (opt_fastaout_rev != nullptr)
                {
                  fasta_print_general(fp_fastaout_rev,
                                      nullptr,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_abundance(reverse_handle),
                                      kept,
                                      res2.ee,
                                      -1,
                                      -1,
                                      nullptr,
                                      0.0);
                }

              if (opt_fastqout_rev != nullptr)
                {
                  fastq_print_general(fp_fastqout_rev,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_quality(reverse_handle) + res2.start,
                                      fastx_get_abundance(reverse_handle),
                                      kept,
                                      res2.ee);
                }
            }
        }

      progress_update(fastx_get_position(forward_handle));
    }

  progress_done();

  if ((reverse_handle != nullptr) and fastx_next(reverse_handle, false, chrmap_no_change))
    {
      fatal("More reverse reads than forward reads");
    }

  if (not opt_quiet)
    {
      std::fprintf(stderr,
              "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
              kept,
              truncated,
              discarded);
    }

  if (opt_log != nullptr)
    {
      std::fprintf(fp_log,
              "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
              kept,
              truncated,
              discarded);
    }

  if (reverse_handle != nullptr)
    {
      if (opt_fastaout_rev != nullptr)
        {
          std::fclose(fp_fastaout_rev);
        }

      if (opt_fastqout_rev != nullptr)
        {
          std::fclose(fp_fastqout_rev);
        }

      if (opt_fastaout_discarded_rev != nullptr)
        {
          std::fclose(fp_fastaout_discarded_rev);
        }

      if (opt_fastqout_discarded_rev != nullptr)
        {
          std::fclose(fp_fastqout_discarded_rev);
        }

      fastx_close(reverse_handle);
    }

  if (opt_fastaout != nullptr)
    {
      std::fclose(fp_fastaout);
    }

  if (opt_fastqout != nullptr)
    {
      std::fclose(fp_fastqout);
    }

  if (opt_fastaout_discarded != nullptr)
    {
      std::fclose(fp_fastaout_discarded);
    }

  if (opt_fastqout_discarded != nullptr)
    {
      std::fclose(fp_fastqout_discarded);
    }

  fastx_close(forward_handle);
}


auto fastq_filter() -> void
{
  filter(true, opt_fastq_filter);
}


auto fastx_filter() -> void
{
  filter(false, opt_fastx_filter);
}
