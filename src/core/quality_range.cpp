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
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
  OF THE POSSIBILITY OF SUCH DAMAGE.

*/



#include "core/quality_range.hpp"
#include "vsearch.hpp"  // struct Parameters
#include "utils/decimal_digits.hpp"  // decimal::to_text
#include "utils/fatal.hpp"  // fatal
#include <cstdint>  // int64_t
#include <string>  // std::string


namespace vsearch {

auto classify_quality(int64_t const quality_score,
                      struct Parameters const & parameters) noexcept -> QualityBound
{
  if (quality_score < parameters.opt_fastq_qmin) { return QualityBound::below_qmin; }
  if (quality_score > parameters.opt_fastq_qmax) { return QualityBound::above_qmax; }
  return QualityBound::in_range;
}


auto quality_out_of_range_message(QualityBound const bound,
                                  int64_t const quality_score,
                                  struct Parameters const & parameters,
                                  QualityLocation const location) -> std::string
{
  /* decimal::to_text, not std::to_string: on libstdc++ <= 10 the latter is a
     std::vsnprintf call with a format string (see decimal_digits.hpp). With
     no location the text is byte-for-byte what the merged copies produced. */

  /* Appended after the bound, so a test grepping for the bound alone still
     matches; only a whole-line match notices. */
  auto const where = location.known()
    ? " in entry no " + decimal::to_text(location.record)
      + " starting on line " + decimal::to_text(location.line)
    : std::string{};

  if (bound == QualityBound::below_qmin)
    {
      return "FASTQ quality value (" + decimal::to_text(quality_score)
        + ") below qmin (" + decimal::to_text(parameters.opt_fastq_qmin) + ")" + where;
    }
  return "FASTQ quality value (" + decimal::to_text(quality_score)
    + ") above qmax (" + decimal::to_text(parameters.opt_fastq_qmax) + ")" + where + "\n"
      "By default, quality values range from 0 to 93\n"
      "(0 to 62 with --fastq_ascii 64).\n"
      "To allow higher quality values, "
      "please use the option --fastq_qmax " + decimal::to_text(quality_score);
}


auto check_quality_score(int64_t const quality_score,
                         struct Parameters const & parameters,
                         QualityLocation const location) -> void
{
  auto const bound = classify_quality(quality_score, parameters);
  if (bound == QualityBound::in_range) { return; }
  // route through fatal(), which writes stderr and the --log file, and in a
  // library session throws instead of std::exit()ing
  fatal(quality_out_of_range_message(bound, quality_score, parameters, location));
}


auto check_quality_range(QualitySymbolRange const range,
                         struct Parameters const & parameters) -> void
{
  if (not range.seen()) { return; }
  auto const ascii_offset = parameters.opt_fastq_ascii;
  check_quality_score(range.lowest - ascii_offset, parameters);
  check_quality_score(range.highest - ascii_offset, parameters);
}

}  // namespace vsearch
