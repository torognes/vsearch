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


#pragma once

#include "core/fastx_char_class.hpp"  // vsearch::ascii::is_printable
#include "utils/decimal_digits.hpp"  // decimal::to_text
#include <cstdint>  // uint64_t
#include <string>  // std::string


namespace vsearch {

/* One wording for every illegal input byte, in one place.

   Four sites reject a byte while reading FASTA/FASTQ -- the FASTA sequence,
   the FASTQ sequence, the FASTQ quality line, and the shared header scan --
   and each used to phrase its own message. They disagreed on three
   independent things: whether the location came first ("Invalid line 2 in
   FASTQ file: ...") or last ("... in sequence on line 2 of FASTA file"),
   whether an unprintable byte restructured the sentence ("Illegal
   unprintable ASCII character no 1") or was parenthesised ("Illegal sequence
   character (unprintable, no 1)"), and whether the offending field was named
   at all. Nothing in the manual pins any of the three, so the disagreement
   was accident, not contract.

   The template below is the one shape all four now produce:

     Illegal <field> character 'X' on line N of <FORMAT> file
     Illegal <field> character (unprintable, no D) on line N of <FORMAT> file

   Callers keep their own deferred-error handling: on a worker thread the
   message is recorded and reported from the main thread rather than
   fatal()ed in place (see the deferred-error note in fastx.hpp). */

// Which field of the record the byte was found in; supplies the noun.
enum struct IllegalField { sequence, quality, header };

// Which format the message names; the header scan serves both.
enum struct IllegalFormat { fasta, fastq, fasta_or_fastq };

/* Both spelled as a single return so they can be constexpr under C++11 (no
   if, no switch), and named for the message rather than for their return type:
   plain `format_name` at vsearch scope would be an invitation to collide. */
constexpr auto illegal_field_noun(IllegalField const field) noexcept -> char const * {
  return field == IllegalField::sequence ? "sequence"
    : field == IllegalField::quality ? "quality"
    : "header";
}

constexpr auto illegal_format_name(IllegalFormat const format) noexcept -> char const * {
  return format == IllegalFormat::fasta ? "FASTA"
    : format == IllegalFormat::fastq ? "FASTQ"
    : "FASTA/FASTQ";
}

/* Distinct parameter types throughout, so no two arguments can be swapped
   silently: the two enums name the field and the format, the byte is an
   unsigned char, and the line number is a uint64_t. */
inline auto illegal_character_message(IllegalField const field,
                                      IllegalFormat const format,
                                      unsigned char const symbol,
                                      std::uint64_t const line_number) -> std::string {
  /* decimal::to_text, not std::to_string: on libstdc++ <= 10 the latter is a
     std::vsnprintf call with a format string (see decimal_digits.hpp). */
  std::string const rendered_byte =
    ascii::is_printable(symbol)
    ? "'" + std::string(1, static_cast<char>(symbol)) + "'"
    : "(unprintable, no " + decimal::to_text(symbol) + ")";
  return std::string("Illegal ") + illegal_field_noun(field) + " character "
    + rendered_byte
    + " on line " + decimal::to_text(line_number)
    + " of " + illegal_format_name(format) + " file";
}

}  // namespace vsearch
