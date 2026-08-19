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

#include "core/fastx.hpp"  // byte_range
#include <array>


namespace vsearch {

/* The FASTA and FASTQ parsers each hold a 256-entry table saying what to do
   with an input byte. The three tables cannot be merged -- 193 of 256 bytes
   deliberately differ, and that difference *is* the documented format
   difference -- but they are all projections of the single partition below.
   Grouping the bytes by the triple (fasta action, fastq sequence action, fastq
   quality action) yields eight distinct groups; 'space' and 'del_or_high' share
   a triple and are split anyway, so that no class name needs a footnote. */
enum struct CharClass : unsigned char {
  iupac,            // ACGTURYSWKMDBHVN, either case: the accepted letters
  printable_other,  // 33-126 minus iupac and '-' '.': legal as quality only
  space,            // 32: outside quality's 33-126, so quality rejects it too
  del_or_high,      // 127-255: as above
  control,          // 0-8, 14-31: fatal in every column
  blank_control,    // tab, VT, FF: fasta strips silently, fastq rejects
  dot_dash,         // '-' '.': fatal in a sequence, legal as a quality symbol
  line_feed,        // 10: line and record structure
  carriage_return,  // 13: stripped silently everywhere
};


namespace ascii {

  /* <cctype> is not constexpr, and utils/ascii_case.hpp wraps std::toupper, so
     the classification below folds case itself. ASCII only, which is all the
     IUPAC letters need. */
  constexpr auto fold_to_upper(unsigned char const byte) noexcept -> unsigned char {
    return (byte >= 'a') and (byte <= 'z')
      ? static_cast<unsigned char>(byte - 'a' + 'A')
      : byte;
  }

  constexpr auto is_iupac_upper(unsigned char const byte) noexcept -> bool {
    return (byte == 'A') or (byte == 'C') or (byte == 'G') or (byte == 'T')
        or (byte == 'U') or (byte == 'R') or (byte == 'Y') or (byte == 'S')
        or (byte == 'W') or (byte == 'K') or (byte == 'M') or (byte == 'D')
        or (byte == 'B') or (byte == 'H') or (byte == 'V') or (byte == 'N');
  }

  constexpr auto is_blank_control(unsigned char const byte) noexcept -> bool {
    return (byte == '\t') or (byte == '\v') or (byte == '\f');
  }

  constexpr auto is_dot_or_dash(unsigned char const byte) noexcept -> bool {
    return (byte == '-') or (byte == '.');
  }

}  // namespace ascii


/* The whole per-byte partition of both parsers, in one screen: a wrong cell has
   nowhere to hide here, whereas the byte-0 bug survived since 2015 as one cell
   of one table disagreeing with its own class, in a sea of 256 identical-looking
   tokens. Order matters: the three named blank controls are peeled off before
   the "below space is a control" test that would otherwise swallow them.

   C++11 constexpr allows a single return statement -- no if, no switch -- so the
   chain below can only be nested conditional operators. */
// NOLINTBEGIN(readability-avoid-nested-conditional-operator)
constexpr auto class_of(unsigned char const byte) noexcept -> CharClass {
  return byte == '\n' ? CharClass::line_feed
    : byte == '\r' ? CharClass::carriage_return
    : ascii::is_blank_control(byte) ? CharClass::blank_control
    : byte < ' ' ? CharClass::control
    : byte == ' ' ? CharClass::space
    : byte > '~' ? CharClass::del_or_high
    : ascii::is_dot_or_dash(byte) ? CharClass::dot_dash
    : ascii::is_iupac_upper(ascii::fold_to_upper(byte)) ? CharClass::iupac
    : CharClass::printable_other;
}
// NOLINTEND(readability-avoid-nested-conditional-operator)


/* A per-class policy is a function rather than a nine-entry array, because
   std::array's const operator[] is only constexpr from C++14 on a conforming
   library: libstdc++ marks it constexpr ungated (so an array would compile on
   GCC 4.9), but libc++ gates it on _LIBCPP_CONSTEXPR_SINCE_CXX14 and the macOS
   builds use libc++. A function of ternaries needs nothing beyond C++11.

   expanded_policy() materialises the flat 256-entry table the parsers index, at
   compile time, so the table lands in .rodata just as a hand-written one would
   and the lookup stays a single load. The hot loops must keep indexing that flat
   table -- never policy(class_of(byte)) per byte, which would put two dependent
   computations in a loop that is 31% of a dereplication run. */
template <unsigned... Bytes> struct byte_indices {};

template <unsigned Count, unsigned... Bytes>
struct make_byte_indices : make_byte_indices<Count - 1, Count - 1, Bytes...> {};

template <unsigned... Bytes>
struct make_byte_indices<0, Bytes...> { using type = byte_indices<Bytes...>; };

template <typename Action, Action (*policy)(CharClass), unsigned... Bytes>
constexpr auto expand(byte_indices<Bytes...> /* deduces the byte pack */) noexcept
  -> std::array<Action, byte_range> {
  return {{ policy(class_of(static_cast<unsigned char>(Bytes)))... }};
}

template <typename Action, Action (*policy)(CharClass)>
constexpr auto expanded_policy() noexcept -> std::array<Action, byte_range> {
  return expand<Action, policy>(typename make_byte_indices<byte_range>::type{});
}

}  // namespace vsearch
