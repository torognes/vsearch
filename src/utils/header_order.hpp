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

#pragma once


#include "view.hpp"  // View<char>
#include <algorithm>  // std::lexicographical_compare, std::min
#include <cstddef>  // std::size_t


// Ordering for sequence headers, matching std::strcmp byte for byte.
//
// std::strcmp compares bytes as unsigned char. Plain char is signed on x86-64
// and on the Windows target, but unsigned on ARM and PowerPC Linux, so
// comparing two View<char> with View::operator< (a std::lexicographical_compare
// over char) would order any header carrying a high-bit byte differently from
// the historical output, *and* differently from one architecture to the next.
// Headers are arbitrary bytes and may carry UTF-8 or Latin-1, so that case is
// reachable: in the "C" locale a UTF-8 'e' begins with 0xC3, which strcmp sorts
// after 'z' and signed char sorts before 'A'.
//
// Comparing through unsigned char keeps one order everywhere, for the same
// reason utils/ascii_case.hpp casts before calling <cctype>. Use these rather
// than View::operator< whenever the result is user-visible, which for headers
// it always is: they are the tie-break of every --sortby* and --derep* output.
//
// Equality needs no helper: View::operator== is signedness-independent, and it
// compares sizes first, so it also short-circuits where strcmp would walk to
// the terminating null.

inline auto header_less(View<char> const lhs, View<char> const rhs) -> bool {
  return std::lexicographical_compare(
      lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(),
      [](char const left, char const right) -> bool {
        return static_cast<unsigned char>(left) < static_cast<unsigned char>(right);
      });
}

// Three-way form, for the comparators that tie-break on "equal headers" and so
// need to tell "equal" from "greater" in one pass. Returns 0 if the headers are
// identical, -1 if lhs sorts first, +1 if rhs sorts first -- the sign
// convention of std::strcmp. A plain loop rather than std::mismatch: the
// four-iterator overload that would bound both ranges is C++14.
inline auto header_compare(View<char> const lhs, View<char> const rhs) -> int {
  auto const common_length = std::min(lhs.size(), rhs.size());
  for (std::size_t index = 0; index < common_length; ++index) {
    auto const lhs_byte = static_cast<unsigned char>(lhs[index]);
    auto const rhs_byte = static_cast<unsigned char>(rhs[index]);
    if (lhs_byte != rhs_byte) {
      return (lhs_byte < rhs_byte) ? -1 : +1;
    }
  }
  if (lhs.size() == rhs.size()) {
    return 0;
  }
  return (lhs.size() < rhs.size()) ? -1 : +1;
}
