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


// How View and Span order their elements.
//
// The primary template defers to the element type's own operator<, which is
// what a bare std::lexicographical_compare would have done anyway.
//
// The specialization for char is the reason this trait exists: it orders bytes
// as unsigned char. That is what std::strcmp does, and what std::string does
// through std::char_traits<char>::lt -- but it is *not* what comparing char
// with operator< does, because char is signed on x86-64 and on the Windows
// target while ARM and PowerPC Linux default it to unsigned. Without the
// specialization, View<char>{} < View<char>{} would order any byte with its
// high bit set differently from std::string, differently from std::strcmp, and
// differently from one architecture to the next.
//
// That divergence is reachable and user-visible in vsearch: sequence headers
// are arbitrary bytes and routinely carry UTF-8 or Latin-1 (an accented author
// name, a non-ASCII locality), and the header is the tie-break of every
// --sortby* and --derep* ordering. Measured over 2692881 pairs from a byte
// alphabet straddling 0x80, the unspecialized comparison disagrees with
// std::strcmp on 1378400 of them where char is signed, and on none where it is
// unsigned.
//
// Only char is specialized. signed char and unsigned char say what they mean
// and are left to the primary template, the same line std::char_traits draws.

template <typename Type>
struct element_order {
  static auto less(Type const & lhs, Type const & rhs) -> bool {
    return lhs < rhs;
  }

  // Three-way, with the sign convention of std::strcmp: negative if lhs sorts
  // first, positive if rhs does, zero if the two are equivalent. Expressed with
  // two less() calls because the primary template can assume nothing beyond a
  // strict weak ordering.
  static auto compare(Type const & lhs, Type const & rhs) -> int {
    if (less(lhs, rhs)) { return -1; }
    if (less(rhs, lhs)) { return +1; }
    return 0;
  }
};

template <>
struct element_order<char> {
  static auto less(char const lhs, char const rhs) -> bool {
    return static_cast<unsigned char>(lhs) < static_cast<unsigned char>(rhs);
  }

  // One byte comparison rather than the primary template's two: char is
  // totally ordered once it is read as unsigned char.
  static auto compare(char const lhs, char const rhs) -> int {
    auto const lhs_byte = static_cast<unsigned char>(lhs);
    auto const rhs_byte = static_cast<unsigned char>(rhs);
    if (lhs_byte == rhs_byte) { return 0; }
    return (lhs_byte < rhs_byte) ? -1 : +1;
  }
};


// Functor wrapper, for the std algorithms that take a comparison object.
template <typename Type>
struct element_less {
  auto operator()(Type const & lhs, Type const & rhs) const -> bool {
    return element_order<Type>::less(lhs, rhs);
  }
};
