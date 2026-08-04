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


#include "utils/view.hpp"  // View<char>
#include <array>
#include <cassert>
#include <cstddef>  // std::size_t, std::ptrdiff_t
#include <cstdint>  // uint64_t
#include <iterator>  // std::next, std::prev, std::distance
#include <limits>  // std::numeric_limits
#include <type_traits>  // std::integral_constant, std::is_integral, std::is_same,
                        // std::is_signed, std::remove_cv


/* The decimal digits of an integer, written into a caller-supplied stack
   buffer: no format string parsed at run time, no allocation, no varargs.

   This is the shared half of what replaces vsearch's integer std::fprintf
   calls. It is deliberately sink-agnostic -- it hands back a View<char> over
   the digits and knows nothing about where they go. print_view.hpp adds the
   std::FILE * sink on top; a caller appending to a std::string or to a
   character buffer (the CIGAR builders in utils/cigar.cpp and
   core/linmemalign.cpp) uses this header directly and keeps <cstdio> out of
   its include list, the same split print_view.hpp already makes for
   View<char>.

   Why not std::to_string: before libstdc++ 11 its integral overloads are
   __gnu_cxx::__to_xstring(&std::vsnprintf, "%ld", ...) -- the run-time format
   string this exists to remove -- and vsearch still supports GCC 4.9. The
   digits are also produced without the temporary std::string.

   Why not the PRI macros and a format string: those macros exist only because
   a format string has to name the width of its argument. With the digits
   produced here there is nothing to name, which is what lets <cinttypes>
   leave the tree.

   C++17 refactoring: std::to_chars */

namespace decimal {

  /* Widest decimal form of any type accepted below, in characters:
       uint64_t   18446744073709551615  -> 20 digits
       int64_t   -9223372036854775808   -> 19 digits plus a sign, also 20
     digits10 is the number of decimal digits that survive a round trip, one
     fewer than the widest representation, hence the + 1. */
  constexpr std::size_t max_width = std::numeric_limits<uint64_t>::digits10 + 1;

  /* Storage for one number. It lives in the caller's frame, so the View
     returned by to_decimal() borrows it: the view must not outlive the buffer
     it was filled from. Every caller consumes the view in the same statement
     or the next one. */
  using Buffer = std::array<char, max_width>;


  namespace detail {

    /* A value split into the parts the digit loop needs. Kept together so that
       one dispatch answers both questions. */
    struct Signed_magnitude {
      uint64_t magnitude;
      bool     negative;
    };

    /* Tag dispatch on signedness rather than a plain 'if (value < 0)': for an
       unsigned Integer that comparison is always false, which GCC reports
       under -Wtype-limits, and vsearch's debug build enables it. */
    template <typename Integer>
    constexpr auto split_sign(Integer const value, std::true_type /* is_signed */) -> Signed_magnitude {
      /* The conversion to uint64_t comes before the negation on purpose:
         negating the most negative value of a signed type overflows in that
         type. Unsigned arithmetic is modulo 2^64, so
         0 - static_cast<uint64_t>(value) is the magnitude for every negative
         input, the minimum included -- that is the case a hand-written -value
         misses. */
      return value < 0
        ? Signed_magnitude{uint64_t{0} - static_cast<uint64_t>(value), true}
        : Signed_magnitude{static_cast<uint64_t>(value), false};
    }

    template <typename Integer>
    constexpr auto split_sign(Integer const value, std::false_type /* is_signed */) -> Signed_magnitude {
      return Signed_magnitude{static_cast<uint64_t>(value), false};
    }

  }  // namespace detail


  /* The decimal form of 'value', as a view into 'buffer'.

     Written back to front, because that is the order repeated division yields
     the digits in; the view therefore starts somewhere inside the buffer
     rather than at its beginning, which is why the width is returned with the
     pointer instead of being re-derived by the caller.

     noexcept: it writes only into the caller's array and allocates nothing. */
  template <typename Integer>
  auto to_decimal(Buffer & buffer, Integer const value) noexcept -> View<char> {
    static_assert(std::is_integral<Integer>::value,
                  "to_decimal() formats integers; pass an integral type");
    /* A char argument would print its numeric code -- 'A' as 65 -- which is
       never what a caller means, and the mistake compiles silently because
       char is an integral type. Rejecting it costs one line. */
    static_assert(not std::is_same<typename std::remove_cv<Integer>::type, char>::value,
                  "to_decimal() of a char would print its code, not the character");
    static_assert(sizeof(Integer) <= sizeof(uint64_t),
                  "max_width is derived from uint64_t; a wider type would overflow the buffer");

    auto const split =
      detail::split_sign(value,
                         std::integral_constant<bool, std::is_signed<Integer>::value>{});

    auto * const buffer_end = std::next(buffer.data(), static_cast<std::ptrdiff_t>(buffer.size()));
    auto * cursor = buffer_end;
    auto rest = split.magnitude;
    /* do-while, not while: zero has one digit and must be written. */
    do {
      cursor = std::prev(cursor);
      *cursor = static_cast<char>('0' + (rest % 10));
      rest /= 10;
    } while (rest != 0);

    if (split.negative) {
      cursor = std::prev(cursor);
      *cursor = '-';
    }

    /* The buffer cannot be overrun: max_width above is the widest decimal form
       of any type the static_asserts admit, sign included, and the loop writes
       exactly one character per decimal digit. This states that invariant
       rather than enforcing it -- by the time it could fail the writes have
       already happened -- and is here as a tripwire for a future change to
       max_width or to the accepted types. */
    assert(cursor >= buffer.data());

    auto const width = static_cast<std::size_t>(std::distance(cursor, buffer_end));
    return View<char>{cursor, width};
  }

}  // namespace decimal


/* tests:

   Checked against std::snprintf as the reference, under -D_GLIBCXX_DEBUG and
   the debug warning set, and instantiated for the four cross-compilation
   targets. The cases below are the ones a hand-rolled digit loop gets wrong.

   #include <cassert>
   #include <cstdint>
   #include <string>

   template <typename Integer>
   auto as_string(Integer const value) -> std::string {
     decimal::Buffer buffer {};
     auto const digits = decimal::to_decimal(buffer, value);
     return std::string(digits.data(), digits.size());
   }

   auto main() -> int {
     assert(decimal::max_width == 20);
     assert(as_string(uint64_t{0}) == "0");                    // do-while, not while
     assert(as_string(UINT64_MAX) == "18446744073709551615");  // 20 digits, fills the buffer
     assert(as_string(INT64_MAX) == "9223372036854775807");
     assert(as_string(INT64_MIN) == "-9223372036854775808");   // -value would overflow
     assert(as_string(-1) == "-1");
     assert(as_string(4294967295U) == "4294967295");           // the "%u" case
     assert(as_string(static_cast<short>(-7)) == "-7");
     for (uint64_t i = 0; i <= 100000; ++i) { assert(not as_string(i).empty()); }

     // rejected at compile time, deliberately:
     // as_string('A');   // char is integral, so this would print 65
     // as_string(1.5);   // not an integral type
   }
*/
