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

#include <cassert>
#include <cstdint>  // uint16_t, uint32_t
#include <type_traits>

#ifndef NDEBUG
#include <limits>
#endif

// round up to the next multiple of x
//
// common idiom:
//   ((input + ( x - 1 )) / x ) * x
//
// when x is a power of two (1, 2, 4, 8, 16, 32, 64, ...):
//   (input + ( x - 1 )) & ~( x - 1 )
//
// for instance: let's round up 11 to the next multiple of 8, i.e. 16
//
//  - multiple (8):    b0000'1000
//  -  1:              b0000'0001 
//                     ----------  (minus)
//  - stub:            b0000'0111
//
//
//  - stub:            b0000'0111
//                     ----------  (bitwise not (~))
//  - bitmask          b1111'1000
//
//
//  - 11:              b0000'1011
//  - stub:            b0000'0111
//                     ----------  (add)
//                     b0001'0010  => temporary result is 18
//
// 
//  - 18:              b0001'0010
//  - bitmask:         b1111'1000
//                     ----------  (bitwise and (&))
//                     b0001'0000  => result is 16
//
// note: if input > max + x, then round up is > max => Error!


// C++14 refactoring: constexpr
// C++17 refactoring: [[nodiscard]]
template <typename Unsigned = std::uint16_t>
auto round_up_to_8(Unsigned const input) -> Unsigned {
  static_assert(std::is_same<Unsigned, std::uint16_t>::value            \
                or std::is_same<Unsigned, std::uint32_t>::value,
                "Invalid type! Only uint16_t, or uint32_t can be used.");
  // add stub to guarantee overflow into the next bucket, then zero
  // out the remainder bits
  static constexpr Unsigned step = 8;
  assert(step >= 1);
  static constexpr Unsigned stub = step - 1;
  static constexpr auto bitmask = static_cast<Unsigned>(~stub);
  assert(input <= std::numeric_limits<Unsigned>::max() - stub);
  return (input + stub) & bitmask;
}


// refactoring: C++14 tests

/*

// test return types:
static_assert(std::is_same<decltype(round_up_to_8<std::uint16_t>(0)),
std::uint16_t>::value, "");
static_assert(std::is_same<decltype(round_up_to_8<std::uint32_t>(0)),
std::uint32_t>::value, "");

// test default type parameter:
static_assert(round_up_to_8(0U) == round_up_to_8<std::uint16_t>(0), "");
static_assert(round_up_to_8<>(0U) == round_up_to_8<std::uint16_t>(0), "");

// test return values (std::uint16_t)
static_assert(round_up_to_8<std::uint16_t>(0) == 0, "");
static_assert(round_up_to_8<std::uint16_t>(1) == 8, "");
static_assert(round_up_to_8<std::uint16_t>(7) == 8, "");
static_assert(round_up_to_8<std::uint16_t>(8) == 8, "");
static_assert(round_up_to_8<std::uint16_t>(9) == 16, "");
static_assert(round_up_to_8<std::uint16_t>(15) == 16, "");
static_assert(round_up_to_8<std::uint16_t>(65'528) == 65'528, "");  // max uint16_t - 7

// test return values (std::uint32_t)
static_assert(round_up_to_8<std::uint32_t>(0) == 0, "");
static_assert(round_up_to_8<std::uint32_t>(1) == 8, "");
static_assert(round_up_to_8<std::uint32_t>(7) == 8, "");
static_assert(round_up_to_8<std::uint32_t>(8) == 8, "");
static_assert(round_up_to_8<std::uint32_t>(9) == 16, "");
static_assert(round_up_to_8<std::uint32_t>(15) == 16, "");
static_assert(round_up_to_8<std::uint32_t>(4'294'967'288) == 4'294'967'288, ""); // max uint32_t - 7

*/


// refactoring:

// class with template member functions?

// usage:
//
// auto const y = Round_up(x).to_next_multiple<8>()
// auto const y = Round_up<uint32_t>(x).to_next_multiple<8>() 
// auto const y = Round_up(x).to_next_power<2>()
// ...
