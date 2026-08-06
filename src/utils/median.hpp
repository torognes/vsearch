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

#include "view.hpp"
#include <cassert>
#include <cstddef>


// Median of a range already sorted by *decreasing* value, reading the value
// out of each element through value_of.
//
// Shared by --sortbysize (median amplicon abundance, over its deck) and the
// --derep_* family (median cluster size, over the live prefix of its hash
// table), which computed the same median separately.
//
// The range is a View rather than an iterator pair because the two callers
// disagree about where it ends: sortbysize means its whole deck, while derep
// means only the buckets that ended up used, the empty ones having sorted to
// the tail. .first(clusters) says that, and asserts it, where a hand-built
// end iterator would only have implied it.
//
// The order must be decreasing, not merely sorted: the even-count branch
// subtracts in the projected type, which is unsigned in both callers, and
// would wrap on an increasing range. The assert below states that.
//
// refactoring C++17 [[nodiscard]]
template <typename Type, typename Projection>
auto median_of_descending(View<Type> const values, Projection value_of) noexcept -> double
{
  // the result is a round value, or a value with a remainder of 0.5
  static constexpr auto half = 0.5;

  if (values.empty()) {
    return 0.0;
  }

  // plain division on the view's std::size_t length. ldiv would have needed a
  // narrowing cast to long (32-bit on the Windows target), its remainder is
  // recomputed with % just below anyway, and its quotient came back signed and
  // so had to be converted again at each subscript -- size() / 2 is already
  // unsigned.
  auto const midpoint = values.size() / 2;

  // odd number of elements: index is zero-based, so if size == 3, midpoint == 1
  if ((values.size() % 2) != 0) {
    return static_cast<double>(value_of(values[midpoint]));
  }

  // even number of elements
  // index is zero-based, so if size == 4, midpoint == 2, lhs index == 1
  auto const rhs_size = value_of(values[midpoint]);
  auto const lhs_size = value_of(values[midpoint - 1]);
  // sorted by decreasing abundance: lhs size >= rhs size
  assert(lhs_size >= rhs_size);
  // limit risk of integer additon overflow:
  // a >= b ; (a + b) / 2 == b + (a - b) / 2
  return static_cast<double>(rhs_size)
    + (static_cast<double>(lhs_size - rhs_size) * half);
}

