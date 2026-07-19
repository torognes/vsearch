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


#include <algorithm>  // std::equal, std::lexicographical_compare, std::min
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator> // std::prev, std::next
#include <limits>


// A read-only counterpart of Span (see span.hpp), inspired by std::span
// (C++20): a non-owning view over a contiguous sequence of elements of any type
// Type. Unlike Span, its single constructor stores the pointer as
// `Type const *`, so a pointer to const data is carried as const rather than
// laundered to a mutable pointer through a const_cast. The interface mirrors
// the read-only half of Span (same member names, e.g. subspan) so that a
// read-only Span<T> parameter can become a View<T> without editing the body.

template <typename Type = char>
class View {
public:
  // Empty view (null pointer, zero length) via the in-class member initializers;
  // the explicit constructor below otherwise suppresses the implicit default.
  View() noexcept = default;

  explicit View(Type const * start, std::size_t const length) noexcept
    : start_ {start},
      length_ {length} {
    assert((start != nullptr) or (length == 0));
    assert(length <= std::numeric_limits<std::ptrdiff_t>::max());
  }

  // Operators
  auto operator==(View<Type> const & other) const -> bool {
    return size() == other.size()
      and std::equal(cbegin(), cend(), other.cbegin());
  }
  auto operator<(View<Type> const & other) const -> bool {
    return std::lexicographical_compare(cbegin(), cend(),
                                        other.cbegin(), other.cend());
  }

  // Iterators
  auto begin() const noexcept -> Type const * { return data(); }
  auto cbegin() const noexcept -> Type const * { return data(); }
  auto end() const noexcept -> Type const * {
    auto const distance = static_cast<std::ptrdiff_t>(size());
    return std::next(data(), distance);
  }
  auto cend() const noexcept -> Type const * { return end(); }
  auto rbegin() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(end());
  }
  auto crbegin() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(cend());
  }
  auto rend() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(begin());
  }
  auto crend() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(cbegin());
  }

  // Element access
  auto front() const noexcept -> Type const & {
    assert(not empty());
    return *data();
  }
  auto back() const noexcept -> Type const & {
    assert(not empty());
    return *std::prev(end());
  }
  auto data() const noexcept -> Type const * { return start_; }
  auto operator[](std::size_t const index) const noexcept -> Type const & {
    assert(index < size());
    auto const distance = static_cast<std::ptrdiff_t>(index);
    return *std::next(data(), distance);
  }

  // Observers
  auto size() const noexcept -> std::size_t { return length_; }
  auto size_bytes() const noexcept -> std::size_t {
    assert(size() <= (std::numeric_limits<std::size_t>::max() / sizeof(Type)));
    return size() * sizeof(Type);
  }
  auto empty() const noexcept -> bool { return size() == 0; }

  // Subviews
  auto subspan(std::size_t const offset, std::size_t const count) const noexcept -> View {
    assert(offset <= size());
    assert(count <= size() - offset);
    auto const distance = static_cast<std::ptrdiff_t>(offset);
    auto const * new_start = std::next(data(), distance);
    return View{new_start, count};
  }
  auto first(std::size_t const count) const noexcept -> View {
    return subspan(0, count);
  }
  auto last(std::size_t const count) const noexcept -> View {
    assert(count <= size());
    return subspan(size() - count, count);
  }
  auto drop(std::size_t const count) const noexcept -> View {
    // drop n first items, return empty if n is >= size()
    auto const offset = std::min(count, size());
    assert(offset <= size());
    return subspan(offset, size() - offset);
  }

private:
  Type const * start_ {};
  std::size_t  length_ {};
};
