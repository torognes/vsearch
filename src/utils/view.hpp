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


#include "element_order.hpp"  // element_order, element_less
#include <algorithm>  // std::equal, std::lexicographical_compare, std::min
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator> // std::prev, std::next
#include <limits>
#include <type_traits>  // std::is_arithmetic, std::remove_cv


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

  explicit View(Type const * const start, std::size_t const length) noexcept
    : start_ {start},
      length_ {length} {
    assert((start != nullptr) or (length == 0));
    assert(length <= std::numeric_limits<std::ptrdiff_t>::max());
  }

  // Operators
  //
  // The four comparison members below are restricted to arithmetic element
  // types, and are noexcept because of it: for an arithmetic Type the ordering
  // bottoms out in a built-in comparison (see element_order.hpp), which cannot
  // throw, whereas an arbitrary Type's operator< can. The restriction is what
  // makes the promise honest.
  //
  // Deliberately checked per-member rather than at class scope: a member
  // function of a class template is instantiated only when used, so a View
  // over a non-arithmetic Type stays perfectly legal to declare, iterate and
  // index, and only an attempt to *compare* such a view is an error.
  auto operator==(View<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a View requires an arithmetic element type");
    return size() == other.size()
      and std::equal(cbegin(), cend(), other.cbegin());
  }
  auto operator!=(View<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a View requires an arithmetic element type");
    return not (*this == other);
  }
  // Ordering goes through element_order (see element_order.hpp), so that a
  // View<char> orders its bytes as unsigned char, like std::strcmp and
  // std::string, rather than as a possibly-signed char.
  auto operator<(View<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a View requires an arithmetic element type");
    return std::lexicographical_compare(cbegin(), cend(),
                                        other.cbegin(), other.cend(),
                                        element_less<Type>{});
  }

  // Three-way form, with the sign convention of std::strcmp and
  // std::string::compare, for callers that must tell "equal" from
  // "greater" in a single pass. C++11 offers no three-way algorithm;
  // std::lexicographical_compare_three_way arrives in C++20.
  auto compare(View<Type> const & other) const noexcept -> int {
    static_assert(comparable, "comparing a View requires an arithmetic element type");
    auto const common_length = std::min(size(), other.size());
    for (std::size_t index = 0; index < common_length; ++index) {
      auto const order = element_order<Type>::compare((*this)[index], other[index]);
      if (order != 0) { return order; }
    }
    if (size() == other.size()) { return 0; }
    return (size() < other.size()) ? -1 : +1;
  }

  // Iterators
  constexpr auto begin() const noexcept -> Type const * { return data(); }
  constexpr auto cbegin() const noexcept -> Type const * { return data(); }
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
  constexpr auto data() const noexcept -> Type const * { return start_; }
  auto operator[](std::size_t const index) const noexcept -> Type const & {
    assert(index < size());
    auto const distance = static_cast<std::ptrdiff_t>(index);
    return *std::next(data(), distance);
  }

  // Observers
  constexpr auto size() const noexcept -> std::size_t { return length_; }
  auto size_bytes() const noexcept -> std::size_t {
    assert(size() <= (std::numeric_limits<std::size_t>::max() / sizeof(Type)));
    return size() * sizeof(Type);
  }
  constexpr auto empty() const noexcept -> bool { return size() == 0; }

  // Subviews
  auto subspan(std::size_t const offset, std::size_t const count) const noexcept -> View {
    assert(offset <= size());
    assert(count <= size() - offset);
    auto const distance = static_cast<std::ptrdiff_t>(offset);
    auto const * const new_start = std::next(data(), distance);
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
  // Predicate behind the comparison members' static_assert above. remove_cv is
  // needed because std::is_arithmetic<char const> is false, and View<Type const>
  // is an ordinary instantiation that must stay comparable.
  static constexpr bool comparable =
    std::is_arithmetic<typename std::remove_cv<Type>::type>::value;

  Type const * start_ {};
  std::size_t  length_ {};
};


// A View over a whole container, with the element type deduced so that
// call sites neither spell it out nor reach for data():
//
//   make_view(cigar_string_)
//
// rather than
//
//   View<char>{cigar_string_.data(), cigar_string_.size()}
//
// A prefix or an interior slice comes from composing with the members
// that already exist -- make_view(vec).first(count), or
// make_view(vec).subspan(offset, count), which also carries the bounds
// assertions that an open-coded &vec[offset] cannot. That composition is
// why there is no offset/count overload here.
//
// noexcept: the containers this is called with (std::vector, std::array,
// std::string) all have noexcept data() and size(), so the only other
// operation left is View's noexcept constructor. A container whose
// accessors can throw is not a supported argument.
//
// constexpr: what C++11 requires of a constexpr function is the shape of
// the declaration, not that every call can be folded -- one return
// statement, a literal return type (View has an implicitly constexpr
// defaulted default constructor and two literal members), and literal
// parameter types (a reference type is one). No standard container has a
// constexpr data() before C++17, so under C++11 the specializations used
// here are simply not constant expressions; the keyword costs nothing and
// starts working the day the standard level moves.
template <typename Container>
constexpr auto make_view(Container const & container) noexcept
  -> View<typename Container::value_type> {
  return View<typename Container::value_type>{container.data(), container.size()};
}
