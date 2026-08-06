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
#include "view.hpp"  // View, the read-only counterpart returned by operator View
#include <algorithm>  // std::equal, std::lexicographical_compare, std::min
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdlib>  // std::size_t
#include <iterator> // std::prev, std::next
#include <type_traits>  // std::is_arithmetic, std::remove_cv

#ifndef NDEBUG
#include <limits>
#endif


// TODO:
//  - add friend function for custom hashing,
//  - goal is to be able to build std::map and std::set of Span<char>


// simple version of std::span (C++20)
//
// only valid for contiguous sequences of elements (vectors or arrays)
// of any type Type (except std::vector<bool>?)

template <typename Type = char>
class Span {
public:
  // Empty span (null pointer, zero length) via the in-class member initializers;
  // the explicit constructors below otherwise suppress the implicit default.
  Span() noexcept = default;

  explicit Span(Type * const start, std::size_t const length) noexcept
    : start_ {start},
      length_ {length} {
    assert((start != nullptr) or (length == 0));
    assert(length <= max_ptrdiff);
  }

  // Conversion to the read-only View<Type> over the same extent. This is the
  // only direction that is sound: it adds a restriction (drops the permission
  // to write) and can never invent access that the Span did not already hold.
  // There is deliberately no View-to-Span conversion, because that direction
  // would have to launder a `Type const *` into a `Type *`; see view.hpp.
  //
  // Kept explicit, as in swarm: handing a mutable Span to a View-consuming API
  // is a narrowing of intent that the call site should spell out, via
  // static_cast<View<Type>>(span) or View<Type>{span}. An implicit conversion
  // would also make an overload pair on Span<T>/View<T> silently ambiguous.
  explicit operator View<Type>() const noexcept {
    return View<Type>{start_, length_};
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
  // function of a class template is instantiated only when used, so
  // Span<Type> over a non-arithmetic Type stays perfectly legal to declare,
  // iterate and index -- Span<struct hit> in core/searchcore.cpp does exactly
  // that -- and only an attempt to *compare* such a span is an error.
  auto operator==(Span<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a Span requires an arithmetic element type");
    return size() == other.size()
      and std::equal(cbegin(), cend(), other.cbegin());
  }
  auto operator!=(Span<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a Span requires an arithmetic element type");
    return not (*this == other);
  }
  // Ordering goes through element_order (see element_order.hpp), so that a
  // Span<char> orders its bytes as unsigned char, like std::strcmp and
  // std::string, rather than as a possibly-signed char.
  auto operator<(Span<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a Span requires an arithmetic element type");
    return std::lexicographical_compare(cbegin(), cend(),
                                        other.cbegin(), other.cend(),
                                        element_less<Type>{});
  }

  // Three-way form, with the sign convention of std::strcmp and
  // std::string::compare, for callers that must tell "equal" from
  // "greater" in a single pass. C++11 offers no three-way algorithm;
  // std::lexicographical_compare_three_way arrives in C++20.
  auto compare(Span<Type> const & other) const noexcept -> int {
    static_assert(comparable, "comparing a Span requires an arithmetic element type");
    auto const common_length = std::min(size(), other.size());
    for (std::size_t index = 0; index < common_length; ++index) {
      auto const order = element_order<Type>::compare((*this)[index], other[index]);
      if (order != 0) { return order; }
    }
    if (size() == other.size()) { return 0; }
    return (size() < other.size()) ? -1 : +1;
  }

  // Iterators
  auto begin() const noexcept -> Type * { return data(); }
  auto cbegin() const noexcept -> Type const * { return data(); }
  auto end() const noexcept -> Type * {
    auto const distance = static_cast<std::ptrdiff_t>(size());
    return std::next(data(), distance);
  }
  auto cend() const noexcept -> Type const * {
    return end();
  }
  auto rbegin() const noexcept -> std::reverse_iterator<Type *> {
    return std::reverse_iterator<Type *>(end());
  }
  auto crbegin() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(cend());
  }
  auto rend() const noexcept -> std::reverse_iterator<Type *> {
    return std::reverse_iterator<Type *>(begin());
  }
  auto crend() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(cbegin());
  }

  // Element access
  // C++17 refactoring: [[nodiscard]]
  auto front() const noexcept -> Type & {
    assert(not empty());
    return *data();
  }
  auto back() const noexcept -> Type & {
    assert(not empty());
    return *std::prev(end());
  }
  auto data() const noexcept -> Type * { return start_; }
  auto operator[](std::size_t const index) const noexcept -> Type & {
    assert(index < size());
    auto const distance = static_cast<std::ptrdiff_t>(index);
    return *std::next(data(), distance);
  }

  // Observers
  auto size() const noexcept -> std::size_t { return length_; }
  auto size_bytes() const noexcept -> std::size_t {
    assert(size() <= (max_size / sizeof(Type)));
    return size() * sizeof(Type);
  }
  auto empty() const noexcept -> bool { return size() == 0; }

  // Subviews
  auto subspan(std::size_t const offset, std::size_t const count) const noexcept -> Span {
    assert(offset <= size());
    assert(count <= size() - offset);
    auto const distance = static_cast<std::ptrdiff_t>(offset);
    auto * const new_start = std::next(data(), distance);
    return Span{new_start, count};
  }
  auto first(std::size_t const count) const noexcept -> Span {
    return subspan(0, count);
  }
  auto last(std::size_t const count) const noexcept -> Span {
    assert(count <= size());
    return subspan(size() - count, count);
  }
  auto drop(std::size_t const count) const noexcept -> Span {
    // drop n first items, return empty if n is >= size()
    auto const offset = std::min(count, size());
    assert(offset <= size());
    return subspan(offset, size() - offset);
  }

private:
  // Predicate behind the comparison members' static_assert above. remove_cv is
  // needed because std::is_arithmetic<char const> is false, and Span<Type const>
  // is an ordinary read-only instantiation that must stay comparable.
  static constexpr bool comparable =
    std::is_arithmetic<typename std::remove_cv<Type>::type>::value;

  Type * start_ {};
  std::size_t length_ {};

#ifndef NDEBUG
  // Upper bounds for the debug-build assertions above. Kept private so that
  // including this header does not export these names into the global
  // namespace (where they could shadow, or be shadowed by, an unrelated
  // max_size / max_ptrdiff elsewhere).
  // C++17 refactoring: [[maybe_unused]]
  static constexpr std::ptrdiff_t max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
  static constexpr std::size_t max_size = std::numeric_limits<std::size_t>::max();
#endif
};


// A Span over a whole container; see make_view() in view.hpp for the
// rationale, and for why a prefix or a slice composes from the existing
// members (first(), subspan()) instead of an overload here.
//
// The container is taken by non-const reference on purpose: a const
// container then fails to compile rather than quietly yielding a mutable
// span over data it does not own the right to modify.
// constexpr for the same reason as make_view(): see the note there.
template <typename Container>
constexpr auto make_span(Container & container) noexcept
  -> Span<typename Container::value_type> {
  return Span<typename Container::value_type>{container.data(), container.size()};
}


// tests:

// #include <algorithm>
// #include <vector>

// auto main() -> int {
//   std::vector<char> v = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
//   auto s = Span<char>{v.data(), 5};
//   assert(s.size() == 5);
//   assert(s.size_bytes() == 5);
//   assert(! s.empty());
//   assert(s.front() == 'a');
//   assert(s.back() == 'e');
//   assert(s[1] == 'b');
//   for (auto c: s) {
//     printf("%c\n", c);
//   }
//   printf("\n");
//   auto report = [](char const &c) -> void { printf("%c\n", c); };
//   std::for_each(s.begin(), s.end(), report);
//   printf("\n");
//   // assert(s[5] == 'f');
//   auto s1 = Span<char>{v.data(), 0};
//   std::for_each(s1.begin(), s1.end(), report);
//   printf("\n");
//   auto s2 = Span<char>{v.data(), 10};
//   auto s3 = s2.first(2);
//   std::for_each(s3.begin(), s3.end(), report);
//   printf("\n");
//   auto s4 = s2.first(10);
//   std::for_each(s4.begin(), s4.end(), report);
//   printf("\n");
//   auto s5 = s2.first(0);
//   std::for_each(s5.begin(), s5.end(), report);
//   printf("\n");
//   auto s6 = s2.last(0);
//   std::for_each(s6.begin(), s6.end(), report);
//   printf("\n");
//   auto s7 = s2.last(2);
//   std::for_each(s7.begin(), s7.end(), report);
//   printf("\n");
//   std::for_each(s7.rbegin(), s7.rend(), report);
//   printf("\n");
//   std::for_each(s7.crbegin(), s7.crend(), report);
//   printf("\n");
// }
