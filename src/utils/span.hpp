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
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdlib>  // std::size_t
#include <iterator> // std::next

#ifndef NDEBUG
#include <limits>
// C++17 refactoring: [[maybe_unused]]
constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
#endif


// simple version of std::span
// only valid for vectors or arrays of chars

class Span {
public:
  Span(char * start, std::size_t length)
    : m_start {start},
      m_length {length} {
    assert(start != nullptr);
    assert(length <= max_ptrdiff);
  }

  // Iterators
  auto begin() const -> char * { return m_start; }
  auto cbegin() const -> char const * { return m_start; }
  auto end() const -> char * {
    auto const distance = static_cast<std::ptrdiff_t>(m_length);
    return std::next(m_start, distance);
  }
  auto cend() const -> char const * {
    auto const distance = static_cast<std::ptrdiff_t>(m_length);
    return std::next(m_start, distance);
  }

  // Element access
  auto front() const -> char const & {
    assert(not empty());
    return *m_start;
  }
  auto back() const -> char const & {
    assert(not empty());
    return *std::prev(end());
  }
  auto data() const -> char * { return m_start; }
  auto operator[](std::size_t index) -> char & {
    assert(index < m_length);
    auto const distance = static_cast<std::ptrdiff_t>(index);
    return *std::next(m_start, distance);
  }

  // Observers
  auto size() const -> std::size_t { return m_length; }
  auto size_bytes() const -> std::size_t { return m_length * sizeof(char); }
  auto empty() const -> bool { return m_length == 0; }

private:
  char * m_start {};
  std::size_t m_length {};
};


// tests:

// #include <algorithm>
// #include <vector>

// auto main() -> int {
//   std::vector<char> v = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
//   auto s = Span{v.data(), 5};
//   assert(s.size() == 5);
//   assert(s.size_bytes() == 5);
//   assert(! s.empty());
//   assert(s.front() == 'a');
//   assert(s.back() == 'e');
//   assert(s[1] == 'b');
//   for (auto c: s) {
//     printf("%c\n", c);
//   }
//   std::for_each(s.begin(), s.end(), [](char &c) -> void { printf("%c\n", c); });
//   auto s1 = Span{v.data(), 0};
//   assert(s[5] == 'f');
// }
