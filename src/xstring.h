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

#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdio>  // std::size_t, std::snprintf
#include <cstring>  // std::strlen, std::strcpy
#include <iterator> // std::prev, std::next


static std::array<char, 1> empty_string = {""};

class xstring
{
private:
  char * string_ {};
  std::size_t length_ {};
  std::size_t alloc_ {};

 public:
  xstring() = default;

  ~xstring() {
    if (capacity() != 0) {
      xfree(string_);
    }
    alloc_ = 0;
    string_ = nullptr;
    length_ = 0;
  }

  // Iterators
  auto end() const -> char * {
    auto const distance = static_cast<std::ptrdiff_t>(size());
    return std::next(data(), distance);
  }

  // Modifiers
  auto clear() -> void {
    length_ = 0;
  }

  auto add_c(char a_char) -> void {
    static constexpr std::size_t needed = 1;
    auto const new_capacity = size() + needed + 1;
    if (new_capacity > capacity()) {
      reserve(new_capacity);
    }
    back() = a_char;
    length_ += 1;
    back() = '\0';
  }

  auto add_d(int a_number) -> void {
    auto const needed = std::snprintf(nullptr, 0, "%d", a_number);
    if (needed < 0) {
      fatal("snprintf failed");
    }

    auto const new_capacity = size() + needed + 1;
    if (new_capacity > capacity()) {
      reserve(new_capacity);
    }
    static_cast<void>(std::snprintf(end(), needed + 1, "%d", a_number));
    length_ += needed;
  }

  // Element access
  auto data() const -> char * {
    if (empty()) {
      return empty_string.data();
    }
    return string_;
  }
  auto back() const -> char & {
    assert(not empty());
    return *std::prev(end());
  }

  // Capacity
  auto capacity() const -> std::size_t { return alloc_; }
  auto empty() const -> bool { return size() == 0; }
  auto reserve(std::size_t const new_capacity) -> void {
    assert(new_capacity > capacity());
    alloc_ = new_capacity;
    string_ = static_cast<char *>(xrealloc(string_, alloc_));
  }
  auto size() const -> std::size_t { return length_; }
};
