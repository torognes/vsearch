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

#include "cigar_operations.hpp"
#include "fatal.hpp"
#include "span.hpp"
#include <algorithm>  // std::max
#include <cassert>
#include <cstdio>  // std::size_t
#include <cstdlib>  // std::strtoll
#include <iterator>  // std::next
#include <vector>
#include <utility>  // std::pair

#ifndef NDEBUG
#include <limits>
#endif


// CIGAR string example: 3M2I3MD
// document the format here, and in vsearch-cigar(5)
// also 'X' (mismatch) and 'N'


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto convert_to_operation(char const operation) -> Operation {
    assert(operation == 'M' or operation == 'I' or operation == 'D');
    if (operation == 'I') {
      return Operation::insertion;
    }
    if (operation == 'D') {
      return Operation::deletion;
    }
    return Operation::match;
  }


  auto convert_from_operation(Operation const operation) -> char {
    switch (operation) {
    case Operation::match:
      return 'M';
      break;
    case Operation::deletion:
      return 'D';
      break;
    case Operation::insertion:
      return 'I';
      break;
    }
    // C++23 refactoring: default: std::unreachable();
    __builtin_unreachable();
  }

}  // end of anonymous namespace


// duplicate: msa.cc
auto find_runlength_of_leftmost_operation(char const * first_character,
                                          char ** first_non_digit) -> long long {
  // std::strtoll:
  // - start from the 'first_character' pointed to,
  // - consume as many characters as possible to form a valid integer,
  // - advance pointer to the first non-digit character,
  // - return the valid integer
  // - if there is no valid integer: pointer is not advanced and strtoll() returns zero
  static constexpr auto decimal_base = 10;
  auto const runlength = std::strtoll(first_character,
                                      first_non_digit,
                                      decimal_base);
  assert(runlength <= std::numeric_limits<int>::max());

  // in cigar strings, runlength of 1 are implicit (no digit)
  return std::max(runlength, 1LL);  // is in [1, INT_MAX]
}


// refactoring C++23: std::pair generator
auto parse_cigar_string(Span<char> const cigar_string) -> std::vector<std::pair<Operation, long long>> {
  std::vector<std::pair<Operation, long long>> parsed_cigar;

  auto * position = cigar_string.begin();
  auto * cigar_end = cigar_string.end();

  while (position < cigar_end)
    {
      // Consume digits (if any), return the position of the
      // first char (M, D, or I), store it, move cursor to the next byte.
      auto ** next_operation = &position;
      auto const run = find_runlength_of_leftmost_operation(position, next_operation);
      // do not dereference if outside of cigar_string! (= missing operation!)
      if (*next_operation >= cigar_end) {
        // fail if ill-formed (ex: '12M1'), we could also silently skip over
        fatal("ill-formed CIGAR string");
      }
      // operations: match (M), insertion (I), or deletion (D)
      auto const operation = **next_operation;
      position = std::next(*next_operation);
      parsed_cigar.emplace_back(convert_to_operation(operation), run);
    }
  return parsed_cigar;
}


auto parse_cigar_string_char(Span<char> const cigar_string) -> std::vector<std::pair<char, long long>> {
  std::vector<std::pair<char, long long>> parsed_cigar;
  auto const cigar_pairs = parse_cigar_string(cigar_string);
  parsed_cigar.reserve(cigar_pairs.size());
  for (auto const & a_pair: cigar_pairs) {
    auto const operation = convert_from_operation(a_pair.first);
    auto const runlength = a_pair.second;
    parsed_cigar.emplace_back(operation, runlength);
  }
  return parsed_cigar;
}


auto print_uncompressed_cigar(std::FILE * output_handle, Span<char> const cigar_string) -> void {
  auto const cigar_pairs = parse_cigar_string_char(cigar_string);
  for (auto const & a_pair: cigar_pairs) {
    auto const operation = a_pair.first;
    auto const runlength = a_pair.second;
    // refactoring? std::fprintf("%s", std::string(runlength, operation).c_str());
    for (auto i = 0LL; i < runlength; ++i) {
      static_cast<void>(std::fprintf(output_handle, "%c", operation));
    }
  }
}


// compiler explorer tests:

// auto main() -> int {
//     // basic example
//     std::vector<char> ex1 = {'2', 'I'};
//     auto res1 = parse_cigar_string(ex1);
//     for (auto const & pair: res1) {
//         std::cout << pair.first << '\t' << pair.second << '\n';
//     }
//     std::cout << '\n';

//     // empty cigar
//     std::vector<char> ex2 = {};
//     auto res2 = parse_cigar_string(ex2);
//     if (res2.empty()) {
//         std::cout << "empty vector\n";
//     }
//     std::cout << '\n';

//     // ommitted run
//     std::vector<char> ex3 = {'M'};
//     auto res3 = parse_cigar_string(ex3);
//     for (auto const & pair: res3) {
//         std::cout << pair.first << '\t' << pair.second << '\n';
//     }
//     std::cout << '\n';

//     // long run
//     std::vector<char> ex4 = {'1', '2', '3', '4', 'M'};
//     auto res4 = parse_cigar_string(ex4);
//     for (auto const & pair: res4) {
//         std::cout << pair.first << '\t' << pair.second << '\n';
//     }
//     std::cout << '\n';

//     // null run
//     std::vector<char> ex5 = {'0', 'M'};
//     auto res5 = parse_cigar_string(ex5);
//     for (auto const & pair: res5) {
//         std::cout << pair.first << '\t' << pair.second << '\n';
//     }
//     std::cout << '\n';

//    // chained operations
//     std::vector<char> ex6 = {'M', 'I', 'D'};
//     auto res6 = parse_cigar_string(ex6);
//     for (auto const & pair: res6) {
//         std::cout << pair.first << '\t' << pair.second << '\n';
//     }
//     std::cout << '\n';

//    // no operation?
//     std::vector<char> ex7 = {'1'};
//     auto res7 = parse_cigar_string(ex7);
//     if (res7.empty()) {
//         std::cout << "empty vector\n";
//     }
//     std::cout << '\n';
// }
