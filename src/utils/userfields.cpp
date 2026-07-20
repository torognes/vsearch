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

#include "userfields.hpp"
#include "vsearch.hpp"  // struct Parameters
#include <algorithm>  // std::equal, std::find
#include <cstdint>  // uint64_t
#include <cstring>  // std::strlen
#include <vector>  // std::vector::clear, push_back


// refactoring: C++11 std::array does not allow conversion from litteral string to char *
static const char * userfields_names[] =
  {
    "query",  // 0
    "target", // 1
    "evalue", // 2
    "id",     // 3
    "pctpv",
    "pctgaps",
    "pairs",
    "gaps",
    "qlo",
    "qhi",
    "tlo",
    "thi",
    "pv",
    "ql",
    "tl",
    "qs",
    "ts",
    "alnlen",
    "opens",
    "exts",
    "raw",
    "bits",
    "aln",
    "caln",
    "qstrand",
    "tstrand",
    "qrow",
    "trow",
    "qframe",
    "tframe",
    "mism",
    "ids",
    "qcov",
    "tcov",   // 33
    "id0",
    "id1",
    "id2",
    "id3",
    "id4",    // 38
    "qilo",   // 39
    "qihi",
    "tilo",
    "tihi",   // 42
    nullptr,
  };

auto parse_userfields_arg(char const * arg, struct Parameters & parameters) -> bool
{
  // Parses the userfields option argument, e.g. query+target+id+alnlen+mism
  // and returns true if it is ok or false if not.
  static constexpr auto separator = '+';
  char const * ptr = arg;
  char const * end_of_string = ptr + std::strlen(ptr); // pointer to end of string

  /* Discard any fields left by a previous --userfields (a repeated option on
     the CLI, or a second library session using the same Parameters). */
  parameters.opt_userfields.clear();

  char const * next_separator = nullptr;

  while (true)
    {
      next_separator = std::find(ptr, end_of_string, separator);

      auto const field_length = static_cast<uint64_t>(next_separator - ptr);

      if (field_length == 0)
        {  // empty token (e.g. "a++b", "+a", "a+") -> bad argument. Previously
           // rejected only incidentally by the name-lookup falling through; the
           // explicit check makes the intent clear (L2c).
          return false;
        }

      char const * const * valid_userfield = userfields_names;

      while (*valid_userfield != nullptr)
        {
          if ((std::strlen(*valid_userfield) == field_length) and std::equal(ptr, ptr + field_length, *valid_userfield))
            {
              break;
            }
          ++valid_userfield;
        }

      if (*valid_userfield == nullptr)
        {    // reached end of list -> unrecognized field
          return false; // bad argument
        }

      auto const nth_valid_userfield = static_cast<int>(valid_userfield - userfields_names);
      parameters.opt_userfields.push_back(nth_valid_userfield);

      ptr = next_separator;

      if (ptr == end_of_string)
        {  // reached end of argument
          return true;
        }

      ++ptr;
    }
}
