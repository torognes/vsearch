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
#include <algorithm>  // std::find, std::find_if
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstring>  // std::strlen
#include <string>  // std::string
#include <utility>  // std::pair
#include <vector>  // std::vector::clear, push_back


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

/* Name and field, side by side. Order is presentation only: the parser
   matches on the name and returns the enumerator sitting next to it, so
   entries may be inserted or reordered freely. */
constexpr std::size_t userfield_count = 53;

/* Adding a Userfield enumerator without adding its name here would leave the
   new field unparseable from the command line; tie the two counts together
   (thir is the last enumerator). */
static_assert(static_cast<int>(Userfield::thir) + 1 ==
              static_cast<int>(userfield_count),
              "one name per Userfield enumerator");

/* Held by a function-local static rather than a namespace-scope object. The
   hazard this shape guarded against -- the names were std::string, whose
   dynamic initialization could throw where -fno-exceptions cannot catch --
   lifted when the entries became pointers to string literals; the shape is
   kept unchanged, one change at a time. */
auto valid_userfields()
  -> std::array<std::pair<char const *, Userfield>, userfield_count> const &
{
  static std::array<std::pair<char const *, Userfield>, userfield_count> const names =
  {{
    {"query",    Userfield::query},
    {"target",   Userfield::target},
    {"evalue",   Userfield::evalue},
    {"id",       Userfield::id},
    {"pctpv",    Userfield::pctpv},
    {"pctgaps",  Userfield::pctgaps},
    {"pairs",    Userfield::pairs},
    {"gaps",     Userfield::gaps},
    {"qlo",      Userfield::qlo},
    {"qhi",      Userfield::qhi},
    {"tlo",      Userfield::tlo},
    {"thi",      Userfield::thi},
    {"pv",       Userfield::pv},
    {"ql",       Userfield::ql},
    {"tl",       Userfield::tl},
    {"qs",       Userfield::qs},
    {"ts",       Userfield::ts},
    {"alnlen",   Userfield::alnlen},
    {"opens",    Userfield::opens},
    {"exts",     Userfield::exts},
    {"raw",      Userfield::raw},
    {"bits",     Userfield::bits},
    {"aln",      Userfield::aln},
    {"caln",     Userfield::caln},
    {"qstrand",  Userfield::qstrand},
    {"tstrand",  Userfield::tstrand},
    {"qrow",     Userfield::qrow},
    {"trow",     Userfield::trow},
    {"qframe",   Userfield::qframe},
    {"tframe",   Userfield::tframe},
    {"mism",     Userfield::mism},
    {"ids",      Userfield::ids},
    {"qcov",     Userfield::qcov},
    {"tcov",     Userfield::tcov},
    {"id0",      Userfield::id0},
    {"id1",      Userfield::id1},
    {"id2",      Userfield::id2},
    {"id3",      Userfield::id3},
    {"id4",      Userfield::id4},
    {"qilo",     Userfield::qilo},
    {"qihi",     Userfield::qihi},
    {"tilo",     Userfield::tilo},
    {"tihi",     Userfield::tihi},
    {"diffs",    Userfield::diffs},
    {"mid",      Userfield::mid},
    {"qseq",     Userfield::qseq},
    {"tseq",     Userfield::tseq},
    {"qrowdots", Userfield::qrowdots},
    {"trowdots", Userfield::trowdots},
    {"qlor",     Userfield::qlor},
    {"qhir",     Userfield::qhir},
    {"tlor",     Userfield::tlor},
    {"thir",     Userfield::thir},
  },};
  return names;
}

}  // anonymous namespace

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
  auto const & valid_names = valid_userfields();

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

      std::string const field(ptr, static_cast<std::size_t>(field_length));
      auto const * const found =
        std::find_if(valid_names.cbegin(), valid_names.cend(),
                     [&field](std::pair<char const *, Userfield> const & entry) -> bool {
                       return field == entry.first;
                     });

      if (found == valid_names.cend())
        {    // reached end of list -> unrecognized field
          return false; // bad argument
        }

      parameters.opt_userfields.push_back(found->second);

      ptr = next_separator;

      if (ptr == end_of_string)
        {  // reached end of argument
          return true;
        }

      ++ptr;
    }
}
