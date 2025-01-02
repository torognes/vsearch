/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "vsearch.h"
#include <cstdint>  // uint64_t
#include <cstring>  // std::strcmp, std::strchr, std::strlen


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
    "tcov",  // 33
    "id0",
    "id1",
    "id2",
    "id3",
    "id4", // 38
    "qilo", // 39
    "qihi",
    "tilo",
    "tihi", // 42
    nullptr
  };

int * userfields_requested = nullptr;
int userfields_requested_count = 0;


auto parse_userfields_arg(char * arg) -> int
{
  // Parses the userfields option argument, e.g. query+target+id+alnlen+mism
  // and returns 1 if it is ok or 0 if not.
  static constexpr auto separator = '+';
  char * ptr = arg;
  char * end_of_string = ptr + std::strlen(ptr); // pointer to end of string

  // refactoring:
  // auto const userfields_requested_count = std::count(v.cbegin(), v.cend(), '+');
  userfields_requested_count = 1;
  while (ptr < end_of_string)
    {
      if (*ptr == separator)
        {
          ++userfields_requested_count;
        }
      ++ptr;
    }

  userfields_requested = static_cast<int *>(xmalloc(sizeof(int) * (uint64_t) userfields_requested_count));

  ptr = arg;  // reset to the start of the string

  char * next_separator = nullptr;

  auto nth_field = 0;

  while (true)
    {
      next_separator = std::strchr(ptr, separator);
      if (next_separator == nullptr)
        {
          next_separator = end_of_string;
        }

      auto const field_length = static_cast<uint64_t>(next_separator - ptr);

      char ** valid_userfield = (char **) userfields_names;

      while (*valid_userfield != nullptr)
        {
          if ((std::strncmp(ptr, *valid_userfield, field_length) == 0) and (std::strlen(*valid_userfield) == field_length))
            {
              break;
            }
          ++valid_userfield;
        }

      if (*valid_userfield == nullptr)
        {    // reached end of list -> unrecognized field
          return 0; // bad argument
        }

      int const i = (int) (((const char **) valid_userfield) - userfields_names);
      userfields_requested[nth_field] = i;
      ++nth_field;

      ptr = next_separator;

      if (ptr == end_of_string)
        {  // reached end of argument
          return 1;
        }

      ++ptr;
    }
}
