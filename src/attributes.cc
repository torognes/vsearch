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
#include <algorithm>  // std::swap
#include <array>
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstdlib>  // std::atol
#include <cstring>  // std::strlen, std::strstr, std::strspn

constexpr auto n_expected_attributes = std::size_t{3};  // 3 attributes: size, ee, length

auto header_find_attribute(const char * header,
                           int header_length,
                           const char * attribute,
                           int * start,
                           int * end,
                           bool allow_decimal) -> bool
{
  /*
    Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
    in the header string, where "size=" is the specified attribute.
    If allow_decimal is true, a dot (.) is allowed within the digits.
  */

  const char * digit_chars = "0123456789";
  const char * digit_chars_decimal = "0123456789.";

  if ((header == nullptr) or (attribute == nullptr))
    {
      return false;
    }

  auto const hlen = header_length;
  auto const alen = static_cast<int>(std::strlen(attribute));

  auto i = 0;

  while (i < hlen - alen)
    {
      auto * r = (char *) std::strstr(header + i, attribute);

      /* no match */
      if (r == nullptr)
        {
          break;
        }

      i = r - header;

      /* check for ';' in front */
      if ((i > 0) and (header[i - 1] != ';'))
        {
          i += alen + 1;
          continue;
        }

      auto const digits
        = (int) std::strspn(header + i + alen,
                       (allow_decimal ? digit_chars_decimal : digit_chars));

      /* check for at least one digit */
      if (digits == 0)
        {
          i += alen + 1;
          continue;
        }

      /* check for ';' after */
      if ((i + alen + digits < hlen) and (header[i + alen + digits] != ';'))
        {
          i += alen + digits + 2;
          continue;
        }

      /* ok */
      *start = i;
      *end = i + alen + digits;
      return true;
    }
  return false;
}


auto header_get_size(char * header, int header_length) -> int64_t {
  /* read size/abundance annotation */
  auto start = 0;
  auto end = 0;
  auto const attribute_is_present = header_find_attribute(header,
                                                          header_length,
                                                          "size=",
                                                          &start,
                                                          &end,
                                                          false);
  if (not attribute_is_present) {
    return 0;
  }

  auto const abundance = std::atol(header + start + 5);

  if (abundance == 0) {
    fatal("Invalid (zero) abundance annotation in FASTA file header");
  }

  return abundance;
}


auto look_for_attribute(char const *header, int const header_length,
                        int &nth_attribute, std::array<int, n_expected_attributes> &attribute_start,
                        std::array<int, n_expected_attributes> &attribute_end,
                        char const* attribute_text,
                        bool const strip_attribute) -> void {
  auto start = 0;
  auto end = 0;
  if (not strip_attribute) { return; }

  auto const attribute_is_present = header_find_attribute(header,
                                                          header_length,
                                                          attribute_text,
                                                          & start,
                                                          & end,
                                                          false);
  if (not attribute_is_present) { return; }
  attribute_start[nth_attribute] = start;
  attribute_end[nth_attribute] = end;
  ++nth_attribute;
}


auto header_fprint_strip(FILE * output_handle,
                         char * header,
                         int header_length,
                         bool strip_size,
                         bool strip_ee,
                         bool strip_length) -> void
{
  auto nth_attribute = 0;
  std::array<int, n_expected_attributes> attribute_start {{}};
  std::array<int, n_expected_attributes> attribute_end {{}};

  /* look for size attribute */
  look_for_attribute(header, header_length,
                     nth_attribute, attribute_start,
                     attribute_end,
                     "size=", strip_size);

  /* look for ee attribute */
  look_for_attribute(header, header_length,
                     nth_attribute, attribute_start,
                     attribute_end,
                     "ee=", strip_ee);

  /* look for length attribute */
  look_for_attribute(header, header_length,
                     nth_attribute, attribute_start,
                     attribute_end,
                     "length=", strip_length);

  /* sort */

  auto last_swap = 0;
  auto limit = nth_attribute - 1;
  while (limit > 0)
    {
      for(auto i = 0; i < limit; i++)
        {
          if (attribute_start[i] > attribute_start[i + 1])
            {
              std::swap(attribute_start[i], attribute_start[i + 1]);
              std::swap(attribute_end[i], attribute_end[i + 1]);
              last_swap = i;
            }
        }
      limit = last_swap;
    }

  /* print */

  if (nth_attribute == 0)
    {
      std::fprintf(output_handle, "%.*s", header_length, header);
    }
  else
    {
      auto prev_end = 0;
      for (auto i = 0; i < nth_attribute; i++)
        {
          /* print part of header in front of this attribute */
          if (attribute_start[i] > prev_end + 1)
            {
              std::fprintf(output_handle, "%.*s",
                      attribute_start[i] - prev_end - 1,
                      header + prev_end);
            }
          prev_end = attribute_end[i];
        }

      /* print the rest, if any */
      if (header_length > prev_end + 1)
        {
          std::fprintf(output_handle, "%.*s",
                  header_length - prev_end,
                  header + prev_end);
        }
    }
}
