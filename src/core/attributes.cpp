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

#include "utils/view.hpp"
#include "core/attributes.hpp"  // View<char>
#include "utils/fatal.hpp"
#include <algorithm>  // std::swap
#include <array>
#include <cerrno>  // errno
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstdlib>  // std::strtoll


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct Attribute {
    char const * text = nullptr;
    int length = 0;  // length of the text field
    bool allow_decimal = false;  // integer or float

    constexpr Attribute(char const * new_text, int new_length, bool new_allow_decimal)
      : text(new_text), length(new_length), allow_decimal(new_allow_decimal) {}

  };


  struct Attributes {
    Attribute ee {"ee=", 3, true};
    Attribute length {"length=", 7, false};
    Attribute size {"size=", 5, false};
  };

  constexpr Attributes attributes;


  constexpr auto n_expected_attributes = std::size_t{3};  // 3 attributes: size, ee, length


  auto header_find_attribute(char const * header,
                             int const header_length,
                             Attribute const attribute,
                             int * start,
                             int * end) -> bool
  {
    /*
      Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
      in the header string, where "size=" is the specified attribute.
      If allow_decimal is true, a dot (.) is allowed within the digits.
    */

    if ((header == nullptr) or (attribute.text == nullptr))
      {
        return false;
      }

    auto const * const header_end = header + header_length;
    auto const * const attribute_text_end = attribute.text + attribute.length;

    auto offset = 0;

    while (offset < header_length - attribute.length)
      {
        /* find the next occurrence of the attribute text, bounded by
           header_length (no dependence on a trailing '\0') */
        auto const * const first_occurence
          = std::search(header + offset, header_end, attribute.text, attribute_text_end);

        /* no match */
        if (first_occurence == header_end)
          {
            break;
          }

        offset = static_cast<int>(first_occurence - header);

        /* check for ';' in front */
        if ((offset > 0) and (header[offset - 1] != ';'))
          {
            offset += attribute.length + 1;
            continue;
          }

        /* count the value's digits, likewise bounded by header_length */
        auto const * value_it = header + offset + attribute.length;
        while ((value_it < header_end) and
               (((*value_it >= '0') and (*value_it <= '9')) or
                (attribute.allow_decimal and (*value_it == '.'))))
          {
            ++value_it;
          }
        auto const digits = static_cast<int>(value_it - (header + offset + attribute.length));

        /* check for at least one digit */
        if (digits == 0)
          {
            offset += attribute.length + 1;
            continue;
          }

        /* check for ';' after */
        if ((offset + attribute.length + digits < header_length) and (header[offset + attribute.length + digits] != ';'))
          {
            offset += attribute.length + digits + 2;
            continue;
          }

        /* ok */
        *start = offset;
        *end = offset + attribute.length + digits;
        return true;
      }
    return false;
  }


  auto look_for_attribute(char const * header, int const header_length,
                          int & nth_attribute, std::array<int, n_expected_attributes> &attribute_start,
                          std::array<int, n_expected_attributes> &attribute_end,
                          Attribute const attribute) -> void {
    auto start = 0;
    auto end = 0;

    auto const attribute_is_present = header_find_attribute(header,
                                                            header_length,
                                                            attribute,
                                                            & start,
                                                            & end);
    if (not attribute_is_present) { return; }
    attribute_start[static_cast<std::size_t>(nth_attribute)] = start;
    attribute_end[static_cast<std::size_t>(nth_attribute)] = end;
    ++nth_attribute;
  }


}  // end of anonymous namespace


auto header_get_size(char const * header, int const header_length) -> int64_t {
  /* read size/abundance annotation */
  static constexpr auto decimal_base = 10;
  auto start = 0;
  auto end = 0;
  auto const attribute_is_present = header_find_attribute(header,
                                                          header_length,
                                                          attributes.size,
                                                          &start,
                                                          &end);
  if (not attribute_is_present) {
    return 0;  // refactoring: return 1 by default?
  }

  char * next_character = nullptr;
  // C++17 refactoring: replace strtoll with std::from_chars
  auto const abundance = std::strtoll(header + start + attributes.size.length, &next_character, decimal_base);
  auto const range_error = (errno == ERANGE);

  if (range_error) {
    fatal("Invalid (range error) abundance annotation in FASTA file header");
  }

  if (abundance == 0) {
    fatal("Invalid (zero) abundance annotation in FASTA file header");
  }

  return abundance;
}


auto annotation_separator(bool & trailing_separator) -> char const * {
  /*
    Return the separator to prepend to the next annotation. When the text
    printed so far already ends with the separator ';' (e.g. a header or a
    label suffix ending with ';'), reuse it rather than emit a second one,
    so that annotations are merged with a single ';' (see issue #271).
  */
  if (trailing_separator) {
    trailing_separator = false;
    return "";
  }
  return ";";
}


auto header_fprint_strip(std::FILE * output_handle,
                         View<char> const header_view,
                         bool const strip_size,
                         bool const strip_ee,
                         bool const strip_length) -> bool
{
  auto const * const header = header_view.data();
  auto const header_length = static_cast<int>(header_view.size());

  auto nth_attribute = 0;
  std::array<int, n_expected_attributes> attribute_start {{}};
  std::array<int, n_expected_attributes> attribute_end {{}};

  /* look for size attribute */
  if (strip_size) {
    look_for_attribute(header, header_length,
                       nth_attribute, attribute_start,
                       attribute_end,
                       attributes.size);
  }

  /* look for ee attribute */
  if (strip_ee) {
    look_for_attribute(header, header_length,
                       nth_attribute, attribute_start,
                       attribute_end,
                       attributes.ee);
  }

  /* look for length attribute */
  if (strip_length) {
    look_for_attribute(header, header_length,
                       nth_attribute, attribute_start,
                       attribute_end,
                       attributes.length);
  }

  /* sort */

  auto last_swap = 0;
  auto limit = nth_attribute - 1;
  while (limit > 0)
    {
      for (auto i = 0; i < limit; ++i)
        {
          if (attribute_start[static_cast<std::size_t>(i)] > attribute_start[static_cast<std::size_t>(i + 1)])
            {
              std::swap(attribute_start[static_cast<std::size_t>(i)], attribute_start[static_cast<std::size_t>(i + 1)]);
              std::swap(attribute_end[static_cast<std::size_t>(i)], attribute_end[static_cast<std::size_t>(i + 1)]);
              last_swap = i;
            }
        }
      limit = last_swap;
    }

  /* print */

  auto last_index = -1;  // index in 'header' of the last emitted character

  if (nth_attribute == 0)
    {
      std::fprintf(output_handle, "%.*s", header_length, header);
      if (header_length > 0) { last_index = header_length - 1; }
    }
  else
    {
      auto prev_end = 0;
      for (auto i = 0; i < nth_attribute; ++i)
        {
          /* print part of header in front of this attribute */
          if (attribute_start[static_cast<std::size_t>(i)] > prev_end + 1)
            {
              std::fprintf(output_handle, "%.*s",
                      attribute_start[static_cast<std::size_t>(i)] - prev_end - 1,
                      header + prev_end);
              last_index = attribute_start[static_cast<std::size_t>(i)] - 2;
            }
          prev_end = attribute_end[static_cast<std::size_t>(i)];
        }

      /* print the rest, if any */
      if (header_length > prev_end + 1)
        {
          std::fprintf(output_handle, "%.*s",
                  header_length - prev_end,
                  header + prev_end);
          last_index = header_length - 1;
        }
    }

  /* report whether the last emitted character is the annotation separator */
  return (last_index >= 0) and (header[last_index] == ';');
}
