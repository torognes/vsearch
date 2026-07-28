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
#include "utils/print_view.hpp"  // fprint
#include <algorithm>  // std::find_if, std::search, std::sort
#include <array>
#include <cerrno>  // errno
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstdlib>  // std::strtoll
#include <iterator>  // std::next, std::distance


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct Attribute {
    char const * text = nullptr;
    int length = 0;  // length of the text field
    bool allow_decimal = false;  // integer or float

    constexpr Attribute(char const * new_text, int new_length, bool new_allow_decimal)
      : text(new_text), length(new_length), allow_decimal(new_allow_decimal) {}

    /* the name as a window, for the scans below; the members stay a pointer and
       a length so that the table of attributes remains constexpr (View's
       constructor asserts, so it cannot be constexpr in C++11) */
    auto view() const noexcept -> View<char>
    {
      return View<char>{text, static_cast<std::size_t>(length)};
    }
  };


  struct Attributes {
    Attribute ee {"ee=", 3, true};
    Attribute length {"length=", 7, false};
    Attribute size {"size=", 5, false};
  };

  constexpr Attributes attributes;


  constexpr auto n_expected_attributes = std::size_t{3};  // 3 attributes: size, ee, length


  /* where one attribute sits inside a header: [start, end) in bytes, with
     'start' on the first byte of the name and 'end' one past the last digit of
     its value. Members carry no default initializer, so the struct stays an
     aggregate in C++11 (which does not allow them here) and '{}' zero-fills it
     into a 'not present' span. */
  struct Attribute_span
  {
    std::size_t start;
    std::size_t end;
    bool present;
  };


  /* a byte belonging to an attribute's value: a digit, or the decimal point for
     the attributes that allow one */
  auto is_value_character(char const symbol, bool const allow_decimal) -> bool
  {
    return ((symbol >= '0') and (symbol <= '9'))
      or (allow_decimal and (symbol == '.'));
  }


  auto header_find_attribute(View<char> const header,
                             Attribute const & attribute) -> Attribute_span
  {
    /*
      Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
      in the header string, where "size=" is the specified attribute.
      If allow_decimal is true, a dot (.) is allowed within the digits.
    */

    if ((header.data() == nullptr) or (attribute.text == nullptr))
      {
        return {};
      }

    auto const name = attribute.view();

    auto offset = std::size_t{0};

    /* the bound is written as an addition, not as
       'offset < header.size() - name.size()': in unsigned arithmetic that
       subtraction wraps to a huge value whenever the header is shorter than the
       attribute name, and the loop below would run on out-of-range offsets */
    while (offset + name.size() < header.size())
      {
        /* find the next occurrence of the attribute text, bounded by the
           header's size (no dependence on a trailing '\0') */
        auto const * const first_occurence
          = std::search(std::next(header.cbegin(), static_cast<std::ptrdiff_t>(offset)),
                        header.cend(), name.cbegin(), name.cend());

        /* no match */
        if (first_occurence == header.cend())
          {
            break;
          }

        offset = static_cast<std::size_t>(std::distance(header.cbegin(), first_occurence));

        /* check for ';' in front */
        if ((offset > 0) and (header[offset - 1] != ';'))
          {
            offset += name.size() + 1;
            continue;
          }

        /* count the value's digits, likewise bounded by the header's size */
        auto const value = header.drop(offset + name.size());
        auto const * const value_end =
          std::find_if(value.cbegin(), value.cend(),
                       [&attribute](char const symbol) -> bool
                       { return not is_value_character(symbol, attribute.allow_decimal); });
        auto const digits = static_cast<std::size_t>(std::distance(value.cbegin(), value_end));

        /* check for at least one digit */
        if (digits == 0)
          {
            offset += name.size() + 1;
            continue;
          }

        /* check for ';' after */
        auto const value_end_offset = offset + name.size() + digits;
        if ((value_end_offset < header.size()) and (header[value_end_offset] != ';'))
          {
            offset += name.size() + digits + 2;
            continue;
          }

        /* ok */
        return Attribute_span{offset, value_end_offset, true};
      }
    return {};
  }


}  // end of anonymous namespace


auto header_get_size(View<char> const header) -> int64_t {
  /* read size/abundance annotation */
  static constexpr auto decimal_base = 10;
  auto const annotation = header_find_attribute(header, attributes.size);
  if (not annotation.present) {
    return 0;  // refactoring: return 1 by default?
  }

  char * next_character = nullptr;
  // C++17 refactoring: replace strtoll with std::from_chars
  auto const value_offset = annotation.start + attributes.size.view().size();
  auto const * const value = std::next(header.data(),
                                       static_cast<std::ptrdiff_t>(value_offset));
  auto const abundance = std::strtoll(value, &next_character, decimal_base);
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
  /* the attributes found, ordered by position in the header by the sort below;
     one array of spans, where two parallel start/end arrays used to be kept
     side by side */
  std::array<Attribute_span, n_expected_attributes> found {{}};
  auto nth_attribute = std::size_t{0};

  auto collect = [&](bool const wanted, Attribute const & attribute) -> void {
    if (not wanted) { return; }
    auto const span = header_find_attribute(header_view, attribute);
    if (span.present) {
      found[nth_attribute] = span;
      ++nth_attribute;
    }
  };

  /* look for size attribute */
  collect(strip_size, attributes.size);

  /* look for ee attribute */
  collect(strip_ee, attributes.ee);

  /* look for length attribute */
  collect(strip_length, attributes.length);

  /* sort */

  /* by position in the header; the attribute names differ, so no two spans
     share a start and the order among equals never comes up */
  std::sort(found.begin(),
            std::next(found.begin(), static_cast<std::ptrdiff_t>(nth_attribute)),
            [](Attribute_span const & lhs, Attribute_span const & rhs) -> bool
            { return lhs.start < rhs.start; });

  /* print */

  auto last_emitted = View<char>{};  // the last chunk printed, empty until one is

  auto emit = [output_handle, &last_emitted](View<char> const chunk) -> void {
    if (chunk.empty()) { return; }
    fprint(output_handle, chunk);
    last_emitted = chunk;
  };

  if (nth_attribute == 0)
    {
      emit(header_view);
    }
  else
    {
      auto prev_end = std::size_t{0};
      for (std::size_t i = 0; i < nth_attribute; ++i)
        {
          /* print part of header in front of this attribute */
          if (found[i].start > prev_end + 1)
            {
              emit(header_view.subspan(prev_end, found[i].start - prev_end - 1));
            }
          prev_end = found[i].end;
        }

      /* print the rest, if any */
      if (header_view.size() > prev_end + 1)
        {
          emit(header_view.drop(prev_end));
        }
    }

  /* report whether the last emitted character is the annotation separator */
  return (not last_emitted.empty()) and (last_emitted.back() == ';');
}
