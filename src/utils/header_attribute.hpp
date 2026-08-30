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

#include "utils/view.hpp"  // View
#include <algorithm>  // std::find_if, std::search
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::distance, std::next


/* One scanner for the annotation grammar

     (^|;)name=value(;|$)

   which vsearch reads from FASTA/FASTQ headers under seven names. Before this
   header the same grammar had four implementations: the table-driven scan in
   core/attributes.cpp (size=, ee=, length=), a near-copy in core/tax.cpp whose
   own comment said "very similar to header_find_attribute()" (tax=), and two
   regular-expression paths in core/otutable.cpp -- POSIX regcomp/regexec where
   <regex.h> exists, std::regex where it does not, which on the targets we
   cross-compile means everywhere except Windows (sample=, barcodelabel=, otu=,
   tax=).

   Four implementations drifted, as four implementations do: core/tax.cpp could
   not see a "tax=" whose '=' was the header's last byte, where the regular
   expression could, so the two disagreed on ";tax=" with an empty value. See
   TBD_20260827_attribute_reading_unification.md for the census, the divergence
   table, and the measurement showing which side of it was observable.

   Only the value's character set was ever an intended difference between the
   four, and that is what Value_chars names.

   Header-only because two of the consumers are per-record: header_get_size()
   runs once per input record under --sizein, and header_fprint_strip() scans up
   to three times per printed record. */
namespace vsearch
{

  /* which bytes a value may be made of -- the one deliberate difference
     between the callers. 'digits' and 'digits_and_dot' stop at anything else,
     so they also stop at a ';'; 'not_semicolon' takes everything up to the
     separator or the header's end. */
  enum struct Value_chars { digits, digits_and_dot, not_semicolon };


  /* One attribute's name and value policy.

     No default member initializers and no constructor, both deliberately: an
     NSDMI makes a class a non-aggregate in C++11, which would force a
     hand-written constexpr constructor before a table of these could be brace
     initialized. As a plain aggregate it needs neither, and the tables below
     stay constexpr on GCC 4.9.

     The name is kept as a pointer and a length rather than a View for the same
     reason: View's constructor asserts, so it cannot be constexpr in C++11. */
  struct Attribute
  {
    char const * text;
    int length;  // length of the text field
    Value_chars value_chars;
    bool allow_empty;  // is "name=" with nothing after it a match?

    /* the name as a window, for the scans below */
    auto view() const noexcept -> View<char>
    {
      return View<char>{text, static_cast<std::size_t>(length)};
    }
  };


  /* where one attribute sits inside a header: [start, end) in bytes, with
     'start' on the first byte of the name and 'end' one past the last byte of
     its value. Members carry no default initializer, so the struct stays an
     aggregate in C++11 (which does not allow them here) and '{}' zero-fills it
     into a 'not present' span. */
  struct Attribute_span
  {
    std::size_t start;
    std::size_t end;
    bool present;
  };


  /* a byte belonging to an attribute's value, under that attribute's policy */
  inline auto is_value_character(char const symbol,
                                 Value_chars const policy) noexcept -> bool
  {
    switch (policy)
      {
      case Value_chars::digits:
        return (symbol >= '0') and (symbol <= '9');
      case Value_chars::digits_and_dot:
        return ((symbol >= '0') and (symbol <= '9')) or (symbol == '.');
      case Value_chars::not_semicolon:
        return symbol != ';';
      }
    /* not a default case on purpose: an exhaustive switch over an enum struct
       is what makes -Wswitch report a policy added later. Unreachable, but the
       return is what makes the function's type honest to -Wreturn-type. */
    return false;
  }


  inline auto header_find_attribute(View<char> const header,
                                    Attribute const & attribute) -> Attribute_span
  {
    /*
      Identify the first occurrence of the pattern (^|;)size=([0-9]+)(;|$)
      in the header string, where "size=" is the specified attribute and the
      value's character set is the attribute's own (see Value_chars).
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
       attribute name, and the loop below would run on out-of-range offsets.

       '<=', not '<': a name whose '=' is the header's last byte leaves no room
       for a value, and the '<' bound made that position unreachable -- both
       for a header that is exactly "tax=" and for one where a rejected
       candidate has already pushed 'offset' there ("xtax=;tax=" pushes it to
       6, which is also header.size() - name.size()). That is the empty-value
       case, and whether it matches is 'allow_empty''s decision to make, not
       the bound's. For the three numeric annotations, which forbid an empty
       value, the length == 0 test below still rejects it and nothing
       changes. */
    while (offset + name.size() <= header.size())
      {
        /* find the next occurrence of the attribute text, bounded by the
           header's size (no dependence on a trailing '\0') */
        auto const * const first_occurrence
          = std::search(std::next(header.cbegin(), static_cast<std::ptrdiff_t>(offset)),
                        header.cend(), name.cbegin(), name.cend());

        /* no match */
        if (first_occurrence == header.cend())
          {
            break;
          }

        offset = static_cast<std::size_t>(std::distance(header.cbegin(), first_occurrence));

        /* check for ';' in front */
        if ((offset > 0) and (header[offset - 1] != ';'))
          {
            offset += name.size() + 1;
            continue;
          }

        /* measure the value, likewise bounded by the header's size */
        auto const value = header.drop(offset + name.size());
        auto const * const value_end =
          std::find_if(value.cbegin(), value.cend(),
                       [&attribute](char const symbol) -> bool
                       { return not is_value_character(symbol, attribute.value_chars); });
        auto const length = static_cast<std::size_t>(std::distance(value.cbegin(), value_end));

        /* check for at least one value byte, unless the attribute allows none */
        if ((length == 0) and (not attribute.allow_empty))
          {
            offset += name.size() + 1;
            continue;
          }

        /* check for ';' after; vacuous under Value_chars::not_semicolon, whose
           value stops only at a ';' or at the header's end */
        auto const value_end_offset = offset + name.size() + length;
        if ((value_end_offset < header.size()) and (header[value_end_offset] != ';'))
          {
            offset += name.size() + length + 2;
            continue;
          }

        /* ok */
        return Attribute_span{offset, value_end_offset, true};
      }
    return {};
  }


  /* The leftmost of several names that fill one slot, by position in the
     header rather than by position in the table.

     This exists for --otutabout's "(sample|barcodelabel)=" alternation: POSIX
     alternation is leftmost-first over the subject, not first-listed over the
     pattern, so a header carrying both annotations is answered by whichever
     comes first in it. Probing the names in order and taking the first one
     found would answer with 'sample=' every time, which is a different
     function -- and a regression against every released vsearch. */
  inline auto header_find_first_attribute(View<char> const header,
                                          View<Attribute> const alternatives) -> Attribute_span
  {
    auto leftmost = Attribute_span{};
    for (auto const & alternative : alternatives)
      {
        auto const span = header_find_attribute(header, alternative);
        if (not span.present)
          {
            continue;
          }
        if ((not leftmost.present) or (span.start < leftmost.start))
          {
            leftmost = span;
          }
      }
    return leftmost;
  }


  /* the value alone, for the callers that want it rather than the whole cut.
     Empty both when the attribute is absent and when it is present with an
     empty value; the callers that need to tell those apart read span.present
     (see core/otutable.cpp, where the difference decides a column). */
  inline auto attribute_value(View<char> const header,
                              Attribute_span const span,
                              Attribute const & attribute) noexcept -> View<char>
  {
    if (not span.present)
      {
        return {};
      }
    auto const name_length = static_cast<std::size_t>(attribute.length);
    auto const value_start = span.start + name_length;
    assert(value_start <= span.end);
    assert(span.end <= header.size());
    return header.subspan(value_start, span.end - value_start);
  }

}  // namespace vsearch
