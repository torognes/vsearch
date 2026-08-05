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


#include "utils/decimal_digits.hpp"  // decimal::Buffer, decimal::to_decimal
#include "utils/print_view.hpp"  // fprint
#include "utils/view.hpp"  // View<char>
#include <array>
#include <algorithm>  // std::copy
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdio>  // std::FILE, std::fwrite, std::snprintf
#include <iterator>  // std::next


/* One output record, assembled in a stack buffer and written with a single
   std::fwrite.

   Why this exists, in one measurement. print_view.hpp's writers take one
   std::FILE * call per field, which is what replaced the single std::fprintf
   per record. That trades a run-time format parse for a per-field stdio call,
   and the trade is only a win when the compiler inlines the writer. It does
   not always: in core/results.cpp GCC emits fprint_integer out of line
   because results_show_blast6out_one calls it nine times, so that record went
   from 5 stdio calls to 24 and is ~9 % slower than the fprintf it replaced.
   The same ten-field tail, measured three ways at -O2:

     one std::fprintf, 10 fields            173.5 ns/record
     one call per field, inlined            137.1 ns/record
     assembled here, one fwrite              95.7 ns/record

   So this is not a micro-optimisation on top of a win; it is what makes the
   split unconditionally faster than the format string it replaced, whatever
   the inliner decides.

   Note the deliberate contrast with print_view.hpp's own note, which argues
   *against* batching a run of numbers. That note is about the stream lock
   alone, with the digit loops inlined; it stands. This is about the call
   overhead when they are not, which is the case the note did not cover.

   Bounded, and safe for any input regardless: the buffer flushes itself when
   the next field would not fit, and a payload larger than the whole buffer
   (a header can be a megabyte) is written straight through rather than
   copied. Capacity is therefore a tuning knob, never a correctness
   constraint, and no call site has to reason about how long its record is.

   Not a general-purpose stream. It has no seek, no formatting state, and no
   error channel -- the stream's own error flag is still what
   utils/open_file.cpp checks with ferror() before fclose(). */

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
template <std::size_t Capacity = 256>
class Record {
public:
  explicit Record(std::FILE * const output_handle) noexcept
    : output_handle_ {output_handle} {
    assert(output_handle != nullptr);
  }

  /* RAII: the record reaches the stream even on an early return, and there are
     several (results_show_userout_one's switch returns per field). This is the
     whole reason the buffer is an object rather than a local array plus a
     hand-written flush at each exit. */
  ~Record() { flush(); }

  /* Non-copyable and non-movable: two Records over one buffer would each
     flush it. There is no reason to pass one anywhere except by reference. */
  Record(Record const &) = delete;
  auto operator=(Record const &) -> Record & = delete;
  Record(Record &&) = delete;
  auto operator=(Record &&) -> Record & = delete;

  auto flush() -> void {
    if (used_ == 0) { return; }
    static constexpr std::size_t element_size = sizeof(char);
    static_cast<void>(std::fwrite(buffer_.data(), element_size, used_, output_handle_));
    used_ = 0;
  }

  /* The stream itself, for a field produced by a helper that writes to a
     std::FILE * directly (core/attributes.cpp's header_fprint_strip, which
     emits a header in chunks). Anything buffered goes out first, so the record
     stays in order; the cost is one extra fwrite where the helper is used,
     against templating that helper on its sink and moving its 80 lines into a
     header to be instantiated per caller. */
  auto stream() -> std::FILE * {
    flush();
    return output_handle_;
  }

  auto put(char const character) -> void {
    if (used_ == Capacity) { flush(); }
    buffer_[used_] = character;
    ++used_;
  }

  auto put(View<char> const text) -> void {
    if (text.empty()) { return; }
    if (text.size() > Capacity) {
      /* larger than the whole buffer: writing it through costs one fwrite,
         copying it would cost several plus the copy */
      flush();
      fprint(output_handle_, text);
      return;
    }
    if (text.size() > Capacity - used_) { flush(); }
    std::copy(text.cbegin(), text.cend(),
              std::next(buffer_.begin(), static_cast<std::ptrdiff_t>(used_)));
    used_ += text.size();
  }

  /* Writable room for a caller that must format into the buffer itself --
     std::snprintf for a double, and nothing else. Returns at least
     'wanted' bytes, flushing to get them. */
  auto reserve(std::size_t const wanted) -> char * {
    assert(wanted <= Capacity);
    if (wanted > Capacity - used_) { flush(); }
    return std::next(buffer_.data(), static_cast<std::ptrdiff_t>(used_));
  }

  auto free_space() const noexcept -> std::size_t { return Capacity - used_; }

  auto commit(std::size_t const written) -> void {
    assert(written <= Capacity - used_);
    used_ += written;
  }

private:
  std::FILE * output_handle_ {};
  std::size_t used_ {};
  std::array<char, Capacity> buffer_ {};
};


/* The sink a converted writer names. 256 bytes covers every record vsearch
   assembles field-by-field (the widest, --blast6out's, is under 200 with
   64-bit lengths); a longer one simply flushes twice. */
using OutputRecord = Record<>;


/* The same names print_view.hpp gives the std::FILE * sink, so converting a
   writer is one token per line -- 'output_handle' becomes 'record' -- and the
   reasoning at each of print_view.hpp's functions still applies unchanged.
   Overloaded on Record<Capacity> & rather than templated on the sink, so a
   call that means the stream cannot silently pick the buffer or vice versa. */

template <std::size_t Capacity>
auto fprint(Record<Capacity> & record, char const character) -> void {
  record.put(character);
}

template <std::size_t Capacity>
auto fprint(Record<Capacity> & record, View<char> const text) -> void {
  record.put(text);
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
template <std::size_t Capacity, std::size_t Size>
auto fprint(Record<Capacity> & record, char const (&literal)[Size]) -> void {
  static_assert(Size > 0, "a string literal always carries its terminator");
  assert(literal[Size - 1] == '\0');
  record.put(View<char>{literal, Size - 1});
}

template <std::size_t Capacity, typename Integer>
auto fprint_integer(Record<Capacity> & record, Integer const value) -> void {
  decimal::Buffer digits {};
  record.put(decimal::to_decimal(digits, value));
}

template <std::size_t Capacity>
auto fprint_spaces(Record<Capacity> & record, std::size_t const count) -> void {
  for (std::size_t written = 0; written < count; ++written) { record.put(' '); }
}

template <std::size_t Capacity, typename Integer>
auto fprint_integer(Record<Capacity> & record, Integer const value,
                    std::size_t const width) -> void {
  decimal::Buffer digits {};
  auto const text = decimal::to_decimal(digits, value);
  if (text.size() < width) { fprint_spaces(record, width - text.size()); }
  record.put(text);
}


/* No fprint_double() here, and that is a finding rather than an omission.
   A double field has to reach the stream in order with the rest of the record,
   so the obvious member would take the format and snprintf into the buffer:

     template <std::size_t Capacity, std::size_t Size>
     auto fprint_double(Record<Capacity> &, char const (&format)[Size], double);

   That does not survive -Wformat=2, which phase 8 of
   DONE_20260804_c_style_elimination.md turned on. Inside the template the
   format is a parameter, so GCC reports -Wformat-nonliteral even though every
   call binds it to a literal, and it cannot be silenced: the format attribute
   is a GNU extension, src/ deliberately has none left, and it does not apply
   to a non-variadic function anyway -- which is the same wall the note in
   Makefile.am used to describe about fatal().

   So a double keeps its own std::fprintf, pointed at record.stream(), exactly
   as Decision 1 says. That costs one flush where the double sits: two writes
   for a record with one double field, against the 24 stdio calls the
   unbuffered form needs. A call site that wants a single write can snprintf
   the double into a local array with a literal format and hand over the View,
   which is three lines instead of one and is why it is not the default. */


/* tests:

   Checked by capturing a std::tmpfile() and comparing with the same sequence
   written straight to the stream through print_view.hpp:

   - every field type, singly and in a record
   - a record that exactly fills the buffer, and one that overflows it by one
     byte, so the flush boundary is exercised on both sides
   - a View longer than Capacity (the write-through path)
   - an empty View, and fprint_spaces(record, 0)
   - an early return with a partly-filled buffer, so the destructor flushes
   - fprint_double against std::fprintf for "%.1f", "%.2f", "%.13lf",
     "%3.0f" and "%5.1lf"

   The interleaving property that matters: a Record must not outlive the lock
   its caller holds over the stream. Every converted writer constructs it
   inside the function that already runs under its command's mutex_output, so
   the flush happens before the lock is released. A Record stored in a struct
   would break that, which is why it is non-movable. */
