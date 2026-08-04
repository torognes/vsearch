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
#include "utils/view.hpp"  // View<char>
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdio>  // std::FILE, std::fputc, std::fwrite


/* Emit the bytes of a View to a stream.

   This is what a "%.*s" conversion was used for before: pass a counted run of
   characters that is not NUL-terminated (a slice of a sequence, a header inside
   a shared buffer). std::fwrite takes both its element size and its count as a
   std::size_t, and the count is what View::size() already returns, so nothing is
   narrowed to the int that "%.*s" requires -- and there is no format string to
   parse at run time.

   Not exactly interchangeable with "%.*s": that conversion stops at an embedded
   NUL as well as at the precision, while fwrite always emits the full byte
   count. vsearch's headers, sequences and quality strings contain no embedded
   NUL (the readers reject them), so the two agree on every input reaching these
   calls.

   It lives here rather than in view.hpp so that View stays a pure data type
   with no <cstdio> dependency. */

inline auto fprint(std::FILE * output_handle, View<char> const text) -> void
{
  /* An empty view may carry a null pointer, and passing one to fwrite is
     undefined even with a zero count. This is a real case, not a defensive
     guard: the --profile output (see msa.cpp) describes a cluster with no
     centroid sequence and hands over an empty sequence view. */
  if (text.empty()) { return; }
  /* fwrite's two size arguments are an element size and a count, in that order,
     and transposing them compiles. Naming the element size says which is which;
     sizeof(char) is 1 by definition, so this is about the reader, not about a
     platform where it could differ. Same shape as sff_convert.cpp's byte_size,
     kept local so a widely-included header exports no such name (see the note
     on span.hpp's private bounds constants). */
  static constexpr std::size_t element_size = sizeof(char);
  static_cast<void>(std::fwrite(text.data(), element_size, text.size(), output_handle));
}


/* Emit one character: what std::fputc was used for.

   The only thing this adds over the call it wraps is that the discarded return
   value is dealt with once, here, instead of at every call site -- vsearch had
   115 'static_cast<void>(std::fputc(...))' spellings, every one of them with a
   character literal.

   Beware of an int argument: 'fprint(handle, count)' compiles, converts to char
   and writes one byte rather than the digits of 'count'. The debug build's
   -Wconversion catches it for a variable, but not for a constant that fits;
   fprint_integer() below is what a number wants. */
inline auto fprint(std::FILE * output_handle, char const character) -> void
{
  static_cast<void>(std::fputc(character, output_handle));
}


/* Emit a string literal: what std::fputs was used for.

   std::fwrite with the array's own bound, not std::fputs, so the length is a
   compile-time constant by construction rather than by optimisation. GCC does
   fold fputs of a literal into exactly this fwrite, and does so even under
   _FORTIFY_SOURCE, but that is a property of one compiler and spelling it out
   costs nothing. (It is fprintf that fortification stops GCC from folding,
   which is why the fprintf calls elsewhere are worth replacing and these are
   not.)

   The parameter is a reference to an array so that Size arrives with it; a
   char const * would have to be walked at run time. The C array is therefore
   deliberate, as in fatal.hpp's explicit_decay.

   Intended for literals, and the contract is the array's bound, not a
   terminator. vsearch does declare partially-filled char arrays --
   cluster.hpp's centroid_label[1024] and cigar[4096], chimera.hpp's four
   label[1024] buffers -- and passing one here would emit all 1023 bytes, not
   the string inside it. Those buffers reach a stream through the View overload
   above, with the length their writer already knows. The assert catches an
   unterminated array; a filled-then-truncated one is the gap it cannot close,
   hence this note. */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
template <std::size_t Size>
auto fprint(std::FILE * output_handle, char const (&literal)[Size]) -> void
{
  static_assert(Size > 0, "a string literal always carries its terminator");
  assert(literal[Size - 1] == '\0');
  /* Size counts the terminating NUL, which is not part of the output. An empty
     literal leaves count 0, which fwrite accepts: unlike the View overload
     above, the pointer here cannot be null. */
  static constexpr std::size_t element_size = sizeof(char);
  static_cast<void>(std::fwrite(literal, element_size, Size - 1, output_handle));
}


/* Emit one integer, in decimal, to a stream: what a "%d" or a "%" PRIu64
   conversion was used for. The digits come from decimal_digits.hpp, so this is
   the same fwrite as above over a locally-produced view.

   Deliberately not a batching writer. A tab-separated run of numbers could be
   accumulated in a fixed buffer and emitted with a single fwrite, which
   measures faster in isolation -- but std::FILE * already buffers, so the only
   thing such a buffer saves is the stdio call count, i.e. the stream lock,
   against a second buffer in front of stdio's own and a lifetime to manage at
   every call site. (putc_unlocked would beat both, and is a POSIX-ism MinGW
   spells differently, so it is out.) */
template <typename Integer>
auto fprint_integer(std::FILE * output_handle, Integer const value) -> void
{
  decimal::Buffer buffer {};
  fprint(output_handle, decimal::to_decimal(buffer, value));
}


/* Emit a run of blanks: what a "%*s" conversion with an empty string was used
   for, and the padding half of the writer below.

   One fputc per blank rather than one fwrite over a static run of spaces,
   because every width in vsearch is small and known: the constant-width table
   fields go up to "%10d", and the two run-time widths (showalign.cpp,
   results.cpp) are the decimal width of a sequence length. A static blank
   array in a header this widely included would buy nothing measurable and
   export a name (see the note on span.hpp's private bounds constants). */
inline auto fprint_spaces(std::FILE * output_handle, std::size_t const count) -> void
{
  for (std::size_t written = 0; written < count; ++written) {
    fprint(output_handle, ' ');
  }
}


/* Emit one integer right-aligned in a field of 'width' characters: what a
   "%10d" or a "%*" PRId64 conversion was used for.

   Matches printf's behaviour on both edges: a value wider than the field is
   written in full rather than truncated, and the blanks go to the left of a
   minus sign. It deliberately does not reproduce "%-10d" (left-aligned; no
   site needs it) or "%02x" (zero-padded hex; two sites, which keep their
   fprintf).

   The two trailing arguments are both integers and are therefore swappable at
   a call site; nothing here can catch that, so a converted table line is worth
   reading against the format string it replaced. */
template <typename Integer>
auto fprint_integer(std::FILE * output_handle, Integer const value, std::size_t const width) -> void
{
  decimal::Buffer buffer {};
  auto const digits = decimal::to_decimal(buffer, value);
  if (digits.size() < width) {
    fprint_spaces(output_handle, width - digits.size());
  }
  fprint(output_handle, digits);
}


/* tests:

   Three overloads share the name 'fprint', so which one a call picks is the
   part worth pinning down. Checked by capturing the bytes written to a
   std::tmpfile():

   fprint(stream, '\t')                -> "\t"          the char overload
   fprint(stream, "\t*\n")             -> "\t*\n"       the literal overload
   fprint(stream, "")                  -> ""            count 0, pointer not null
   fprint(stream, View<char>{d, 3})    -> "abc"         the View overload
   fprint(stream, View<char>{})        -> ""            returns before fwrite
   fprint_integer(stream, 4294967295U) -> "4294967295"
   fprint_integer(stream, -5, 6)       -> "    -5"      blanks before the sign
   fprint_integer(stream, 1234567, 3)  -> "1234567"     never truncates
   fprint_spaces(stream, 0)            -> ""

   No call is ambiguous: a char and a char-array reference are each an exact
   match for their own overload, and View's two-argument constructor is
   explicit, so a literal cannot reach the View overload instead. The
   non-template View overload also outranks the array template for a View
   argument.

   There is deliberately no fprint(std::FILE *, char const *) overload. Against
   the array template, a string-literal argument gives both candidates an Exact
   Match sequence -- array-to-pointer is an lvalue transformation, excluded from
   the subsequence comparison -- and a non-template beats a template
   specialization on the tiebreaker, so adding one would silently route every
   literal through strlen instead of its compile-time bound. A run-time
   char const * therefore goes through View. */
