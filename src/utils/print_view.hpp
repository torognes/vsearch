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


#include "utils/view.hpp"  // View<char>
#include <cstddef>  // std::size_t
#include <cstdio>  // std::FILE, std::fwrite


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
