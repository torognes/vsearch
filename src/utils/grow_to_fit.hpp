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

#include <cstddef>  // std::size_t


/* Grow a reused scratch buffer so it can hold at least min_size elements,
   never shrinking it.

   The buffers this serves are high-water marks kept across records: a short
   record after a long one reuses the allocation and does no work at all. That
   is the whole idiom, hand-written at more than a dozen sites, and the point
   of naming it is that the name says "never shrinks" where an unguarded
   resize() would not.

   min_size is the caller's, deliberately. Six of the call sites ask for
   "length + 1" because they hand the buffer to a helper that terminates it,
   and folding that "+ 1" in here would hide the reason for the extra byte at
   the exact sites where a reader needs to see it -- the audit in
   TBD_20260901_grow_to_fit.md found one buffer whose terminator was being read
   by accident. Passing min_size in also lets the buffers that are not sized
   from a record length (the kmer tables in core/unique.cpp) use the same
   helper.

   Deliberately unconstrained: a container without size()/resize() already
   produces a better diagnostic from the body ("no member named 'resize'; did
   you mean 'size'?", naming the type and the line) than a static_assert would.
   A type trait becomes worth its ~20 lines of C++11 void_t machinery only if
   an overload appears that needs SFINAE to disambiguate -- a fill-value form,
   say -- at which point the constraint has to move to enable_if anyway,
   because a static_assert fires instead of removing a candidate.

   The std::max spelling was measured against this guarded one: byte-identical
   codegen at -O2. This form is kept because it reads as what it is. */
namespace vsearch
{

  template <typename Container>
  auto grow_to_fit(Container & container, std::size_t const min_size) -> void {
    if (container.size() < min_size) {
      container.resize(min_size);
    }
  }

}  // namespace vsearch
