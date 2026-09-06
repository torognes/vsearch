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

#include <cassert>  // assert


/* The one place a nucleotide becomes a table index.

   It lives one level above the per-map namespaces, so each map reaches it
   by unqualified lookup and no map owns a copy -- centralizing this cast
   is the whole reason the character maps are accessors rather than bare
   arrays.

   Why the maps are shaped the way they are. vsearch is built without LTO,
   so a mapper defined in a .cpp costs a real call per base at every call
   site; before this refactoring that was 106 emitted calls, map_4bit()
   alone accounting for 1.41 % of an --allpairs_global run. Each map is
   therefore a header, and its table a function-local `static constexpr`:

   - `static constexpr` is constant-initialized straight into .rodata, so
     it carries no thread-safe-static guard. A `static const` initialized
     by a call would cost a load and a branch on that guard at every call,
     inside the hot loop -- verified down to GCC 4.9;
   - the enclosing function is `inline`, so its local static has vague
     linkage: every translation unit emits a COMDAT copy and the linker
     keeps exactly one, which is what C++17 spells `inline constexpr`. */
namespace vsearch
{
  namespace maps
  {

    inline auto to_uchar(char const nucleotide) noexcept -> unsigned char {
      assert(nucleotide >= 0);
      return static_cast<unsigned char>(nucleotide);
    }

  }  // namespace maps
}  // namespace vsearch
