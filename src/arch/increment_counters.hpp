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

#include "utils/span.hpp"  // Span
#include "utils/view.hpp"  // View

/* increment_counters_from_bitmap: increment the 16-bit k-mer counters
   selected by the 1-bits of a bitmap (the search/sintax/cluster hot path).
   The SIMD implementation is architecture-specific and lives in one backend
   per ISA under arch/ (selected by the build, no preprocessor ladder):
     arch/x86_64/SSE2/    - SSE2 intrinsics   -> _sse2 variant
     arch/x86_64/SSSE3/   - SSSE3 intrinsics  -> _ssse3 variant
     arch/aarch64/        - NEON
     arch/ppc64le/        - AltiVec
     arch/simde/          - SSE2 intrinsics via SIMDE (portable fallback)
   On x86_64 both the _sse2 and _ssse3 variants are built and chosen at runtime
   by cpu_features_detect(); every other target builds a single plain-named
   variant. */

using count_t = unsigned short;

/* Every backend consumes the bitmap sixteen bits at a time: two bytes in, and
   the sixteen counters those bits select out (two 128-bit vectors of eight
   16-bit counters each). counters.size() counters therefore take
   (counters.size() + counters_per_round - 1) / counters_per_round rounds. */
constexpr auto counters_per_round = 16U;
constexpr auto bytes_per_round = 2U;  // the same sixteen bits, as whole bytes

/* counters.size() is the old 'totalbits' parameter: one counter per bitmap
   bit, so the span's own length is the bound and there is nothing left to pass
   beside it -- nor two same-shaped pointers left to transpose at a call site.

   Both ends are rounded up to a whole round, so the last one writes up to
   counters_per_round - 1 counters past counters.size() and reads up to fifteen
   bytes past the bitmap bits it needs. Both are deliberate and already paid
   for: the callers reserve 16 counters of headroom (overflow_padding in
   core/search.cpp, core/cluster.cpp and core/chimera.cpp) and the index pads
   every bitmap by 127 bits (simd_padding in Dbindex::set_bitmap_width). */

#ifdef __x86_64__
auto increment_counters_from_bitmap_sse2(Span<count_t> counters,
                                         View<unsigned char> bitmap) -> void;
auto increment_counters_from_bitmap_ssse3(Span<count_t> counters,
                                          View<unsigned char> bitmap) -> void;
#else
auto increment_counters_from_bitmap(Span<count_t> counters,
                                    View<unsigned char> bitmap) -> void;
#endif
