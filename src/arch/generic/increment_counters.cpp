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

#include "arch/increment_counters.hpp"
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // int16_t, uint16_t
#include <cstring>  // std::memcpy
#include <iterator>  // std::next


/*
  Generic backend: plain C++11, no intrinsics and no SIMDE. Used on every
  target that has no hand-written backend of its own (RISC-V, MIPS, s390x, ...);
  aarch64, ppc64le and the two x86 ISA levels have their own.

  This replaced a SIMDE backend that compiled the x86 SSE intrinsics through
  SIMDE's translation layer. On a target with no vector unit for SIMDE to map
  onto, that emulates 128-bit lane semantics in scalar integer operations: it
  measured 305 instructions against the 53 below, on mips64el at -O2, with
  neither form emitting a single vector instruction. The loop is deliberately
  left rolled rather than unrolled sixteen ways, because on the ISAs that do
  have a per-lane variable shift (NEON, AltiVec) GCC vectorizes the rolled form
  of its own accord and unrolling defeats that.

  Kept out of arch/x86_64/: below AVX2 the x86 ISAs have no per-lane variable
  shift, so this loop cannot vectorize there and PSHUFB really is the algorithm
  -- see arch/x86_64/SSSE3/.
*/
auto increment_counters_from_bitmap(Span<count_t> const counters,
                                    View<unsigned char> const bitmap) -> void
{
  /* the ceiling the SIMD backends saturate at: _mm_subs_epi16, vqsubq_s16 and
     vec_subs are signed saturating, so they cap at INT16_MAX rather than
     wrapping the unsigned-short counter at 65536 */
  static constexpr auto counter_max = int16_t{32767};

  auto const rounds = (counters.size() + counters_per_round - 1) / counters_per_round;

  for (auto round = std::size_t{0}; round < rounds; ++round)
    {
      /* std::next takes a signed difference_type, and these products are
         std::size_t; the casts keep -Wsign-conversion quiet without widening
         anything, since a round index cannot reach PTRDIFF_MAX */
      auto const bitmap_offset =
        static_cast<std::ptrdiff_t>(round * bytes_per_round);
      auto const counter_offset =
        static_cast<std::ptrdiff_t>(round * counters_per_round);

      uint16_t bits = 0;
      std::memcpy(&bits, std::next(bitmap.data(), bitmap_offset), bytes_per_round);
      auto * const slice = std::next(counters.data(), counter_offset);

      /* Written as a fold rather than a branch: the increment is zeroed at the
         ceiling instead of jumping over the add, which is the shape a
         vectorizer can take. Saturation costs nothing expressed this way. */
      for (auto lane = 0U; lane < counters_per_round; ++lane)
        {
          auto const value = static_cast<int16_t>(slice[lane]);
          auto const selected = static_cast<int16_t>((bits >> lane) & 1U);
          slice[lane] = static_cast<count_t>(
            value + ((value < counter_max) ? selected : int16_t{0}));
        }
    }
}
