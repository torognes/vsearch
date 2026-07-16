/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

// x86-only translation unit: runtime SIMD detection via CPUID is an x86
// concern (on ARM/POWER the ISA level is fixed at compile time). The build
// compiles this file only under TARGET_X86_64, so no __x86_64__ guard is
// needed here -- the Makefile owns the architecture selection.

#include "arch/cpu_features.hpp"
#include "vsearch.hpp"  // struct Parameters
#include "utils/fatal.hpp"  // fatal
#include <cstdint>  // int64_t
#include <cpuid.h>  // __cpuid_count, bit_* feature masks


namespace {
  struct cpuid_registers {
    unsigned int eax {0};
    unsigned int ebx {0};
    unsigned int ecx {0};
    unsigned int edx {0};
  };

  // All call sites query sub-leaf 0, so the sub-leaf is fixed here rather
  // than passed in (avoids two adjacent same-type parameters).
  auto get_cpuid(unsigned int const leaf) noexcept -> cpuid_registers {
    cpuid_registers registers {};
    __cpuid_count(leaf, 0U, registers.eax, registers.ebx, registers.ecx, registers.edx);
    return registers;
  }

  // Read the low 32 bits of XCR0 via XGETBV. Must only be called when CPUID
  // reports OSXSAVE, otherwise XGETBV raises #UD (illegal instruction).
  // cpu_features.cpp is compiled with the baseline target, so the _xgetbv
  // intrinsic is unavailable; use the equivalent one-instruction asm.
  auto read_xcr0() noexcept -> unsigned int {
    unsigned int xcr0_lo {0};
    unsigned int xcr0_hi {0};
    __asm__ __volatile__("xgetbv" : "=a"(xcr0_lo), "=d"(xcr0_hi) : "c"(0U));
    static_cast<void>(xcr0_hi);
    return xcr0_lo;
  }
}  // anonymous namespace


auto cpu_features_detect(struct Parameters & parameters) -> void
{
  // Feature masks (bit_MMX, bit_SSE, ...) come from <cpuid.h>. bit_OSXSAVE
  // is not defined by older <cpuid.h> versions (GCC 4.x), so spell it out.
  static constexpr unsigned int basic_leaf_mask = 0xffU;  // CPUID.0:EAX low byte
  static constexpr unsigned int extended_features_leaf = 7U;
  static constexpr unsigned int bit_osxsave = 0x08000000U;  // CPUID.1:ECX bit 27
  static constexpr unsigned int xcr0_avx_state = 0x6U;  // XMM | YMM

  cpuid_registers const leaf0 = get_cpuid(0U);
  unsigned int const maxlevel = leaf0.eax & basic_leaf_mask;

  if (maxlevel >= 1U)
    {
      cpuid_registers const leaf1 = get_cpuid(1U);
      parameters.mmx_present    = static_cast<int64_t>((leaf1.edx & bit_MMX)    != 0U);
      parameters.sse_present    = static_cast<int64_t>((leaf1.edx & bit_SSE)    != 0U);
      parameters.sse2_present   = static_cast<int64_t>((leaf1.edx & bit_SSE2)   != 0U);
      parameters.sse3_present   = static_cast<int64_t>((leaf1.ecx & bit_SSE3)   != 0U);
      parameters.ssse3_present  = static_cast<int64_t>((leaf1.ecx & bit_SSSE3)  != 0U);
      parameters.sse41_present  = static_cast<int64_t>((leaf1.ecx & bit_SSE4_1) != 0U);
      parameters.sse42_present  = static_cast<int64_t>((leaf1.ecx & bit_SSE4_2) != 0U);
      parameters.popcnt_present = static_cast<int64_t>((leaf1.ecx & bit_POPCNT) != 0U);

      // AVX/AVX2 are only usable if the OS has enabled saving of the YMM
      // register state: CPUID must report OSXSAVE (leaf-1 ECX bit 27) and
      // XCR0 (read via XGETBV) must have both the SSE (bit 1) and AVX
      // (bit 2) state-enable bits set. Without this check an AVX-capable
      // CPU running on an old OS would be over-reported.
      bool const osxsave_present = (leaf1.ecx & bit_osxsave) != 0U;
      bool const avx_os_enabled =
        osxsave_present and ((read_xcr0() & xcr0_avx_state) == xcr0_avx_state);
      bool const avx_supported = (leaf1.ecx & bit_AVX) != 0U;
      parameters.avx_present = static_cast<int64_t>(avx_supported and avx_os_enabled);

      if (maxlevel >= extended_features_leaf)
        {
          cpuid_registers const leaf7 = get_cpuid(extended_features_leaf);
          bool const avx2_supported = (leaf7.ebx & bit_AVX2) != 0U;
          parameters.avx2_present = static_cast<int64_t>(avx2_supported and avx_os_enabled);
        }
    }
}


auto cpu_features_test(struct Parameters const & parameters) -> void
{
  if (parameters.sse2_present == 0)
    {
      fatal("Sorry, this program requires a cpu with SSE2.");
    }
}
