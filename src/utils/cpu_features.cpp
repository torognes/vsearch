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

#include "cpu_features.hpp"
#include <cstdint>  // int64_t
#include <cstdio>  // fprintf, stderr
#ifdef __x86_64__
#include <cpuid.h>  // __cpuid_count
#endif


/* cpu features available */

int64_t altivec_present = 0;
int64_t neon_present = 0;
int64_t mmx_present = 0;
int64_t sse_present = 0;
int64_t sse2_present = 0;
int64_t sse3_present = 0;
int64_t ssse3_present = 0;
int64_t sse41_present = 0;
int64_t sse42_present = 0;
int64_t popcnt_present = 0;
int64_t avx_present = 0;
int64_t avx2_present = 0;


#ifdef __x86_64__
namespace {
  struct cpuid_registers {
    unsigned int eax {0};
    unsigned int ebx {0};
    unsigned int ecx {0};
    unsigned int edx {0};
  };

  auto get_cpuid(unsigned int const leaf, unsigned int const subleaf) noexcept -> cpuid_registers {
    cpuid_registers registers {};
    __cpuid_count(leaf, subleaf, registers.eax, registers.ebx, registers.ecx, registers.edx);
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
#endif


auto cpu_features_detect() -> void
{
#ifdef __aarch64__
#ifdef __ARM_NEON
  /* may check /proc/cpuinfo for asimd or neon */
  neon_present = 1;
#else
#error ARM Neon not present
#endif
#elif __PPC__
  altivec_present = 1;
#elif __x86_64__
  cpuid_registers const leaf0 = get_cpuid(0, 0);
  unsigned int const maxlevel = leaf0.eax & 0xff;

  if (maxlevel >= 1)
    {
      cpuid_registers const leaf1 = get_cpuid(1, 0);
      mmx_present    = (leaf1.edx >> 23U) & 1U;
      sse_present    = (leaf1.edx >> 25U) & 1U;
      sse2_present   = (leaf1.edx >> 26U) & 1U;
      sse3_present   = (leaf1.ecx >>  0U) & 1U;
      ssse3_present  = (leaf1.ecx >>  9U) & 1U;
      sse41_present  = (leaf1.ecx >> 19U) & 1U;
      sse42_present  = (leaf1.ecx >> 20U) & 1U;
      popcnt_present = (leaf1.ecx >> 23U) & 1U;

      // AVX/AVX2 are only usable if the OS has enabled saving of the YMM
      // register state: CPUID must report OSXSAVE (leaf-1 ECX bit 27) and
      // XCR0 (read via XGETBV) must have both the SSE (bit 1) and AVX
      // (bit 2) state-enable bits set. Without this check an AVX-capable
      // CPU running on an old OS would be over-reported.
      static constexpr unsigned int xcr0_avx_state = 0x6U;  // XMM | YMM
      bool const osxsave_present = ((leaf1.ecx >> 27U) & 1U) != 0U;
      bool const avx_os_enabled =
        osxsave_present and ((read_xcr0() & xcr0_avx_state) == xcr0_avx_state);
      bool const avx_supported = ((leaf1.ecx >> 28U) & 1U) != 0U;
      avx_present = (avx_supported and avx_os_enabled) ? 1 : 0;

      if (maxlevel >= 7)
        {
          cpuid_registers const leaf7 = get_cpuid(7, 0);
          bool const avx2_supported = ((leaf7.ebx >> 5U) & 1U) != 0U;
          avx2_present = (avx2_supported and avx_os_enabled) ? 1 : 0;
        }
    }
#else
    // simde
#endif
}


auto cpu_features_show() -> void
{
  fprintf(stderr, "CPU features:");
  if (neon_present != 0)
    {
      fprintf(stderr, " neon");
    }
  if (altivec_present != 0)
    {
      fprintf(stderr, " altivec");
    }
  if (mmx_present != 0)
    {
      fprintf(stderr, " mmx");
    }
  if (sse_present != 0)
    {
      fprintf(stderr, " sse");
    }
  if (sse2_present != 0)
    {
      fprintf(stderr, " sse2");
    }
  if (sse3_present != 0)
    {
      fprintf(stderr, " sse3");
    }
  if (ssse3_present != 0)
    {
      fprintf(stderr, " ssse3");
    }
  if (sse41_present != 0)
    {
      fprintf(stderr, " sse4.1");
    }
  if (sse42_present != 0)
    {
      fprintf(stderr, " sse4.2");
    }
  if (popcnt_present != 0)
    {
      fprintf(stderr, " popcnt");
    }
  if (avx_present != 0)
    {
      fprintf(stderr, " avx");
    }
  if (avx2_present != 0)
    {
      fprintf(stderr, " avx2");
    }
  fprintf(stderr, "\n");
}
