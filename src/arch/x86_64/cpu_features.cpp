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
}  // anonymous namespace


auto cpu_features_detect(struct Parameters & parameters) -> void
{
  // Feature masks (bit_SSE2, bit_SSSE3) come from <cpuid.h>. Only the two
  // features the code actually branches on are probed; the ten flags this
  // function also used to set were never read anywhere in the tree.
  //
  // Should AVX or AVX2 ever be needed: reporting bit_AVX / bit_AVX2 is not
  // sufficient on its own. The OS must also have enabled saving of the YMM
  // register state, i.e. CPUID must report OSXSAVE (leaf-1 ECX bit 27, not
  // defined by older <cpuid.h> versions such as GCC 4.x, so spell it out)
  // and XCR0 -- read via XGETBV, which raises #UD unless OSXSAVE is set --
  // must have both the SSE (bit 1) and AVX (bit 2) state-enable bits set.
  // Without that check an AVX-capable CPU running on an old OS is
  // over-reported. AVX2 additionally lives in leaf 7 (EBX), which is only
  // available when CPUID reports a maximum basic leaf of at least 7.
  static constexpr unsigned int basic_leaf_mask = 0xffU;  // CPUID.0:EAX low byte

  cpuid_registers const leaf0 = get_cpuid(0U);
  unsigned int const maxlevel = leaf0.eax & basic_leaf_mask;

  if (maxlevel >= 1U)
    {
      cpuid_registers const leaf1 = get_cpuid(1U);
      parameters.sse2_present  = static_cast<int64_t>((leaf1.edx & bit_SSE2)  != 0U);
      parameters.ssse3_present = static_cast<int64_t>((leaf1.ecx & bit_SSSE3) != 0U);
    }
}


auto cpu_features_test(struct Parameters const & parameters) -> void
{
  if (parameters.sse2_present == 0)
    {
      fatal("Sorry, this program requires a cpu with SSE2.");
    }
}
