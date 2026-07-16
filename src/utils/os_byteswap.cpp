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

// Operating System specific commands to swap bytes
// C++23 refactoring: replace with std::byteswap()


#if defined(_MSC_VER) || defined(_WIN32)

#include <cstdint>  // uint16_t, uint32_t, uint64_t
#include <stdlib.h>

auto bswap_16(uint16_t bsx) noexcept -> uint16_t {
  return _byteswap_ushort(bsx);
}

auto bswap_32(uint32_t bsx) noexcept -> uint32_t {
  return _byteswap_ulong(bsx);
}

auto bswap_64(uint64_t bsx) noexcept -> uint64_t {
  return _byteswap_uint64(bsx);
}


#elif defined(__APPLE__)

// Mac OS X / Darwin features
#include <cstdint>  // uint16_t, uint32_t, uint64_t
#include <libkern/OSByteOrder.h>

auto bswap_16(uint16_t bsx) noexcept -> uint16_t {
  return OSSwapInt16(bsx);
}

auto bswap_32(uint32_t bsx) noexcept -> uint32_t {
  return OSSwapInt32(bsx);
}

auto bswap_64(uint64_t bsx) noexcept -> uint64_t {
  return OSSwapInt64(bsx);
}


#elif defined(__FreeBSD__)

#include <cstdint>  // uint16_t, uint32_t, uint64_t
#include <sys/endian.h>

auto bswap_16(uint16_t bsx) noexcept -> uint16_t {
  return bswap16(bsx);
}

auto bswap_32(uint32_t bsx) noexcept -> uint32_t {
  return bswap32(bsx);
}

auto bswap_64(uint64_t bsx) noexcept -> uint64_t {
  return bswap64(bsx);
}


#elif defined(__NetBSD__)

#include <cstdint>  // uint16_t, uint32_t, uint64_t
#include <sys/types.h>
#include <machine/bswap.h>

auto bswap_16(uint16_t bsx) noexcept -> uint16_t {
  return bswap16(bsx);
}

auto bswap_32(uint32_t bsx) noexcept -> uint32_t {
  return bswap32(bsx);
}

auto bswap_64(uint64_t bsx) noexcept -> uint64_t {
  return bswap64(bsx);
}


#else

// Linux and other operating systems. Use the compiler byteswap builtins
// (GCC 4.8+ / Clang) when available — identical to what glibc's <byteswap.h>
// macros expand to — and fall back to a portable shift-based implementation
// for any other compiler, so the build no longer depends on <byteswap.h>
// being present.
#include <cstdint>  // uint16_t, uint32_t, uint64_t

auto bswap_16(uint16_t bsx) noexcept -> uint16_t {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap16(bsx);
#else
  return static_cast<uint16_t>((bsx >> 8U) | (bsx << 8U));
#endif
}

auto bswap_32(uint32_t bsx) noexcept -> uint32_t {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap32(bsx);
#else
  return ((bsx & UINT32_C(0x000000FF)) << 24U) |
         ((bsx & UINT32_C(0x0000FF00)) << 8U)  |
         ((bsx & UINT32_C(0x00FF0000)) >> 8U)  |
         ((bsx & UINT32_C(0xFF000000)) >> 24U);
#endif
}

auto bswap_64(uint64_t bsx) noexcept -> uint64_t {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap64(bsx);
#else
  return ((bsx & UINT64_C(0x00000000000000FF)) << 56U) |
         ((bsx & UINT64_C(0x000000000000FF00)) << 40U) |
         ((bsx & UINT64_C(0x0000000000FF0000)) << 24U) |
         ((bsx & UINT64_C(0x00000000FF000000)) << 8U)  |
         ((bsx & UINT64_C(0x000000FF00000000)) >> 8U)  |
         ((bsx & UINT64_C(0x0000FF0000000000)) >> 24U) |
         ((bsx & UINT64_C(0x00FF000000000000)) >> 40U) |
         ((bsx & UINT64_C(0xFF00000000000000)) >> 56U);
#endif
}

#endif
