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

/*
  DynamicLibraries loads the optional compression libraries (zlib, bzip2)
  at runtime and resolves the handful of symbols vsearch needs. The
  constructor opens the libraries and the destructor closes them (RAII),
  so a single instance owned by main() replaces the former file-scope
  globals and the manual dynlibs_open()/dynlibs_close() pair. The
  accessors forward to the resolved function pointers. They take and
  return type-erased void* handles and are defined out of line in
  dynlibs.cpp, so this header pulls in neither zlib.h nor bzlib.h: the
  real gzFile/BZFILE types and the HAVE_*_H guards are confined to that
  single translation unit. Callers test compression::gzip_supported /
  compression::bzip2_supported (below) at run time instead of guarding
  with #ifdef.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"  // HAVE_ZLIB_H, HAVE_BZLIB_H
#endif
#include <cstdio>  // std::FILE

namespace compression
{
  /* Whether support for each optional compression library was compiled in.
     These constants are the single source of truth for the HAVE_*_H feature
     tests: callers branch on them at run time (the compiler elides the dead
     branch) instead of scattering #ifdef guards throughout the code base. */
#ifdef HAVE_ZLIB_H
  constexpr bool gzip_supported = true;
#else
  constexpr bool gzip_supported = false;
#endif
#ifdef HAVE_BZLIB_H
  constexpr bool bzip2_supported = true;
#else
  constexpr bool bzip2_supported = false;
#endif
}

class DynamicLibraries
{
public:
  DynamicLibraries() noexcept;
  ~DynamicLibraries();

  // owns OS library handles: neither copyable nor movable
  DynamicLibraries(DynamicLibraries const &) = delete;
  auto operator=(DynamicLibraries const &) -> DynamicLibraries & = delete;
  DynamicLibraries(DynamicLibraries &&) = delete;
  auto operator=(DynamicLibraries &&) -> DynamicLibraries & = delete;

  // gzip: handles are the opaque zlib gzFile, exposed here as void*
  auto gzip_available() const noexcept -> bool { return gz_lib != nullptr; }
  auto gz_open(int file_descriptor) const noexcept -> void *;
  auto gz_close(void * stream) const noexcept -> int;
  auto gz_read(void * stream, void * buffer, unsigned length) const noexcept -> int;
  auto gzip_version() const noexcept -> char const *;
  auto gzip_compile_flags() const noexcept -> unsigned long;

  // bzip2: handles are the opaque BZFILE, exposed here as void*. bz_read
  // returns the byte count, or a negative value on a genuine read error
  // (the BZ_* status codes are interpreted in dynlibs.cpp, not leaked out).
  auto bzip2_available() const noexcept -> bool { return bz2_lib != nullptr; }
  auto bz_open(std::FILE * file) const noexcept -> void *;
  auto bz_close(void * stream) const noexcept -> void;
  auto bz_read(void * stream, void * buffer, int length) const noexcept -> int;

private:
  void * gz_lib = nullptr;
  void * bz2_lib = nullptr;

  /* zlib/bzip2 entry points, resolved at construction and stored type-erased
     as a generic function-pointer type; dynlibs.cpp casts each back to its
     real signature at the call site, exactly as dynlib::symbol() is cast. */
  void (*gzdopen_p)() = nullptr;
  void (*gzclose_p)() = nullptr;
  void (*gzread_p)() = nullptr;
  void (*bz_read_open_p)() = nullptr;
  void (*bz_read_close_p)() = nullptr;
  void (*bz_read_p)() = nullptr;
};
