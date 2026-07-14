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
  globals and the manual dynlibs_open()/dynlibs_close() pair. The thin
  accessors forward to the resolved function pointers; keeping them in the
  header lets the compiler inline the calls in the read hot loop. The
  gzFile/BZFILE types come from zlib.h/bzlib.h, which this header includes
  directly (under the same HAVE_* guards).
*/

#ifdef HAVE_CONFIG_H
#include "config.h"  // HAVE_ZLIB_H, HAVE_BZLIB_H
#endif
#ifdef HAVE_ZLIB_H
#include <zlib.h>  // gzFile
#endif
#ifdef HAVE_BZLIB_H
#include <bzlib.h>  // BZFILE
#endif

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

#ifdef HAVE_ZLIB_H
  auto gzip_available() const noexcept -> bool { return gz_lib != nullptr; }
  auto gzdopen(int file_descriptor, char const * mode) const noexcept -> gzFile
  { return gzdopen_p(file_descriptor, mode); }
  auto gzclose(gzFile file) const noexcept -> int { return gzclose_p(file); }
  auto gzread(gzFile file, void * buffer, unsigned length) const noexcept -> int
  { return gzread_p(file, buffer, length); }
  auto gzip_version() const noexcept -> char const *;
  auto gzip_compile_flags() const noexcept -> unsigned long;
#endif

#ifdef HAVE_BZLIB_H
  auto bzip2_available() const noexcept -> bool { return bz2_lib != nullptr; }
  auto bz_read_open(int * bzerror, FILE * file, int verbosity, int use_small,
                    void * unused, int unused_length) const noexcept -> BZFILE *
  { return BZ2_bzReadOpen_p(bzerror, file, verbosity, use_small, unused, unused_length); }
  auto bz_read_close(int * bzerror, BZFILE * file) const noexcept -> void
  { BZ2_bzReadClose_p(bzerror, file); }
  auto bz_read(int * bzerror, BZFILE * file, void * buffer, int length) const noexcept -> int
  { return BZ2_bzRead_p(bzerror, file, buffer, length); }
#endif

private:
#ifdef HAVE_ZLIB_H
  void * gz_lib = nullptr;
  gzFile ZEXPORT (*gzdopen_p) OF((int, char const *)) = nullptr;
  int ZEXPORT (*gzclose_p) OF((gzFile)) = nullptr;
  int ZEXPORT (*gzread_p) OF((gzFile, void *, unsigned)) = nullptr;
#endif
#ifdef HAVE_BZLIB_H
  void * bz2_lib = nullptr;
  BZFILE * (*BZ2_bzReadOpen_p)(int *, FILE *, int, int, void *, int) = nullptr;
  void (*BZ2_bzReadClose_p)(int *, BZFILE *) = nullptr;
  int (*BZ2_bzRead_p)(int *, BZFILE *, void *, int) = nullptr;
#endif
};
