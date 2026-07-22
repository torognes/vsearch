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

#include "vsearch.hpp"
#include "os/dynlibs.hpp"
#include "utils/dynlib_loader.hpp"
#include "utils/fatal.hpp"
#ifdef HAVE_BZLIB_H
#include <bzlib.h>  // BZFILE
#endif
#include <cstdio>  // std::FILE
#include <string>
#ifdef HAVE_ZLIB_H
#include <zlib.h>  // gzFile
#endif


#ifdef HAVE_ZLIB_H
# ifdef _WIN32
const std::string gz_libname = "zlib1.dll";
# elif defined(__APPLE__)
const std::string gz_libname = "libz.dylib";
# elif defined(__FreeBSD__)
const std::string gz_libname = "libz.so.6";
# else
const std::string gz_libname = "libz.so.1";
# endif
#endif

#ifdef HAVE_BZLIB_H
# ifdef _WIN32
const std::string bz2_libname = "libbz2.dll";
# elif defined(__APPLE__)
const std::string bz2_libname = "libbz2.dylib";
# elif defined(__FreeBSD__)
const std::string bz2_libname = "libbz2.so.4";
# else
const std::string bz2_libname = "libbz2.so.1";
# endif
#endif

#ifdef HAVE_BZLIB_H
/* fixed arguments passed to BZ2_bzReadOpen (see bz_open below) */
constexpr auto BZ_VERBOSE_0 = 0;
// constexpr auto BZ_VERBOSE_1 = 1;
// constexpr auto BZ_VERBOSE_2 = 2;
// constexpr auto BZ_VERBOSE_3 = 3;
// constexpr auto BZ_VERBOSE_4 = 4;
constexpr auto BZ_MORE_MEM = 0;  /* faster decompression using more memory */
// constexpr auto BZ_LESS_MEM = 1;  /* slower decompression but requires less memory */
#endif


DynamicLibraries::DynamicLibraries() noexcept
{
#ifdef HAVE_ZLIB_H
  gz_lib = dynlib::open(gz_libname.data());
  if (gz_lib != nullptr)
    {
      // dynlib::symbol returns a generic void(*)(); each accessor casts it
      // back to the real zlib signature at the point of call.
      gzdopen_p = dynlib::symbol(gz_lib, "gzdopen");
      gzclose_p = dynlib::symbol(gz_lib, "gzclose");
      gzread_p = dynlib::symbol(gz_lib, "gzread");
      if (not ((gzdopen_p != nullptr) && (gzclose_p != nullptr) && (gzread_p != nullptr)))
        {
          fatal("Invalid compression library (zlib)");
        }
    }
#endif

#ifdef HAVE_BZLIB_H
  bz2_lib = dynlib::open(bz2_libname.data());
  if (bz2_lib != nullptr)
    {
      bz_read_open_p = dynlib::symbol(bz2_lib, "BZ2_bzReadOpen");
      bz_read_close_p = dynlib::symbol(bz2_lib, "BZ2_bzReadClose");
      bz_read_p = dynlib::symbol(bz2_lib, "BZ2_bzRead");
      if (not ((bz_read_open_p != nullptr) && (bz_read_close_p != nullptr) && (bz_read_p != nullptr)))
        {
          fatal("Invalid compression library (bz2)");
        }
    }
#endif
}


DynamicLibraries::~DynamicLibraries()
{
  // handles are plain void*, so the teardown needs no HAVE_* guards: an
  // unloaded library simply left its handle nullptr.
  if (gz_lib != nullptr)
    {
      dynlib::close(gz_lib);
    }
  if (bz2_lib != nullptr)
    {
      dynlib::close(bz2_lib);
    }
}


// gzip_version / gzip_compile_flags name no zlib type, so they need no guard;
// they are only ever called when gzip_available() is true (gz_lib loaded).
auto DynamicLibraries::gzip_version() const noexcept -> char const *
{
  auto * const version_fn = reinterpret_cast<char const * (*)()>(
    dynlib::symbol(gz_lib, "zlibVersion"));
  return version_fn();
}


auto DynamicLibraries::gzip_compile_flags() const noexcept -> unsigned long
{
  auto * const flags_fn = reinterpret_cast<unsigned long (*)()>(
    dynlib::symbol(gz_lib, "zlibCompileFlags"));
  return flags_fn();
}


/* The accessors below are the only place the real gzFile/BZFILE types and the
   BZ_* status codes appear. Each casts a type-erased symbol back to its real
   signature and casts the void* handle back to the library type. When the
   library is not compiled in, stub definitions keep the (unconditional)
   declarations link-complete; they are never reached because gzip_available()
   / bzip2_available() stay false. */

#ifdef HAVE_ZLIB_H

auto DynamicLibraries::gz_open(int const file_descriptor) const noexcept -> void *
{
  auto * const open_fn = reinterpret_cast<gzFile (*)(int, char const *)>(gzdopen_p);
  return open_fn(file_descriptor, "rb");
}


auto DynamicLibraries::gz_close(void * const stream) const noexcept -> int
{
  auto * const close_fn = reinterpret_cast<int (*)(gzFile)>(gzclose_p);
  return close_fn(static_cast<gzFile>(stream));
}


auto DynamicLibraries::gz_read(void * const stream, void * const buffer, unsigned const length) const noexcept -> int
{
  auto * const read_fn = reinterpret_cast<int (*)(gzFile, void *, unsigned)>(gzread_p);
  return read_fn(static_cast<gzFile>(stream), buffer, length);
}

#else

auto DynamicLibraries::gz_open(int) const noexcept -> void * { return nullptr; }
auto DynamicLibraries::gz_close(void *) const noexcept -> int { return -1; }
auto DynamicLibraries::gz_read(void *, void *, unsigned) const noexcept -> int { return -1; }

#endif


#ifdef HAVE_BZLIB_H

auto DynamicLibraries::bz_open(std::FILE * const file) const noexcept -> void *
{
  auto * const open_fn = reinterpret_cast<BZFILE * (*)(int *, FILE *, int, int, void *, int)>(bz_read_open_p);
  int bz_error = 0;
  return open_fn(& bz_error, file, BZ_VERBOSE_0, BZ_MORE_MEM, nullptr, 0);
}


auto DynamicLibraries::bz_close(void * const stream) const noexcept -> void
{
  auto * const close_fn = reinterpret_cast<void (*)(int *, BZFILE *)>(bz_read_close_p);
  int bz_error = 0;
  close_fn(& bz_error, static_cast<BZFILE *>(stream));
}


auto DynamicLibraries::bz_read(void * const stream, void * const buffer, int const length) const noexcept -> int
{
  auto * const read_fn = reinterpret_cast<int (*)(int *, BZFILE *, void *, int)>(bz_read_p);
  int bz_error = 0;
  int const bytes_read = read_fn(& bz_error, static_cast<BZFILE *>(stream), buffer, length);
  if ((bytes_read < 0) or
      not ((bz_error == BZ_OK) or
           (bz_error == BZ_STREAM_END) or
           (bz_error == BZ_SEQUENCE_ERROR)))
    {
      return -1;
    }
  return bytes_read;
}

#else

auto DynamicLibraries::bz_open(std::FILE *) const noexcept -> void * { return nullptr; }
auto DynamicLibraries::bz_close(void *) const noexcept -> void { }
auto DynamicLibraries::bz_read(void *, void *, int) const noexcept -> int { return -1; }

#endif
