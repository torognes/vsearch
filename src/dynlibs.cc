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

#include "vsearch.h"
#include "dynlibs.h"
#include "utils/dynlib_loader.hpp"
#include "utils/fatal.hpp"
#include <cstdio>  // std::FILE
#include <string>


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


DynamicLibraries::DynamicLibraries() noexcept
{
#ifdef HAVE_ZLIB_H
  gz_lib = dynlib::open(gz_libname.data());
  if (gz_lib != nullptr)
    {
      gzdopen_p = reinterpret_cast<gzFile (*)(int, const char*)>(
        dynlib::symbol(gz_lib, "gzdopen"));
      gzclose_p = reinterpret_cast<int (*)(gzFile)>(
        dynlib::symbol(gz_lib, "gzclose"));
      gzread_p = reinterpret_cast<int (*)(gzFile, void*, unsigned)>(
        dynlib::symbol(gz_lib, "gzread"));
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
      BZ2_bzReadOpen_p = reinterpret_cast<BZFILE* (*)(int*, FILE*, int, int, void*, int)>(
        dynlib::symbol(bz2_lib, "BZ2_bzReadOpen"));
      BZ2_bzReadClose_p = reinterpret_cast<void (*)(int*, BZFILE*)>(
        dynlib::symbol(bz2_lib, "BZ2_bzReadClose"));
      BZ2_bzRead_p = reinterpret_cast<int (*)(int*, BZFILE*, void*, int)>(
        dynlib::symbol(bz2_lib, "BZ2_bzRead"));
      if (not ((BZ2_bzReadOpen_p != nullptr) && (BZ2_bzReadClose_p != nullptr) && (BZ2_bzRead_p != nullptr)))
        {
          fatal("Invalid compression library (bz2)");
        }
    }
#endif
}


DynamicLibraries::~DynamicLibraries()
{
#ifdef HAVE_ZLIB_H
  if (gz_lib != nullptr)
    {
      dynlib::close(gz_lib);
    }
#endif

#ifdef HAVE_BZLIB_H
  if (bz2_lib != nullptr)
    {
      dynlib::close(bz2_lib);
    }
#endif
}


#ifdef HAVE_ZLIB_H
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
#endif
