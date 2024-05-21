/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#ifdef HAVE_ZLIB_H
# ifdef _WIN32
const char gz_libname[] = "zlib1.dll";
HMODULE gz_lib;
# else
#  ifdef __APPLE__
const char gz_libname[] = "libz.dylib";
#  else
const char gz_libname[] = "libz.so.1";
#  endif
void * gz_lib;
# endif

gzFile ZEXPORT (*gzdopen_p) OF((int, const char *));
int ZEXPORT (*gzclose_p) OF((gzFile));
int ZEXPORT (*gzread_p) OF((gzFile, void *, unsigned));

#endif

#ifdef HAVE_BZLIB_H
# ifdef _WIN32
const char bz2_libname[] = "libbz2.dll";
HMODULE bz2_lib;
# else
#  ifdef __APPLE__
const char bz2_libname[] = "libbz2.dylib";
#  else
const char bz2_libname[] = "libbz2.so.1";
#  endif
void * bz2_lib;
# endif

BZFILE* (*BZ2_bzReadOpen_p)(int*, FILE*, int, int, void*, int);
void (*BZ2_bzReadClose_p)(int*, BZFILE*);
int (*BZ2_bzRead_p)(int*, BZFILE*, void*, int);

#endif

void dynlibs_open()
{
#ifdef HAVE_ZLIB_H
#ifdef _WIN32
  gz_lib = LoadLibraryA(gz_libname);
#else
  gz_lib = dlopen(gz_libname, RTLD_LAZY);
#endif
  if (gz_lib)
    {
      gzdopen_p = (gzFile (*)(int, const char*))
        arch_dlsym(gz_lib, "gzdopen");
      gzclose_p = (int (*)(gzFile))
        arch_dlsym(gz_lib, "gzclose");
      gzread_p = (int (*)(gzFile, void*, unsigned))
        arch_dlsym(gz_lib, "gzread");
      if (not (gzdopen_p && gzclose_p && gzread_p))
        {
          fatal("Invalid compression library (zlib)");
        }
    }
#endif

#ifdef HAVE_BZLIB_H
#ifdef _WIN32
  bz2_lib = LoadLibraryA(bz2_libname);
#else
  bz2_lib = dlopen(bz2_libname, RTLD_LAZY);
#endif
  if (bz2_lib)
    {
      BZ2_bzReadOpen_p = (BZFILE* (*)(int*, FILE*, int, int, void*, int))
        arch_dlsym(bz2_lib, "BZ2_bzReadOpen");
      BZ2_bzReadClose_p = (void (*)(int*, BZFILE*))
        arch_dlsym(bz2_lib, "BZ2_bzReadClose");
      BZ2_bzRead_p = (int (*)(int*, BZFILE*, void*, int))
        arch_dlsym(bz2_lib, "BZ2_bzRead");
      if (not (BZ2_bzReadOpen_p && BZ2_bzReadClose_p && BZ2_bzRead_p))
        {
          fatal("Invalid compression library (bz2)");
        }
    }
#endif
}

void dynlibs_close()
{
#ifdef HAVE_ZLIB_H
  if (gz_lib)
    {
#ifdef _WIN32
      FreeLibrary(gz_lib);
#else
      dlclose(gz_lib);
#endif
    }
  gz_lib = nullptr;
#endif

#ifdef HAVE_BZLIB_H
  if (bz2_lib)
    {
#ifdef _WIN32
      FreeLibrary(bz2_lib);
#else
      dlclose(bz2_lib);
#endif
    }
  bz2_lib = nullptr;
#endif
}
