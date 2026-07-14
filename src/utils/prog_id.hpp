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

#ifdef HAVE_CONFIG_H
#include "config.h"  // PACKAGE, PACKAGE_VERSION
#endif

/* Program identity strings reported in the --version banner and in the
   headers of several output formats. Split out of the umbrella vsearch header
   so the per-architecture / per-OS preprocessor ladders that build the
   "os_cpu" tag are confined to the few TUs that print them (vsearch.cc,
   core/otutable.cpp, core/results.cpp). Only the identity *string* is defined
   here; the SIMD intrinsics and the OS runtime headers live elsewhere. */

#define PROG_NAME PACKAGE
#define PROG_VERSION PACKAGE_VERSION

#ifdef __x86_64__
#define PROG_CPU "x86_64"
#elif __PPC__
#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#else
#error Big endian ppc64 CPUs not supported
#endif
#elif __aarch64__
#define PROG_CPU "aarch64"
#else
#define PROG_CPU "simde"
#endif

#ifdef _WIN32
#define PROG_OS "win"
#elif __APPLE__
#define PROG_OS "macos"
#elif __linux__
#define PROG_OS "linux"
#elif __FreeBSD__
#define PROG_OS "freebsd"
#elif __NetBSD__
#define PROG_OS "netbsd"
#else
#define PROG_OS "unknown"
#endif

#define PROG_ARCH PROG_OS "_" PROG_CPU
