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

#include <cstdint>  // uint64_t
#include <cstdio>   // std::FILE, std::size_t


#ifdef _WIN32
using xstat_t = struct __stat64;
#else
using xstat_t = struct stat;
#endif

auto arch_get_memused() -> uint64_t;
auto arch_get_memtotal() -> uint64_t;
auto arch_get_cores() -> long;
auto arch_get_user_system_time(double * user_time, double * system_time) -> void;
auto arch_srandom() -> void;
auto arch_random() -> uint64_t;
auto xmalloc(std::size_t size) -> void *;
auto xrealloc(void * ptr, std::size_t size) -> void *;
auto xfree(void * ptr) -> void;

auto xfstat(int file_descriptor, xstat_t * buf) -> int;
auto xstat(const char * path, xstat_t  * buf) -> int;
auto xlseek(int file_descriptor, uint64_t offset, int whence) -> uint64_t;
auto xftello(std::FILE * stream) -> uint64_t;

auto xopen_read(const char * path) -> int;
auto xopen_write(const char * path) -> int;

auto xstrcasestr(const char * haystack, const char * needle) -> const char *;

#ifdef _WIN32
auto arch_dlsym(HMODULE handle, const char * symbol) -> FARPROC;
#else
auto arch_dlsym(void * handle, const char * symbol) -> void *;
#endif
