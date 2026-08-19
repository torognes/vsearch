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

#include <cstdint>  // uint64_t
#include <cstdio>   // std::FILE, std::size_t

// <sys/stat.h> provides struct stat and, on mingw, the __stat64 -> _stat64
// macro that xstat_t depends on. It must precede the alias below so this
// header is self-contained rather than relying on the includer to pull in
// <sys/stat.h> in the right order first.
#include <sys/stat.h>


#ifdef _WIN32
using xstat_t = struct __stat64;
#else
using xstat_t = struct stat;
#endif

auto system_get_memused() -> uint64_t;

/* The machine's memory. Prefer system_get_memlimit() below for anything that
   budgets or reports what this run may use: this figure is the host's, and a
   container's cgroup limit is invisible in it, so inside one it can be orders
   of magnitude too large. */
auto system_get_memtotal() -> uint64_t;

/* The memory this process may actually use: system_get_memtotal(), reduced to
   the smallest limit that applies to it. On Linux that is the cgroup memory
   limit, which is how Slurm, Docker, podman and Kubernetes cap a job -- they
   share the host's kernel, so there is no smaller "total" for
   system_get_memtotal() to report (issue #584). Platforms with no such
   mechanism visible return system_get_memtotal() unchanged.

   A virtual machine needs nothing here: it has its own kernel, and
   system_get_memtotal() already reports what the hypervisor gave the guest.

   Never zero, and never larger than system_get_memtotal(). */
auto system_get_memlimit() -> uint64_t;

/* The machine's online cores. As with system_get_memtotal(), prefer
   system_get_available_cores() below wherever the question is what this run
   may use. */
auto system_get_cores() -> long;

/* The cores this process may actually run on: system_get_cores(), reduced by
   every restriction that applies to it. On Linux that is the CPU affinity mask
   (taskset, a cpuset cgroup, Slurm's --cpu-bind, "docker --cpuset-cpus") and
   the cgroup CPU quota ("docker --cpus", a Kubernetes CPU limit); on Windows
   the process affinity mask.

   Deliberately reads no environment variable: OMP_NUM_THREADS belongs to an
   OpenMP runtime, which vsearch is not, and the kernel interfaces above
   already cover the schedulers that set it.

   Never zero, and never larger than system_get_cores(), except that a host
   where sysconf() cannot answer at all sees that failure passed through
   unchanged by both. */
auto system_get_available_cores() -> long;
auto xmalloc(std::size_t size) -> void *;
auto xfree(void * ptr) -> void;

auto xfstat(int file_descriptor, xstat_t * buf) -> int;
auto xstat(const char * path, xstat_t  * buf) -> int;
/* Where the descriptor currently is. POSIX has no tell() for a descriptor,
   so this is the lseek(fd, 0, SEEK_CUR) idiom; the offset and whence used to
   be parameters, but the one caller only ever asked for the position, so the
   name now says what the call does. The stream counterpart is xftello(). */
auto xtell_fd(int file_descriptor) -> uint64_t;
auto xftello(std::FILE * stream) -> uint64_t;
