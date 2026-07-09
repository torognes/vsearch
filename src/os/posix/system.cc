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

#include "system.h"
#include "utils/fatal.hpp"
#include <algorithm>  // std::max
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::size_t
#include <cstdlib>  // posix_memalign, std::realloc, std::free
#include <fcntl.h>  // open, O_RDONLY, O_WRONLY, O_CREAT, O_TRUNC
#include <sys/resource.h>  // getrusage, RUSAGE_SELF, struct rusage
#include <sys/stat.h>  // fstat, stat, struct stat, S_IRUSR, S_IWUSR
#include <unistd.h>  // sysconf, _SC_*, lseek, off_t

/* total physical RAM is queried through different interfaces depending on
   the Unix flavour; include the matching header (mirrors the per-OS choice
   made in vsearch.h). <unistd.h> above must precede this test because it
   defines the _SC_PHYS_PAGES macro. */
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <sys/sysctl.h>  // sysctl / sysctlbyname
#elif ! (defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE))
#include <sys/sysinfo.h>  // sysinfo
#endif


constexpr auto vsearch_memalignment = 16;


auto system_get_memused() -> uint64_t
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);

#ifdef __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  return r_usage.ru_maxrss;
#else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  return static_cast<uint64_t>(r_usage.ru_maxrss) * 1024;
#endif
}


auto system_get_memtotal() -> uint64_t
{
#if defined(__APPLE__)

  int mib [] = { CTL_HW, HW_MEMSIZE };
  int64_t ram = 0;
  std::size_t length = sizeof(ram);
  if(sysctl(mib, 2, &ram, &length, NULL, 0) == -1)
    fatal("Cannot determine amount of RAM");
  return ram;

#elif defined(__FreeBSD__)

  /* sysctlbyname("hw.physmem") writes a uint64_t directly, avoiding the
     32-bit overflow of the older sysctl({CTL_HW, HW_PHYSMEM}) interface
     on hosts with >= 4 GB */
  uint64_t ram = 0;
  std::size_t length = sizeof(ram);
  if (sysctlbyname("hw.physmem", &ram, &length, nullptr, 0) != 0)
    fatal("Cannot determine amount of RAM");
  return ram;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  int64_t const phys_pages = sysconf(_SC_PHYS_PAGES);
  int64_t const pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) or (pagesize == -1))
    {
      fatal("Cannot determine amount of RAM");
    }
  return static_cast<uint64_t>(pagesize * phys_pages);

#else

  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;

#endif
}


auto system_get_cores() -> long
{
  return sysconf(_SC_NPROCESSORS_ONLN);
}


auto xmalloc(std::size_t size) -> void *
{
  static constexpr auto minimal_allocation = std::size_t{1};
  size = std::max(size, minimal_allocation);
  void * ptr = nullptr;
  if (posix_memalign(&ptr, vsearch_memalignment, size) != 0)
    {
      ptr = nullptr;
    }
  if (ptr == nullptr)
    {
      fatal("Unable to allocate enough memory.");
    }
  return ptr;
}


auto xrealloc(void * ptr, std::size_t size) -> void *
{
  /* NOTE: unlike xmalloc (posix_memalign), the POSIX branch here uses plain
     realloc, which only guarantees max_align_t alignment, not xmalloc's
     vsearch_memalignment (16 bytes). Buffers that require 16-byte alignment
     for aligned SIMD loads/stores (the align_simd.cc VECTOR_SHORT arrays) must
     therefore be allocated with xmalloc and never grown through xrealloc. As
     audited, no such buffer currently passes through xrealloc — every caller
     resizes byte/scalar data (input buffers, sequence storage, query strings,
     k-mer lists). If that ever changes, emulate an aligned realloc here. */
  static constexpr auto minimal_allocation = std::size_t{1};
  size = std::max(size, minimal_allocation);
  void * new_ptr = std::realloc(ptr, size);
  if (new_ptr == nullptr)
    {
      fatal("Unable to reallocate enough memory.");
    }
  return new_ptr;
}


auto xfree(void * ptr) -> void
{
  if (ptr != nullptr)
    {
      std::free(ptr);
    }
  else
    {
      fatal("Trying to free a null pointer");
    }
}


auto xfstat(int file_descriptor, xstat_t * buf) -> int
{
  return fstat(file_descriptor, buf);  // return zero if success
}


auto xstat(const char * path, xstat_t * buf) -> int
{
  return stat(path, buf);
}


auto xlseek(int file_descriptor, uint64_t offset, int whence) -> uint64_t
{
  return static_cast<uint64_t>(lseek(file_descriptor, static_cast<off_t>(offset), whence));  // libC or linuxism: replace with std::fseek()?
}


// refactoring: only used in fastx.cc
auto xftello(std::FILE * stream) -> uint64_t
{
  return static_cast<uint64_t>(ftello(stream));
}


// refactoring: only used in udb.cc
auto xopen_read(const char * path) -> int
{
  return open(path, O_RDONLY);
}


// refactoring: only used in udb.cc
auto xopen_write(const char * path) -> int
{
  return open(path,
              O_WRONLY | O_CREAT | O_TRUNC,
              S_IRUSR | S_IWUSR);
}
