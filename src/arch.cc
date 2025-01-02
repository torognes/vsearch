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
#include "dynlibs.h"
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t
#include <cstdlib>  // std::realloc, std::free
#include <string.h>  // strcasestr


constexpr auto memalignment = 16;


auto arch_get_memused() -> uint64_t
{
#ifdef _WIN32

  PROCESS_MEMORY_COUNTERS pmc;
  GetProcessMemoryInfo(GetCurrentProcess(),
                       &pmc,
                       sizeof(PROCESS_MEMORY_COUNTERS));
  return pmc.PeakWorkingSetSize;

#else

  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);

# ifdef __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  return r_usage.ru_maxrss;
# else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  return r_usage.ru_maxrss * 1024;
# endif

#endif
}


auto arch_get_memtotal() -> uint64_t
{
#ifdef _WIN32

  MEMORYSTATUSEX ms;
  ms.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&ms);
  return ms.ullTotalPhys;

#elif defined(__APPLE__)

  int mib [] = { CTL_HW, HW_MEMSIZE };
  int64_t ram = 0;
  size_t length = sizeof(ram);
  if(sysctl(mib, 2, &ram, &length, NULL, 0) == -1)
    fatal("Cannot determine amount of RAM");
  return ram;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  int64_t const phys_pages = sysconf(_SC_PHYS_PAGES);
  int64_t const pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) || (pagesize == -1))
    {
      fatal("Cannot determine amount of RAM");
    }
  return pagesize * phys_pages;

#else

  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;

#endif
}


auto arch_get_cores() -> long
{
#ifdef _WIN32
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwNumberOfProcessors;
#else
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}


auto arch_get_user_system_time(double * user_time, double * system_time) -> void
{
  *user_time = 0;
  *system_time = 0;
#ifdef _WIN32
  HANDLE hProcess = GetCurrentProcess();
  FILETIME ftCreation, ftExit, ftKernel, ftUser;
  ULARGE_INTEGER ul;
  GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser);
  ul.u.HighPart = ftUser.dwHighDateTime;
  ul.u.LowPart = ftUser.dwLowDateTime;
  *user_time = ul.QuadPart * 100.0e-9;
  ul.u.HighPart = ftKernel.dwHighDateTime;
  ul.u.LowPart = ftKernel.dwLowDateTime;
  *system_time = ul.QuadPart * 100.0e-9;
#else
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);
  * user_time = r_usage.ru_utime.tv_sec * 1.0
    + r_usage.ru_utime.tv_usec * 1.0e-6;
  * system_time = r_usage.ru_stime.tv_sec * 1.0
    + r_usage.ru_stime.tv_usec * 1.0e-6;
#endif
}


auto arch_srandom() -> void
{
  /* initialize pseudo-random number generator */
  unsigned int seed = opt_randseed;
  if (seed == 0)
    {
#ifdef _WIN32
      srand(GetTickCount());
#else
      int const fd = open("/dev/urandom", O_RDONLY);
      if (fd < 0)
        {
          fatal("Unable to open /dev/urandom");
        }
      if (read(fd, & seed, sizeof(seed)) < 0)
        {
          fatal("Unable to read from /dev/urandom");
        }
      close(fd);
      srandom(seed);
#endif
    }
  else
    {
#ifdef _WIN32
      srand(seed);
#else
      srandom(seed);
#endif
    }
}


auto arch_random() -> uint64_t
{
#ifdef _WIN32
  return rand();
#else
  return random();
#endif
}


auto xmalloc(size_t size) -> void *
{
  if (size == 0)
    {
      size = 1;
    }
  void * ptr = nullptr;
#ifdef _WIN32
  ptr = _aligned_malloc(size, memalignment);
#else
  if (posix_memalign(& ptr, memalignment, size) != 0)
    {
      ptr = nullptr;
    }
#endif
  if (ptr == nullptr)
    {
      fatal("Unable to allocate enough memory.");
    }
  return ptr;
}


auto xrealloc(void * ptr, size_t size) -> void *
{
  if (size == 0)
    {
      size = 1;
    }
#ifdef _WIN32
  void * t = _aligned_realloc(ptr, size, memalignment);
#else
  void * t = realloc(ptr, size);
#endif
  if (not t)
    {
      fatal("Unable to reallocate enough memory.");
    }
  return t;
}


auto xfree(void * ptr) -> void
{
  if (ptr)
    {
#ifdef _WIN32
      _aligned_free(ptr);
#else
      free(ptr);
#endif
    }
  else
    {
      fatal("Trying to free a null pointer");
    }
}


auto xfstat(int file_descriptor, xstat_t * buf) -> int
{
#ifdef _WIN32
  return _fstat64(file_descriptor, buf);
#else
  return fstat(file_descriptor, buf);
#endif
}


auto xstat(const char * path, xstat_t * buf) -> int
{
#ifdef _WIN32
  return _stat64(path, buf);
#else
  return stat(path, buf);
#endif
}


auto xlseek(int file_descriptor, uint64_t offset, int whence) -> uint64_t
{
#ifdef _WIN32
  return _lseeki64(file_descriptor, offset, whence);
#else
  return lseek(file_descriptor, offset, whence);
#endif
}


auto xftello(std::FILE * stream) -> uint64_t
{
#ifdef _WIN32
  return _ftelli64(stream);
#else
  return ftello(stream);
#endif
}


auto xopen_read(const char * path) -> int
{
#ifdef _WIN32
  return _open(path, _O_RDONLY | _O_BINARY);
#else
  return open(path, O_RDONLY);
#endif
}


auto xopen_write(const char * path) -> int
{
#ifdef _WIN32
  return _open(path,
               _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY,
               _S_IREAD | _S_IWRITE);
#else
  return open(path,
              O_WRONLY | O_CREAT | O_TRUNC,
              S_IRUSR | S_IWUSR);
#endif
}


auto xstrcasestr(const char * haystack, const char * needle) -> const char *
{
#ifdef _WIN32
  return StrStrIA(haystack, needle);
#else
  return strcasestr(haystack, needle);
#endif
}


#ifdef _WIN32
auto arch_dlsym(HMODULE handle, const char * symbol) -> FARPROC
{
  return GetProcAddress(handle, symbol);
}
#else
auto arch_dlsym(void * handle, const char * symbol) -> void *
{
  return dlsym(handle, symbol);
}
#endif
