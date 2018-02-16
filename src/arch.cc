/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

const int memalignment = 16;

uint64_t arch_get_memused()
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

uint64_t arch_get_memtotal()
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

  int64_t phys_pages = sysconf(_SC_PHYS_PAGES);
  int64_t pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) || (pagesize == -1))
    fatal("Cannot determine amount of RAM");
  return pagesize * phys_pages;

#else

  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;

#endif
}

long arch_get_cores()
{
#ifdef _WIN32
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwNumberOfProcessors;
#else
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

void arch_get_user_system_time(double * user_time, double * system_time)
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

void arch_srandom()
{
  /* initialize pseudo-random number generator */
  unsigned int seed = opt_randseed;
  if (seed == 0)
    {
#ifdef _WIN32
      srand(GetTickCount());
#else
      int fd = open("/dev/urandom", O_RDONLY);
      if (fd < 0)
        fatal("Unable to open /dev/urandom");
      if (read(fd, & seed, sizeof(seed)) < 0)
        fatal("Unable to read from /dev/urandom");
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

uint64_t arch_random()
{
#ifdef _WIN32
  return rand();
#else
  return random();
#endif
}

void * xmalloc(size_t size)
{
  if (size == 0)
    size = 1;
  void * t = 0;
#ifdef _WIN32
  t = _aligned_malloc(size, memalignment);
#else
  if (posix_memalign(& t, memalignment, size))
    t = 0;
#endif
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  if (size == 0)
    size = 1;
#ifdef _WIN32
  void * t = _aligned_realloc(ptr, size, memalignment);
#else
  void * t = realloc(ptr, size);
#endif
  if (!t)
    fatal("Unable to reallocate enough memory.");
  return t;
}

void xfree(void * ptr)
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
    fatal("Trying to free a null pointer");
}

int xfstat(int fd, xstat_t * buf)
{
#ifdef _WIN32
  return _fstat64(fd, buf);
#else
  return fstat(fd, buf);
#endif
}

int xstat(const char * path, xstat_t  * buf)
{
#ifdef _WIN32
  return _stat64(path, buf);
#else
  return stat(path, buf);
#endif
}

uint64_t xlseek(int fd, uint64_t offset, int whence)
{
#ifdef _WIN32
  return _lseeki64(fd, offset, whence);
#else
  return lseek(fd, offset, whence);
#endif
}

uint64_t xftello(FILE * stream)
{
#ifdef _WIN32
  return _ftelli64(stream);
#else
  return ftello(stream);
#endif
}

int xopen_read(const char * path)
{
#ifdef _WIN32
  return _open(path, _O_RDONLY | _O_BINARY);
#else
  return open(path, O_RDONLY);
#endif
}

int xopen_write(const char * path)
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
