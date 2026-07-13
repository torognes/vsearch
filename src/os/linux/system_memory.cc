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

#include "os/system.hpp"
#include "utils/fatal.hpp"
#include <cstdint>  // int64_t, uint64_t
#include <sys/resource.h>  // getrusage, RUSAGE_SELF, struct rusage
#include <unistd.h>  // sysconf, _SC_PHYS_PAGES, _SC_PAGESIZE

/* <unistd.h> above defines _SC_PHYS_PAGES; sysconf() is preferred and the
   sysinfo() branch is the fallback for the rare Unix that lacks it. This file
   also serves as the generic-Unix backend (the build's default when the host
   is neither Windows, macOS nor FreeBSD). */
#if ! (defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE))
#include <sys/sysinfo.h>  // sysinfo
#endif


auto system_get_memused() -> uint64_t
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);
  /* Linux: ru_maxrss gives the size in kilobytes  */
  return static_cast<uint64_t>(r_usage.ru_maxrss) * 1024;
}


auto system_get_memtotal() -> uint64_t
{
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

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
