/*
    Copyright (C) 2014 Torbjorn Rognes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "vsearch.h"

unsigned long arch_get_memused()
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);
  
#if defined __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  return r_usage.ru_maxrss;
#else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  return r_usage.ru_maxrss * 1024;
#endif
}

unsigned long arch_get_memtotal()
{
#ifdef __APPLE__
  int mib [] = { CTL_HW, HW_MEMSIZE };
  int64_t ram = 0;
  size_t length = sizeof(ram);
  if(-1 == sysctl(mib, 2, &ram, &length, NULL, 0))
    fatal("Cannot determine amount of RAM");
  return ram;
#else
  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;
#endif
}

void arch_srandom_init()
{
  /* init random generator with some random data */
#ifdef __APPLE__
  srandomdev();
#else
  unsigned int seed;
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0)
    fatal("Unable to open /dev/urandom");
  if (read(fd, & seed, sizeof(seed)) < 0)
    fatal("Unable to read from /dev/urandom");
#endif
}
