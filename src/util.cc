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


long gcd(long a, long b)
{
  if (b == 0)
  {
    return a;
  }
  else
  {
    return gcd(b, a % b);
  }
}

void fatal(const char * msg)
{
  fprintf(stderr, "Error: %s\n", msg);
  exit(1);
}

void fatal(const char * format, const char * message)
{
  fprintf(stderr, format, message);
  fprintf(stderr, "\n");
  exit(1);
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t;
  posix_memalign(& t, alignment, size);

  if (t==NULL)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

char * xstrchrnul(char *s, int c)
{
  char * r = strchr(s, c);

  if (r)
    return r;
  else
    return (char *)s + strlen(s);
}

unsigned long hash_fnv_1a_64(unsigned char * s, unsigned long n)
{
  const unsigned long fnv_offset = 14695981039346656037UL;
  const unsigned long fnv_prime = 1099511628211;

  unsigned long hash = fnv_offset;

  for(unsigned long i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = (hash ^ c) * fnv_prime;
    }

  return hash;
}

long getusec(void)
{
  struct timeval tv;
  if(gettimeofday(&tv,0) != 0) return 0;
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

void show_rusage()
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);
  
  fprintf(stderr, "Time: %.3fs (user)", r_usage.ru_utime.tv_sec * 1.0 + (double) r_usage.ru_utime.tv_usec * 1.0e-6);
  fprintf(stderr, " %.3fs (sys)", r_usage.ru_stime.tv_sec * 1.0 + r_usage.ru_stime.tv_usec * 1.0e-6);

#if defined __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  fprintf(stderr, " Memory: %.0fMB\n", r_usage.ru_maxrss * 1.0e-6);
#else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  fprintf(stderr, " Memory: %.0fMB\n", r_usage.ru_maxrss * 1.0e-3);
#endif
}

