/*
    Copyright (C) 2014-2015 Torbjorn Rognes

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
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  long phys_pages = sysconf(_SC_PHYS_PAGES);
  long pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) || (pagesize == -1))
    fatal("Cannot determine amount of RAM");
  return pagesize * phys_pages;

#elif defined(__APPLE__)

  int mib [] = { CTL_HW, HW_MEMSIZE };
  int64_t ram = 0;
  size_t length = sizeof(ram);
  if(sysctl(mib, 2, &ram, &length, NULL, 0) == -1)
    fatal("Cannot determine amount of RAM");
  return ram;

#else

  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;

#endif
}

void fprint_seq_digest_sha1(FILE * fp, char * seq, int seqlen)
{
  /* Print a hexadecimal representation of the sha1 hash of the sequence to fp.
     First normalize string by uppercasing it and replacing U's with T's. */

  char * normalized = (char*) xmalloc(seqlen+1);
  string_normalize(normalized, seq, seqlen);

#if __APPLE__
  #define DIG_LEN_SHA1 CC_SHA1_DIGEST_LENGTH
  unsigned char digest[DIG_LEN_SHA1];
  CC_SHA1(normalized, (CC_LONG) seqlen, digest);
#else
  #define DIG_LEN_SHA1 SHA_DIGEST_LENGTH
  unsigned char digest[DIG_LEN_SHA1];
  SHA1((const unsigned char*)normalized, (size_t) seqlen, digest);
#endif

  free(normalized);
  fprint_hex(fp, digest, DIG_LEN_SHA1);
}

void fprint_seq_digest_md5(FILE * fp, char * seq, int seqlen)
{
  /* Print a hexadecimal representation of the md5 hash of the sequence to fp.
     First normalize string by uppercasing it and replacing U's with T's. */

  char * normalized = (char*) xmalloc(seqlen+1);
  string_normalize(normalized, seq, seqlen);

#ifdef __APPLE__
  #define DIG_LEN_MD5 CC_MD5_DIGEST_LENGTH
  unsigned char digest[DIG_LEN_MD5];
  CC_MD5(normalized, (CC_LONG) seqlen, digest);
#else
  #define DIG_LEN_MD5 MD5_DIGEST_LENGTH
  unsigned char digest[DIG_LEN_MD5];
  MD5((const unsigned char*)normalized, (size_t) seqlen, digest);
#endif

  free(normalized);
  fprint_hex(fp, digest, DIG_LEN_MD5);
}
