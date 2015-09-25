/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

static const char hexdigits[] = "0123456789abcdef";

void get_hex_seq_digest_sha1(char * hex, char * seq, int seqlen)
{
  /* Save hexadecimal representation of the SHA1 hash of the sequence.
     The string array digest must be large enough (LEN_HEX_DIG_SHA1).
     First normalize string by uppercasing it and replacing U's with T's. */

  char * normalized = (char*) xmalloc(seqlen+1);
  string_normalize(normalized, seq, seqlen);

  unsigned char digest[LEN_DIG_SHA1];

#if __APPLE__
  CC_SHA1(normalized, (CC_LONG) seqlen, digest);
#else
  SHA1((const unsigned char*)normalized, (size_t) seqlen, digest);
#endif

  free(normalized);

  for(int i=0; i<LEN_DIG_SHA1; i++)
    {
      hex[2*i+0] = hexdigits[digest[i] >> 4];
      hex[2*i+1] = hexdigits[digest[i] & 15];
    }
  hex[2*LEN_DIG_SHA1] = 0;
}

void get_hex_seq_digest_md5(char * hex, char * seq, int seqlen)
{
  /* Save hexadecimal representation of the MD5 hash of the sequence.
     The string array digest must be large enough (LEN_HEX_DIG_MD5).
     First normalize string by uppercasing it and replacing U's with T's. */

  char * normalized = (char*) xmalloc(seqlen+1);
  string_normalize(normalized, seq, seqlen);

  unsigned char digest[LEN_DIG_MD5];

#if __APPLE__
  CC_MD5(normalized, (CC_LONG) seqlen, digest);
#else
  MD5((const unsigned char*)normalized, (size_t) seqlen, digest);
#endif

  free(normalized);

  for(int i=0; i<LEN_DIG_MD5; i++)
    {
      hex[2*i+0] = hexdigits[digest[i] >> 4];
      hex[2*i+1] = hexdigits[digest[i] & 15];
    }
  hex[2*LEN_DIG_MD5] = 0;
}


void fprint_seq_digest_sha1(FILE * fp, char * seq, int seqlen)
{
  char digest[LEN_HEX_DIG_SHA1];
  get_hex_seq_digest_sha1(digest, seq, seqlen);
  fprintf(fp, "%s", digest);
}

void fprint_seq_digest_md5(FILE * fp, char * seq, int seqlen)
{
  char digest[LEN_HEX_DIG_MD5];
  get_hex_seq_digest_md5(digest, seq, seqlen);
  fprintf(fp, "%s", digest);
}
