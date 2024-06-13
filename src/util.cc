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
#include "md5.h"


//#define SHOW_RUSAGE

static const char * progress_prompt;
static uint64_t progress_next;
static uint64_t progress_size;
static uint64_t progress_pct;
static bool progress_show;

auto progress_init(const char * prompt, uint64_t size) -> void
{
  progress_show = isatty(fileno(stderr)) && (! opt_quiet) && (! opt_no_progress);
  progress_prompt = prompt;
  progress_size = size;
  progress_pct = 0;
  progress_next = ((progress_pct + 1) * progress_size + 99) / 100;

  if (! opt_quiet)
    {
      fprintf(stderr, "%s", prompt);
      if (progress_show)
        {
          fprintf(stderr, " %d%%", 0);
        }
    }
}

auto progress_update(uint64_t progress) -> void
{
  if ((progress >= progress_next) && progress_show)
    {
      if (progress_size > 0)
        {
          progress_pct = 100 * progress / progress_size;
          fprintf(stderr,
                  "  \r%s %" PRIu64 "%%",
                  progress_prompt,
                  progress_pct);
          progress_next = ((progress_pct + 1) * progress_size + 99) / 100;
        }
      else
        {
          fprintf(stderr, "  \r%s 0%%", progress_prompt);
        }
    }
}

auto progress_done() -> void
{
  if (! opt_quiet)
    {
      if (progress_show)
        {
          fprintf(stderr, "  \r%s", progress_prompt);
        }
      fprintf(stderr, " %d%%\n", 100);
    }
}

__attribute__((noreturn))
auto fatal(const char * msg) -> void
{
  fprintf(stderr, "\n\n");
  fprintf(stderr, "Fatal error: %s\n", msg);

  if (fp_log)
    {
      fprintf(fp_log, "\n\n");
      fprintf(fp_log, "Fatal error: %s\n", msg);
    }

  exit(EXIT_FAILURE);
}

__attribute__((noreturn))
auto fatal(const char * format,
           const char * message) -> void
{
  fprintf(stderr, "\n\nFatal error: ");
  fprintf(stderr, format, message);
  fprintf(stderr, "\n");

  if (opt_log)
    {
      fprintf(fp_log, "\n\nFatal error: ");
      fprintf(fp_log, format, message);
      fprintf(fp_log, "\n");
    }

  exit(EXIT_FAILURE);
}

auto xstrdup(const char *s) -> char *
{
  size_t const len = strlen(s);
  char * p = (char *) xmalloc(len + 1);
  return strcpy(p, s);
}


auto xstrchrnul(char *s, int c) -> char *
{
  char * r = strchr(s, c);

  if (r)
    {
      return r;
    }
  else
    {
      return s + strlen(s);
    }
}

auto xsprintf(char * * ret, const char * format, ...) -> int
{
  va_list ap;
  va_start(ap, format);
  int len = vsnprintf(nullptr, 0, format, ap);
  va_end(ap);
  if (len < 0)
    {
      fatal("Error with vsnprintf in xsprintf");
    }
  char * p = (char *) xmalloc(len + 1);
  va_start(ap, format);
  len = vsnprintf(p, len + 1, format, ap);
  va_end(ap);
  *ret = p;
  return len;
}

auto hash_cityhash64(char * s, uint64_t n) -> uint64_t
{
  return CityHash64((const char *) s, n);
}

auto hash_cityhash128(char * s, uint64_t n) -> uint128
{
  return CityHash128((const char *) s, n);
}

auto getusec() -> int64_t
{
  struct timeval tv;
  if (gettimeofday(&tv,nullptr) != 0)
    {
      return 0;
    }
  return (tv.tv_sec * 1000000) + tv.tv_usec;
}

auto show_rusage() -> void
{
#ifdef SHOW_RUSAGE
  double user_time = 0.0;
  double system_time = 0.0;

  arch_get_user_system_time(&user_time, &system_time);

  double megabytes = arch_get_memused() / 1024.0 / 1024.0;

  fprintf(stderr, "Time: %.3fs (user) %.3fs (sys) Memory: %.0lfMB\n",
          user_time, system_time, megabytes);

  if (opt_log)
    fprintf(fp_log, "Time: %.3fs (user) %.3fs (sys) Memory: %.0lfMB\n",
            user_time, system_time, megabytes);
#endif
}

auto reverse_complement(char * rc, char * seq, int64_t len) -> void
{
  /* Write the reverse complementary sequence to rc.
     The memory for rc must be long enough for the rc of the sequence
     (identical to the length of seq + 1. */

  for (int64_t i = 0; i < len; i++)
    {
      rc[i] = chrmap_complement[(int) (seq[len - 1 - i])];
    }
  rc[len] = 0;
}

auto random_init() -> void
{
  arch_srandom();
}

auto random_int(int64_t n) -> int64_t
{
  /*
    Generate a random integer in the range 0 to n-1, inclusive.
    n must be > 0
    The random() function returns a random number in the range
    0 to 2147483647 (=2^31-1=RAND_MAX), inclusive.
    We should avoid some of the upper generated numbers to
    avoid modulo bias.
  */

  int64_t const random_max = RAND_MAX;
  int64_t const limit = random_max - ((random_max + 1) % n);
  int64_t r = arch_random();
  while (r > limit)
    {
      r = arch_random();
    }
  return r % n;
}

auto random_ulong(uint64_t n) -> uint64_t
{
  /*
    Generate a random integer in the range 0 to n-1, inclusive,
    n must be > 0
  */

  uint64_t const random_max = ULONG_MAX;
  uint64_t const limit = random_max - ((random_max - n + 1) % n);
  uint64_t r = ((arch_random() << 48) ^ (arch_random() << 32) ^
                (arch_random() << 16) ^ (arch_random()));
  while (r > limit)
    {
      r = ((arch_random() << 48) ^ (arch_random() << 32) ^
           (arch_random() << 16) ^ (arch_random()));
    }
  return r % n;
}

auto string_normalize(char * normalized, char * s, unsigned int len) -> void
{
  /* convert string to upper case and replace U by T */
  char * p = s;
  char * q = normalized;
  for (unsigned int i = 0; i < len; i++)
    {
      *q++ = chrmap_normalize[(int) (*p++)];
    }
  *q = 0;
}

auto fprint_hex(FILE * fp, unsigned char * data, int len) -> void
{
  for (int i = 0; i < len; i++)
    {
      fprintf(fp, "%02x", data[i]);
    }
}

auto SHA1(const unsigned char * d, unsigned long n, unsigned char * md) -> void
{
  if (! md)
    {
      fatal("Error in computing SHA1 digest");
    }
  SHA1_CTX c;
  SHA1_Init(&c);
  SHA1_Update(&c, d, n);
  SHA1_Final(&c, md);
}

auto MD5(void * d, unsigned long n, unsigned char * md) -> void
{
  if (! md)
    {
      fatal("Error in computing MD5 digest");
    }
  MD5_CTX c;
  MD5_Init(&c);
  MD5_Update(&c, d, n);
  MD5_Final(md, &c);
}

static const char hexdigits[] = "0123456789abcdef";

auto get_hex_seq_digest_sha1(char * hex, char * seq, int seqlen) -> void
{
  /* Save hexadecimal representation of the SHA1 hash of the sequence.
     The string array digest must be large enough (LEN_HEX_DIG_SHA1).
     First normalize string by uppercasing it and replacing U's with T's. */

  char * normalized = (char *) xmalloc(seqlen+1);
  string_normalize(normalized, seq, seqlen);

  unsigned char digest[LEN_DIG_SHA1];

  SHA1((const unsigned char*) normalized, (size_t) seqlen, digest);

  xfree(normalized);

  for (int i = 0; i < LEN_DIG_SHA1; i++)
    {
      hex[(2 * i) + 0] = hexdigits[digest[i] >> 4];
      hex[(2 * i) + 1] = hexdigits[digest[i] & 15];
    }
  hex[2 * LEN_DIG_SHA1] = 0;
}

auto get_hex_seq_digest_md5(char * hex, char * seq, int seqlen) -> void
{
  /* Save hexadecimal representation of the MD5 hash of the sequence.
     The string array digest must be large enough (LEN_HEX_DIG_MD5).
     First normalize string by uppercasing it and replacing U's with T's. */

  char * normalized = (char *) xmalloc(seqlen + 1);
  string_normalize(normalized, seq, seqlen);

  unsigned char digest[MD5_DIGEST_LENGTH];

  MD5(normalized, (size_t) seqlen, digest);

  xfree(normalized);

  for (int i = 0; i < MD5_DIGEST_LENGTH; i++)
    {
      hex[(2 * i) + 0] = hexdigits[digest[i] >> 4];
      hex[(2 * i) + 1] = hexdigits[digest[i] & 15];
    }
  hex[2 * MD5_DIGEST_LENGTH] = 0;
}


auto fprint_seq_digest_sha1(FILE * fp, char * seq, int seqlen) -> void
{
  char digest[LEN_HEX_DIG_SHA1];
  get_hex_seq_digest_sha1(digest, seq, seqlen);
  fprintf(fp, "%s", digest);
}

auto fprint_seq_digest_md5(FILE * fp, char * seq, int seqlen) -> void
{
  char digest[LEN_HEX_DIG_MD5];
  get_hex_seq_digest_md5(digest, seq, seqlen);
  fprintf(fp, "%s", digest);
}

auto fopen_input(const char * filename) -> FILE *
{
  /* open the input stream given by filename, but use stdin if name is - */
  if (strcmp(filename, "-") == 0)
    {
      int const fd = dup(STDIN_FILENO);
      if (fd < 0)
        {
          return nullptr;
        }
      else
        {
          return fdopen(fd, "rb");
        }
    }
  else
    {
      return fopen(filename, "rb");
    }
}

auto fopen_output(const char * filename) -> FILE *
{
  /* open the output stream given by filename, but use stdout if name is - */
  if (strcmp(filename, "-") == 0)
    {
      int const fd = dup(STDOUT_FILENO);
      if (fd < 0)
        {
          return nullptr;
        }
      else
        {
          return fdopen(fd, "w");
        }
    }
  else
    {
      return fopen(filename, "w");
    }
}
