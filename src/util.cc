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
#include "city.h"
#include "md5.h"
#include "utils/maps.hpp"
#include <array>
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <climits>  // ULONG_MAX, RAND_MAX
#include <cstdarg>  // va_list
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t, std::vsnprintf, std::fopen
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::strlen, std::strcmp, std::strcpy, std::strchr
#include <ctime>  // timeval, gettimeofday
#include <iterator>  // std::next
#include <limits>
#include <vector>


//#define SHOW_RUSAGE
constexpr auto one_hundred_percent = 100UL;
constexpr auto nighty_nine_percent = 99UL;
static const char * progress_prompt;
static uint64_t progress_next;
static uint64_t progress_size;
static uint64_t progress_pct;
static bool progress_show;


auto progress_init(const char * prompt, uint64_t size) -> void
{
  progress_show = isatty(fileno(stderr)) and (not opt_quiet) and (not opt_no_progress);
  progress_prompt = prompt;
  progress_size = size;
  progress_pct = 0;
  progress_next = ((progress_pct + 1) * progress_size + nighty_nine_percent) / one_hundred_percent;

  if (opt_quiet) { return; }
  std::fprintf(stderr, "%s", prompt);
  if (not progress_show) { return; }
  std::fprintf(stderr, " %d%%", 0);
}


auto progress_update(uint64_t progress) -> void
{
  if ((progress < progress_next) or not progress_show) { return; }
  if (progress_size == 0) {
    std::fprintf(stderr, "  \r%s 0%%", progress_prompt);
    return;
  }
  progress_pct = one_hundred_percent * progress / progress_size;
  std::fprintf(stderr,
          "  \r%s %" PRIu64 "%%",
          progress_prompt,
          progress_pct);
  progress_next = ((progress_pct + 1) * progress_size + nighty_nine_percent) / one_hundred_percent;
}


auto progress_done() -> void
{
  if (opt_quiet) { return; }
  if (progress_show)
    {
      std::fprintf(stderr, "  \r%s", progress_prompt);
    }
  std::fprintf(stderr, " %ld%%\n", one_hundred_percent);
}


__attribute__((noreturn))
auto fatal(const char * msg) -> void
{
  std::fprintf(stderr, "\n\n");
  std::fprintf(stderr, "Fatal error: %s\n", msg);

  if (fp_log != nullptr)
    {
      std::fprintf(fp_log, "\n\n");
      std::fprintf(fp_log, "Fatal error: %s\n", msg);
    }

  std::exit(EXIT_FAILURE);
}


__attribute__((noreturn))
auto fatal(const char * format,
           const char * message) -> void
{
  std::fprintf(stderr, "\n\nFatal error: ");
  std::fprintf(stderr, format, message);
  std::fprintf(stderr, "\n");

  if (opt_log != nullptr)
    {
      std::fprintf(fp_log, "\n\nFatal error: ");
      std::fprintf(fp_log, format, message);
      std::fprintf(fp_log, "\n");
    }

  std::exit(EXIT_FAILURE);
}


auto xstrdup(char const * src) -> char *
{
  auto const len = std::strlen(src);
  auto * dest = (char *) xmalloc(len + 1);
  return std::strcpy(dest, src);
}


auto xstrchrnul(char * str, int target) -> char *
{
  // find the first occurrence to static_cast<char>(target)
  auto * first_occurrence = std::strchr(str, target);

  if (first_occurrence != nullptr) {
    return first_occurrence;
  }
  return std::next(str, static_cast<long>(std::strlen(str)));
}


auto xsprintf(char * * ret, const char * format, ...) -> int
{
  std::va_list args;
  va_start(args, format);
  auto len = std::vsnprintf(nullptr, 0, format, args);
  va_end(args);
  if (len < 0)
    {
      fatal("Error with vsnprintf in xsprintf");
    }
  auto * buffer = (char *) xmalloc(len + 1);
  va_start(args, format);
  len = std::vsnprintf(buffer, len + 1, format, args);
  va_end(args);
  *ret = buffer;
  return len;
}


auto hash_cityhash64(char * sequence, uint64_t length) -> uint64_t
{
  return CityHash64((const char *) sequence, length);
}


auto hash_cityhash128(char * sequence, uint64_t length) -> uint128
{
  return CityHash128((const char *) sequence, length);
}


auto show_rusage() -> void
{
#ifdef SHOW_RUSAGE
  static constexpr auto a_megabyte = 1024.0 * 1024.0;
  double user_time = 0.0;
  double system_time = 0.0;

  arch_get_user_system_time(&user_time, &system_time);

  double const megabytes = arch_get_memused() / a_megabyte;

  std::fprintf(stderr, "Time: %.3fs (user) %.3fs (sys) Memory: %.0lfMB\n",
          user_time, system_time, megabytes);

  if (opt_log)
    std::fprintf(fp_log, "Time: %.3fs (user) %.3fs (sys) Memory: %.0lfMB\n",
            user_time, system_time, megabytes);
#endif
}


auto reverse_complement(char * rc_seq, char * seq, int64_t len) -> void
{
  /* Write the reverse complementary sequence to rc_seq.
     The memory for rc_seq must be long enough for the rc_seq of the sequence
     (identical to the length of seq + 1). */

  for (auto i = 0LL; i < len; ++i) {
    auto const unsigned_char = static_cast<unsigned char>(*std::next(seq, len - 1 - i));
    auto const complement_char = static_cast<char>(chrmap_complement_vector[unsigned_char]);
    *std::next(rc_seq, i) = complement_char;
  }
  *std::next(rc_seq, len) = '\0';
}


auto random_init() -> void
{
  arch_srandom();
}


auto random_int(int64_t upper_limit) -> int64_t
{
  /*
    Generate a random integer in the range 0 to n-1, inclusive.
    n must be > 0
    The random() function returns a random number in the range
    0 to 2147483647 (=2^31-1=RAND_MAX), inclusive.
    We should avoid some of the upper generated numbers to
    avoid modulo bias.
  */
  assert(upper_limit != 0);
  int64_t const random_max = RAND_MAX;
  int64_t const limit = random_max - ((random_max + 1) % upper_limit);
  auto random_value = static_cast<int64_t>(arch_random());
  while (random_value > limit)
    {
      random_value = static_cast<int64_t>(arch_random());
    }
  return random_value % upper_limit;
}


auto random_ulong(uint64_t upper_limit) -> uint64_t
{
  /*
    Generate a random integer in the range 0 to n-1, inclusive,
    n must be > 0
  */
  assert(upper_limit != 0U);
  static constexpr auto shift_16_bits = 16U;
  static constexpr auto shift_32_bits = 32U;
  static constexpr auto shift_48_bits = 48U;
  auto const random_max = std::numeric_limits<uint64_t>::max();
  auto const limit = random_max - ((random_max - upper_limit + 1) % upper_limit);
  auto random_value = ((arch_random() << shift_48_bits) ^ (arch_random() << shift_32_bits) ^
                       (arch_random() << shift_16_bits) ^ (arch_random()));
  while (random_value > limit)
    {
      random_value = ((arch_random() << shift_48_bits) ^ (arch_random() << shift_32_bits) ^
                      (arch_random() << shift_16_bits) ^ (arch_random()));
    }
  return random_value % upper_limit;
}


auto string_normalize(char * normalized, char * raw_seq, unsigned int len) -> void
{
  /* convert string to upper case and replace U by T */
  for (auto i = 0U; i < len; ++i)
    {
      auto const unsigned_char = static_cast<unsigned char>(*raw_seq);
      auto const normalized_char = chrmap_normalize_vector[unsigned_char];
      *normalized = static_cast<unsigned char>(normalized_char);
      std::advance(normalized, 1);
      std::advance(raw_seq, 1);
    }
  *normalized = '\0';
}


auto fprint_hex(std::FILE * output_handle, unsigned char * data, int len) -> void
{
  for (auto i = 0; i < len; ++i)
    {
      std::fprintf(output_handle, "%02x", *std::next(data, i));
    }
}


auto SHA1(const unsigned char * data, unsigned long len, unsigned char * digest) -> void
{
  if (digest == nullptr)
    {
      fatal("Error in computing SHA1 digest");
    }
  SHA1_CTX a_context;
  SHA1_Init(&a_context);
  SHA1_Update(&a_context, data, len);
  SHA1_Final(&a_context, digest);
}


auto MD5(void * data, unsigned long len, unsigned char * digest) -> void
{
  if (digest == nullptr)
    {
      fatal("Error in computing MD5 digest");
    }
  MD5_CTX a_context;
  MD5_Init(&a_context);
  MD5_Update(&a_context, data, len);
  MD5_Final(digest, &a_context);
}


static const std::array<char, 16> hexdigits = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};

auto get_hex_seq_digest_sha1(char * hex, char * seq, int seqlen) -> void
{
  /* Save hexadecimal representation of the SHA1 hash of the sequence.
     The string array digest must be large enough (len_hex_dig_sha1).
     First normalize string by uppercasing it and replacing U's with T's. */

  std::vector<char> normalized(seqlen + 1);
  string_normalize(normalized.data(), seq, seqlen);

  std::vector<unsigned char> digest(sha1_digest_length);

  SHA1((const unsigned char *) normalized.data(), (size_t) seqlen, digest.data());

  for (int i = 0; i < sha1_digest_length; ++i)
    {
      hex[(2 * i) + 0] = hexdigits[digest[i] >> 4U];
      hex[(2 * i) + 1] = hexdigits[digest[i] & 15U];
    }
  hex[2 * sha1_digest_length] = '\0';
}


auto get_hex_seq_digest_md5(char * hex, char * seq, int seqlen) -> void
{
  /* Save hexadecimal representation of the MD5 hash of the sequence.
     The string array digest must be large enough (len_hex_dig_md5).
     First normalize string by uppercasing it and replacing U's with T's. */

  char * normalized = (char *) xmalloc(seqlen + 1);
  string_normalize(normalized, seq, seqlen);

  unsigned char digest[md5_digest_length];

  MD5(normalized, (size_t) seqlen, digest);

  xfree(normalized);

  for (int i = 0; i < md5_digest_length; i++)
    {
      hex[(2 * i) + 0] = hexdigits[digest[i] >> 4U];
      hex[(2 * i) + 1] = hexdigits[digest[i] & 15U];
    }
  hex[2 * md5_digest_length] = 0;
}


auto fprint_seq_digest_sha1(FILE * fp, char * seq, int seqlen) -> void
{
  char digest[len_hex_dig_sha1];
  get_hex_seq_digest_sha1(digest, seq, seqlen);
  std::fprintf(fp, "%s", digest);
}


auto fprint_seq_digest_md5(FILE * fp, char * seq, int seqlen) -> void
{
  char digest[len_hex_dig_md5];
  get_hex_seq_digest_md5(digest, seq, seqlen);
  std::fprintf(fp, "%s", digest);
}


auto fopen_input(const char * filename) -> FILE *
{
  /* open the input stream given by filename, but use stdin if name is - */
  if (std::strcmp(filename, "-") == 0)
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
  if (std::strcmp(filename, "-") == 0)
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
