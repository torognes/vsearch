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

#include "sequence_digest.hpp"
#include "fatal.hpp"  // fatal
#include "string_normalize.hpp"  // string_normalize
#include "vendored/md5.h"  // MD5_CTX, MD5_Init, MD5_Update, MD5_Final
#include "vendored/sha1.h"  // SHA1_CTX, SHA1_Init, SHA1_Update, SHA1_Final
#include <array>  // std::array
#include <cstddef>  // std::size_t
#include <cstdio>  // std::FILE, std::fprintf
#include <iterator>  // std::advance
#include <vector>  // std::vector


namespace {

constexpr auto drop_lower_nibble = 4U;
constexpr auto mask_upper_nibble = 15U;
constexpr std::array<char, 16> hexdigits =
  {{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'}};


auto SHA1(unsigned char const * data, unsigned long const len, unsigned char * digest) -> void
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


auto MD5(void * data, unsigned long const len, unsigned char * digest) -> void
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

}  // namespace


auto get_hex_seq_digest_sha1(char * hex, char const * seq, int const seqlen) -> void
{
  /* Save hexadecimal representation of the SHA1 hash of the sequence.
     The string array digest must be large enough (len_hex_dig_sha1).
     First normalize string by uppercasing it and replacing U's with T's.

     The hash always runs on this private, per-call `normalized` copy, never
     on `seq` directly: the vendored SHA-1 may reinterpret/byte-swap its input
     in place, so it must only ever touch a throwaway, suitably aligned buffer.
     Being per-call (no shared state), this is also what makes it safe to call
     from the search worker threads (e.g. --relabel_sha1). Do not hash `seq`
     directly, and do not enable SHA1HANDSOFF (a shared static workspace) -- see
     S26 in CODE_REVIEW.md. */

  std::vector<char> normalized(static_cast<std::size_t>(seqlen) + 1);
  string_normalize(normalized.data(), seq, static_cast<unsigned int>(seqlen));

  std::vector<unsigned char> digest(sha1_digest_length);

  SHA1(reinterpret_cast<unsigned char const *>(normalized.data()),
       static_cast<std::size_t>(seqlen),
       digest.data());

  for (auto const & element: digest) {
    *hex = hexdigits[element >> drop_lower_nibble];
    std::advance(hex, 1);
    *hex = hexdigits[element & mask_upper_nibble];
    std::advance(hex, 1);
  }
  *hex = '\0';
}


auto get_hex_seq_digest_md5(char * hex, char const * seq, int const seqlen) -> void
{
  /* Save hexadecimal representation of the MD5 hash of the sequence.
     The string array digest must be large enough (len_hex_dig_md5).
     First normalize string by uppercasing it and replacing U's with T's.

     As with get_hex_seq_digest_sha1 above, the hash runs on this private,
     per-call `normalized` copy (never `seq`): the vendored MD5 may read it
     via an aligned reinterpret, and being per-call it is thread-safe. */

  std::vector<char> normalized(static_cast<std::size_t>(seqlen) + 1);
  string_normalize(normalized.data(), seq, static_cast<unsigned int>(seqlen));

  std::vector<unsigned char> digest(md5_digest_length);

  MD5(normalized.data(), static_cast<std::size_t>(seqlen), digest.data());

  for (auto const & element: digest) {
    *hex = hexdigits[element >> drop_lower_nibble];
    std::advance(hex, 1);
    *hex = hexdigits[element & mask_upper_nibble];
    std::advance(hex, 1);
  }
  *hex = '\0';
}


auto fprint_seq_digest_sha1(std::FILE * output_handle, char const * seq, int const seqlen) -> void
{
  std::vector<char> hex_digest(len_hex_dig_sha1);
  get_hex_seq_digest_sha1(hex_digest.data(), seq, seqlen);
  std::fprintf(output_handle, "%s", hex_digest.data());
}


auto fprint_seq_digest_md5(std::FILE * output_handle, char const * seq, int const seqlen) -> void
{
  std::vector<char> hex_digest(len_hex_dig_md5);
  get_hex_seq_digest_md5(hex_digest.data(), seq, seqlen);
  std::fprintf(output_handle, "%s", hex_digest.data());
}
