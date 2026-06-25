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

#pragma once

#include <cstdint> // uint64_t
#include <cstdio>  // std::FILE, std::size_t
#include <cstddef> // std::size_t
#include <limits>  // std::numeric_limits
#include <random>  // std::mt19937_64
#include <utility> // std::swap
#include "utils/fatal.hpp" // fatal() (used by the templates below)


constexpr auto md5_digest_length = 16;
constexpr auto sha1_digest_length = 20;
constexpr auto len_hex_dig_md5 = (2 * md5_digest_length) + 1;
constexpr auto len_hex_dig_sha1 = (2 * sha1_digest_length) + 1;

auto xstrdup(char const * src) -> char *;
auto xsprintf(char * * ret, char const * format, ...) -> int;
auto hash_cityhash64(char const * sequence, uint64_t length) -> uint64_t;
auto hash_cityhash128(char const * sequence, uint64_t length) -> uint128;
auto show_rusage() -> void;

auto progress_init(char const * prompt, uint64_t size) -> void;
auto progress_update(uint64_t progress) -> void;
auto progress_done() -> void;

auto random_init() -> void;

/* ---- Cross-platform reproducible RNG ------------------------------------

   All of the pieces below produce bit-identical results on every platform
   for a given seed: SplitMix64 is fixed integer arithmetic, std::mt19937_64
   has a standard-specified sequence and seeding, and random_bounded() /
   random_shuffle() are implemented here rather than via
   std::uniform_int_distribution / std::shuffle (both implementation-defined,
   so not portable). This is what makes --randseed reproducible across builds
   and operating systems. */

/* SplitMix64: tiny 64-bit PRNG, one word of state. Models a
   UniformRandomBitGenerator so it works with the helpers below exactly like
   std::mt19937_64. It is cheap to (re)seed, which is why sintax uses one per
   query (seeded from the base seed and the query number) to make --randseed
   reproducible independently of the number of threads. */
class SplitMix64 {
public:
  using result_type = uint64_t;
  explicit SplitMix64(uint64_t seed_value) : state_(seed_value) {}
  auto seed(uint64_t seed_value) -> void { state_ = seed_value; }
  auto operator()() -> uint64_t;
  static constexpr auto min() -> uint64_t { return 0; }
  static constexpr auto max() -> uint64_t { return std::numeric_limits<uint64_t>::max(); }
private:
  uint64_t state_;
};

/* The process-wide 64-bit base seed, established once by random_init():
   opt_randseed if non-zero (full 64 bits, no truncation), otherwise a value
   from the operating system. */
auto random_base_seed() -> uint64_t;

/* Derive a well-separated sub-stream seed from a base seed and an index
   (e.g. base = random_base_seed(), index = query number). */
auto random_substream_seed(uint64_t base, uint64_t index) -> uint64_t;

/* Unbiased integer in [0, range), range > 0 (fatal otherwise). Lemire's
   multiply-shift method with rejection; falls back to rejection on the
   remainder where a 128-bit product is unavailable. Templated on any
   UniformRandomBitGenerator (SplitMix64 or std::mt19937_64) whose output
   spans the full 64 bits. */
template <typename URBG>
auto random_bounded(URBG & generator, uint64_t range) -> uint64_t
{
  if (range == 0) { fatal("Internal error: random_bounded() called with range 0"); }
#ifdef __SIZEOF_INT128__
  __uint128_t product = static_cast<__uint128_t>(generator()) * range;
  auto low = static_cast<uint64_t>(product);
  if (low < range)
    {
      uint64_t const threshold = (-range) % range;  /* == 2^64 mod range */
      while (low < threshold)
        {
          product = static_cast<__uint128_t>(generator()) * range;
          low = static_cast<uint64_t>(product);
        }
    }
  return static_cast<uint64_t>(product >> 64U);
#else
  /* portable fallback: rejection sampling to avoid modulo bias */
  uint64_t const limit = std::numeric_limits<uint64_t>::max()
                       - (std::numeric_limits<uint64_t>::max() % range);
  uint64_t value = generator();
  while (value >= limit) { value = generator(); }
  return value % range;
#endif
}

/* Portable in-place Fisher-Yates shuffle over data[0 .. count), replacing
   std::shuffle (which is implementation-defined). */
template <typename T, typename URBG>
auto random_shuffle(T * data, std::size_t count, URBG & generator) -> void
{
  for (std::size_t i = count; i > 1; --i)
    {
      auto const j = static_cast<std::size_t>(random_bounded(generator, i));  /* [0, i) */
      std::swap(data[i - 1], data[j]);
    }
}

auto string_normalize(char * normalized, char const * raw_seq, unsigned int len) -> void;

auto reverse_complement(char * rc_seq, char const * seq, int64_t len) -> void;

auto get_hex_seq_digest_sha1(char * hex, char const * seq, int seqlen) -> void;
auto get_hex_seq_digest_md5(char * hex, char const * seq, int seqlen) -> void;

auto fprint_seq_digest_sha1(std::FILE * output_handle, char const * seq, int seqlen) -> void;
auto fprint_seq_digest_md5(std::FILE * output_handle, char const * seq, int seqlen) -> void;

auto fopen_output(char const * filename) -> std::FILE *;
