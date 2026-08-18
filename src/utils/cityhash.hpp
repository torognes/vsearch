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

#include "utils/view.hpp"
#include "vendored/city.h"  // CityHash64, CityHash128, uint128
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t


/* The project's hashing entry points. Defined inline rather than in a
   translation unit of their own: each is a single forwarding call, they sit
   in the inner loop of dereplication, of the k-mer index and of read
   merging, and vsearch is not built with link-time optimisation -- an
   out-of-line one-line forwarder would add a call per hash that the compiler
   could not see through. */

inline auto hash_cityhash64(View<char> const sequence) -> uint64_t
{
  return CityHash64(sequence.data(), sequence.size());
}


inline auto hash_cityhash128(View<char> const sequence) -> uint128
{
  return CityHash128(sequence.data(), sequence.size());
}


/* Hash a 2-bit-packed k-mer held in an unsigned int, reading only the bytes
   its 2 * kmer_length significant bits occupy -- ceil(2 * kmer_length / 8),
   which is (kmer_length + 3) / 4.

   Reading fewer bytes than the unsigned int holds makes the result depend on
   the host's byte order (see View::as_bytes() in utils/view.hpp, which views
   the object representation): on a little-endian host these are the bytes
   that carry the significant bits, on a big-endian one they would be the
   zero padding, and every k-mer of a given length would hash alike. Every
   target vsearch supports is little-endian, and core/udb.cpp asserts as much
   for the UDB format. Both callers (core/kmerhash.cpp, core/unique.cpp)
   probe linearly and compare the full k-mer on the way, so a big-endian port
   would lose the hash's spread, not its correctness.

   The k-mer is taken by reference because the byte view must name a live
   object; the view never escapes this function. */
inline auto hash_packed_kmer(unsigned int const & kmer, int const kmer_length) -> uint64_t
{
  assert(kmer_length > 0);
  auto const significant_bytes = static_cast<std::size_t>((kmer_length + 3) / 4);
  assert(significant_bytes <= sizeof(kmer));
  return hash_cityhash64(View<unsigned int>{&kmer, 1}.as_bytes().first(significant_bytes));
}
