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

#include "core/mask.hpp"  // Masking
#include "utils/fatal_allocator.hpp"  // FatalAllocator
#include <cstdint>  // uint64_t
#include <vector>  // std::vector

/* Finds the unique k-mers (words) in a sequence, reusing its scratch buffers
   across calls (see unique.cpp for the bitmap and hash variants). Owns those
   buffers via std::vector (RAII); one instance per thread. */
class Uniquer
{
public:
  auto count(int wordlength,
             int seqlen,
             char const * seq,
             unsigned int * listlen,
             unsigned int const * * list,
             Masking seqmask) -> void;

  // noexcept: reads the already-built hash/bitmap only (no allocation), so it
  // cannot fatal()/throw, unlike count() which grows its buffers.
  auto count_shared(int wordlength,
                    int listlen,
                    unsigned int const * list) const noexcept -> unsigned int;

private:
  struct bucket
  {
    unsigned int kmer;
    unsigned int count;
  };

  auto count_bitmap(int wordlength,
                    int seqlen,
                    char const * seq,
                    unsigned int * listlen,
                    unsigned int const * * list,
                    Masking seqmask) -> void;

  auto count_hash(int wordlength,
                  int seqlen,
                  char const * seq,
                  unsigned int * listlen,
                  unsigned int const * * list,
                  Masking seqmask) -> void;

  std::vector<bucket, FatalAllocator<bucket>> hash_ {};
  std::vector<unsigned int, FatalAllocator<unsigned int>> list_ {};
  std::vector<uint64_t, FatalAllocator<uint64_t>> bitmap_ {};
  unsigned int hash_mask_ {0};
};
