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

#include <cassert>
#include <cstdint>  // uint64_t


/* Sizing arithmetic for the hand-rolled open-addressing tables.

   Five copies of two expressions lived at the call sites before this header:
   the 2/3-fill-rate loop in commands/derep_prefix.cpp and core/dbhash.cpp,
   character for character, and the load-factor-1/2 loop in core/unique.cpp and
   twice in core/kmerhash.cpp. They are not worth a helper for their length --
   they are worth one because the reasoning that protects them was recorded at
   some copies and not others, exactly the shape of the dbindex bitmap mincount
   bug, which had to be fixed in two places because its sizing arithmetic was
   likewise duplicated.

   Why two named functions rather than one taking a fill rate: a signature of
   the form (elements, numerator, denominator) is three adjacent parameters of
   the same type, so transposing two of them still compiles and still returns a
   plausible size. Naming the rate removes the failure mode entirely and reads
   better at the call site.

   The width is load-bearing, and the reasoning recorded at core/unique.cpp and
   again on kh_handle_s (utils/kmer_hash_struct.hpp) applies to every caller:
   "the hash grows to 2 * sequence length, which exceeds INT_MAX for sequences
   above ~1.07 Gnt (the doubling would otherwise overflow int before reaching
   the target)". Taking and returning uint64_t also makes the intermediate
   multiplication wrap rather than overflow a signed type, which is what
   3 * dbsequencecount did.

   Neither function allocates or throws, hence noexcept.

   // C++14 refactoring: both can become constexpr -- a C++11 constexpr
   // function body is limited to a single return statement, so expressing
   // either loop that way today would mean writing it recursively for no gain.

   What is deliberately not here: commands/derep_smallmem.cpp sizes its table
   by x1.5 growth to a non-power-of-two and steps it with a modulo, so none of
   this applies to it. See TBD_20260825_flatmap_helpers.md. */
namespace vsearch
{

  /* Smallest power of two that holds 'elements' entries at a fill rate of at
     most two thirds -- the persistent tables, sized once from a known count
     and never grown (commands/derep_prefix.cpp, core/dbhash.cpp).

     Returns 1 rather than 0 for an empty table, which is what the loops this
     replaced did: they started at 1 and simply never ran. Callers derive a
     mask as (result - 1), so a zero would be a wraparound. */
  inline auto table_size_two_thirds(uint64_t const elements) noexcept -> uint64_t
  {
    uint64_t size = 1;
    while (3 * elements > 2 * size)
      {
        size <<= 1U;
      }
    assert((size & (size - 1)) == 0);  /* the mask below depends on it */
    return size;
  }


  /* Smallest power of two that holds 'elements' entries at a fill rate of at
     most one half -- the per-sequence scratch tables, sized to the sequence so
     that the load factor is picked by construction and no rehash is ever
     needed (core/unique.cpp, core/kmerhash.cpp).

     Same zero case, and for the same reason, as above. */
  inline auto table_size_half(uint64_t const elements) noexcept -> uint64_t
  {
    uint64_t size = 1;
    while (size < 2 * elements)
      {
        size <<= 1U;
      }
    assert((size & (size - 1)) == 0);  /* the mask below depends on it */
    return size;
  }

}  // namespace vsearch
