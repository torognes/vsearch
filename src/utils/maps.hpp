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


#include <cassert>  // assert


/* Which of the two character maps a parser applies to the bases it accepts.
   Passed as a template argument rather than as the table pointer it stands
   for, because every call site knows the answer at compile time: 22 of the 30
   fastx_s::next() calls name chrmap_no_change(), 7 name chrmap_upcase(), and
   Database::read() chooses between them from its own literal argument. */
enum struct Mapping : unsigned char {
  none,    /* the accepted byte is written through unchanged */
  upcase,  /* ... after ASCII case folding */
};


/* Map one byte that a fastx parser has already classified as Action::accept.
   Both cases replace a 256-byte table lookup, and are equivalences rather
   than approximations because of what can reach them:

   - the FASTA and FASTQ sequence accept sets are the same 32 IUPAC letters
     (ABCDGHKMNRSTUVWY and their lowercase, verified byte-for-byte), and on
     that set chrmap_no_change() is the identity and chrmap_upcase() is plain
     ASCII case folding;
   - the FASTQ quality accept set is bytes 33 to 126, and the map it used was
     the identity on all 256 values, so Mapping::none covers it too.

   Only Mapping::upcase constrains its input, and the assert states that. */
template <Mapping mapping>
inline auto map_accepted_base(char const base) -> char {
  /* refactoring (C++17): `if constexpr`. The plain `if` already costs nothing
     -- mapping is part of the type, so the dead half is folded away and both
     spellings emit the same instructions -- but `if constexpr` would also stop
     type-checking the dead branch, which would matter if a case ever held code
     valid for one mapping alone. */
  if (mapping == Mapping::none) {
    return base;
  }
  assert(((base >= 'A') and (base <= 'Z')) or ((base >= 'a') and (base <= 'z')));
  /* clearing the 'a' - 'A' bit upper-cases an ASCII letter and leaves an
     already upper-case one alone; unsigned throughout, because char is signed
     on x86-64 and Windows and a bitwise operator on a signed operand is a
     trap even where the value cannot reach the sign bit */
  static constexpr auto lowercase_bit = 0x20U;
  return static_cast<char>(static_cast<unsigned char>(base) & ~lowercase_bit);
}


auto chrmap_no_change() -> unsigned char const *;

auto chrmap_normalize() -> unsigned char const *;

auto chrmap_upcase() -> unsigned char const *;

auto chrmap_complement() -> unsigned char const *;

auto chrmap_2bit() -> unsigned int const *;

auto chrmap_4bit() -> unsigned char const *;

auto chrmap_mask_ambig() -> unsigned int const *;

auto chrmap_mask_lower() -> unsigned int const *;

auto map_uppercase(char nucleotide) -> char;

auto map_2bit(char nucleotide) -> unsigned int;

auto map_4bit(char nucleotide) -> unsigned char;

auto map_complement(char nucleotide) -> char;

auto map_mask_ambig(char nucleotide) -> unsigned int;

auto map_mask_lower(char nucleotide) -> unsigned int;

auto is_equivalent_4bit(char lhs, char rhs) -> bool;

auto is_equivalent_4bit_rhs(char lhs, char rhs) -> bool;

auto is_ambiguous_4bit(unsigned char nucleotide) -> bool;

auto is_same_4bit(char lhs, char rhs) -> bool;
