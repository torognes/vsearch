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


#include <cctype>  // std::isalnum, std::isdigit


// Guarded wrappers around the single-character <cctype> functions.
//
// std::toupper and its neighbours take an int whose value must be
// representable as unsigned char, or be EOF. A plain char argument is
// therefore undefined behaviour wherever char is signed -- x86-64 and the
// Windows target, while ARM and PowerPC Linux default char to unsigned
// instead -- and the byte has its high bit set. vsearch reaches that case
// through sequence headers, which are arbitrary bytes and may carry UTF-8 or
// Latin-1 (an accented author name, a non-ASCII locality).
//
// Converting through unsigned char first is the fix. Keeping it here states
// the contract once rather than at every call site, and the is_* wrappers
// return bool so that callers need no "!= 0" comparison.
//
// vsearch never calls std::setlocale, so these run in the "C" locale, where
// they are plain ASCII case folding and classification.

// The three case helpers spell the "C" locale mapping out rather than calling
// libc, which reaches it through __ctype_toupper_loc(): a thread-local load and
// an indirect table lookup, per character, in a function the compiler cannot
// inline or vectorise. dust_core() folds case over every base of every sequence
// and --fastx_mask classifies every base again, which made libc's toupper 1.6 %
// and isupper 1.0 % of a masking run. The arithmetic below is the same mapping
// by the paragraph above: only one case's letters move, everything else is
// identity. is_alnum and is_digit keep calling libc: they are not on a hot
// path, and is_alnum in particular classifies far more than these letters do.
inline auto to_upper(char const character) -> char {
  static constexpr auto letters_in_alphabet = 26U;
  static constexpr auto case_bit = 32U;  // 0x20, the distance between cases
  auto const byte = static_cast<unsigned char>(character);
  auto const is_lower_case = static_cast<unsigned int>(byte - 'a') < letters_in_alphabet;
  return static_cast<char>(byte - (is_lower_case ? case_bit : 0U));
}

inline auto to_lower(char const character) -> char {
  static constexpr auto letters_in_alphabet = 26U;
  static constexpr auto case_bit = 32U;  // 0x20, the distance between cases
  auto const byte = static_cast<unsigned char>(character);
  auto const is_upper_case = static_cast<unsigned int>(byte - 'A') < letters_in_alphabet;
  return static_cast<char>(byte + (is_upper_case ? case_bit : 0U));
}

inline auto is_alnum(char const character) -> bool {
  return std::isalnum(static_cast<unsigned char>(character)) != 0;
}

inline auto is_digit(char const character) -> bool {
  return std::isdigit(static_cast<unsigned char>(character)) != 0;
}

inline auto is_upper(char const character) -> bool {
  static constexpr auto letters_in_alphabet = 26U;
  return static_cast<unsigned int>(static_cast<unsigned char>(character) - 'A') < letters_in_alphabet;
}
