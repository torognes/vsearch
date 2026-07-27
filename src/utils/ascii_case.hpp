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


#include <cctype>  // std::isalnum, std::isdigit, std::isupper, std::tolower, std::toupper


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

inline auto to_upper(char const character) -> char {
  return static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
}

inline auto to_lower(char const character) -> char {
  return static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
}

inline auto is_alnum(char const character) -> bool {
  return std::isalnum(static_cast<unsigned char>(character)) != 0;
}

inline auto is_digit(char const character) -> bool {
  return std::isdigit(static_cast<unsigned char>(character)) != 0;
}

inline auto is_upper(char const character) -> bool {
  return std::isupper(static_cast<unsigned char>(character)) != 0;
}
