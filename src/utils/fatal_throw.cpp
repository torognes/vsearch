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

#include "fatal.hpp"  // fatal_detail::exit_or_throw, VsearchError, throw_on_fatal
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::snprintf
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <string>  // std::string (VsearchError payload)

// fatal_detail::throw_on_fatal() is defined in vsearch_api.cpp, next to the
// VsearchSession ctor/dtor that are its only mutators. Only this translation
// unit reads it.


namespace {
  // Render a printf-style (format, args...) pair into a std::string, so the
  // thrown VsearchError carries the same text fatal() printed. Two-pass
  // std::snprintf: measure, then format into a right-sized buffer (C++11
  // guarantees std::string storage is contiguous and null-terminated, so
  // writing length+1 bytes into &text[0] is well defined).
  auto format_message(char const * format, char const * argument) -> std::string {
    int const length = std::snprintf(nullptr, 0, format, argument);
    if (length <= 0) { return std::string(); }
    std::string text(static_cast<std::size_t>(length), '\0');
    std::snprintf(&text[0], static_cast<std::size_t>(length) + 1, format, argument);
    return text;
  }

  auto format_message(char const * format, char const symbol,
                      uint64_t const line_number) -> std::string {
    int const length = std::snprintf(nullptr, 0, format, symbol, line_number);
    if (length <= 0) { return std::string(); }
    std::string text(static_cast<std::size_t>(length), '\0');
    std::snprintf(&text[0], static_cast<std::size_t>(length) + 1, format, symbol, line_number);
    return text;
  }
}


// Library build (-fexceptions): when a session is active, unwind back to the
// consumer (who can catch VsearchError, skip the bad input and continue)
// instead of killing the process; with no session active it exits exactly as
// the CLI does. Selected over fatal_exit.cpp by libvsearch_core's source list
// in Makefile.am, not the preprocessor.
namespace fatal_detail {
  __attribute__((noreturn))
  auto exit_or_throw(char const * message) -> void {
    if (throw_on_fatal()) { throw VsearchError{message}; }
    std::exit(EXIT_FAILURE);
  }

  __attribute__((noreturn))
  auto exit_or_throw(char const * format, char const * message) -> void {
    if (throw_on_fatal()) { throw VsearchError{format_message(format, message)}; }
    std::exit(EXIT_FAILURE);
  }

  __attribute__((noreturn))
  auto exit_or_throw(char const * format, char const symbol,
                     uint64_t const line_number) -> void {
    if (throw_on_fatal()) { throw VsearchError{format_message(format, symbol, line_number)}; }
    std::exit(EXIT_FAILURE);
  }
}
