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

#include "fatal.hpp"  // fatal, fatal_detail::exit_or_throw
#include "logfile.hpp"  // log_file::handle
#include <cstdint> // uint64_t
#include <cstdio>  // std::fprintf, std::fputc, std::fputs


[[noreturn]]
auto fatal(char const * message) -> void {
  std::fputs("\n\n", stderr);
  std::fprintf(stderr, "Fatal error: %s\n", message);

  auto * const log = log_file::handle();
  if (log != nullptr) {
    std::fputs("\n\n", log);
    std::fprintf(log, "Fatal error: %s\n", message);
  }

  fatal_detail::exit_or_throw(message);  // CLI: exit; library: throw when in a session
}


[[noreturn]]
auto fatal(char const * format,
           char const * message) -> void {
  std::fputs("\n\nFatal error: ", stderr);
  std::fprintf(stderr, format, message);
  std::fputc('\n', stderr);

  auto * const log = log_file::handle();
  if (log != nullptr) {
    std::fputs("\n\nFatal error: ", log);
    std::fprintf(log, format, message);
    std::fputc('\n', log);
  }

  fatal_detail::exit_or_throw(format, message);
}


// used in fastx.cc
[[noreturn]]
auto fatal(char const * format,
           char const symbol,
           uint64_t const line_number) -> void {
  std::fputs("\n\nFatal error: ", stderr);
  std::fprintf(stderr, format, symbol, line_number);
  std::fputc('\n', stderr);

  auto * const log = log_file::handle();
  if (log != nullptr) {
    std::fputs("\n\nFatal error: ", log);
    std::fprintf(log, format, symbol, line_number);
    std::fputc('\n', log);
  }

  fatal_detail::exit_or_throw(format, symbol, line_number);
}
