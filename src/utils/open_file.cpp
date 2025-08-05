/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "fatal.hpp"
#include "open_file.hpp"
#include <unistd.h>  // dup, STDIN_FILENO, STDOUT_FILENO
#include <cassert>
#include <cerrno>  // errno
#include <cstdio>  // std::fopen, fdopen
#include <cstring>  // std::strcmp


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // C++17 refactoring:
  // constexpr std::string_view a_dash = "-";
  // simpler string comparisons: if (filename == a_dash) {

  // type-safe string mode wrapper
  struct ModeString {
    explicit constexpr ModeString(char const * str) noexcept
      : mode(str) {
    }
    char const * mode;
  };


  // Safely wrapping fopen()
  auto open_file(char const * filename,
                 ModeString const & mode) -> FileHandle {
    assert(filename != nullptr);
    assert(mode.mode != nullptr);
    return FileHandle{std::fopen(filename, mode.mode)};
  }


  auto open_file_descriptor(int const file_descriptor,
                            ModeString const & mode) -> FileHandle {
    assert(file_descriptor >= 0);
    assert(mode.mode != nullptr);
    return FileHandle{fdopen(file_descriptor, mode.mode)};
  }


  auto check_file_descriptor(int const file_descriptor) -> void {
    assert(file_descriptor >= -1);
    if (file_descriptor != -1) {
      return;
    }
    if (errno == EBADF) {
      fatal("original fd is not an open file descriptor.");
    }
    if (errno == EMFILE) {
      fatal("too many open file descriptors.");
    }
    fatal("cannot duplicate input or output stream.");
  }

}  // end of anonymous namespace


// read_file, file to read, open_input_file, open_istream
auto open_input_file(char const * filename) -> FileHandle {
  if (filename == nullptr) {
    std::FILE * empty = nullptr;
    return FileHandle{empty};
  }
  auto const mode = ModeString{"rb"};  // r: reading, b: non-UNIX environments
  /* open the input stream given by filename, but if name is '-' then
     use a duplicate of stdin (fd = STDIN_FILENO = 0) */
  if (std::strcmp(filename, "-") == 0) {
    auto const file_descriptor = dup(STDIN_FILENO);
    check_file_descriptor(file_descriptor);
    return open_file_descriptor(file_descriptor, mode);
  }
  return open_file(filename, mode);
}


// write_file, file to write, open_output_file, open_ostream
auto open_output_file(char const * filename) -> FileHandle {
  if (filename == nullptr) {
    std::FILE * empty = nullptr;
    return FileHandle{empty};
  }
  auto const mode = ModeString{"w"};  // w: writing, no b?
  /* open the output stream given by filename, but if name is '-' then
     use a duplicate of stdout (fd = STDOUT_FILENO = 1) */
  if (std::strcmp(filename, "-") == 0) {
    auto const file_descriptor = dup(STDOUT_FILENO);
    check_file_descriptor(file_descriptor);
    return open_file_descriptor(file_descriptor, mode);
  }
  return open_file(filename, mode);
}
