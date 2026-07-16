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

#include "fatal.hpp"
#include "open_file.hpp"
#include <unistd.h>  // dup, STDIN_FILENO, STDOUT_FILENO
#include <cassert>
#include <cerrno>  // errno
#include <cstdio>  // std::fopen, fdopen
#include <cstring>  // std::strcmp
#include <string>  // std::string (building the open-failure message)


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


  auto is_dash(char const * filename) -> bool {
    assert(filename != nullptr);
    return std::strcmp(filename, "-") == 0;
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


  // Open a stdio stream and return a raw FILE * for the caller to wrap in the
  // matching RAII handle. A filename of '-' is served by a duplicate of the
  // given standard stream (STDIN_FILENO for input, STDOUT_FILENO for output).
  auto open_stream(char const * filename,
                   ModeString const & mode,
                   int const standard_fileno) -> std::FILE * {
    assert(filename != nullptr);
    assert(mode.mode != nullptr);
    if (is_dash(filename)) {
      auto const file_descriptor = dup(standard_fileno);
      check_file_descriptor(file_descriptor);
      return fdopen(file_descriptor, mode.mode);
    }
    return std::fopen(filename, mode.mode);
  }


  // Both output openers report an open failure the same way; only the
  // null-filename policy differs between them.
  auto fatal_output_open_failed(char const * filename,
                                OutputOption const option) -> void {
    assert(filename != nullptr);
    std::string const message = std::string("unable to open output file for ")
      + option.name + " (" + filename + ")";
    fatal(message.c_str());
  }

}  // end of anonymous namespace


// read_file, file to read, open_input_file, open_istream
auto open_input_file(char const * filename) -> FileHandle {
  if (filename == nullptr) {
    return FileHandle{nullptr};
  }
  auto const mode = ModeString{"rb"};  // r: reading, b: non-UNIX environments
  /* open the input stream given by filename, but if name is '-' then
     use a duplicate of stdin (fd = STDIN_FILENO = 0) */
  return FileHandle{open_stream(filename, mode, STDIN_FILENO)};
}


auto CheckedCloseOutputHandle::operator()(std::FILE * file_handle) noexcept -> void {
  if (file_handle == nullptr) {
    return;
  }
  /* A write error (full disk, quota, broken pipe) is often deferred by stdio
     until the buffer is flushed, so check fflush and the error flag before
     closing; fclose also flushes and can report the same error. Fail loudly
     rather than leave a silently truncated output file. */
  if ((std::fflush(file_handle) != 0) or (std::ferror(file_handle) != 0)) {
    fatal("Unable to write to output file (disk full, quota exceeded, or broken pipe?)");
  }
  if (std::fclose(file_handle) != 0) {
    fatal("Unable to close output file (disk full or quota exceeded?)");
  }
}


// write_file, file to write, open_output_file, open_ostream
auto open_output_file(char const * filename) -> OutputFileHandle {
  if (filename == nullptr) {
    return OutputFileHandle{nullptr};
  }
  auto const mode = ModeString{"wb"};  // w: writing, b: binary (no \n->\r\n on non-UNIX), matches input "rb"
  /* open the output stream given by filename, but if name is '-' then
     use a duplicate of stdout (fd = STDOUT_FILENO = 1) */
  return OutputFileHandle{open_stream(filename, mode, STDOUT_FILENO)};
}


auto open_mandatory_output_file(char const * filename,
                                OutputOption const option) -> OutputFileHandle {
  if (filename == nullptr) {
    fatal("output file must be specified with %s", option.name);
  }
  auto output_handle = open_output_file(filename);
  if (not output_handle) {
    fatal_output_open_failed(filename, option);
  }
  return output_handle;
}


auto open_optional_output_file(char const * filename,
                               OutputOption const option) -> OutputFileHandle {
  if (filename == nullptr) {
    return OutputFileHandle{nullptr};
  }
  auto output_handle = open_output_file(filename);
  if (not output_handle) {
    fatal_output_open_failed(filename, option);
  }
  return output_handle;
}
