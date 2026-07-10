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

#include <cstdio>  // std::FILE, std:fclose, std::fopen
#include <memory>  // std::unique_ptr


// Make sure std::FILE * cannot leak

// for reference, see:
// https://stackoverflow.com/questions/26360916/using-custom-deleter-with-unique-ptr

// note: taking a reference to a standard library function (such as
// &std::fclose) is unspecified behavior (or ill-formed):
//
// using FileHandle = std::unique_ptr<std::FILE, decltype(&std::fclose)>;  UB!

// prefer a deleter struct with an operator() to call std::fclose
struct CloseFileHandle {
  auto operator()(std::FILE * file_handle) noexcept -> void {
    static_cast<void>(std::fclose(file_handle));
  }
};

using FileHandle = std::unique_ptr<std::FILE, CloseFileHandle>;


// Output streams get a checked close: a deferred write error (full disk, quota
// exceeded, broken pipe) is surfaced as a fatal error rather than left as a
// silently truncated output file. Input streams keep the plain CloseFileHandle
// above (a stale read-error flag must not turn into a write failure).
struct CheckedCloseOutputHandle {
  auto operator()(std::FILE * file_handle) noexcept -> void;  // defined in open_file.cpp
};

using OutputFileHandle = std::unique_ptr<std::FILE, CheckedCloseOutputHandle>;


// A command-line option name (for example "--output"), wrapped in its own type
// so it cannot be transposed with the filename argument of the output openers
// (both would otherwise be char const *).
struct OutputOption {
  explicit constexpr OutputOption(char const * str) noexcept
    : name(str) {
  }
  char const * name;
};


/* Open a named input stream in binary read mode. A null filename yields an
   empty handle; "-" reads a duplicate of stdin. On failure the handle is empty
   (the caller checks it). */
auto open_input_file(char const * filename) -> FileHandle;

/* Low-level named output opener in binary write mode: a null filename yields an
   empty handle, "-" writes a duplicate of stdout, and an open failure yields an
   empty handle. Prefer the mandatory/optional wrappers below, which also report
   the error; this primitive is used where the caller does its own checking. */
auto open_output_file(char const * filename) -> OutputFileHandle;

/* Open a mandatory named output stream (binary write mode). Fatal if the option
   was not given (filename == nullptr), naming <option>, or if the file cannot
   be opened. "-" writes a duplicate of stdout. */
auto open_mandatory_output_file(char const * filename, OutputOption option) -> OutputFileHandle;

/* Open an optional named output stream (binary write mode). A null filename
   (option not given) yields an empty handle; a named file that cannot be opened
   is fatal, naming <option>. "-" writes a duplicate of stdout. */
auto open_optional_output_file(char const * filename, OutputOption option) -> OutputFileHandle;
