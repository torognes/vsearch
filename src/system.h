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

#include <cstdint>  // uint64_t
#include <cstdio>   // std::FILE, std::size_t


#ifdef _WIN32
using xstat_t = struct __stat64;
#else
using xstat_t = struct stat;
#endif

auto system_get_memused() -> uint64_t;
auto system_get_memtotal() -> uint64_t;
auto system_get_cores() -> long;
auto xmalloc(std::size_t size) -> void *;
auto xrealloc(void * ptr, std::size_t size) -> void *;
auto xfree(void * ptr) -> void;

auto xfstat(int file_descriptor, xstat_t * buf) -> int;
auto xstat(const char * path, xstat_t  * buf) -> int;
auto xlseek(int file_descriptor, uint64_t offset, int whence) -> uint64_t;
auto xftello(std::FILE * stream) -> uint64_t;

auto xopen_read(const char * path) -> int;
auto xopen_write(const char * path) -> int;
auto xclose(int file_descriptor) -> int;


/* RAII owner of a raw file descriptor (as returned by xopen_read/xopen_write):
   closes it with xclose() on destruction, so the fd-based read/write/lseek I/O
   in udb.cc no longer relies on manual close() discipline (the leak and
   double-close bug class the FILE * handles removed, applied to the descriptor
   side). Move-only; a negative descriptor is the empty sentinel. A close error
   is not surfaced here, matching the plain close of an input FILE *; the one
   caller that must check it (the UDB writer, where a failed close can mean a
   truncated file) release()es the descriptor and xclose()es it explicitly. */
class FileDescriptor {
public:
  FileDescriptor() noexcept = default;
  explicit FileDescriptor(int descriptor) noexcept : file_descriptor(descriptor) {}

  FileDescriptor(FileDescriptor const &) = delete;
  auto operator=(FileDescriptor const &) -> FileDescriptor & = delete;

  FileDescriptor(FileDescriptor && other) noexcept
    : file_descriptor(other.file_descriptor) {
    other.file_descriptor = -1;
  }
  auto operator=(FileDescriptor && other) noexcept -> FileDescriptor & {
    if (this != &other) {
      reset();
      file_descriptor = other.file_descriptor;
      other.file_descriptor = -1;
    }
    return *this;
  }

  ~FileDescriptor() { reset(); }

  auto get() const noexcept -> int { return file_descriptor; }

  /* relinquish ownership: hand back the descriptor and forget it (no close), so
     the caller can run a checked xclose() itself */
  auto release() noexcept -> int {
    int const released = file_descriptor;
    file_descriptor = -1;
    return released;
  }

  /* close now (a no-op on an empty handle); the destructor is the fallback */
  auto reset() noexcept -> void {
    if (file_descriptor >= 0) {
      static_cast<void>(xclose(file_descriptor));
    }
    file_descriptor = -1;
  }

private:
  int file_descriptor = -1;
};
