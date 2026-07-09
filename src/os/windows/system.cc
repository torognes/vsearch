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

#include <windows.h>  // SYSTEM_INFO, GetSystemInfo
#include <algorithm>  // std::max
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, _ftelli64
#include <fcntl.h>  // _O_RDONLY, _O_WRONLY, _O_CREAT, _O_TRUNC, _O_BINARY
#include <io.h>  // _open, _lseeki64
#include <malloc.h>  // _aligned_malloc, _aligned_realloc, _aligned_free
#include <sys/stat.h>  // _fstat64, _stat64, _S_IREAD, _S_IWRITE
// <sys/stat.h> must precede "system.h": it defines the "#define __stat64
// _stat64" alias (via _mingw_stat64.h) that system.h's xstat_t = struct
// __stat64 relies on, so that xstat_t resolves to the real struct _stat64
// expected by _fstat64/_stat64 rather than a distinct forward-declared type.
#include "system.h"
#include "utils/fatal.hpp"

/* system_get_memused()/system_get_memtotal() are genuinely per-OS and live in
   the os/<os>/system_memory.cc backends (mirrors swarm); everything below is
   Win32-specific but not memory-related. */


constexpr auto vsearch_memalignment = 16;


auto system_get_cores() -> long
{
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwNumberOfProcessors;
}


auto xmalloc(std::size_t size) -> void *
{
  static constexpr auto minimal_allocation = std::size_t{1};
  size = std::max(size, minimal_allocation);
  void * ptr = _aligned_malloc(size, vsearch_memalignment);
  if (ptr == nullptr)
    {
      fatal("Unable to allocate enough memory.");
    }
  return ptr;
}


auto xrealloc(void * ptr, std::size_t size) -> void *
{
  /* NOTE: _aligned_realloc preserves the vsearch_memalignment (16-byte)
     alignment that _aligned_malloc gave the block, so this Windows branch,
     unlike the POSIX xrealloc (plain realloc, only max_align_t), is safe for
     the 16-byte-aligned SIMD buffers (the align_simd.cc VECTOR_SHORT arrays).
     The portable contract nevertheless still holds: those buffers are
     allocated with xmalloc and, as audited, never grown through xrealloc —
     every caller resizes byte/scalar data (input buffers, sequence storage,
     query strings, k-mer lists). Keep honouring it so the POSIX build, which
     cannot rely on this, stays correct. */
  static constexpr auto minimal_allocation = std::size_t{1};
  size = std::max(size, minimal_allocation);
  void * new_ptr = _aligned_realloc(ptr, size, vsearch_memalignment);
  if (new_ptr == nullptr)
    {
      fatal("Unable to reallocate enough memory.");
    }
  return new_ptr;
}


auto xfree(void * ptr) -> void
{
  if (ptr != nullptr)
    {
      _aligned_free(ptr);
    }
  else
    {
      fatal("Trying to free a null pointer");
    }
}


auto xfstat(int file_descriptor, xstat_t * buf) -> int
{
  return _fstat64(file_descriptor, buf);
}


auto xstat(const char * path, xstat_t * buf) -> int
{
  return _stat64(path, buf);
}


auto xlseek(int file_descriptor, uint64_t offset, int whence) -> uint64_t
{
  return _lseeki64(file_descriptor, offset, whence);
}


// refactoring: only used in fastx.cc
auto xftello(std::FILE * stream) -> uint64_t
{
  return _ftelli64(stream);
}


// refactoring: only used in udb.cc
auto xopen_read(const char * path) -> int
{
  return _open(path, _O_RDONLY | _O_BINARY);
}


// refactoring: only used in udb.cc
auto xopen_write(const char * path) -> int
{
  return _open(path,
               _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY,
               _S_IREAD | _S_IWRITE);
}
