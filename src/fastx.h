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

#include <array>
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t


constexpr auto byte_range = 256U;

struct fastx_buffer_s
{
  char * data = nullptr;
  uint64_t length = 0;
  uint64_t alloc = 0;
  uint64_t position = 0;
};

auto buffer_init(struct fastx_buffer_s * buffer) -> void;
auto buffer_free(struct fastx_buffer_s * buffer) -> void;
auto buffer_extend(struct fastx_buffer_s * dest_buffer,
                   char * source_buf,
                   uint64_t len) -> void;
auto buffer_makespace(struct fastx_buffer_s * buffer, uint64_t size) -> void;

struct fastx_s
{
  bool is_pipe = false;
  bool is_fastq = false;
  bool is_empty = false;

  std::FILE * fp = nullptr;

#ifdef HAVE_ZLIB_H
  gzFile fp_gz = nullptr;
#endif

#ifdef HAVE_BZLIB_H
  BZFILE * fp_bz = nullptr;
#endif

  struct fastx_buffer_s file_buffer;

  struct fastx_buffer_s header_buffer;
  struct fastx_buffer_s sequence_buffer;
  struct fastx_buffer_s plusline_buffer;
  struct fastx_buffer_s quality_buffer;

  uint64_t file_size = 0;
  uint64_t file_position = 0;

  uint64_t lineno = 0;
  uint64_t lineno_start = 0;
  int64_t seqno = 0;

  uint64_t stripped_all = 0;
  std::array<uint64_t, byte_range> stripped {{}};

  int format = 0;
};

using fastx_handle = struct fastx_s *;


/* fastx input */

auto fastx_is_fastq(fastx_handle input_handle) -> bool;
auto fastx_is_empty(fastx_handle input_handle) -> bool;
auto fastx_is_pipe(fastx_handle input_handle) -> bool;
auto fastx_filter_header(fastx_handle input_handle, bool truncateatspace) -> void;
auto fastx_open(const char * filename) -> fastx_handle;
auto fastx_close(fastx_handle input_handle) -> void;
auto fastx_next(fastx_handle input_handle,
                bool truncateatspace,
                const unsigned char * char_mapping) -> bool;
auto fastx_get_position(fastx_handle input_handle) -> uint64_t;
auto fastx_get_size(fastx_handle input_handle) -> uint64_t;
auto fastx_get_lineno(fastx_handle input_handle) -> uint64_t;
auto fastx_get_seqno(fastx_handle input_handle) -> uint64_t;
auto fastx_get_header(fastx_handle input_handle) -> char *;
auto fastx_get_sequence(fastx_handle input_handle) -> char *;
auto fastx_get_header_length(fastx_handle input_handle) -> uint64_t;
auto fastx_get_sequence_length(fastx_handle input_handle) -> uint64_t;

auto fastx_get_quality(fastx_handle input_handle) -> char *;
auto fastx_get_abundance(fastx_handle input_handle) -> int64_t;

auto fastx_file_fill_buffer(fastx_handle input_handle) -> uint64_t;
