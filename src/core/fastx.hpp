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

#ifdef HAVE_CONFIG_H
#include "config.h"  // HAVE_ZLIB_H, HAVE_BZLIB_H
#endif

#include "utils/view.hpp"  // View
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t
#include <cstring>  // std::memchr
#include <iterator>  // std::next, std::distance

#ifdef HAVE_ZLIB_H
#include <zlib.h>  // gzFile
#endif
#ifdef HAVE_BZLIB_H
#include <bzlib.h>  // BZFILE
#endif


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
                   char const * source_buf,
                   uint64_t len) -> void;
auto buffer_makespace(struct fastx_buffer_s * buffer, uint64_t size) -> void;

enum struct Format : unsigned char { undefined, plain, bzip, gzip };

class DynamicLibraries;  // set from parameters.dyn_libs in fastx_open()

struct fastx_s
{
  bool is_pipe = false;
  bool is_fastq = false;
  bool is_empty = false;

  std::FILE * fp = nullptr;

  /* runtime-loaded compression libraries, borrowed (non-owning) from the
     DynamicLibraries instance in main(); nullptr in library-only builds */
  DynamicLibraries const * libraries = nullptr;

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

  Format format = Format::undefined;

  /* Deferred error reporting (prototype for CC3). When defer_errors is
     set, a parse error records its message here and makes fastx_next()
     return false, instead of calling fatal() (std::exit()) on the spot.
     A command that reads this handle from worker threads enables this so
     the worker can stop cooperatively; the error is then reported from
     the main thread after the worker pool has joined, avoiding a
     std::exit() that races sibling threads. Default false → behavior is
     unchanged (immediate fatal) for every other caller. */
  bool defer_errors = false;
  bool error = false;
  std::array<char, 512> errmsg {{}};
};

using fastx_handle = struct fastx_s *;


/* fastx input */

auto fastx_is_fastq(struct fastx_s const * input_handle) -> bool;
auto fastx_is_empty(struct fastx_s const * input_handle) -> bool;
auto fastx_is_pipe(struct fastx_s const * input_handle) -> bool;
auto fastx_filter_header(fastx_handle input_handle, bool truncateatspace) -> void;
auto fastx_get_error(struct fastx_s const * input_handle) -> bool;
auto fastx_get_errmsg(struct fastx_s const * input_handle) -> char const *;
auto fastx_set_deferred_error(fastx_handle input_handle, char const * message) -> void;
auto fastx_open(const char * filename, struct Parameters const & parameters) -> fastx_handle;
auto fastx_close(fastx_handle input_handle, struct Parameters const & parameters) -> void;
auto fastx_next(fastx_handle input_handle,
                bool truncateatspace,
                const unsigned char * char_mapping) -> bool;
auto fastx_get_position(struct fastx_s const * input_handle) -> uint64_t;
auto fastx_get_size(struct fastx_s const * input_handle) -> uint64_t;
auto fastx_get_lineno(struct fastx_s const * input_handle) -> uint64_t;
auto fastx_get_seqno(struct fastx_s const * input_handle) -> uint64_t;
auto fastx_get_header(struct fastx_s const * input_handle) -> char const *;
auto fastx_get_sequence(struct fastx_s const * input_handle) -> char const *;
auto fastx_get_header_length(struct fastx_s const * input_handle) -> uint64_t;
auto fastx_get_sequence_length(struct fastx_s const * input_handle) -> uint64_t;

// Reject a sequence too long for the int length bookkeeping used downstream.
// Called from fasta_next/fastx_next so every read is bounded at one choke
// point, symmetric with fastx_filter_header. On a worker thread an over-long
// sequence records a deferred error (reported from the main thread), otherwise
// it is fatal.
auto fastx_filter_sequence_length(fastx_handle input_handle) -> void;

auto fastx_get_quality(struct fastx_s const * input_handle) -> char const *;
auto fastx_get_abundance(struct fastx_s const * input_handle) -> int64_t;

auto fastx_file_fill_buffer(fastx_handle input_handle) -> uint64_t;


/* Line-reading primitives shared by the FASTA and FASTQ record parsers.

   A read-only view of one input line sitting in the file buffer:
     view          bytes up to and including the LF when the line is complete
     has_newline   true when the fragment ends at an LF (the line is complete)

   The three fasta/fastq readers all scanned for the next '\n' by hand, each
   repeating the same std::memchr + length + pointer arithmetic. The scan is
   deliberately kept separate from the "is there more input?" check: a caller
   first calls the existing fastx_file_fill_buffer() (0 == end of input) and, in
   the loops that stop at a record-boundary sentinel, tests the first buffered
   byte BEFORE scanning, so no memchr is spent on a line that is about to be
   handed to the next record. scan_line_fragment() then locates the LF and
   consume_fragment() advances the read position once the caller has copied the
   fragment out. Policy that genuinely differs between the loops (raw vs filtered
   copy, lineno accounting, EOF handling, sentinels) stays with the caller. */
struct Line_fragment
{
  View<char> view;
  bool has_newline;
};

// Locate the next LF in the file buffer and return the fragment starting at the
// current read position, WITHOUT refilling. PRECONDITION: fastx_file_fill_buffer()
// has just reported at least one unconsumed byte. Kept inline in the header so
// the hot parser loops in fasta.cpp / fastq.cpp keep their fully-inlined scan (no
// whole-program optimisation is assumed at build time).
inline auto scan_line_fragment(fastx_handle input_handle) -> Line_fragment
{
  auto & file_buffer = input_handle->file_buffer;
  auto const rest = file_buffer.length - file_buffer.position;
  auto * const start = std::next(file_buffer.data,
                                 static_cast<std::ptrdiff_t>(file_buffer.position));
  auto * const line_end = static_cast<char *>(std::memchr(start, '\n', rest));
  auto const has_newline = (line_end != nullptr);
  auto const length = has_newline
    ? static_cast<std::size_t>(std::distance(start, line_end)) + 1
    : static_cast<std::size_t>(rest);
  return Line_fragment{View<char>{start, length}, has_newline};
}

// Advance the file-buffer read position past a fragment already copied out.
inline auto consume_fragment(fastx_handle input_handle,
                             Line_fragment const & fragment) -> void
{
  input_handle->file_buffer.position += fragment.view.size();
}
