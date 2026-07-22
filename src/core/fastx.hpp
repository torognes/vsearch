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

#include "core/seq_record.hpp"  // SeqRecord (returned by fastx_record)
#include "utils/fatal_allocator.hpp"  // FatalAllocator
#include "utils/view.hpp"  // View
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t
#include <iterator>  // std::next, std::distance
#include <memory>  // std::unique_ptr
#include <string>  // std::char_traits
#include <vector>


constexpr auto byte_range = 256U;

/* One growable byte buffer used by the streaming FASTA/FASTQ reader: the file
   refill buffer and the four per-record buffers (header, sequence, plusline,
   quality). It owns its storage through a std::vector using FatalAllocator, so
   an out-of-memory condition ends the program through fatal() -- exactly as the
   former xmalloc/xrealloc buffers did -- and the storage is released
   automatically when the owning fastx_s is destroyed, so no explicit free is
   needed (RAII). This mirrors how Database owns its packed buffers.

   'length' is the logical number of bytes in use (extend() always keeps them
   NUL-terminated); 'position' is the read cursor into the file refill buffer
   (the per-record buffers leave it at zero). Both are public because the parser
   hot loops in fasta.cpp / fastq.cpp advance and test them directly. data()
   hands out the owned storage and alloc() its current capacity. */
class FastxBuffer
{
public:
  uint64_t length = 0;
  uint64_t position = 0;

  auto data() noexcept -> char * { return storage_.data(); }
  auto data() const noexcept -> char const * { return storage_.data(); }
  auto alloc() const noexcept -> uint64_t { return storage_.size(); }

  /* Reset to a single empty, NUL-terminated block (former buffer_init). */
  auto init() -> void;
  /* Ensure at least 'size' more bytes fit after 'length', rounding the
     allocation up to the nearest block (former buffer_makespace). */
  auto makespace(uint64_t size) -> void;
  /* Append 'len' bytes from 'source', then a terminating NUL (former
     buffer_extend). */
  auto extend(char const * source, uint64_t len) -> void;

private:
  std::vector<char, FatalAllocator<char>> storage_;
};

enum struct Format : unsigned char { undefined, plain, bzip, gzip };

class DynamicLibraries;  // set from parameters.dyn_libs in fastx_open()

/* Deleter that closes an open gzip/bzip2 stream through the borrowed
   DynamicLibraries facade, routing to the matching close function by Format.
   operator() is defined in fastx.cpp, where DynamicLibraries is complete. */
struct CompressedStreamDeleter
{
  DynamicLibraries const * libraries = nullptr;
  Format format = Format::undefined;

  CompressedStreamDeleter() = default;  // empty handle: never invoked
  CompressedStreamDeleter(DynamicLibraries const * libs, Format fmt) noexcept
    : libraries(libs), format(fmt) {}

  auto operator()(void * stream) const noexcept -> void;
};

/* Owning handle for an open compressed stream: reset()/destruction closes it
   through the deleter. Empty (and a no-op to destroy) for a plain file. The
   handle is type-erased to void* so this header needs neither zlib.h nor
   bzlib.h; the real gzFile/BZFILE types stay inside DynamicLibraries. */
using CompressedStream = std::unique_ptr<void, CompressedStreamDeleter>;

struct Line_fragment;  // defined below; named by the friend declarations

struct fastx_s
{
private:
  /* The data members are private: commands reach the reader only through the
     member accessors below (mirroring Database). The free functions that make
     up the reader's implementation -- the opener, the FASTA/FASTQ record
     parsers and the shared line/buffer primitives, all split across fastx.cpp,
     fasta.cpp and fastq.cpp -- are granted access as friends. */
  friend auto fastx_open(char const * filename, struct Parameters const & parameters) -> std::unique_ptr<fastx_s>;
  friend auto fastx_file_fill_buffer(fastx_s * input_handle) -> uint64_t;
  friend auto fastx_filter_header(fastx_s * input_handle, bool truncateatspace) -> void;
  friend auto fastx_filter_sequence_length(fastx_s * input_handle) -> void;
  friend auto fasta_next(fastx_s * input_handle, bool truncateatspace, unsigned char const * char_mapping) -> bool;
  friend auto fasta_filter_sequence(fastx_s * input_handle, unsigned char const * char_mapping) -> void;
  friend auto fastq_next(fastx_s * input_handle, bool truncateatspace, unsigned char const * char_mapping) -> bool;
  friend auto scan_line_fragment(fastx_s * input_handle) -> Line_fragment;
  friend auto consume_fragment(fastx_s * input_handle, Line_fragment const & fragment) -> void;

  bool is_pipe = false;
  bool is_fastq = false;
  bool is_empty = false;

  std::FILE * fp = nullptr;

  /* runtime-loaded compression libraries, borrowed (non-owning) from the
     DynamicLibraries instance in main(); nullptr in library-only builds */
  DynamicLibraries const * libraries = nullptr;

  /* the active compressed stream (gzip or bzip2), or empty for a plain file;
     an RAII handle that closes the stream on reset()/destruction (see
     CompressedStream / CompressedStreamDeleter above). */
  CompressedStream compressed_stream;

  FastxBuffer file_buffer;

  FastxBuffer header_buffer;
  FastxBuffer sequence_buffer;
  FastxBuffer plusline_buffer;
  FastxBuffer quality_buffer;

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

public:
  /* Read API, mirroring Database's accessors. The former fastx_get_, fasta_get_
     and fastq_get_ free functions were three near-identical families dispatching
     on the format; they collapse into this single member set (the FASTA/FASTQ
     difference is confined to get_quality/quality_view). The trivial accessors
     are inline here, exactly as Database inlines its getters, so returning a
     record field stays as cheap as the former free function. */

  // Format of the input. An empty input is accepted as FASTQ, preserving the
  // historical fastx_is_fastq() behaviour.
  auto is_fastq_input() const noexcept -> bool { return is_fastq or is_empty; }
  auto is_empty_input() const noexcept -> bool { return is_empty; }
  auto is_pipe_input() const noexcept -> bool { return is_pipe; }
  // The detected format is FASTQ (first byte '@'); false for FASTA and for an
  // empty input. Unlike is_fastq_input(), an empty input reports false here.
  auto is_fastq_format() const noexcept -> bool { return is_fastq; }

  // Deferred-error mode: a caller reading this handle from worker threads turns
  // it on so a parse error is recorded (see set_deferred_error) instead of
  // calling fatal() from a worker; defers_errors() reads the flag.
  auto enable_deferred_errors() noexcept -> void { defer_errors = true; }
  auto defers_errors() const noexcept -> bool { return defer_errors; }

  // Count one invalid character stripped from the input (for the end-of-input
  // report_stripped_warning()). Used by the FASTA/FASTQ sequence/quality
  // filters, which run over every input byte.
  auto record_stripped(unsigned char symbol) noexcept -> void
  {
    ++stripped_all;
    ++stripped[symbol];
  }

  // Current record; the returned pointers/views stay valid only until the next
  // next()/close() call on this handle.
  auto get_header() const noexcept -> char const * { return header_buffer.data(); }
  auto get_sequence() const noexcept -> char const * { return sequence_buffer.data(); }
  auto get_header_length() const noexcept -> uint64_t { return header_buffer.length; }
  auto get_sequence_length() const noexcept -> uint64_t { return sequence_buffer.length; }
  // Quality is meaningful only for FASTQ; a FASTA record reports none (nullptr),
  // matching the former fastx_get_quality().
  auto get_quality() const noexcept -> char const *
  {
    return is_fastq ? quality_buffer.data() : nullptr;
  }
  auto get_quality_length() const noexcept -> uint64_t { return quality_buffer.length; }
  auto get_abundance() const -> int64_t;               // 1 when ;size= is absent
  auto get_abundance_and_presence() const -> int64_t;  // 0 when ;size= is absent

  // View/SeqRecord companions, mirroring Database::sequence_view()/record().
  auto header_view() const noexcept -> View<char>
  {
    return View<char>{get_header(), get_header_length()};
  }
  auto sequence_view() const noexcept -> View<char>
  {
    return View<char>{get_sequence(), get_sequence_length()};
  }
  auto quality_view() const noexcept -> View<char>
  {
    return View<char>{is_fastq ? quality_buffer.data() : nullptr,
                      is_fastq ? get_sequence_length() : uint64_t{0}};
  }
  auto record() const -> SeqRecord
  {
    return SeqRecord{header_view(), sequence_view(), quality_view()};
  }

  // Stream position, total size and running counters.
  auto get_position() const noexcept -> uint64_t { return file_position; }
  auto get_size() const noexcept -> uint64_t { return file_size; }
  auto get_lineno() const noexcept -> uint64_t { return lineno_start; }
  auto get_seqno() const noexcept -> uint64_t { return static_cast<uint64_t>(seqno); }

  // Deferred-error protocol (see the defer_errors note above).
  auto get_error() const noexcept -> bool { return error; }
  auto get_errmsg() const noexcept -> char const * { return errmsg.data(); }
  auto set_deferred_error(char const * message) -> void;

  // Advance to the next record, dispatching to the FASTA or FASTQ parser by
  // format. Returns false at end of input or on a deferred parse error.
  auto next(bool truncateatspace, unsigned char const * char_mapping) -> bool;

  // Emit the end-of-input warning about invalid characters stripped from the
  // input (to stderr and, when open, the log file). This is the user-facing
  // half of the former fastx_close(); the handle's storage is now released by
  // its owning std::unique_ptr, so there is no separate close/delete step.
  auto report_stripped_warning(struct Parameters const & parameters) const -> void;

  /* Frees the owned resources (open files and buffers). Having it here means a
     fastx_s held in a std::unique_ptr is cleaned up automatically when the
     stack unwinds — e.g. when fatal() throws in a library session part-way
     through fastx_open() or a read loop. It must not throw (destructors run
     during unwinding), so it never calls fatal(); fastx_close() keeps the
     user-facing stripped-character warning and then deletes the handle. */
  ~fastx_s();
};

using fastx_handle = struct fastx_s *;


/* fastx input */

/* The record read API (get_header/get_sequence/get_quality/..., record(),
   is_fastq_input(), next(), report_stripped_warning(), the deferred-error
   protocol) is now a member set on fastx_s above, mirroring Database. fastx_open
   returns a std::unique_ptr<fastx_s> that owns the handle (RAII: closed and
   freed when it goes out of scope), so there is no fastx_close free function.
   These remaining free functions are not simple accessors: the opener and the
   two in-parser filters. */
auto fastx_filter_header(fastx_handle input_handle, bool truncateatspace) -> void;
auto fastx_open(const char * filename, struct Parameters const & parameters) -> std::unique_ptr<fastx_s>;

// Reject a sequence too long for the int length bookkeeping used downstream.
// Called from fasta_next/fastq_next so every read is bounded at one choke
// point, symmetric with fastx_filter_header. On a worker thread an over-long
// sequence records a deferred error (reported from the main thread), otherwise
// it is fatal.
auto fastx_filter_sequence_length(fastx_handle input_handle) -> void;

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
  auto * const start = std::next(file_buffer.data(),
                                 static_cast<std::ptrdiff_t>(file_buffer.position));
  auto const * const line_end = std::char_traits<char>::find(start, rest, '\n');
  auto const has_newline = (line_end != nullptr);
  auto const length = has_newline
    ? static_cast<std::size_t>(line_end - start) + 1
    : static_cast<std::size_t>(rest);
  return Line_fragment{View<char>{start, length}, has_newline};
}

// Advance the file-buffer read position past a fragment already copied out.
inline auto consume_fragment(fastx_handle input_handle,
                             Line_fragment const & fragment) -> void
{
  input_handle->file_buffer.position += fragment.view.size();
}
