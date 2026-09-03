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

#include "core/quality_range.hpp"  // QualityLocation
#include "core/seq_record.hpp"  // SeqRecord (returned by fastx_record)
#include "utils/fatal_allocator.hpp"  // FatalAllocator
#include "utils/maps.hpp"  // Mapping
#include "utils/quality_encoding.hpp"  // QualitySymbolRange, sanger_ascii_offset
#include "utils/span.hpp"  // Span
#include "utils/view.hpp"  // View
#include <array>
#include <cassert>  // assert
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
  /* the bytes in use as a single read-only window, for callers that scan or
     compare the whole buffer instead of indexing into it */
  auto view() const noexcept -> View<char> { return View<char>{data(), length}; }
  /* the bytes the read cursor has not passed yet: the window the line scanner
     works on. Only the file refill buffer advances 'position', so for the
     per-record buffers this is the same window as view(). */
  auto unread() const noexcept -> View<char>
  {
    assert(position <= length);
    auto const distance = static_cast<std::ptrdiff_t>(position);
    return View<char>{std::next(data(), distance), length - position};
  }
  /* the byte at the read cursor, for the record-boundary sentinel tests in the
     FASTA/FASTQ parsers. Deliberately a direct read and not unread().front():
     building a window to look at one byte cost 7% of the runtime on an input of
     one-nucleotide records, where the test runs once per line. */
  auto peek() const noexcept -> char
  {
    assert(position < length);
    return data()[position];
  }
  /* writable companion to view(), for the in-place filters that rewrite the
     bytes in use (see fastx_filter_header); mirrors Database::mutable_sequence()
     in that only a non-const buffer hands out a mutable window */
  auto span() noexcept -> Span<char> { return Span<char>{data(), length}; }

  /* Reset to a single empty, NUL-terminated block (former buffer_init). */
  auto init() -> void;
  /* Ensure at least 'size' more bytes fit after 'length', rounding the
     allocation up to the nearest block (former buffer_makespace). */
  auto makespace(uint64_t size) -> void;
  /* Append the bytes of 'source', then a terminating NUL (former
     buffer_extend). */
  auto extend(View<char> source) -> void;

private:
  std::vector<char, FatalAllocator<char>> storage_;
};

enum struct Format : unsigned char { undefined, plain, bzip, gzip };

class DynamicLibraries;  // set from parameters.runtime.dyn_libs in fastx_open()

/* Deleter that closes an open gzip/bzip2 stream through the borrowed
   DynamicLibraries facade, routing to the matching close function by Format.
   operator() is defined in fastx.cpp, where DynamicLibraries is complete. */
struct CompressedStreamDeleter
{
  DynamicLibraries const * libraries = nullptr;
  Format format = Format::undefined;

  CompressedStreamDeleter() = default;  // empty handle: never invoked
  CompressedStreamDeleter(DynamicLibraries const * libs, Format const fmt) noexcept
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
  template <Mapping mapping>
  friend auto fasta_filter_sequence(fastx_s * input_handle) -> void;
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

  /* The lowest and highest quality symbol seen anywhere in this file, and the
     offset the caller asked for. fastq_next() compares the two once, when it
     reaches end of file, and warns if the symbols contradict the offset --
     which is the only detection left now that --fastq_qmax accepts the whole
     representable range (see DONE_20260825_quality_range.md). Filled by the
     quality branch of buffer_filter_extend(), which already reads every one
     of those bytes. */
  QualitySymbolRange quality_range;
  int quality_offset = sanger_ascii_offset;
  bool warn_on_suspicious_offset = true;

  /* The encoding is a property of the file, not of a record, so the range is
     sampled from the first records only and the per-byte updates stop after
     that. Tracking every byte of a 400k-record file cost 3.9% of
     --fastq_stats, the most parser-bound command; bounding it makes the cost
     O(1) in file size (0.1% of that same run) while still reading tens of
     thousands of symbols, far more than any encoding heuristic needs. The
     choice is made once per record, never per byte. */
  static constexpr int64_t offset_sample_records = 10000;
  auto samples_quality_range() const noexcept -> bool {
    return seqno < offset_sample_records;
  }

  /* Below this many records the observed range is not evidence of anything:
     a file holding a single 'K' is a legal Sanger Q42 and a legal Illumina
     1.5+ Q11, and warning about it would fire on every small hand-made test
     file. --fastq_chars still prints its guess for such a file, because a
     guess the user asked for may be speculative where an unsolicited warning
     may not. */
  static constexpr int64_t minimum_records_for_offset_guess = 100;

  /* The same two bytes for the current record alone, reset by fastq_next().
     A command that checks the --fastq_qmin/--fastq_qmax window per record
     asks for this with track_quality_range() and reads it back with
     quality_symbol_range(), instead of walking the quality string a second
     time of its own. Opt-in because the per-byte updates are not free (see
     offset_sample_records above): a command that does not check the window
     pays one branch per line fragment and nothing else. */
  QualitySymbolRange record_quality_range;
  bool track_record_quality_range = false;

  auto tracks_quality_range() const noexcept -> bool {
    return track_record_quality_range or samples_quality_range();
  }

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

  // The shared body of the two set_deferred_error() overloads below: the
  // message arrives already measured, so neither of them needs std::strlen.
  auto record_deferred_error(View<char> message) -> void;

  /* Compare the observed quality symbols against quality_offset and warn if
     they disagree; called once by fastq_next() at end of file. Defined in
     fastq.cpp, where the rationale lives. */
  auto warn_if_offset_looks_wrong() -> void;

public:
  /* Read API, mirroring Database's accessors. The former fastx_get_, fasta_get_
     and fastq_get_ free functions were three near-identical families dispatching
     on the format; they collapsed into a single member set, since narrowed to
     the View/SeqRecord accessors below once every caller consumed views (the
     FASTA/FASTQ difference is confined to quality_view). The trivial accessors
     are inline here, exactly as Database inlines its getters, so returning a
     record field stays as cheap as the former free function. */

  /* --fastq_chars exists to diagnose the offset and prints its own guess, so
     it silences the reader's warning rather than saying the same thing twice. */
  auto silence_offset_warning() noexcept -> void { warn_on_suspicious_offset = false; }

  /* Ask the parser to record the quality-symbol range of every record, and
     read back the range of the record next() just returned. Call the setter
     once, after opening: it costs the reader two comparisons per quality
     byte, and saves the caller a second pass over the same string. The range
     is empty (seen() false) for a FASTA record or an empty quality line. */
  auto track_quality_range() noexcept -> void { track_record_quality_range = true; }
  auto quality_symbol_range() const noexcept -> QualitySymbolRange { return record_quality_range; }

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
  auto record_stripped(unsigned char const symbol) noexcept -> void
  {
    ++stripped_all;
    ++stripped[symbol];
  }

  // Current record; the returned views stay valid only until the next
  // next()/close() call on this handle.
  auto get_abundance() const -> int64_t;               // 1 when ;size= is absent
  auto get_abundance_and_presence() const -> int64_t;  // 0 when ;size= is absent

  // View/SeqRecord accessors, mirroring Database::sequence_view()/record().
  // The header and the sequence are exactly the bytes their buffer holds, so
  // they hand out the buffer's own window; quality does not, see below.
  auto header_view() const noexcept -> View<char>
  {
    return header_buffer.view();
  }
  auto sequence_view() const noexcept -> View<char>
  {
    return sequence_buffer.view();
  }
  // Deliberately not quality_buffer.view(): a FASTA record has no quality at
  // all (nullptr, empty), and for FASTQ the window is sized from the sequence,
  // the length fastq_next() has just checked the quality string against.
  auto quality_view() const noexcept -> View<char>
  {
    return View<char>{is_fastq ? quality_buffer.data() : nullptr,
                      is_fastq ? sequence_buffer.length : uint64_t{0}};
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

  /* Where the current record sits in the input, for a message that names it.
     seqno counts from -1 (set by fastx_open), so it is 0-based once a record
     has been read and the ordinal is one more; keeping that arithmetic here
     means no caller repeats it, and no caller can fill the two members of
     QualityLocation in the wrong order. */
  auto quality_location() const noexcept -> vsearch::QualityLocation {
    return vsearch::QualityLocation{get_seqno() + 1, get_lineno()};
  }

  // Deferred-error protocol (see the defer_errors note above).
  auto get_error() const noexcept -> bool { return error; }
  auto get_errmsg() const noexcept -> char const * { return errmsg.data(); }
  auto set_deferred_error(char const * message) -> void;
  // The assembled-message form, for the callers that build the text with
  // std::string concatenation: the string knows its own length, so the copy
  // below takes it from size() instead of measuring the same bytes again with
  // std::strlen -- and stays correct if a message ever carries an embedded
  // '\0'. A string literal still selects the overload above (array-to-pointer
  // is an exact match, the std::string conversion is user-defined), exactly as
  // described for fatal() in utils/fatal.hpp.
  auto set_deferred_error(std::string const & message) -> void;

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
  // the bytes not yet consumed, so that neither the remaining length nor the
  // start of the scan is computed here, and the fragment below is carved out of
  // a bounded window
  auto const unread = input_handle->file_buffer.unread();
  auto const * const line_end =
    std::char_traits<char>::find(unread.data(), unread.size(), '\n');
  auto const has_newline = (line_end != nullptr);
  auto const length = has_newline
    ? static_cast<std::size_t>(std::distance(unread.data(), line_end)) + 1
    : unread.size();
  return Line_fragment{unread.first(length), has_newline};
}

// Advance the file-buffer read position past a fragment already copied out.
inline auto consume_fragment(fastx_handle input_handle,
                             Line_fragment const & fragment) -> void
{
  input_handle->file_buffer.position += fragment.view.size();
}
