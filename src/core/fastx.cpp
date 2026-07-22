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

#include "vsearch.hpp"
#include "core/attributes.hpp"  // header_get_size
#include "core/buffer_headroom.hpp"  // buffer_headroom
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "os/dynlibs.hpp"
#include "os/system.hpp"  // xstat_t, xfstat, xlseek, xftello, S_ISREG
#include "utils/fatal.hpp"
#include "utils/logfile.hpp"  // log_file::handle
#include "utils/make_unique.hpp"  // make_unique
#include "utils/open_file.hpp"  // open_input_file
#include "utils/span.hpp"
#include <sys/stat.h>
#include <unistd.h>  // dup, STDOUT_FILENO
#include <algorithm>  // std::copy_n, std::equal, std::find_first_of
#include <array>
#include <cassert>  // assert
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t, std::fread, std::fileno
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::memcmp, std::strcmp
#include <iterator> // std::distance
#include <limits>  // std::numeric_limits
#include <memory>  // std::unique_ptr
#include <vector>


/* file compression and format detector */
/* basic file buffering function for fastq and fastx parsers */

constexpr uint64_t fastx_buffer_alloc = 8192;

constexpr std::array<unsigned char, 2> magic_gzip = {0x1f, 0x8b};
constexpr std::array<unsigned char, 2> magic_bzip = {'B', 'Z'};


auto FastxBuffer::init() -> void
{
  /* one block, holding an empty NUL-terminated string; resize() (not reserve())
     so data() is valid across the whole block, matching the former xmalloc. */
  storage_.resize(fastx_buffer_alloc);
  storage_[0] = 0;
  length = 0;
  position = 0;
}


auto FastxBuffer::makespace(uint64_t const size) -> void
{
  /* make sure there is space for 'size' more chars in the buffer */

  if (length + size > alloc())
    {
      /* grow to fit 'size' more characters, rounding the allocation up to the
         nearest block, as the former xrealloc-based growth did. resize()
         preserves the bytes already in use and zero-fills the new tail (which
         extend() overwrites); it never shrinks here because makespace() only
         runs when the requested total exceeds the current capacity. */
      auto const new_alloc =
        ((length + size + fastx_buffer_alloc - 1) / fastx_buffer_alloc)
        * fastx_buffer_alloc;
      storage_.resize(static_cast<std::size_t>(new_alloc));
    }
}


auto FastxBuffer::extend(char const * const source, uint64_t const len) -> void
{
  makespace(len + 1);
  std::copy_n(source, len, data() + length);
  length += len;
  data()[length] = 0;
}


/* fastx_s member functions that are too large to inline in the header (they
   pull in header_get_size, the format dispatch, or the error buffer). The
   trivial record accessors are inline in fastx.hpp. */

auto fastx_s::get_abundance() const -> int64_t
{
  // return 1 if the ;size= annotation is not present
  auto const size = header_get_size(header_buffer.data(),
                                    static_cast<int>(header_buffer.length));
  if (size > 0)
    {
      return size;
    }
  return 1;
}


auto fastx_s::get_abundance_and_presence() const -> int64_t
{
  // return 0 if the ;size= annotation is not present
  return header_get_size(header_buffer.data(), static_cast<int>(header_buffer.length));
}


auto fastx_s::set_deferred_error(char const * const message) -> void
{
  /* record the first deferred parse error and flag the handle (see the
     deferred-error note in fastx.hpp). First error wins; later ones are
     ignored so the message reflects the earliest failure. */
  if (not error)
    {
      std::snprintf(errmsg.data(), errmsg.size(), "%s", message);
      error = true;
    }
}


auto fastx_s::next(bool const truncateatspace, unsigned char const * char_mapping) -> bool
{
  /* deferred-error mode (see fastx.hpp): if a previous call already recorded a
     parse error, report no further records so every worker stops; and if the
     current record triggered a deferred error, treat it as unusable (return
     false) rather than handing back a bogus record. */
  if (error)
    {
      return false;
    }
  bool const got_record = is_fastq
    ? fastq_next(this, truncateatspace, char_mapping)
    : fasta_next(this, truncateatspace, char_mapping);
  if (error)
    {
      return false;
    }
  return got_record;
}


auto find_header_end_first_blank(Span<char> raw_header) -> std::size_t {
  static const std::vector<char> blanks {' ', '\t', '\0', '\r', '\n'};
  auto * result = std::find_first_of(raw_header.begin(), raw_header.end(),
                                     blanks.begin(), blanks.end());
  if (result != raw_header.end()) {
    *result = '\0';
  }
  return static_cast<std::size_t>(std::distance(raw_header.begin(), result));
}


auto find_header_end(Span<char> raw_header) -> std::size_t {
  static const std::vector<char> blanks {'\0', '\r', '\n'};
  auto * result = std::find_first_of(raw_header.begin(), raw_header.end(),
                                     blanks.begin(), blanks.end());
  if (result != raw_header.end()) {
    *result = '\0';
  }
  return static_cast<std::size_t>(std::distance(raw_header.begin(), result));
}


// Emit a pre-formatted warning line to stderr and, when open, the log file.
// The caller formats the message itself (mirroring the snprintf-into-a-buffer
// idiom used by the fatal()/deferred paths below), so this stays a simple
// reporter and the message's format string is validated by -Wformat at the
// call site instead of being forwarded as a runtime argument to fprintf.
auto warn(char const * const message) -> void {
  std::fprintf(stderr, "\nWARNING: %s\n", message);

  auto * const log = log_file::handle();
  if (log != nullptr) {
    std::fprintf(log, "\nWARNING: %s\n", message);
  }
}


auto fastx_filter_header(fastx_handle input_handle, bool const truncateatspace) -> void {
  // truncate header (in-place)
  auto raw_header = Span<char>{input_handle->header_buffer.data(), input_handle->header_buffer.length};
  auto const count = truncateatspace ? find_header_end_first_blank(raw_header) : find_header_end(raw_header);
  input_handle->header_buffer.length = count;

  /* Reject a header too long for the int header-length bookkeeping used
     downstream (searchinfo_s::query_head_len). Such a header would be narrowed
     to a negative int, breaking the per-query header allocation. The
     sequence length is bounded by --maxseqlength, but the header length is not,
     so guard it here, at the single point all FASTA/FASTQ (and DB) reads pass
     through. Mirrors the deferred/fatal handling of the illegal-character check
     below: on a worker thread the error is recorded and reported from the main
     thread, never fatal()ed here. */
  static constexpr auto max_header_length =
    static_cast<std::size_t>(std::numeric_limits<int>::max() - buffer_headroom);
  if (count > max_header_length) {
    std::array<char, 256> message {{}};
    std::snprintf(message.data(), message.size(),
                  "FASTA/FASTQ header too long (%" PRIu64 " bytes) on line %"
                  PRIu64 ".\nHeaders longer than %" PRIu64 " bytes are not supported.",
                  static_cast<uint64_t>(count), input_handle->lineno_start,
                  static_cast<uint64_t>(max_header_length));
    if (input_handle->defer_errors) {
      input_handle->set_deferred_error(message.data());
      return;
    }
    fatal(message.data());
  }

  // scan for unusual symbols
  auto const trimmed_header = raw_header.first(count);
  for (auto const symbol: trimmed_header) {
    auto const is_illegal = ((symbol == 127) or
                             ((symbol > '\0') and (symbol < ' ') and not (symbol == '\t')));
    if (is_illegal) {
      if (input_handle->defer_errors) {
        /* Record the error and stop scanning instead of exiting here:
           this may run on a worker thread (see the deferred-error note in
           fastx.h). The caller reports it from the main thread. */
        std::array<char, 256> message {{}};
        std::snprintf(message.data(), message.size(),
                      "Illegal character encountered in FASTA/FASTQ header.\n"
                      "Unprintable ASCII character no %d on line %" PRIu64 ".",
                      symbol, input_handle->lineno_start);
        input_handle->set_deferred_error(message.data());
        return;
      }
      fatal("Illegal character encountered in FASTA/FASTQ header.\n"
            "Unprintable ASCII character no %d on line %" PRIu64 ".",
            symbol, input_handle->lineno_start);
    }
    auto const symbol_unsigned = static_cast<unsigned char>(symbol);
    auto const is_not_ascii = (symbol_unsigned > 127);
    if (is_not_ascii) {
      std::array<char, 256> message {{}};
      std::snprintf(message.data(), message.size(),
                    "Non-ASCII character encountered in FASTA/FASTQ header.\n"
                    "Character no %d (0x%2x) on line %" PRIu64 ".",
                    symbol_unsigned, symbol_unsigned, input_handle->lineno_start);
      warn(message.data());
    }
  }
}


auto fastx_open(char const * filename, struct Parameters const & parameters) -> std::unique_ptr<fastx_s>
{
  // Held in a unique_ptr so any fatal() below frees the partially-built handle
  // when the stack unwinds (library session), and so the caller owns the handle
  // by RAII: it is closed and freed when the returned unique_ptr goes out of
  // scope, with no explicit close/delete step.
  auto input_handle = make_unique<fastx_s>();

  input_handle->fp = nullptr;
  input_handle->libraries = parameters.dyn_libs;

  input_handle->fp = open_input_file(filename).release();
  if (input_handle->fp == nullptr)
    {
      fatal("Unable to open file for reading (%s)", filename);
    }

  /* Get mode and size of original (uncompressed) file */

  xstat_t fs;
  if (xfstat(fileno(input_handle->fp), & fs) != 0)
    {
      fatal("Unable to get status for input file (%s)", filename);
    }

  /* Treat anything that is not a regular file as a non-rewindable stream:
     named pipes, sockets, and character devices (such as FreeBSD's
     /dev/stdin and the /dev/fd/N entries used by shell process
     substitution) cannot be seeked, sized, or reopened from the start.
     The compression autodetection below consumes bytes and then closes
     and reopens the file to rewind, which silently corrupts such streams,
     so it must be skipped for them. */
  input_handle->is_pipe = not S_ISREG(fs.st_mode);

  if (input_handle->is_pipe)
    {
      input_handle->file_size = 0;
    }
  else
    {
      input_handle->file_size = static_cast<uint64_t>(fs.st_size);
    }

  bool const is_stdin = (std::strcmp(filename, "-") == 0);

  if (input_handle->is_pipe and is_stdin)
    {
      /* stdin cannot be rewound or peeked without consuming bytes,
         so rely on the user-provided flags */
      if (parameters.opt_gzip_decompress)
        {
          input_handle->format = Format::gzip;
        }
      else if (parameters.opt_bzip2_decompress)
        {
          input_handle->format = Format::bzip;
        }
      else
        {
          input_handle->format = Format::plain;
        }
    }
  else if (input_handle->is_pipe)
    {
      /* non-stdin, non-rewindable stream (e.g. bash process
         substitution, a named FIFO, or a character device such as
         /dev/stdin on FreeBSD): the decompress flags were meant for
         stdin, not for arbitrary streams wired in by the shell. Assume
         plain; compressed inputs should be passed as a regular file or
         pre-decompressed upstream of the pipe. */
      input_handle->format = Format::plain;
    }
  else
    {
      /* autodetect compression (plain, gzipped or bzipped) */

      /* read two characters and compare with magic */

      std::array<unsigned char, 2> magic {{}};

      input_handle->format = Format::plain;

      // refactoring: fread() see C++ Weekly - Ep 482 - Safely Wrapping C APIs
      auto const bytes_read = std::fread(magic.data(), 1, 2, input_handle->fp);

      if (bytes_read >= 2)
        {
          if (std::equal(magic.begin(), magic.end(), magic_gzip.begin()))
            {
              input_handle->format = Format::gzip;
            }
          else if (std::equal(magic.begin(), magic.end(), magic_bzip.begin()))
            {
              input_handle->format = Format::bzip;
            }
        }
      else
        {
          /* consider it an empty file or a tiny fasta file, uncompressed */
        }

      /* close and reopen to avoid problems with gzip library */
      /* rewind was not enough */

      std::fclose(input_handle->fp);
      input_handle->fp = open_input_file(filename).release();
      if (input_handle->fp == nullptr)
        {
          fatal("Unable to open file for reading (%s)", filename);
        }
    }

  if (input_handle->format == Format::gzip)
    {
      /* GZIP: Keep original file open, then open as gzipped file as well */
      if ((not compression::gzip_supported) or (input_handle->libraries == nullptr)
          or not input_handle->libraries->gzip_available())
        {
          fatal("Files compressed with gzip are not supported");
        }
      input_handle->compressed_stream = CompressedStream(
          input_handle->libraries->gz_open(fileno(input_handle->fp)),
          CompressedStreamDeleter{input_handle->libraries, Format::gzip});
      if (input_handle->compressed_stream == nullptr)
        { // dup?
          fatal("Unable to open gzip compressed file (%s)", filename);
        }
    }

  if (input_handle->format == Format::bzip)
    {
      /* BZIP2: Keep original file open, then open as bzipped file as well */
      if ((not compression::bzip2_supported) or (input_handle->libraries == nullptr)
          or not input_handle->libraries->bzip2_available())
        {
          fatal("Files compressed with bzip2 are not supported");
        }
      input_handle->compressed_stream = CompressedStream(
          input_handle->libraries->bz_open(input_handle->fp),
          CompressedStreamDeleter{input_handle->libraries, Format::bzip});
      if (input_handle->compressed_stream == nullptr)
        {
          fatal("Unable to open bzip2 compressed file (%s)", filename);
        }
    }

  /* init buffers */

  input_handle->file_position = 0;

  input_handle->file_buffer.init();

  /* start filling up file buffer */

  auto const rest = fastx_file_fill_buffer(input_handle.get());

  /* examine first char and see if it starts with > or @ */

  input_handle->is_empty = true;
  input_handle->is_fastq = false;

  if (rest > 0)
    {
      input_handle->is_empty = false;

      auto filetype = 0;
      auto const * first = input_handle->file_buffer.data();

      if (*first == '>')
        {
          filetype = 1;
        }
      else if (*first == '@')
        {
          filetype = 2;
          input_handle->is_fastq = true;
        }

      if (filetype == 0)
        {
          /* close files if unrecognized file type */

          input_handle->compressed_stream.reset();

          std::fclose(input_handle->fp);
          input_handle->fp = nullptr;

          if (rest >= 2)
            {
              if (std::memcmp(first, magic_gzip.data(), 2) == 0)
                {
                  fatal("File appears to be gzip compressed. Please use --gzip_decompress");
                }

              if (std::memcmp(first, magic_bzip.data(), 2) == 0)
                {
                  fatal("File appears to be bzip2 compressed. Please use --bzip2_decompress");
                }
            }

          fatal("File type not recognized.");
        }
    }

  /* more initialization */

  input_handle->header_buffer.init();
  input_handle->sequence_buffer.init();
  input_handle->plusline_buffer.init();
  input_handle->quality_buffer.init();

  input_handle->stripped_all = 0;

  input_handle->lineno = 1;
  input_handle->lineno_start = 1;
  input_handle->seqno = -1;

  assert(input_handle != nullptr);
  return input_handle;
}


auto CompressedStreamDeleter::operator()(void * const stream) const noexcept -> void
{
  if ((libraries == nullptr) or (stream == nullptr))
    {
      return;
    }
  if (format == Format::gzip)
    {
      libraries->gz_close(stream);
    }
  else if (format == Format::bzip)
    {
      libraries->bz_close(stream);
    }
}


/* Release the owned resources: the compressed stream (if any) and the
   underlying file. compressed_stream is an RAII handle, so reset() closes any
   open gzip/bzip2 stream through its deleter; the five FastxBuffer members
   free their own storage. The stream is closed before fclose(fp) because the
   library stream wraps fp's descriptor. Runs on delete, including during stack
   unwinding when fatal() throws in a library session, so it must not throw: it
   calls no fatal() and null-guards fp, and the deleter no-ops on the empty
   handle of a reader abandoned part-way through fastx_open(). */
fastx_s::~fastx_s()
{
  compressed_stream.reset();

  if (fp != nullptr)
    {
      std::fclose(fp);
    }
}


auto fastx_s::report_stripped_warning(struct Parameters const & parameters) const -> void
{
  /* Warn about stripped chars */

  if (stripped_all != 0U)
    {
      std::fprintf(stderr, "WARNING: %" PRIu64 " invalid characters stripped from %s file:", stripped_all, (is_fastq ? "FASTQ" : "FASTA"));
      for (int i = 0; i < 256; i++)
        {
          if (stripped[static_cast<std::size_t>(i)] != 0U)
            {
              std::fprintf(stderr, " %c(%" PRIu64 ")", i, stripped[static_cast<std::size_t>(i)]);
            }
        }
      std::fprintf(stderr, "\n");
      std::fprintf(stderr, "REMINDER: vsearch does not support amino acid sequences\n");

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log, "WARNING: %" PRIu64 " invalid characters stripped from %s file:", stripped_all, (is_fastq ? "FASTQ" : "FASTA"));
          for (int i = 0; i < 256; i++)
            {
              if (stripped[static_cast<std::size_t>(i)] != 0U)
                {
                  std::fprintf(parameters.fp_log, " %c(%" PRIu64 ")", i, stripped[static_cast<std::size_t>(i)]);
                }
            }
          std::fprintf(parameters.fp_log, "\n");
          std::fprintf(parameters.fp_log, "REMINDER: vsearch does not support amino acid sequences\n");
        }
    }
}


auto fastx_file_fill_buffer(fastx_handle input_handle) -> uint64_t
{
  /* read more data if necessary */
  uint64_t const rest = input_handle->file_buffer.length - input_handle->file_buffer.position;

  if (rest > 0)
    {
      return rest;
    }

  uint64_t space = input_handle->file_buffer.alloc() - input_handle->file_buffer.length;

  if (space == 0)
    {
      /* back to beginning of buffer */
      input_handle->file_buffer.position = 0;
      input_handle->file_buffer.length = 0;
      space = input_handle->file_buffer.alloc();
    }

  int bytes_read = 0;

  switch (input_handle->format)
    {
    case Format::plain:
      bytes_read = static_cast<int>(std::fread(input_handle->file_buffer.data()
                         + input_handle->file_buffer.position,
                         1,
                         space,
                         input_handle->fp));
      break;

    case Format::gzip:
      bytes_read = input_handle->libraries->gz_read(input_handle->compressed_stream.get(),
                               input_handle->file_buffer.data()
                               + input_handle->file_buffer.position,
                               static_cast<unsigned int>(space));
      if (bytes_read < 0)
        {
          fatal("Unable to read gzip compressed file");
        }
      break;

    case Format::bzip:
      bytes_read = input_handle->libraries->bz_read(input_handle->compressed_stream.get(),
                                   input_handle->file_buffer.data()
                                   + input_handle->file_buffer.position,
                                   static_cast<int>(space));
      if (bytes_read < 0)
        {
          fatal("Unable to read from bzip2 compressed file");
        }
      break;

    default:
      fatal("Internal error");
    }

  if (not input_handle->is_pipe)
    {
      if (input_handle->format == Format::gzip)
        {
          /* Circumvent the missing gzoffset function in zlib 1.2.3 and earlier */
          int const fd = dup(fileno(input_handle->fp));
          input_handle->file_position = xlseek(fd, 0, SEEK_CUR);
          close(fd);
        }
      else
        {
          input_handle->file_position = xftello(input_handle->fp);
        }
    }

  input_handle->file_buffer.length += static_cast<uint64_t>(bytes_read);
  return static_cast<uint64_t>(bytes_read);
}


auto fastx_filter_sequence_length(fastx_handle input_handle) -> void
{
  /* Reject a sequence too long for the int sequence-length bookkeeping used
     downstream (e.g. searchinfo_s::qseqlen, sized as qseqlen + buffer_headroom;
     the chimera query buffers; cut's rc_buffer). Such a sequence would be
     narrowed to a negative int and overflow the per-record buffer. Database
     sequences over --maxseqlength are discarded by Database::read, but many commands
     read records directly via next() with no length filter, so
     guard it here, at the single point all FASTA/FASTQ (and DB) reads pass
     through -- symmetric with the header guard in fastx_filter_header. Mirrors
     its deferred/fatal handling: on a worker thread the error is recorded and
     reported from the main thread, never fatal()ed here. */
  static constexpr auto max_sequence_length =
    static_cast<uint64_t>(std::numeric_limits<int>::max() - buffer_headroom);
  auto const length = input_handle->sequence_buffer.length;
  if (length <= max_sequence_length)
    {
      return;
    }
  std::array<char, 256> message {{}};
  std::snprintf(message.data(), message.size(),
                "FASTA/FASTQ sequence too long (%" PRIu64 " nt) on line %"
                PRIu64 ".\nSequences longer than %" PRIu64 " nt are not supported.",
                length, input_handle->lineno_start, max_sequence_length);
  if (input_handle->defer_errors)
    {
      input_handle->set_deferred_error(message.data());
      return;
    }
  fatal(message.data());
}


