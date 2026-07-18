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
#include "core/buffer_headroom.hpp"  // buffer_headroom
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "os/dynlibs.hpp"
#include "os/system.hpp"  // xmalloc, xrealloc, xfree
#include "utils/fatal.hpp"
#include "utils/logfile.hpp"  // log_file::handle
#include "utils/make_unique.hpp"  // make_unique
#include "utils/open_file.hpp"  // open_input_file
#include "utils/span.hpp"
#include <unistd.h>  // dup, STDOUT_FILENO
#include <algorithm>  // std::find_first_of
#include <array>
#include <cassert>  // assert
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t, std::fread, std::fileno
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::memcpy, std::memcmp, std::strcmp
#include <iterator> // std::distance
#include <limits>  // std::numeric_limits
#include <memory>  // std::unique_ptr
#include <vector>


/* file compression and format detector */
/* basic file buffering function for fastq and fastx parsers */

constexpr uint64_t fastx_buffer_alloc = 8192;

#ifdef HAVE_BZLIB_H
constexpr auto BZ_VERBOSE_0 = 0;
// constexpr auto BZ_VERBOSE_1 = 1;
// constexpr auto BZ_VERBOSE_2 = 2;
// constexpr auto BZ_VERBOSE_3 = 3;
// constexpr auto BZ_VERBOSE_4 = 4;
constexpr auto BZ_MORE_MEM = 0;  /* faster decompression using more memory */
// constexpr auto BZ_LESS_MEM = 1;  /* slower decompression but requires less memory */
#endif

constexpr std::array<unsigned char, 2> magic_gzip = {0x1f, 0x8b};
constexpr std::array<unsigned char, 2> magic_bzip = {'B', 'Z'};


auto buffer_init(struct fastx_buffer_s * buffer) -> void
{
  buffer->alloc = fastx_buffer_alloc;
  buffer->data = static_cast<char *>(xmalloc(buffer->alloc));
  buffer->data[0] = 0;
  buffer->length = 0;
  buffer->position = 0;
}


auto buffer_free(struct fastx_buffer_s * buffer) -> void
{
  if (buffer->data != nullptr)
    {
      xfree(buffer->data);
    }
  buffer->data = nullptr;
  buffer->alloc = 0;
  buffer->length = 0;
  buffer->position = 0;
}


auto buffer_makespace(struct fastx_buffer_s * buffer, uint64_t const size) -> void
{
  /* make sure there is space for x more chars in buffer */

  if (buffer->length + size > buffer->alloc)
    {
      /* alloc space for x more characters,
         but round up to nearest block size */
      buffer->alloc =
        ((buffer->length + size + fastx_buffer_alloc - 1) / fastx_buffer_alloc)
        * fastx_buffer_alloc;
      buffer->data = static_cast<char *>(xrealloc(buffer->data, buffer->alloc));
    }
}


auto buffer_extend(struct fastx_buffer_s * dest_buffer,
                   char const * source_buf,
                   uint64_t const len) -> void
{
  buffer_makespace(dest_buffer, len + 1);
  std::memcpy(dest_buffer->data + dest_buffer->length,
         source_buf,
         len);
  dest_buffer->length += len;
  dest_buffer->data[dest_buffer->length] = 0;
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
  auto raw_header = Span<char>{input_handle->header_buffer.data, input_handle->header_buffer.length};
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
      fastx_set_deferred_error(input_handle, message.data());
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
        fastx_set_deferred_error(input_handle, message.data());
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


auto fastx_open(char const * filename, struct Parameters const & parameters) -> fastx_handle
{
  // refactoring: duplicate function to output a struct fastx_s input_handle_s;
  // Held in a unique_ptr so any fatal() below frees the partially-built handle
  // when the stack unwinds (library session); released to the caller on success.
  auto input_handle = make_unique<fastx_s>();

  input_handle->fp = nullptr;
  input_handle->libraries = parameters.dyn_libs;

#ifdef HAVE_ZLIB_H
  input_handle->fp_gz = nullptr;
#endif

#ifdef HAVE_BZLIB_H
  input_handle->fp_bz = nullptr;
  int bzError = 0;
#endif

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
          if (std::memcmp(magic.data(), magic_gzip.data(), 2) == 0)
            {
              input_handle->format = Format::gzip;
            }
          else if (std::memcmp(magic.data(), magic_bzip.data(), 2) == 0)
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
#ifdef HAVE_ZLIB_H
      if ((input_handle->libraries == nullptr) or not input_handle->libraries->gzip_available())
        {
          fatal("Files compressed with gzip are not supported");
        }
      input_handle->fp_gz = input_handle->libraries->gzdopen(fileno(input_handle->fp), "rb");
      if (input_handle->fp_gz == nullptr)
        { // dup?
          fatal("Unable to open gzip compressed file (%s)", filename);
        }
#else
      fatal("Files compressed with gzip are not supported");
#endif
    }

  if (input_handle->format == Format::bzip)
    {
      /* BZIP2: Keep original file open, then open as bzipped file as well */
#ifdef HAVE_BZLIB_H
      if ((input_handle->libraries == nullptr) or not input_handle->libraries->bzip2_available())
        {
          fatal("Files compressed with bzip2 are not supported");
        }
      input_handle->fp_bz = input_handle->libraries->bz_read_open(& bzError, input_handle->fp,
                                     BZ_VERBOSE_0, BZ_MORE_MEM,
                                     nullptr, 0);
      if (input_handle->fp_bz == nullptr)
        {
          fatal("Unable to open bzip2 compressed file (%s)", filename);
        }
#else
      fatal("Files compressed with bzip2 are not supported");
#endif
    }

  /* init buffers */

  input_handle->file_position = 0;

  buffer_init(& input_handle->file_buffer);

  /* start filling up file buffer */

  auto const rest = fastx_file_fill_buffer(input_handle.get());

  /* examine first char and see if it starts with > or @ */

  auto filetype = 0;
  input_handle->is_empty = true;
  input_handle->is_fastq = false;

  if (rest > 0)
    {
      input_handle->is_empty = false;

      auto const * first = input_handle->file_buffer.data;

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

          switch (input_handle->format)
            {
            case Format::plain:
              break;

            case Format::gzip:
#ifdef HAVE_ZLIB_H
              input_handle->libraries->gzclose(input_handle->fp_gz);
              input_handle->fp_gz = nullptr;
              break;
#endif

            case Format::bzip:
#ifdef HAVE_BZLIB_H
              input_handle->libraries->bz_read_close(&bzError, input_handle->fp_bz);
              input_handle->fp_bz = nullptr;
              break;
#endif

            default:
              fatal("Internal error");
            }

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

  buffer_init(& input_handle->header_buffer);
  buffer_init(& input_handle->sequence_buffer);
  buffer_init(& input_handle->plusline_buffer);
  buffer_init(& input_handle->quality_buffer);

  input_handle->stripped_all = 0;

  input_handle->lineno = 1;
  input_handle->lineno_start = 1;
  input_handle->seqno = -1;

  assert(input_handle != nullptr);
  return input_handle.release();
}


auto fastx_is_fastq(struct fastx_s const * input_handle) -> bool
{
  return input_handle->is_fastq or input_handle->is_empty;
}


auto fastx_is_empty(struct fastx_s const * input_handle) -> bool
{
  return input_handle->is_empty;
}


auto fastx_is_pipe(struct fastx_s const * input_handle) -> bool
{
  return input_handle->is_pipe;
}


/* Release the owned resources: the compression handle (if any), the underlying
   file, and the five buffers. Runs on delete, including during stack unwinding
   when fatal() throws in a library session, so it must not throw: it never
   calls fatal() (the "Internal error" default of the old fastx_close switch is
   dropped — format is always plain/gzip/bzip here) and null-guards every member
   so a handle abandoned part-way through fastx_open() is torn down safely. */
fastx_s::~fastx_s()
{
  switch (format)
    {
    case Format::gzip:
#ifdef HAVE_ZLIB_H
      if ((libraries != nullptr) and (fp_gz != nullptr))
        {
          libraries->gzclose(fp_gz);
        }
#endif
      break;

    case Format::bzip:
#ifdef HAVE_BZLIB_H
      if ((libraries != nullptr) and (fp_bz != nullptr))
        {
          int bz_error = 0;
          libraries->bz_read_close(&bz_error, fp_bz);
        }
#endif
      break;

    default:
      break;
    }

  if (fp != nullptr)
    {
      std::fclose(fp);
    }

  buffer_free(& file_buffer);
  buffer_free(& header_buffer);
  buffer_free(& sequence_buffer);
  buffer_free(& plusline_buffer);
  buffer_free(& quality_buffer);
}


auto fastx_close(fastx_handle input_handle, struct Parameters const & parameters) -> void
{
  /* Warn about stripped chars */

  if (input_handle->stripped_all != 0U)
    {
      std::fprintf(stderr, "WARNING: %" PRIu64 " invalid characters stripped from %s file:", input_handle->stripped_all, (input_handle->is_fastq ? "FASTQ" : "FASTA"));
      for (int i = 0; i < 256; i++)
        {
          if (input_handle->stripped[static_cast<std::size_t>(i)] != 0U)
            {
              std::fprintf(stderr, " %c(%" PRIu64 ")", i, input_handle->stripped[static_cast<std::size_t>(i)]);
            }
        }
      std::fprintf(stderr, "\n");
      std::fprintf(stderr, "REMINDER: vsearch does not support amino acid sequences\n");

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log, "WARNING: %" PRIu64 " invalid characters stripped from %s file:", input_handle->stripped_all, (input_handle->is_fastq ? "FASTQ" : "FASTA"));
          for (int i = 0; i < 256; i++)
            {
              if (input_handle->stripped[static_cast<std::size_t>(i)] != 0U)
                {
                  std::fprintf(parameters.fp_log, " %c(%" PRIu64 ")", i, input_handle->stripped[static_cast<std::size_t>(i)]);
                }
            }
          std::fprintf(parameters.fp_log, "\n");
          std::fprintf(parameters.fp_log, "REMINDER: vsearch does not support amino acid sequences\n");
        }
    }

  /* ~fastx_s closes the files and frees the buffers. */
  delete input_handle;
}


auto fastx_file_fill_buffer(fastx_handle input_handle) -> uint64_t
{
  /* read more data if necessary */
  uint64_t const rest = input_handle->file_buffer.length - input_handle->file_buffer.position;

  if (rest > 0)
    {
      return rest;
    }

  uint64_t space = input_handle->file_buffer.alloc - input_handle->file_buffer.length;

  if (space == 0)
    {
      /* back to beginning of buffer */
      input_handle->file_buffer.position = 0;
      input_handle->file_buffer.length = 0;
      space = input_handle->file_buffer.alloc;
    }

  int bytes_read = 0;

#ifdef HAVE_BZLIB_H
  int bzError = 0;
#endif

  switch (input_handle->format)
    {
    case Format::plain:
      bytes_read = static_cast<int>(std::fread(input_handle->file_buffer.data
                         + input_handle->file_buffer.position,
                         1,
                         space,
                         input_handle->fp));
      break;

    case Format::gzip:
#ifdef HAVE_ZLIB_H
      bytes_read = input_handle->libraries->gzread(input_handle->fp_gz,
                               input_handle->file_buffer.data
                               + input_handle->file_buffer.position,
                               static_cast<unsigned int>(space));
      if (bytes_read < 0)
        {
          fatal("Unable to read gzip compressed file");
        }
      break;
#endif

    case Format::bzip:
#ifdef HAVE_BZLIB_H
      bytes_read = input_handle->libraries->bz_read(& bzError,
                                   input_handle->fp_bz,
                                   input_handle->file_buffer.data
                                   + input_handle->file_buffer.position,
                                   static_cast<int>(space));
      if ((bytes_read < 0) or
          not ((bzError == BZ_OK) or
             (bzError == BZ_STREAM_END) or
             (bzError == BZ_SEQUENCE_ERROR)))
        {
          fatal("Unable to read from bzip2 compressed file");
        }
      break;
#endif

    default:
      fatal("Internal error");
    }

  if (not input_handle->is_pipe)
    {
#ifdef HAVE_ZLIB_H
      if (input_handle->format == Format::gzip)
        {
          /* Circumvent the missing gzoffset function in zlib 1.2.3 and earlier */
          int const fd = dup(fileno(input_handle->fp));
          input_handle->file_position = xlseek(fd, 0, SEEK_CUR);
          close(fd);
        }
      else
#endif
        {
          input_handle->file_position = xftello(input_handle->fp);
        }
    }

  input_handle->file_buffer.length += static_cast<uint64_t>(bytes_read);
  return static_cast<uint64_t>(bytes_read);
}


auto fastx_next(fastx_handle input_handle,
                bool const truncateatspace,
                const unsigned char * char_mapping) -> bool
{
  /* deferred-error mode (see fastx.h): if a previous call already
     recorded a parse error, report no further records so every worker
     stops; and if the current record triggered a deferred error, treat
     it as unusable (return false) rather than handing back a bogus record */
  if (input_handle->error)
    {
      return false;
    }
  bool const got_record = input_handle->is_fastq
    ? fastq_next(input_handle, truncateatspace, char_mapping)
    : fasta_next(input_handle, truncateatspace, char_mapping);
  if (input_handle->error)
    {
      return false;
    }
  return got_record;
}


auto fastx_get_error(struct fastx_s const * input_handle) -> bool
{
  return input_handle->error;
}


auto fastx_get_errmsg(struct fastx_s const * input_handle) -> char const *
{
  return input_handle->errmsg.data();
}


auto fastx_set_deferred_error(fastx_handle input_handle, char const * message) -> void
{
  /* record the first deferred parse error and flag the handle (see the
     deferred-error note in fastx.h). First error wins; later ones are
     ignored so the message reflects the earliest failure. */
  if (not input_handle->error)
    {
      std::snprintf(input_handle->errmsg.data(), input_handle->errmsg.size(), "%s", message);
      input_handle->error = true;
    }
}


auto fastx_filter_sequence_length(fastx_handle input_handle) -> void
{
  /* Reject a sequence too long for the int sequence-length bookkeeping used
     downstream (e.g. searchinfo_s::qseqlen, sized as qseqlen + buffer_headroom;
     the chimera query buffers; cut's rc_buffer). Such a sequence would be
     narrowed to a negative int and overflow the per-record buffer. Database
     sequences over --maxseqlength are discarded by Database::read, but many commands
     read records directly via fasta_next/fastx_next with no length filter, so
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
      fastx_set_deferred_error(input_handle, message.data());
      return;
    }
  fatal(message.data());
}


auto fastx_get_position(struct fastx_s const * input_handle) -> uint64_t
{
  if (input_handle->is_fastq)
    {
      return fastq_get_position(input_handle);
    }
  return fasta_get_position(input_handle);
}


auto fastx_get_size(struct fastx_s const * input_handle) -> uint64_t
{
  if (input_handle->is_fastq)
    {
      return fastq_get_size(input_handle);
    }
  return fasta_get_size(input_handle);
}


auto fastx_get_lineno(struct fastx_s const * input_handle) -> uint64_t
{
  if (input_handle->is_fastq)
    {
      return fastq_get_lineno(input_handle);
    }
  return fasta_get_lineno(input_handle);
}


auto fastx_get_seqno(struct fastx_s const * input_handle) -> uint64_t
{
  if (input_handle->is_fastq)
    {
      return fastq_get_seqno(input_handle);
    }
  return fasta_get_seqno(input_handle);
}


auto fastx_get_header(struct fastx_s const * input_handle) -> char const *
{
  if (input_handle->is_fastq)
    {
      return fastq_get_header(input_handle);
    }
  return fasta_get_header(input_handle);
}


auto fastx_get_sequence(struct fastx_s const * input_handle) -> char const *
{
  if (input_handle->is_fastq)
    {
      return fastq_get_sequence(input_handle);
    }
  return fasta_get_sequence(input_handle);
}


auto fastx_get_header_length(struct fastx_s const * input_handle) -> uint64_t
{
  if (input_handle->is_fastq)
    {
      return fastq_get_header_length(input_handle);
    }
  return fasta_get_header_length(input_handle);
}


auto fastx_get_sequence_length(struct fastx_s const * input_handle) -> uint64_t
{
  if (input_handle->is_fastq)
    {
      return fastq_get_sequence_length(input_handle);
    }
  return fasta_get_sequence_length(input_handle);
}


auto fastx_get_quality(struct fastx_s const * input_handle) -> char const *
{
  if (input_handle->is_fastq)
    {
      return fastq_get_quality(input_handle);
    }
  return nullptr;
}


auto fastx_get_abundance(struct fastx_s const * input_handle) -> int64_t
{
  if (input_handle->is_fastq)
    {
      return fastq_get_abundance(input_handle);
    }
  return fasta_get_abundance(input_handle);
}


auto fastx_record(fastx_handle input_handle) -> SeqRecord
{
  auto const sequence_length = fastx_get_sequence_length(input_handle);
  auto const * quality = fastx_get_quality(input_handle);  // nullptr for FASTA
  return SeqRecord{
    View<char>{fastx_get_header(input_handle), fastx_get_header_length(input_handle)},
    View<char>{fastx_get_sequence(input_handle), sequence_length},
    View<char>{quality, quality != nullptr ? sequence_length : uint64_t{0}}};
}
