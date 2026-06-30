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

#include "vsearch.h"
#include "attributes.h"
#include "utils/fatal.hpp"
#include <algorithm>  // std::min
#include <array>
#include <cassert>  // assert
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio> // std::FILE, std::fprintf, std::size_t, std::snprintf
#include <cstring>  // std::memchr, std::strlen
#include <iterator>  // std::next
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {


  enum struct Action : unsigned char {
    warn,    // (0) symbol is stripped, with a warning
    accept,  // (1)
    reject,  // (2) fatal printable symbol ('.', '-')
    show,    // (3) fatal non-printable symbol (0-32, but not 127?)
    skip,    // (4) symbol is stripped, silently
    count    // (5) track the number of lines
  };


  /*
    How to handle input characters for FASTA

    0=warn, 1=accept, 2=reject, 3=show, 4=skip, 5=count

    3,  3,  3,  3,  3,  3,  3,  3,  3,  4,  5,  4,  4,  4,  3,  3,    // 0-15
    3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,    // 16-31
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,  0,    // 32-47
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,    // 48-63
    0,  1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,    // 64-79
    0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,    // 80-95
    0,  1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,    // 96-111
    0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,    // 112-127
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
  */
  const std::vector<Action> char_actions =
    {
      Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::skip,  Action::count,  Action::skip,  Action::skip,  Action::skip,  Action::show,  Action::show,  // 0-15
      Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  Action::show,  // 16-31
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::reject,  Action::reject,  Action::warn,  // 32-47
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  // 48-63
      Action::warn,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::warn,  Action::warn,  Action::accept,  Action::accept,  Action::warn,  Action::warn,  Action::accept,  Action::warn,  Action::accept,  Action::accept,  Action::warn,  // 64-79
      Action::warn,  Action::warn,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::warn,  Action::accept,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  // 80-95
      Action::warn,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::warn,  Action::warn,  Action::accept,  Action::accept,  Action::warn,  Action::warn,  Action::accept,  Action::warn,  Action::accept,  Action::accept,  Action::warn,  // 96-111
      Action::warn,  Action::warn,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::accept,  Action::warn,  Action::accept,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  // 112-127
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,
      Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn,  Action::warn
    };


  auto map_action(char const nucleotide) -> Action {
    auto const current_char = static_cast<unsigned char>(nucleotide);
    return char_actions[current_char];
  }


  auto report_illegal_symbol_and_exit(fastx_handle input_handle, unsigned char symbol, uint64_t line_number) -> void {
    static constexpr std::size_t max_buffer_size = 200;
    std::array<char, max_buffer_size> msg {{}};
    static_cast<void>(std::snprintf(
        msg.data(), max_buffer_size,
        "Illegal character '%c' in sequence on line %" PRIu64 " of FASTA file",
        symbol,
        line_number));
    /* deferred-error mode (see fastx.h): record and return instead of
       exiting, so a worker thread does not std::exit() with siblings live */
    if (input_handle->defer_errors) {
      fastx_set_deferred_error(input_handle, msg.data());
      return;
    }
    fatal(msg.data());
  }


  auto report_unprintable_symbol_and_exit(fastx_handle input_handle, unsigned char symbol, uint64_t line_number) -> void {
    static constexpr std::size_t max_buffer_size = 200;
    std::array<char, max_buffer_size> msg {{}};
    static_cast<void>(std::snprintf(
        msg.data(), max_buffer_size,
        "Illegal unprintable ASCII character no %d in sequence on line %" PRIu64 " of FASTA file",
        symbol,
        line_number));
    if (input_handle->defer_errors) {
      fastx_set_deferred_error(input_handle, msg.data());
      return;
    }
    fatal(msg.data());
  }

}  // end of anonymous namespace


auto fasta_open(const char * filename) -> fastx_handle
{
  auto * input_handle = fastx_open(filename);

  if (fastx_is_fastq(input_handle) and not input_handle->is_empty)
    {
      fatal("FASTA file expected, FASTQ file found (%s)", filename);
    }

  assert(input_handle != nullptr);
  return input_handle;
}


auto fasta_close(fastx_handle input_handle) -> void
{
  fastx_close(input_handle);
}


auto fasta_filter_sequence(fastx_handle input_handle,
                           unsigned char const * char_mapping) -> void
{
  /* Strip unwanted characters from the sequence and raise warnings or
     errors on certain characters. */

  auto * source = input_handle->sequence_buffer.data;
  auto * dest = source;

  while (*source != '\0')
    {
      auto const current_char = static_cast<unsigned char>(*source);

      switch (map_action(*source))
        {
        case Action::warn:
          /* stripped */
          ++input_handle->stripped_all;
          ++input_handle->stripped[current_char];
          break;

        case Action::accept:
          /* legal character */
          *dest = static_cast<char>(char_mapping[current_char]);
          dest = std::next(dest);
          break;

        case Action::reject:
          /* fatal character */
          report_illegal_symbol_and_exit(input_handle, current_char, input_handle->lineno);
          if (input_handle->error) { return; }
          break;

        case Action::show:
          /* fatal unprintable character */
          report_unprintable_symbol_and_exit(input_handle, current_char, input_handle->lineno);
          if (input_handle->error) { return; }
          break;

        case Action::skip:
          /* silently stripped chars (whitespace) */
          break;

        case Action::count:
          /* newline (silently stripped) */
          ++input_handle->lineno;
          break;
        }
      source = std::next(source);
    }

  /* add nullchar after sequence */
  *dest = '\0';
  input_handle->sequence_buffer.length = static_cast<uint64_t>(dest - input_handle->sequence_buffer.data);
}


auto fasta_next(fastx_handle input_handle,
                bool truncateatspace,
                const unsigned char * char_mapping) -> bool
{
  input_handle->lineno_start = input_handle->lineno;

  input_handle->header_buffer.length = 0;
  input_handle->header_buffer.data[0] = 0;
  input_handle->sequence_buffer.length = 0;
  input_handle->sequence_buffer.data[0] = 0;

  std::size_t rest = fastx_file_fill_buffer(input_handle);

  if (rest == 0)
    {
      return false;
    }

  /* read header */

  /* check initial > character */

  if (input_handle->file_buffer.data[input_handle->file_buffer.position] != '>')
    {
      if (input_handle->defer_errors)
        {
          fastx_set_deferred_error(input_handle, "Invalid FASTA - header must start with > character");
          return false;
        }
      std::fprintf(stderr, "Found character %02x\n", static_cast<unsigned char>(input_handle->file_buffer.data[input_handle->file_buffer.position]));
      fatal("Invalid FASTA - header must start with > character");
    }
  ++input_handle->file_buffer.position;
  --rest;

  char const * line_end = nullptr;
  while (line_end == nullptr)
    {
      /* get more data if buffer empty*/
      rest = fastx_file_fill_buffer(input_handle);
      if (rest == 0)
        {
          if (input_handle->defer_errors)
            {
              fastx_set_deferred_error(input_handle, "Invalid FASTA - header must be terminated with newline");
              return false;
            }
          fatal("Invalid FASTA - header must be terminated with newline");
        }

      /* find new line char ('LF') */
      auto * const current_position = std::next(input_handle->file_buffer.data, static_cast<long>(input_handle->file_buffer.position));
      line_end = static_cast<char *>(std::memchr(current_position,
                                                 '\n',
                                                 rest));

      /* copy to header buffer */
      uint64_t len = rest;
      if (line_end != nullptr)
        {
          /* LF found, copy up to and including LF */
          len = static_cast<uint64_t>(line_end - (input_handle->file_buffer.data + input_handle->file_buffer.position) + 1);
          ++input_handle->lineno;
        }
      buffer_extend(& input_handle->header_buffer,
                    input_handle->file_buffer.data + input_handle->file_buffer.position,
                    len);
      input_handle->file_buffer.position += len;
      rest -= len;
    }

  /* read one or more sequence lines */

  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(input_handle);

      /* end if no more data */
      if (rest == 0)
        {
          break;
        }

      /* end if new sequence starts */
      if ((line_end != nullptr) and (input_handle->file_buffer.data[input_handle->file_buffer.position] == '>'))
        {
          break;
        }

      /* find LF */
      auto * const current_position = std::next(input_handle->file_buffer.data, static_cast<std::ptrdiff_t>(input_handle->file_buffer.position));
      line_end = static_cast<char *>(std::memchr(current_position, '\n', rest));

      uint64_t len = rest;
      if (line_end != nullptr)
        {
          /* LF found, copy up to and including LF */
          len = static_cast<uint64_t>(line_end - current_position + 1);
        }
      buffer_extend(& input_handle->sequence_buffer,
                    current_position,
                    len);
      input_handle->file_buffer.position += len;
      rest -= len;
    }

  ++input_handle->seqno;

  fastx_filter_header(input_handle, truncateatspace);
  fasta_filter_sequence(input_handle, char_mapping);

  return true;
}


auto fasta_get_abundance(struct fastx_s const * input_handle) -> int64_t
{
  // return 1 if not present
  auto const size = header_get_size(input_handle->header_buffer.data,
                                 static_cast<int>(input_handle->header_buffer.length));
  if (size > 0)
    {
      return size;
    }
  return 1;
}


auto fasta_get_abundance_and_presence(struct fastx_s const * input_handle) -> int64_t
{
  // return 0 if not present
  return header_get_size(input_handle->header_buffer.data, static_cast<int>(input_handle->header_buffer.length));
}


auto fasta_get_position(struct fastx_s const * input_handle) -> uint64_t
{
  return input_handle->file_position;
}


auto fasta_get_size(struct fastx_s const * input_handle) -> uint64_t
{
  return input_handle->file_size;
}


auto fasta_get_lineno(struct fastx_s const * input_handle) -> uint64_t
{
  return input_handle->lineno_start;
}


auto fasta_get_seqno(struct fastx_s const * input_handle) -> uint64_t
{
  return static_cast<uint64_t>(input_handle->seqno);
}


auto fasta_get_header_length(fastx_handle input_handle) -> uint64_t
{
  return input_handle->header_buffer.length;
}


auto fasta_get_sequence_length(fastx_handle input_handle) -> uint64_t
{
  return input_handle->sequence_buffer.length;
}


auto fasta_get_header(fastx_handle input_handle) -> char const *
{
  return input_handle->header_buffer.data;
}


auto fasta_get_sequence(fastx_handle input_handle) -> char const *
{
  return input_handle->sequence_buffer.data;
}


/* fasta output */

auto fasta_print_sequence(std::FILE * output_handle, char const * seq, uint64_t const len, int const width) -> void
{
  /*
    The actual length of the sequence may be longer than "len", but only
    "len" characters are printed.

    Specify width of lines - zero (or <1) means linearize (all on one line).
  */

  if (width < 1)  // no sequence folding
    {
      /* len may exceed INT_MAX, which the "%.*s" precision (an int) cannot
         express, so write the bytes directly rather than through fprintf. */
      std::fwrite(seq, 1, len, output_handle);
      std::fprintf(output_handle, "\n");
    }
  else  // sequence folding every 'width'
    {
      auto const width_u = static_cast<uint64_t>(width);
      for (uint64_t i = 0; i < len; i += width_u)
        {
          /* each chunk is at most 'width' (an int), so the cast is safe even
             when the whole sequence is longer than INT_MAX */
          auto const chunk = static_cast<int>(std::min(len - i, width_u));
          std::fprintf(output_handle, "%.*s\n", chunk, seq + i);
        }
    }
}


auto fasta_print(std::FILE * output_handle, char const * header,
                 char const * seq, uint64_t const len) -> void
{
  std::fprintf(output_handle, ">%s\n", header);
  fasta_print_sequence(output_handle, seq, len, static_cast<int>(opt_fasta_width));
}


inline auto fprint_seq_label(std::FILE * output_handle, char const * seq, int const len) -> void
{
  /* normalize first? */
  std::fprintf(output_handle, "%.*s", len, seq);
}


auto fasta_print_general(std::FILE * output_handle,
                         char const * prefix,
                         char const * seq,
                         int const len,
                         char const * header,
                         int const header_length,
                         uint64_t const abundance,
                         int64_t const ordinal,
                         double const expected_error,
                         int64_t const clustersize,
                         int const clusterid,
                         char const * score_name,
                         double const score,
                         uint64_t const centroid_size) -> void
{
  std::fprintf(output_handle, ">");

  if (prefix != nullptr)
    {
      std::fprintf(output_handle, "%s", prefix);
    }

  // track whether the text printed so far ends with the annotation
  // separator ';', so that appended annotations are merged with a single
  // separator instead of producing ";;" (see issue #271)
  auto trailing_separator = false;

  if (opt_relabel_self)
    {
      fprint_seq_label(output_handle, seq, len);
    }
  else if (opt_relabel_sha1)
    {
      fprint_seq_digest_sha1(output_handle, seq, len);
    }
  else if (opt_relabel_md5)
    {
      fprint_seq_digest_md5(output_handle, seq, len);
    }
  else if ((opt_relabel != nullptr) and (ordinal > 0))
    {
      std::fprintf(output_handle, "%s%" PRId64, opt_relabel, ordinal);
    }
  else
    {
      bool const strip_size = opt_xsize or (opt_sizeout and (abundance > 0));
      bool const strip_ee = opt_xee or ((opt_eeout or opt_fastq_eeout) and (expected_error >= 0.0));
      bool const strip_length = opt_xlength or opt_lengthout;
      trailing_separator = header_fprint_strip(output_handle,
                                               header,
                                               header_length,
                                               strip_size,
                                               strip_ee,
                                               strip_length);
    }

  if (opt_label_suffix != nullptr)
    {
      std::fprintf(output_handle, "%s", opt_label_suffix);
      if (*opt_label_suffix != '\0')
        {
          trailing_separator = (opt_label_suffix[std::strlen(opt_label_suffix) - 1] == ';');
        }
    }

  if (opt_sample != nullptr)
    {
      std::fprintf(output_handle, "%ssample=%s", annotation_separator(trailing_separator), opt_sample);
    }

  if (clustersize > 0)
    {
      std::fprintf(output_handle, "%sseqs=%" PRId64, annotation_separator(trailing_separator), clustersize);
    }

  if (clusterid >= 0)
    {
      std::fprintf(output_handle, "%sclusterid=%d", annotation_separator(trailing_separator), clusterid);
    }

  if (opt_sizeout and (abundance > 0))
    {
      std::fprintf(output_handle, "%ssize=%" PRIu64, annotation_separator(trailing_separator), abundance);
    }

  if (opt_centroid_sizeout and (centroid_size > 0))
    {
      std::fprintf(output_handle, "%scentroid_size=%" PRIu64, annotation_separator(trailing_separator), centroid_size);
    }

  if ((opt_eeout or opt_fastq_eeout) and (expected_error >= 0.0))
    {
      auto const * separator = annotation_separator(trailing_separator);
      if (expected_error < 0.000000001) {
        std::fprintf(output_handle, "%see=%.13lf", separator, expected_error);
      } else if (expected_error < 0.00000001) {
        std::fprintf(output_handle, "%see=%.12lf", separator, expected_error);
      } else if (expected_error < 0.0000001) {
        std::fprintf(output_handle, "%see=%.11lf", separator, expected_error);
      } else if (expected_error < 0.000001) {
        std::fprintf(output_handle, "%see=%.10lf", separator, expected_error);
      } else if (expected_error < 0.00001) {
        std::fprintf(output_handle, "%see=%.9lf", separator, expected_error);
      } else if (expected_error < 0.0001) {
        std::fprintf(output_handle, "%see=%.8lf", separator, expected_error);
      } else if (expected_error < 0.001) {
        std::fprintf(output_handle, "%see=%.7lf", separator, expected_error);
      } else if (expected_error < 0.01) {
        std::fprintf(output_handle, "%see=%.6lf", separator, expected_error);
      } else if (expected_error < 0.1) {
        std::fprintf(output_handle, "%see=%.5lf", separator, expected_error);
      } else {
        std::fprintf(output_handle, "%see=%.4lf", separator, expected_error);
      }
    }

  if (opt_lengthout)
    {
      std::fprintf(output_handle, "%slength=%d", annotation_separator(trailing_separator), len);
    }

  if (score_name != nullptr)
    {
      std::fprintf(output_handle, "%s%s=%.4lf", annotation_separator(trailing_separator), score_name, score);
    }

  if (opt_relabel_keep and
      (((opt_relabel != nullptr) and (ordinal > 0)) or opt_relabel_sha1 or opt_relabel_md5 or opt_relabel_self))
    {
      std::fprintf(output_handle, " %s", header);
    }

  std::fprintf(output_handle, "\n");

  if (seq != nullptr)
    {
      fasta_print_sequence(output_handle, seq, static_cast<uint64_t>(len), static_cast<int>(opt_fasta_width));
    }
}


// A single uint64_t ordinal parameter: it is the widest unsigned type, so
// every caller (passing int, size_t or uint64_t, all non-negative 1-based
// counters) converts without narrowing or sign-change. Two overloads taking
// int and size_t were ambiguous for a uint64_t argument on platforms where
// uint64_t, size_t and int are all distinct types (e.g. macOS).
auto fasta_print_db_relabel(std::FILE * output_handle,
                            uint64_t seqno,
                            uint64_t ordinal) -> void
{
  fasta_print_general(output_handle,
                      nullptr,
                      db_getsequence(seqno),
                      static_cast<int>(db_getsequencelen(seqno)),
                      db_getheader(seqno),
                      static_cast<int>(db_getheaderlen(seqno)),
                      db_getabundance(seqno),
                      static_cast<int64_t>(ordinal),
                      -1.0,
                      -1, -1,
                      nullptr, 0.0,
                      0);
}


auto fasta_print_db(std::FILE * output_handle, uint64_t seqno) -> void
{
  fasta_print_general(output_handle,
                      nullptr,
                      db_getsequence(seqno),
                      static_cast<int>(db_getsequencelen(seqno)),
                      db_getheader(seqno),
                      static_cast<int>(db_getheaderlen(seqno)),
                      db_getabundance(seqno),
                      0,
                      -1.0,
                      -1, -1,
                      nullptr, 0.0,
                      0);
}
