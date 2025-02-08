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

#include "vsearch.h"
#include "attributes.h"
#include "maps.h"
#include <array>
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::snprintf, std::size_t
#include <cstring>  // std::memcmp, std::memchr, std::strlen


auto fastq_fatal(uint64_t lineno, const char * msg) -> void
{
  char * string = nullptr;
  if (xsprintf(&string,
               "Invalid line %lu in FASTQ file: %s",
               lineno,
               msg) == -1)
    {
      fatal("Out of memory");
    }

  if (string != nullptr)
    {
      fatal(string);
      xfree(string);
    }
  else
    {
      fatal("Out of memory");
    }
}


auto buffer_filter_extend(fastx_handle input_handle,
                          struct fastx_buffer_s * dest_buffer,
                          char * source_buf,
                          uint64_t len,
                          unsigned int * char_action,
                          const unsigned char * char_mapping,
                          bool * ok,
                          char * illegal_char) -> void
{
  buffer_makespace(dest_buffer, len + 1);

  /* Strip unwanted characters from the string and raise warnings or
     errors on certain characters. */

  auto * p = source_buf;
  auto * d = dest_buffer->data + dest_buffer->length;
  auto * q = d;
  *ok = true;

  for (auto i = 0ULL; i < len; i++)
    {
      auto const c = *p++;
      char const m = char_action[(unsigned char) (c)];

      switch(m)
        {
        case 0:
          /* stripped */
          input_handle->stripped_all++;
          input_handle->stripped[(unsigned char) (c)]++;
          break;

        case 1:
          /* legal character */
          *q++ = char_mapping[(unsigned char) (c)];
          break;

        case 2:
          /* fatal character */
          if (*ok)
            {
              *illegal_char = c;
            }
          *ok = false;
          break;

        case 3:
          /* silently stripped chars (whitespace) */
          break;

        case 4:
          /* newline (silently stripped) */
          break;
        }
    }

  /* add zero after sequence */
  *q = 0;
  dest_buffer->length += q - d;
}


auto fastq_open(const char * filename) -> fastx_handle
{
  auto * input_handle = fastx_open(filename);

  if (! fastx_is_fastq(input_handle))
    {
      fatal("FASTQ file expected, FASTA file found (%s)", filename);
    }

  return input_handle;
}


auto fastq_close(fastx_handle input_handle) -> void
{
  fastx_close(input_handle);
}


auto fastq_next(fastx_handle input_handle,
                bool const truncateatspace,
                const unsigned char * char_mapping) -> bool
{
  input_handle->header_buffer.length = 0;
  input_handle->header_buffer.data[0] = 0;
  input_handle->sequence_buffer.length = 0;
  input_handle->sequence_buffer.data[0] = 0;
  input_handle->plusline_buffer.length = 0;
  input_handle->plusline_buffer.data[0] = 0;
  input_handle->quality_buffer.length = 0;
  input_handle->quality_buffer.data[0] = 0;

  input_handle->lineno_start = input_handle->lineno;

  static constexpr auto max_message_length = std::size_t{200};
  std::array<char, max_message_length> message {{}};
  auto ok = true;
  char illegal_char = '\0';

  auto rest = fastx_file_fill_buffer(input_handle);

  /* check end of file */

  if (rest == 0)
    {
      return false;
    }

  /* read header */

  /* check initial @ character */

  if (input_handle->file_buffer.data[input_handle->file_buffer.position] != '@')
    {
      fastq_fatal(input_handle->lineno, "Header line must start with '@' character");
    }
  input_handle->file_buffer.position++;
  rest--;

  char * lf = nullptr;
  while (lf == nullptr)
    {
      /* get more data if buffer empty */
      rest = fastx_file_fill_buffer(input_handle);
      if (rest == 0)
        {
          fastq_fatal(input_handle->lineno, "Unexpected end of file");
        }

      /* find LF */
      lf = (char *) memchr(input_handle->file_buffer.data + input_handle->file_buffer.position,
                           '\n',
                           rest);

      /* copy to header buffer */
      auto len = rest;
      if (lf != nullptr)
        {
          /* LF found, copy up to and including LF */
          len = lf - (input_handle->file_buffer.data + input_handle->file_buffer.position) + 1;
          input_handle->lineno++;
        }
      buffer_extend(&input_handle->header_buffer,
                    input_handle->file_buffer.data + input_handle->file_buffer.position,
                    len);
      input_handle->file_buffer.position += len;
      rest -= len;
    }

  /* read sequence line(s) */
  lf = nullptr;
  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(input_handle);

      /* cannot end here */
      if (rest == 0)
        {
          fastq_fatal(input_handle->lineno, "Unexpected end of file");
        }

      /* end when new line starting with + is seen */
      if ((lf != nullptr) && (input_handle->file_buffer.data[input_handle->file_buffer.position] == '+'))
        {
          break;
        }

      /* find LF */
      lf = (char *) memchr(input_handle->file_buffer.data + input_handle->file_buffer.position,
                           '\n', rest);

      /* copy to sequence buffer */
      auto len = rest;
      if (lf != nullptr)
        {
          /* LF found, copy up to and including LF */
          len = lf - (input_handle->file_buffer.data + input_handle->file_buffer.position) + 1;
          input_handle->lineno++;
        }

      buffer_filter_extend(input_handle,
                           &input_handle->sequence_buffer,
                           input_handle->file_buffer.data + input_handle->file_buffer.position,
                           len,
                           char_fq_action_seq, char_mapping,
                           &ok, &illegal_char);
      input_handle->file_buffer.position += len;
      rest -= len;

      if (! ok)
        {
          if ((illegal_char >= 32) && (illegal_char < 127))
            {
              snprintf(message.data(),
                       max_message_length,
                       "Illegal sequence character '%c'",
                       illegal_char);
            }
          else
            {
              snprintf(message.data(),
                       max_message_length,
                       "Illegal sequence character (unprintable, no %d)",
                       (unsigned char) illegal_char);
            }
          fastq_fatal(input_handle->lineno - ((lf != nullptr) ? 1 : 0), message.data());
        }
    }

  /* read + line */

  /* skip + character */
  input_handle->file_buffer.position++;
  rest--;

  lf = nullptr;
  while (lf == nullptr)
    {
      /* get more data if buffer empty */
      rest = fastx_file_fill_buffer(input_handle);

      /* cannot end here */
      if (rest == 0)
        {
          fastq_fatal(input_handle->lineno, "Unexpected end of file");
        }

      /* find LF */
      lf = (char *) memchr(input_handle->file_buffer.data + input_handle->file_buffer.position,
                           '\n',
                           rest);
      /* copy to plusline buffer */
      auto len = rest;
      if (lf != nullptr)
        {
          /* LF found, copy up to and including LF */
          len = lf - (input_handle->file_buffer.data + input_handle->file_buffer.position) + 1;
          input_handle->lineno++;
        }
      buffer_extend(&input_handle->plusline_buffer,
                    input_handle->file_buffer.data + input_handle->file_buffer.position,
                    len);
      input_handle->file_buffer.position += len;
      rest -= len;
    }

  /* check that the plus line is empty or identical to @ line */

  bool plusline_invalid = false;
  if (input_handle->header_buffer.length == input_handle->plusline_buffer.length)
    {
      if ((memcmp(input_handle->header_buffer.data,
                 input_handle->plusline_buffer.data,
                  input_handle->header_buffer.length) != 0))
        {
          plusline_invalid = true;
        }
    }
  else
    {
      if ((input_handle->plusline_buffer.length > 2) ||
          ((input_handle->plusline_buffer.length == 2) && (input_handle->plusline_buffer.data[0] != '\r')))
        {
          plusline_invalid = true;
        }
    }
  if (plusline_invalid)
    {
      fastq_fatal(input_handle->lineno - ((lf != nullptr) ? 1 : 0),
                  "'+' line must be empty or identical to header");
    }

  /* read quality line(s) */

  lf = nullptr;
  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(input_handle);

      /* end if no more data */
      if (rest == 0)
        {
          break;
        }

      /* end if next entry starts : LF + '@' + correct length */
      if ((lf != nullptr) &&
          (input_handle->file_buffer.data[input_handle->file_buffer.position] == '@') &&
          (input_handle->quality_buffer.length == input_handle->sequence_buffer.length))
        {
          break;
        }

      /* find LF */
      lf = (char *) memchr(input_handle->file_buffer.data + input_handle->file_buffer.position,
                           '\n', rest);

      /* copy to quality buffer */
      auto len = rest;
      if (lf != nullptr)
        {
          /* LF found, copy up to and including LF */
          len = lf - (input_handle->file_buffer.data + input_handle->file_buffer.position) + 1;
          input_handle->lineno++;
        }

      buffer_filter_extend(input_handle,
                           &input_handle->quality_buffer,
                           input_handle->file_buffer.data + input_handle->file_buffer.position,
                           len,
                           char_fq_action_qual, chrmap_identity,
                           &ok, &illegal_char);
      input_handle->file_buffer.position += len;
      rest -= len;

      /* break if quality line already too long */
      if (input_handle->quality_buffer.length > input_handle->sequence_buffer.length)
        {
          break;
        }

      if (! ok)
        {
          if ((illegal_char >= 32) && (illegal_char < 127))
            {
              snprintf(message.data(),
                       max_message_length,
                       "Illegal quality character '%c'",
                       illegal_char);
            }
          else
            {
              snprintf(message.data(),
                       max_message_length,
                       "Illegal quality character (unprintable, no %d)",
                       (unsigned char) illegal_char);
            }
          fastq_fatal(input_handle->lineno - ((lf != nullptr) ? 1 : 0), message.data());
        }
    }

  if (input_handle->sequence_buffer.length != input_handle->quality_buffer.length)
    {
      fastq_fatal(input_handle->lineno - ((lf != nullptr) ? 1 : 0),
                  "Sequence and quality lines must be equally long");
    }

  fastx_filter_header(input_handle, truncateatspace);

  input_handle->seqno++;

  return true;
}


auto fastq_get_quality(fastx_handle input_handle) -> char *
{
  return input_handle->quality_buffer.data;
}


auto fastq_get_quality_length(fastx_handle input_handle) -> uint64_t
{
  return input_handle->quality_buffer.length;
}


auto fastq_get_position(fastx_handle input_handle) -> uint64_t
{
  return input_handle->file_position;
}


auto fastq_get_size(fastx_handle input_handle) -> uint64_t
{
  return input_handle->file_size;
}


auto fastq_get_lineno(fastx_handle input_handle) -> uint64_t
{
  return input_handle->lineno_start;
}


auto fastq_get_seqno(fastx_handle input_handle) -> uint64_t
{
  return input_handle->seqno;
}


auto fastq_get_header_length(fastx_handle input_handle) -> uint64_t
{
  return input_handle->header_buffer.length;
}


auto fastq_get_sequence_length(fastx_handle input_handle) -> uint64_t
{
  return input_handle->sequence_buffer.length;
}


auto fastq_get_header(fastx_handle input_handle) -> char *
{
  return input_handle->header_buffer.data;
}


auto fastq_get_sequence(fastx_handle input_handle) -> char *
{
  return input_handle->sequence_buffer.data;
}


auto fastq_get_abundance(fastx_handle input_handle) -> int64_t
{
  // return 1 if not present
  auto const size = header_get_size(input_handle->header_buffer.data,
                                 input_handle->header_buffer.length);
  if (size > 0)
    {
      return size;
    }
  return 1;
}


auto fastq_get_abundance_and_presence(fastx_handle input_handle) -> int64_t
{
  // return 0 if not present
  return header_get_size(input_handle->header_buffer.data, input_handle->header_buffer.length);
}


inline auto fprint_seq_label(std::FILE * output_handle, char * seq, int len) -> void
{
  /* normalize first? */
  std::fprintf(output_handle, "%.*s", len, seq);
}


auto fastq_print_general(FILE * output_handle,
                         char * seq,
                         int len,
                         char * header,
                         int header_len,
                         char * quality,
                         int abundance,
                         int ordinal,
                         double ee) -> void
{
  std::fprintf(output_handle, "@");

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
  else if ((opt_relabel != nullptr) && (ordinal > 0))
    {
      std::fprintf(output_handle, "%s%d", opt_relabel, ordinal);
    }
  else
    {
      auto const xsize = opt_xsize || (opt_sizeout && (abundance > 0));
      auto const xee = opt_xee || ((opt_eeout || opt_fastq_eeout) && (ee >= 0.0));
      auto const xlength = opt_xlength || opt_lengthout;
      header_fprint_strip(output_handle,
                          header,
                          header_len,
                          xsize,
                          xee,
                          xlength);
    }

  if (opt_label_suffix != nullptr)
    {
      std::fprintf(output_handle, "%s", opt_label_suffix);
    }

  if (opt_sample != nullptr)
    {
      std::fprintf(output_handle, ";sample=%s", opt_sample);
    }

  if (opt_sizeout && (abundance > 0))
    {
      std::fprintf(output_handle, ";size=%u", abundance);
    }

  if ((opt_eeout || opt_fastq_eeout) && (ee >= 0.0))
    {
      if (ee < 0.000000001) {
        std::fprintf(output_handle, ";ee=%.13lf", ee);
      } else if (ee < 0.00000001) {
        std::fprintf(output_handle, ";ee=%.12lf", ee);
      } else if (ee < 0.0000001) {
        std::fprintf(output_handle, ";ee=%.11lf", ee);
      } else if (ee < 0.000001) {
        std::fprintf(output_handle, ";ee=%.10lf", ee);
      } else if (ee < 0.00001) {
        std::fprintf(output_handle, ";ee=%.9lf", ee);
      } else if (ee < 0.0001) {
        std::fprintf(output_handle, ";ee=%.8lf", ee);
      } else if (ee < 0.001) {
        std::fprintf(output_handle, ";ee=%.7lf", ee);
      } else if (ee < 0.01) {
        std::fprintf(output_handle, ";ee=%.6lf", ee);
      } else if (ee < 0.1) {
        std::fprintf(output_handle, ";ee=%.5lf", ee);
      } else {
        std::fprintf(output_handle, ";ee=%.4lf", ee);
      }
    }

  if (opt_lengthout)
    {
      std::fprintf(output_handle, ";length=%d", len);
    }

  if (opt_relabel_keep &&
      (((opt_relabel != nullptr) && (ordinal > 0)) || opt_relabel_sha1 || opt_relabel_md5 || opt_relabel_self))
    {
      std::fprintf(output_handle, " %.*s", header_len, header);
    }

  std::fprintf(output_handle, "\n%.*s\n+\n%.*s\n", len, seq, len, quality);
}


auto fastq_print(std::FILE * output_handle, char * header, char * sequence, char * quality) -> void
{
  auto const slen = static_cast<int>(std::strlen(sequence));
  auto const hlen = static_cast<int>(std::strlen(header));
  fastq_print_general(output_handle, sequence, slen, header, hlen, quality, 0, 0, -1.0);
}
