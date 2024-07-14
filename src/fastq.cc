/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::snprintf
#include <cstring>  // std::memcmp, std::memchr, std::strlen


auto fastq_fatal(uint64_t lineno, const char * msg) -> void
{
  char * string = nullptr;
  if (xsprintf(& string,
               "Invalid line %lu in FASTQ file: %s",
               lineno,
               msg) == -1)
    {
      fatal("Out of memory");
    }

  if (string)
    {
      fatal(string);
      xfree(string);
    }
  else
    {
      fatal("Out of memory");
    }
}


auto buffer_filter_extend(fastx_handle h,
                          struct fastx_buffer_s * dest_buffer,
                          char * source_buf,
                          uint64_t len,
                          unsigned int * char_action,
                          const unsigned char * char_mapping,
                          bool * ok,
                          char * illegal_char) -> void
{
  buffer_makespace(dest_buffer, len+1);

  /* Strip unwanted characters from the string and raise warnings or
     errors on certain characters. */

  char * p = source_buf;
  char * d = dest_buffer->data + dest_buffer->length;
  char * q = d;
  *ok = true;

  for (uint64_t i = 0; i < len; i++)
    {
      char const c = *p++;
      char const m = char_action[(unsigned char) c];

      switch(m)
        {
        case 0:
          /* stripped */
          h->stripped_all++;
          h->stripped[(unsigned char) c]++;
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
  fastx_handle h = fastx_open(filename);

  if (! fastx_is_fastq(h))
    {
      fatal("FASTQ file expected, FASTA file found (%s)", filename);
    }

  return h;
}


auto fastq_close(fastx_handle h) -> void
{
  fastx_close(h);
}


auto fastq_next(fastx_handle h,
                bool truncateatspace,
                const unsigned char * char_mapping) -> bool
{
  h->header_buffer.length = 0;
  h->header_buffer.data[0] = 0;
  h->sequence_buffer.length = 0;
  h->sequence_buffer.data[0] = 0;
  h->plusline_buffer.length = 0;
  h->plusline_buffer.data[0] = 0;
  h->quality_buffer.length = 0;
  h->quality_buffer.data[0] = 0;

  h->lineno_start = h->lineno;

  char msg[200];
  bool ok = true;
  char illegal_char = 0;

  uint64_t rest = fastx_file_fill_buffer(h);

  /* check end of file */

  if (rest == 0)
    {
      return false;
    }

  /* read header */

  /* check initial @ character */

  if (h->file_buffer.data[h->file_buffer.position] != '@')
    {
      fastq_fatal(h->lineno, "Header line must start with '@' character");
    }
  h->file_buffer.position++;
  rest--;

  char * lf = nullptr;
  while (lf == nullptr)
    {
      /* get more data if buffer empty */
      rest = fastx_file_fill_buffer(h);
      if (rest == 0)
        {
          fastq_fatal(h->lineno, "Unexpected end of file");
        }

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);

      /* copy to header buffer */
      uint64_t len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }
      buffer_extend(& h->header_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  /* read sequence line(s) */
  lf = nullptr;
  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(h);

      /* cannot end here */
      if (rest == 0)
        {
          fastq_fatal(h->lineno, "Unexpected end of file");
        }

      /* end when new line starting with + is seen */
      if (lf && (h->file_buffer.data[h->file_buffer.position] == '+'))
        {
          break;
        }

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);

      /* copy to sequence buffer */
      uint64_t len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }

      buffer_filter_extend(h,
                           & h->sequence_buffer,
                           h->file_buffer.data + h->file_buffer.position,
                           len,
                           char_fq_action_seq, char_mapping,
                           & ok, & illegal_char);
      h->file_buffer.position += len;
      rest -= len;

      if (! ok)
        {
          if ((illegal_char >= 32) && (illegal_char < 127))
            {
              snprintf(msg,
                       200,
                       "Illegal sequence character '%c'",
                       illegal_char);
            }
          else
            {
              snprintf(msg,
                       200,
                       "Illegal sequence character (unprintable, no %d)",
                       (unsigned char) illegal_char);
            }
          fastq_fatal(h->lineno - (lf ? 1 : 0), msg);
        }
    }

  /* read + line */

  /* skip + character */
  h->file_buffer.position++;
  rest--;

  lf = nullptr;
  while (lf == nullptr)
    {
      /* get more data if buffer empty */
      rest = fastx_file_fill_buffer(h);

      /* cannot end here */
      if (rest == 0)
        {
          fastq_fatal(h->lineno, "Unexpected end of file");
        }

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);
      /* copy to plusline buffer */
      uint64_t len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }
      buffer_extend(& h->plusline_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  /* check that the plus line is empty or identical to @ line */

  bool plusline_invalid = false;
  if (h->header_buffer.length == h->plusline_buffer.length)
    {
      if (memcmp(h->header_buffer.data,
                 h->plusline_buffer.data,
                 h->header_buffer.length))
        {
          plusline_invalid = true;
        }
    }
  else
    {
      if ((h->plusline_buffer.length > 2) ||
          ((h->plusline_buffer.length == 2) && (h->plusline_buffer.data[0] != '\r')))
        {
          plusline_invalid = true;
        }
    }
  if (plusline_invalid)
    {
      fastq_fatal(h->lineno - (lf ? 1 : 0),
                  "'+' line must be empty or identical to header");
    }

  /* read quality line(s) */

  lf = nullptr;
  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(h);

      /* end if no more data */
      if (rest == 0)
        {
          break;
        }

      /* end if next entry starts : LF + '@' + correct length */
      if (lf &&
          (h->file_buffer.data[h->file_buffer.position] == '@') &&
          (h->quality_buffer.length == h->sequence_buffer.length))
        {
          break;
        }

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);

      /* copy to quality buffer */
      uint64_t len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }

      buffer_filter_extend(h,
                           & h->quality_buffer,
                           h->file_buffer.data + h->file_buffer.position,
                           len,
                           char_fq_action_qual, chrmap_identity,
                           & ok, & illegal_char);
      h->file_buffer.position += len;
      rest -= len;

      /* break if quality line already too long */
      if (h->quality_buffer.length > h->sequence_buffer.length)
        {
          break;
        }

      if (! ok)
        {
          if ((illegal_char >= 32) && (illegal_char < 127))
            {
              snprintf(msg,
                       200,
                       "Illegal quality character '%c'",
                       illegal_char);
            }
          else
            {
              snprintf(msg,
                       200,
                       "Illegal quality character (unprintable, no %d)",
                       (unsigned char) illegal_char);
            }
          fastq_fatal(h->lineno - (lf ? 1 : 0), msg);
        }
    }

  if (h->sequence_buffer.length != h->quality_buffer.length)
    {
      fastq_fatal(h->lineno - (lf ? 1 : 0),
                  "Sequence and quality lines must be equally long");
    }

  fastx_filter_header(h, truncateatspace);

  h->seqno++;

  return true;
}


auto fastq_get_quality(fastx_handle h) -> char *
{
  return h->quality_buffer.data;
}


auto fastq_get_quality_length(fastx_handle h) -> uint64_t
{
  return h->quality_buffer.length;
}


auto fastq_get_position(fastx_handle h) -> uint64_t
{
  return h->file_position;
}


auto fastq_get_size(fastx_handle h) -> uint64_t
{
  return h->file_size;
}


auto fastq_get_lineno(fastx_handle h) -> uint64_t
{
  return h->lineno_start;
}


auto fastq_get_seqno(fastx_handle h) -> uint64_t
{
  return h->seqno;
}


auto fastq_get_header_length(fastx_handle h) -> uint64_t
{
  return h->header_buffer.length;
}


auto fastq_get_sequence_length(fastx_handle h) -> uint64_t
{
  return h->sequence_buffer.length;
}


auto fastq_get_header(fastx_handle h) -> char *
{
  return h->header_buffer.data;
}


auto fastq_get_sequence(fastx_handle h) -> char *
{
  return h->sequence_buffer.data;
}


auto fastq_get_abundance(fastx_handle h) -> int64_t
{
  // return 1 if not present
  int64_t const size = header_get_size(h->header_buffer.data,
                                 h->header_buffer.length);
  if (size > 0)
    {
      return size;
    }
  else
    {
      return 1;
    }
}


auto fastq_get_abundance_and_presence(fastx_handle h) -> int64_t
{
  // return 0 if not present
  return header_get_size(h->header_buffer.data, h->header_buffer.length);
}


inline auto fprint_seq_label(FILE * output_handle, char * seq, int len) -> void
{
  /* normalize first? */
  fprintf(output_handle, "%.*s", len, seq);
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
  fprintf(output_handle, "@");

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
  else if (opt_relabel && (ordinal > 0))
    {
      fprintf(output_handle, "%s%d", opt_relabel, ordinal);
    }
  else
    {
      bool const xsize = opt_xsize || (opt_sizeout && (abundance > 0));
      bool const xee = opt_xee || ((opt_eeout || opt_fastq_eeout) && (ee >= 0.0));
      bool const xlength = opt_xlength || opt_lengthout;
      header_fprint_strip(output_handle,
                          header,
                          header_len,
                          xsize,
                          xee,
                          xlength);
    }

  if (opt_label_suffix)
    {
      fprintf(output_handle, "%s", opt_label_suffix);
    }

  if (opt_sample)
    {
      fprintf(output_handle, ";sample=%s", opt_sample);
    }

  if (opt_sizeout && (abundance > 0))
    {
      fprintf(output_handle, ";size=%u", abundance);
    }

  if ((opt_eeout || opt_fastq_eeout) && (ee >= 0.0))
    {
      if (ee < 0.000000001) {
        fprintf(output_handle, ";ee=%.13lf", ee);
      } else if (ee < 0.00000001) {
        fprintf(output_handle, ";ee=%.12lf", ee);
      } else if (ee < 0.0000001) {
        fprintf(output_handle, ";ee=%.11lf", ee);
      } else if (ee < 0.000001) {
        fprintf(output_handle, ";ee=%.10lf", ee);
      } else if (ee < 0.00001) {
        fprintf(output_handle, ";ee=%.9lf", ee);
      } else if (ee < 0.0001) {
        fprintf(output_handle, ";ee=%.8lf", ee);
      } else if (ee < 0.001) {
        fprintf(output_handle, ";ee=%.7lf", ee);
      } else if (ee < 0.01) {
        fprintf(output_handle, ";ee=%.6lf", ee);
      } else if (ee < 0.1) {
        fprintf(output_handle, ";ee=%.5lf", ee);
      } else {
        fprintf(output_handle, ";ee=%.4lf", ee);
      }
    }

  if (opt_lengthout)
    {
      fprintf(output_handle, ";length=%d", len);
    }

  if (opt_relabel_keep &&
      ((opt_relabel && (ordinal > 0)) || opt_relabel_sha1 || opt_relabel_md5 || opt_relabel_self))
    {
      fprintf(output_handle, " %.*s", header_len, header);
    }

  fprintf(output_handle, "\n%.*s\n+\n%.*s\n", len, seq, len, quality);
}


auto fastq_print(std::FILE * output_handle, char * header, char * sequence, char * quality) -> void
{
  int const slen = std::strlen(sequence);
  int const hlen = std::strlen(header);
  fastq_print_general(output_handle, sequence, slen, header, hlen, quality, 0, 0, -1.0);
}
