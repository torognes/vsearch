/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#define FASTQ_BUFFER_ALLOC 8192

#ifdef HAVE_BZLIB
#define BZ_VERBOSE_0 0
#define BZ_VERBOSE_1 1
#define BZ_VERBOSE_2 2
#define BZ_VERBOSE_3 3
#define BZ_VERBOSE_4 4
#define BZ_MORE_MEM 0  /* faster decompression using more memory */
#define BZ_LESS_MEM 1  /* slower decompression but requires less memory */
#endif

#define FORMAT_PLAIN 1
#define FORMAT_BZIP  2
#define FORMAT_GZIP  3

static unsigned char MAGIC_GZIP[] = "\x1f\x8b";
static unsigned char MAGIC_BZIP[] = "BZ";

static char map_identity[256];

void fastq_fatal(unsigned long lineno, const char * msg)
{
  char * string; 
  if (asprintf(& string,
               "Invalid line %lu in FASTQ file: %s",
               lineno,
               msg) == -1)
    fatal("Out of memory");
  
  if (string)
    {
      fatal(string);
      free(string);
    }
  else
    fatal("Out of memory");
}

void buffer_init(struct fastq_buffer_s * buffer)
{
  buffer->alloc = FASTQ_BUFFER_ALLOC;
  buffer->data = (char*) xmalloc(buffer->alloc);
  buffer->data[0] = 0;
  buffer->length = 0;
  buffer->position = 0;
}

void buffer_free(struct fastq_buffer_s * buffer)
{
  if (buffer->data)
    free(buffer->data);
  buffer->data = 0;
  buffer->alloc = 0;
  buffer->length = 0;
  buffer->position = 0;
}

void buffer_makespace(struct fastq_buffer_s * buffer, unsigned long x)
{
  /* make sure there is space for x more chars in buffer */

  if (buffer->length + x > buffer->alloc)
    {
      /* alloc space for x more characters,
         but round up to nearest block size */
      buffer->alloc = 
        ((buffer->length + x + FASTQ_BUFFER_ALLOC - 1) / FASTQ_BUFFER_ALLOC)
        * FASTQ_BUFFER_ALLOC;
      buffer->data = (char*) xrealloc(buffer->data, buffer->alloc);
    }
}

void buffer_truncate(struct fastq_buffer_s * buffer, bool truncateatspace)
{
  /* Truncate the zero-terminated header string by inserting a new
     terminator (zero byte) at the first space (if truncateatspace)
     or first linefeed. */
  
  if (truncateatspace)
    buffer->length = strcspn(buffer->data, " \n");
  else
    buffer->length = strcspn(buffer->data, "\n");
  
  buffer->data[buffer->length] = 0;
}


fastq_handle fastq_open(const char * filename)
{
  fastq_handle h = (fastq_handle) xmalloc(sizeof(struct fastq_s));
  
  h->fp = NULL;
  h->fp = fopen(filename, "rb");
  if(!h->fp)
    fatal("Error: Unable to open fastq file for reading (%s)", filename);
  
  /* detect compression (plain, gzipped or bzipped) */
  
  unsigned char magic[2];
  h->format = FORMAT_PLAIN;
  if (fread(&magic, 1, 2, h->fp) >= 2)
    {
      if (!memcmp(magic, MAGIC_GZIP, 2))
        h->format = FORMAT_GZIP;
      else if (!memcmp(magic, MAGIC_BZIP, 2))
        h->format = FORMAT_BZIP;
    }

  /* Get size of original (uncompressed) file */

  if (fseek(h->fp, 0, SEEK_END))
    fatal("Error: Unable to seek in fastq file (%s)", filename);
  h->file_size = ftell(h->fp);
  rewind(h->fp);

  if (h->format == FORMAT_GZIP)
    {
      /* GZIP: Close ordinary file and open again as gzipped file */
#ifdef HAVE_ZLIB
      fclose(h->fp);
      if (! (h->fp_gz = gzopen(filename, "rb")))
        fatal("Unable to open gzip compressed fastq file (%s)", filename);
#else
      fatal("Files compressed with gzip are not supported");
#endif
    }

  if (h->format == FORMAT_BZIP)
    {
      /* BZIP2: Keep original file open, then open as bzipped file as well */
#ifdef HAVE_ZLIB
      int bzError;
      if (! (h->fp_bz = BZ2_bzReadOpen(& bzError, h->fp,
                                       BZ_VERBOSE_0, BZ_MORE_MEM, NULL, 0)))
        fatal("Unable to open bzip2 compressed fastq file (%s)", filename);
#else
      fatal("Files compressed with bzip2 are not supported");
#endif
    }

  h->stripped_all = 0;

  for(int i=0; i<256; i++)
    h->stripped[i] = 0;

  h->file_position = 0;

  buffer_init(& h->file_buffer);
  buffer_init(& h->header_buffer);
  buffer_init(& h->sequence_buffer);
  buffer_init(& h->quality_buffer);

  h->lineno = 1;
  h->seqno = -1;

  for(int i=0; i<256; i++)
    map_identity[i] = i;

  return h;
}

void fastq_close(fastq_handle h)
{
  /* Warn about stripped chars */

  if (h->stripped_all)
    {
      fprintf(stderr, "WARNING: invalid characters stripped from fastq file:");
      for (int i=0; i<256;i++)
        if (h->stripped[i])
          fprintf(stderr, " %c(%lu)", i, h->stripped[i]);
      fprintf(stderr, "\n");

      if (opt_log)
        {
          fprintf(fp_log, "WARNING: invalid characters stripped from fastq file:");
          for (int i=0; i<256;i++)
            if (h->stripped[i])
              fprintf(fp_log, " %c(%lu)", i, h->stripped[i]);
          fprintf(fp_log, "\n");
        }
    }

#ifdef HAVE_BZLIB
  int bz_error;
#endif
  
  switch(h->format)
    {
    case FORMAT_PLAIN:
      fclose(h->fp);
      h->fp = 0;
      break;

    case FORMAT_GZIP:
#ifdef HAVE_ZLIB
      gzclose(h->fp_gz);
      h->fp_gz = 0;
      break;
#endif
      
    case FORMAT_BZIP:
#ifdef HAVE_BZLIB
      BZ2_bzReadClose(&bz_error, h->fp_bz);
      h->fp_bz = 0;
      break;
#endif

    default:
      fatal("Internal error");
    }
  
  buffer_free(& h->file_buffer);
  buffer_free(& h->header_buffer);
  buffer_free(& h->sequence_buffer);
  buffer_free(& h->quality_buffer);

  h->file_size = 0;
  h->file_position = 0;

  h->lineno = 0;
  h->seqno = -1;
}


unsigned long fastq_file_fill_buffer(fastq_handle h)
{
  /* read more data if necessary */
  unsigned long rest = h->file_buffer.length - h->file_buffer.position;
  
  if (rest > 0)
    return rest;
  else
    {
      unsigned long space = h->file_buffer.alloc - h->file_buffer.length;

      if (space == 0)
        {
          /* back to beginning of buffer */
          h->file_buffer.position = 0;
          h->file_buffer.length = 0;
          space = h->file_buffer.alloc;
        }
      
      int bytes_read = 0;
      
#ifdef HAVE_BZLIB
      int bzError = 0;
#endif

      switch(h->format)
        {
        case FORMAT_PLAIN:
          bytes_read = fread(h->file_buffer.data
                             + h->file_buffer.position,
                             1,
                             space,
                             h->fp);
          break;

        case FORMAT_GZIP:
#ifdef HAVE_ZLIB
          bytes_read = gzread(h->fp_gz,
                              h->file_buffer.data + h->file_buffer.position,
                              space);
          if (bytes_read < 0)
            fatal("Error reading gzip compressed fastq file");
          break;
#endif
          
        case FORMAT_BZIP:
#ifdef HAVE_BZLIB
          bytes_read = BZ2_bzRead(& bzError,
                                  h->fp_bz,
                                  h->file_buffer.data 
                                  + h->file_buffer.position,
                                  space);
          if ((bytes_read < 0) ||
              ! ((bzError == BZ_OK) ||
                 (bzError == BZ_STREAM_END) ||
                 (bzError == BZ_SEQUENCE_ERROR)))
            fatal("Error reading bzip2 compressed fastq file");
          break;
#endif
          
        default:
          fatal("Internal error");
        }
      
      h->file_buffer.length += bytes_read;
      return bytes_read;
    }
}

void buffer_extend(struct fastq_buffer_s * dest_buffer,
                  char * source_buf,
                  unsigned long len)
{
  buffer_makespace(dest_buffer, len+1);
  memcpy(dest_buffer->data + dest_buffer->length,
         source_buf,
         len);
  dest_buffer->length += len;
  dest_buffer->data[dest_buffer->length] = 0;
}

void buffer_filter_extend(fastq_handle h,
                          struct fastq_buffer_s * dest_buffer,
                          char * source_buf,
                          unsigned long len,
                          unsigned int * char_action,
                          char * char_mapping,
                          unsigned long lineno_start)
{
  buffer_makespace(dest_buffer, len+1);

  /* Strip unwanted characters from the string and raise warnings or
     errors on certain characters. */

  unsigned long lineno = lineno_start;

  char * p = source_buf;
  char * d = dest_buffer->data + dest_buffer->length;
  char * q = d;
  char msg[200];

  for(unsigned long i = 0; i < len; i++)
    {
      char c = *p++;
      char m = char_action[(int)c];

      switch(m)
        {
        case 0:
          /* stripped */
          h->stripped_all++;
          h->stripped[(int)c]++;
          break;
              
        case 1:
          /* legal character */
          *q++ = char_mapping[(int)(c)];
          break;
          
        case 2:
          /* fatal character */
          if ((c>=32) && (c<127))
            snprintf(msg,
                     200,
                     "Illegal character '%c'",
                     c);
          else
            snprintf(msg,
                     200,
                     "Illegal unprintable character %#.2x (hexadecimal)",
                     c);
          fastq_fatal(lineno, msg);
          break;
              
        case 3:
          /* silently stripped chars (whitespace) */
          break;

        case 4:
          /* newline (silently stripped) */
          lineno++;
          break;
        }
    }

  /* add zero after sequence */
  *q = 0;
  dest_buffer->length += q - d;
}

bool fastq_next(fastq_handle h,
                bool truncateatspace,
                char * char_mapping)
{
  h->header_buffer.length = 0;
  h->header_buffer.data[0] = 0;
  h->sequence_buffer.length = 0;
  h->sequence_buffer.data[0] = 0;
  h->quality_buffer.length = 0;
  h->quality_buffer.data[0] = 0;

  unsigned long rest = fastq_file_fill_buffer(h);

  /* check end of file */

  if (rest == 0)
    return 0;

  /* read header */

  /* check initial @ character */
  
  if (h->file_buffer.data[h->file_buffer.position] != '@')
    fastq_fatal(h->lineno, "Header line must start with '@' character");
  h->file_buffer.position++;
  rest--;

  char * lf = 0;
  while (lf == 0)
    {
      /* get more data if buffer empty */
      rest = fastq_file_fill_buffer(h);
      if (rest == 0)
        fastq_fatal(h->lineno, "Unexpected end of file");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);

      /* copy to header buffer */
      unsigned long len = rest;
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

  unsigned long lineno_seq = h->lineno;

  /* read sequence line(s) */
  lf = 0;
  while (1)
    {
      /* get more data, if necessary */
      rest = fastq_file_fill_buffer(h);

      /* cannot end here */
      if (rest == 0)
        fastq_fatal(h->lineno, "Unexpected end of file");
      
      /* end when new line starting with + is seen */
      if (lf && (h->file_buffer.data[h->file_buffer.position] == '+'))
        break;

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);
      
      /* copy to sequence buffer */
      unsigned long len = rest;
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
                           char_fq_action_seq, char_mapping, lineno_seq);
      h->file_buffer.position += len;
      rest -= len;
    }

  if (h->sequence_buffer.length == 0)
    fastq_fatal(lineno_seq, "Empty sequence line");

  unsigned long lineno_plus = h->lineno;

  /* read + line */
  fastq_buffer_s plusline_buffer;
  buffer_init(&plusline_buffer);

  /* skip + character */
  h->file_buffer.position++;
  rest--;
  
  lf = 0;
  while (lf == 0)
    {
      /* get more data if buffer empty */
      rest = fastq_file_fill_buffer(h);

      /* cannot end here */
      if (rest == 0)
        fastq_fatal(h->lineno, "Unexpected end of file");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);
      /* copy to plusline buffer */
      unsigned long len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }
      buffer_extend(& plusline_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  /* check that the plus line is empty or identical to @ line */

  bool plusline_invalid = 0;
  if (h->header_buffer.length == plusline_buffer.length)
    {
      if (memcmp(h->header_buffer.data,
                 plusline_buffer.data,
                 h->header_buffer.length))
        plusline_invalid = 1;
    }
  else
    {
      if ((plusline_buffer.length > 2) ||
          ((plusline_buffer.length == 2) && (plusline_buffer.data[0] != '\r')))
        plusline_invalid = 1;
    }
  if (plusline_invalid)
    fastq_fatal(lineno_plus, 
                "'+' line must be empty or identical to header");

  /* read quality line(s) */

  unsigned long lineno_qual = h->lineno;

  lf = 0;
  while (1)
    {
      /* get more data, if necessary */
      rest = fastq_file_fill_buffer(h);

      /* end if no more data */
      if (rest == 0)
        break;
      
      /* end if next entry starts : LF + '@' + correct length */
      if (lf && 
          (h->file_buffer.data[h->file_buffer.position] == '@') &&
          (h->quality_buffer.length == h->sequence_buffer.length))
        break;
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);
      
      /* copy to quality buffer */
      unsigned long len = rest;
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
                           char_fq_action_qual, map_identity, lineno_qual);
      h->file_buffer.position += len;
      rest -= len;
    }

  if (h->quality_buffer.length == 0)
    fastq_fatal(lineno_seq, "Empty quality line");

  if (h->sequence_buffer.length != h->quality_buffer.length)
    fastq_fatal(lineno_qual,
                "Sequence and quality lines must be equally long");

  buffer_truncate(& h->header_buffer, truncateatspace);

#ifdef HAVE_ZLIB
  if (h->format == FORMAT_GZIP)
    h->file_position = gzoffset(h->fp_gz);
  else
#endif
    h->file_position = ftell(h->fp);
  
  h->seqno++;

  return 1;
}

char * fastq_get_quality(fastq_handle h)
{
  return h->quality_buffer.data;
}

unsigned long fastq_get_quality_length(fastq_handle h)
{
  return h->quality_buffer.length;
}

unsigned long fastq_get_position(fastq_handle h)
{
  return h->file_position;
}

unsigned long fastq_get_size(fastq_handle h)
{
  return h->file_size;
}

unsigned long fastq_get_lineno(fastq_handle h)
{
  return h->lineno;
}

unsigned long fastq_get_seqno(fastq_handle h)
{
  return h->seqno;
}

unsigned long fastq_get_header_length(fastq_handle h)
{
  return h->header_buffer.length;
}

unsigned long fastq_get_sequence_length(fastq_handle h)
{
  return h->sequence_buffer.length;
}

char * fastq_get_header(fastq_handle h)
{
  return h->header_buffer.data;
}

char * fastq_get_sequence(fastq_handle h)
{
  return h->sequence_buffer.data;
}
