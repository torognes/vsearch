/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

/* file compression and format detector */
/* basic file buffering function for fastq and fastx parsers */

#define FASTX_BUFFER_ALLOC 8192

#ifdef HAVE_BZLIB_H
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


void buffer_init(struct fastx_buffer_s * buffer)
{
  buffer->alloc = FASTX_BUFFER_ALLOC;
  buffer->data = (char*) xmalloc(buffer->alloc);
  buffer->data[0] = 0;
  buffer->length = 0;
  buffer->position = 0;
}

void buffer_free(struct fastx_buffer_s * buffer)
{
  if (buffer->data)
    xfree(buffer->data);
  buffer->data = 0;
  buffer->alloc = 0;
  buffer->length = 0;
  buffer->position = 0;
}

void buffer_makespace(struct fastx_buffer_s * buffer, uint64_t x)
{
  /* make sure there is space for x more chars in buffer */

  if (buffer->length + x > buffer->alloc)
    {
      /* alloc space for x more characters,
         but round up to nearest block size */
      buffer->alloc = 
        ((buffer->length + x + FASTX_BUFFER_ALLOC - 1) / FASTX_BUFFER_ALLOC)
        * FASTX_BUFFER_ALLOC;
      buffer->data = (char*) xrealloc(buffer->data, buffer->alloc);
    }
}

void buffer_truncate(struct fastx_buffer_s * buffer, bool truncateatspace)
{
  /* Truncate the zero-terminated header string by inserting a new
     terminator (zero byte) at the first space/tab character
     (if truncateatspace) or first linefeed. */
  
  if (truncateatspace)
    buffer->length = strcspn(buffer->data, " \t\n");
  else
    buffer->length = strcspn(buffer->data, "\n");
  
  buffer->data[buffer->length] = 0;
}

void buffer_extend(struct fastx_buffer_s * dest_buffer,
                  char * source_buf,
                  uint64_t len)
{
  buffer_makespace(dest_buffer, len+1);
  memcpy(dest_buffer->data + dest_buffer->length,
         source_buf,
         len);
  dest_buffer->length += len;
  dest_buffer->data[dest_buffer->length] = 0;
}

fastx_handle fastx_open(const char * filename)
{
  fastx_handle h = (fastx_handle) xmalloc(sizeof(struct fastx_s));

  h->fp = 0;

#ifdef HAVE_ZLIB_H
  h->fp_gz = 0;
#endif

#ifdef HAVE_BZLIB_H
  h->fp_bz = 0;
  int bzError = 0;
#endif
 
  h->fp = fopen_input(filename);
  if (!h->fp)
    fatal("Unable to open file for reading (%s)", filename);

  /* Get mode and size of original (uncompressed) file */

  xstat_t fs;
  if (xfstat(fileno(h->fp), & fs))
    fatal("Unable to get status for input file (%s)", filename);

  h->is_pipe = S_ISFIFO(fs.st_mode);

  if (h->is_pipe)
    h->file_size = 0;
  else
    h->file_size = fs.st_size;

  if (opt_gzip_decompress)
    {
      h->format = FORMAT_GZIP;
    }
  else if (opt_bzip2_decompress)
    {
      h->format = FORMAT_BZIP;
    }
  else if (h->is_pipe)
    {
      h->format = FORMAT_PLAIN;
    }
  else
    {
      /* autodetect compression (plain, gzipped or bzipped) */

      /* read two characters and compare with magic */

      unsigned char magic[2];

      h->format = FORMAT_PLAIN;
      if (fread(&magic, 1, 2, h->fp) < 2)
        fatal("Unable to read from file (%s)", filename);

      if (memcmp(magic, MAGIC_GZIP, 2) == 0)
        h->format = FORMAT_GZIP;
      else if (memcmp(magic, MAGIC_BZIP, 2) == 0)
        h->format = FORMAT_BZIP;

      /* close and reopen to avoid problems with gzip library */
      /* rewind was not enough */

      fclose(h->fp);
      h->fp = fopen_input(filename);
      if (!h->fp)
        fatal("Unable to open file for reading (%s)", filename);
    }

  if (h->format == FORMAT_GZIP)
    {
      /* GZIP: Keep original file open, then open as gzipped file as well */
#ifdef HAVE_ZLIB_H
      if (!gz_lib)
        fatal("Files compressed with gzip are not supported");
      if (! (h->fp_gz = (*gzdopen_p)(fileno(h->fp), "rb"))) // dup?
        fatal("Unable to open gzip compressed file (%s)", filename);
#else
      fatal("Files compressed with gzip are not supported");
#endif
    }

  if (h->format == FORMAT_BZIP)
    {
      /* BZIP2: Keep original file open, then open as bzipped file as well */
#ifdef HAVE_BZLIB_H
      if (!bz2_lib)
        fatal("Files compressed with bzip2 are not supported");
      if (! (h->fp_bz = (*BZ2_bzReadOpen_p)(& bzError, h->fp,
                                         BZ_VERBOSE_0, BZ_MORE_MEM,
                                         NULL, 0)))
        fatal("Unable to open bzip2 compressed file (%s)", filename);
#else
      fatal("Files compressed with bzip2 are not supported");
#endif
    }

  /* init buffers */

  h->file_position = 0;

  buffer_init(& h->file_buffer);

  /* start filling up file buffer */

  uint64_t rest = fastx_file_fill_buffer(h);
  
  if (rest < 2)
    fatal("File too small");


  /* examine first char and see if it starts with > or @ */

  int filetype = 0;
  char * first = h->file_buffer.data;

  if (*first == '>')
    {
      filetype = 1;
      h->is_fastq = 0;
    }
  else if (*first == '@')
    {
      filetype = 2;
      h->is_fastq = 1;
    }

  if (filetype == 0)
    {
      /* close files if unrecognized file type */
      
      switch(h->format)
        {
        case FORMAT_PLAIN:
          break;
          
        case FORMAT_GZIP:
#ifdef HAVE_ZLIB_H
          (*gzclose_p)(h->fp_gz);
          h->fp_gz = 0;
          break;
#endif
          
        case FORMAT_BZIP:
#ifdef HAVE_BZLIB_H
          (*BZ2_bzReadClose_p)(&bzError, h->fp_bz);
          h->fp_bz = 0;
          break;
#endif
          
        default:
          fatal("Internal error");
        }
      
      fclose(h->fp);
      h->fp = 0;

      if (memcmp(first, MAGIC_GZIP, 2) == 0)
        fatal("File appears to be gzip compressed. Please use --gzip_decompress");

      if (memcmp(first, MAGIC_BZIP, 2) == 0)
        fatal("File appears to be bzip2 compressed. Please use --bzip2_decompress");

      fatal("File type not recognized.");

      return 0;
    }

  /* more initialization */

  buffer_init(& h->header_buffer);
  buffer_init(& h->sequence_buffer);
  buffer_init(& h->plusline_buffer);
  buffer_init(& h->quality_buffer);

  h->stripped_all = 0;

  for(int i=0; i<256; i++)
    h->stripped[i] = 0;

  h->lineno = 1;
  h->lineno_start = 1;
  h->seqno = -1;

  return h;
}

bool fastx_is_fastq(fastx_handle h)
{
  return h->is_fastq;
}

void fastx_close(fastx_handle h)
{
  /* Warn about stripped chars */

  if (h->stripped_all)
    {
      fprintf(stderr, "WARNING: %" PRIu64 " invalid characters stripped from %s file:", h->stripped_all, (h->is_fastq ? "FASTQ" : "FASTA"));
      for (int i=0; i<256;i++)
        if (h->stripped[i])
          fprintf(stderr, " %c(%" PRIu64 ")", i, h->stripped[i]);
      fprintf(stderr, "\n");

      if (opt_log)
        {
          fprintf(fp_log, "WARNING: %" PRIu64 " invalid characters stripped from %s file:", h->stripped_all, (h->is_fastq ? "FASTQ" : "FASTA"));
          for (int i=0; i<256;i++)
            if (h->stripped[i])
              fprintf(fp_log, " %c(%" PRIu64 ")", i, h->stripped[i]);
          fprintf(fp_log, "\n");
        }
    }

#ifdef HAVE_BZLIB_H
  int bz_error;
#endif

  switch(h->format)
    {
    case FORMAT_PLAIN:
      break;

    case FORMAT_GZIP:
#ifdef HAVE_ZLIB_H
      (*gzclose_p)(h->fp_gz);
      h->fp_gz = 0;
      break;
#endif

    case FORMAT_BZIP:
#ifdef HAVE_BZLIB_H
      (*BZ2_bzReadClose_p)(&bz_error, h->fp_bz);
      h->fp_bz = 0;
      break;
#endif

    default:
      fatal("Internal error");
    }

  fclose(h->fp);
  h->fp = 0;

  buffer_free(& h->file_buffer);
  buffer_free(& h->header_buffer);
  buffer_free(& h->sequence_buffer);
  buffer_free(& h->plusline_buffer);
  buffer_free(& h->quality_buffer);

  h->file_size = 0;
  h->file_position = 0;

  h->lineno = 0;
  h->seqno = -1;

  xfree(h);
  h=0;
}

uint64_t fastx_file_fill_buffer(fastx_handle h)
{
  /* read more data if necessary */
  uint64_t rest = h->file_buffer.length - h->file_buffer.position;

  if (rest > 0)
    return rest;
  else
    {
      uint64_t space = h->file_buffer.alloc - h->file_buffer.length;

      if (space == 0)
        {
          /* back to beginning of buffer */
          h->file_buffer.position = 0;
          h->file_buffer.length = 0;
          space = h->file_buffer.alloc;
        }

      int bytes_read = 0;

#ifdef HAVE_BZLIB_H
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
#ifdef HAVE_ZLIB_H
          bytes_read = (*gzread_p)(h->fp_gz,
                                   h->file_buffer.data
                                   + h->file_buffer.position,
                                   space);
          if (bytes_read < 0)
            fatal("Unable to read gzip compressed file");
          break;
#endif

        case FORMAT_BZIP:
#ifdef HAVE_BZLIB_H
          bytes_read = (*BZ2_bzRead_p)(& bzError,
                                       h->fp_bz,
                                       h->file_buffer.data
                                       + h->file_buffer.position,
                                       space);
          if ((bytes_read < 0) ||
              ! ((bzError == BZ_OK) ||
                 (bzError == BZ_STREAM_END) ||
                 (bzError == BZ_SEQUENCE_ERROR)))
            fatal("Unable to read from bzip2 compressed file");
          break;
#endif

        default:
          fatal("Internal error");
        }

      if (!h->is_pipe)
        {
#ifdef HAVE_ZLIB_H
          if (h->format == FORMAT_GZIP)
            {
              /* Circumvent the missing gzoffset function in zlib 1.2.3 and earlier */
              int fd = dup(fileno(h->fp));
              h->file_position = xlseek(fd, 0, SEEK_CUR);
              close(fd);
            }
          else
#endif
            h->file_position = xftello(h->fp);
        }

      h->file_buffer.length += bytes_read;
      return bytes_read;
    }
}

bool fastx_next(fastx_handle h,
                bool truncateatspace,
                const unsigned char * char_mapping)
{
  if (h->is_fastq)
    return fastq_next(h, truncateatspace, char_mapping);
  else
    return fasta_next(h, truncateatspace, char_mapping);
}

uint64_t fastx_get_position(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_position(h);
  else
    return fasta_get_position(h);
}


uint64_t fastx_get_size(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_size(h);
  else
    return fasta_get_size(h);
}


uint64_t fastx_get_lineno(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_lineno(h);
  else
    return fasta_get_lineno(h);
}


uint64_t fastx_get_seqno(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_seqno(h);
  else
    return fasta_get_seqno(h);
}

char * fastx_get_header(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_header(h);
  else
    return fasta_get_header(h);
}

char * fastx_get_sequence(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_sequence(h);
  else
    return fasta_get_sequence(h);
}

uint64_t fastx_get_header_length(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_header_length(h);
  else
    return fasta_get_header_length(h);
}

uint64_t fastx_get_sequence_length(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_sequence_length(h);
  else
    return fasta_get_sequence_length(h);
}


char * fastx_get_quality(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_quality(h);
  else
    return 0;
}

int64_t fastx_get_abundance(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_abundance(h);
  else
    return fasta_get_abundance(h);
}

