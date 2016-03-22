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

/* file type detector and wrapper for fastq and fastx parser */

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

int vfastx_detect(const char * filename)
{
#ifdef HAVE_ZLIB_H
  gzFile fp_gz = 0;
#endif

#ifdef HAVE_BZLIB_H
  BZFILE * fp_bz = 0;
#endif

  int format;

  FILE * fp = fopen(filename, "rb");
  if (!fp)
    fatal("Error: Unable to open file for reading (%s)", filename);
  
  /* detect compression (plain, gzipped or bzipped) */
  
  unsigned char magic[2];
  format = FORMAT_PLAIN;
  if (fread(&magic, 1, 2, fp) >= 2)
    {
      if (!memcmp(magic, MAGIC_GZIP, 2))
        format = FORMAT_GZIP;
      else if (!memcmp(magic, MAGIC_BZIP, 2))
        format = FORMAT_BZIP;
    }

  /* close and reopen to avoid problems with gzip library */
  /* rewind was not enough */

  fclose(fp);
  fp = fopen(filename, "rb");
  if (!fp)
    fatal("Error: Unable to open file for reading (%s)", filename);

  if (format == FORMAT_GZIP)
    {
      /* GZIP: Keep original file open, then open as bzipped file as well */
#ifdef HAVE_ZLIB_H
      if (!gz_lib)
        fatal("Files compressed with gzip are not supported");
      if (! (fp_gz = (*gzdopen_p)(fileno(fp), "rb")))
        fatal("Unable to open gzip compressed file (%s)", filename);
#else
      fatal("Files compressed with gzip are not supported");
#endif
    }

  if (format == FORMAT_BZIP)
    {
      /* BZIP2: Keep original file open, then open as bzipped file as well */
#ifdef HAVE_ZLIB_H
      if (!bz2_lib)
        fatal("Files compressed with bzip2 are not supported");
      int bzError;
      if (! (fp_bz = (*BZ2_bzReadOpen_p)(& bzError, fp,
                                       BZ_VERBOSE_0, BZ_MORE_MEM, NULL, 0)))
        fatal("Unable to open bzip2 compressed file (%s)", filename);
#else
      fatal("Files compressed with bzip2 are not supported");
#endif
    }

  /* read one char and see if it starts with > or @ */

  const int BUFFERLEN = 1;
  char buffer[BUFFERLEN];
  
  int bytes_read = 0;
  
#ifdef HAVE_BZLIB_H
  int bzError = 0;
#endif
 
  switch(format)
    {
    case FORMAT_PLAIN:
      bytes_read = fread(buffer,
                         1,
                         BUFFERLEN,
                         fp);
      break;
      
    case FORMAT_GZIP:
#ifdef HAVE_ZLIB_H
      bytes_read = (*gzread_p)(fp_gz,
                             buffer,
                             BUFFERLEN);
      if (bytes_read < 0)
        fatal("Error reading gzip compressed file (%s)", filename);
      break;
#endif
      
    case FORMAT_BZIP:
#ifdef HAVE_BZLIB_H
      bytes_read = (*BZ2_bzRead_p)(& bzError,
                                 fp_bz,
                                 buffer,
                                 BUFFERLEN);
      if ((bytes_read < 0) ||
          ! ((bzError == BZ_OK) ||
             (bzError == BZ_STREAM_END) ||
             (bzError == BZ_SEQUENCE_ERROR)))
        fatal("Error reading bzip2 compressed file (%s)", filename);
      break;
#endif
      
    default:
      fatal("Internal error");
    }

  if (bytes_read < BUFFERLEN)
    fatal("Error reading file (%s)", filename);

  int filetype = 0;
  if (buffer[0] == '>')
    filetype = 1;
  else if (buffer[0] == '@')
    filetype = 2;

  /* close files */

#ifdef HAVE_BZLIB_H
  int bz_error;
#endif
  
  switch(format)
    {
    case FORMAT_PLAIN:
      fclose(fp);
      fp = 0;
      break;

    case FORMAT_GZIP:
#ifdef HAVE_ZLIB_H
      (*gzclose_p)(fp_gz);
      fp_gz = 0;
      break;
#endif
      
    case FORMAT_BZIP:
#ifdef HAVE_BZLIB_H
      (*BZ2_bzReadClose_p)(&bz_error, fp_bz);
      fp_bz = 0;
      break;
#endif

    default:
      fatal("Internal error");
    }

  return filetype;
}

bool fastx_is_fastq(fastx_handle h)
{
  return h->is_fastq;
}

fastx_handle fastx_open(const char * filename)
{
  int filetype = vfastx_detect(filename);

  if (filetype == 0)
    return 0;

  fastx_handle h = (fastx_handle) xmalloc(sizeof(struct fastx_s));

  h->is_fastq = (filetype == 2);

  if (h->is_fastq)
    h->handle.fastq = fastq_open(filename);
  else
    h->handle.fasta = fasta_open(filename);

  return h;
}

void fastx_close(fastx_handle h)
{
  if (h->is_fastq)
    {
      fastq_close(h->handle.fastq);
      h->handle.fastq = 0;
    }
  else
    {
      fasta_close(h->handle.fasta);
      h->handle.fasta = 0;
    }
  free(h);
}

bool fastx_next(fastx_handle h,
                bool truncateatspace,
                char * char_mapping)
{
  if (h->is_fastq)
    return fastq_next(h->handle.fastq, truncateatspace, char_mapping);
  else
    return fasta_next(h->handle.fasta, truncateatspace, char_mapping);
}

unsigned long fastx_get_position(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_position(h->handle.fastq);
  else
    return fasta_get_position(h->handle.fasta);
}


unsigned long fastx_get_size(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_size(h->handle.fastq);
  else
    return fasta_get_size(h->handle.fasta);
}


unsigned long fastx_get_lineno(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_lineno(h->handle.fastq);
  else
    return fasta_get_lineno(h->handle.fasta);
}


unsigned long fastx_get_seqno(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_seqno(h->handle.fastq);
  else
    return fasta_get_seqno(h->handle.fasta);
}

char * fastx_get_header(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_header(h->handle.fastq);
  else
    return fasta_get_header(h->handle.fasta);
}

char * fastx_get_sequence(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_sequence(h->handle.fastq);
  else
    return fasta_get_sequence(h->handle.fasta);
}

unsigned long fastx_get_header_length(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_header_length(h->handle.fastq);
  else
    return fasta_get_header_length(h->handle.fasta);
}

unsigned long fastx_get_sequence_length(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_sequence_length(h->handle.fastq);
  else
    return fasta_get_sequence_length(h->handle.fasta);
}


char * fastx_get_quality(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_quality(h->handle.fastq);
  else
    return 0;
}

long fastx_get_abundance(fastx_handle h)
{
  if (h->is_fastq)
    return fastq_get_abundance(h->handle.fastq);
  else
    return fasta_get_abundance(h->handle.fasta);
}

