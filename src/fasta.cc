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

#define FASTA_BUFFER_ALLOC 8192

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

void buffer_init(struct fasta_buffer_s * buffer)
{
  buffer->alloc = FASTA_BUFFER_ALLOC;
  buffer->data = (char*) xmalloc(buffer->alloc);
  buffer->data[0] = 0;
  buffer->length = 0;
  buffer->position = 0;
}

void buffer_free(struct fasta_buffer_s * buffer)
{
  if (buffer->data)
    free(buffer->data);
  buffer->data = 0;
  buffer->alloc = 0;
  buffer->length = 0;
  buffer->position = 0;
}

fasta_handle fasta_open(const char * filename)
{
  fasta_handle h = (fasta_handle) xmalloc(sizeof(struct fasta_s));
  
  h->fp = NULL;
  h->fp = fopen(filename, "rb");
  if(!h->fp)
    fatal("Error: Unable to open fasta file for reading (%s)", filename);
  
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
    fatal("Error: Unable to seek in fasta file (%s)", filename);
  h->file_size = ftell(h->fp);
  rewind(h->fp);

  if (h->format == FORMAT_GZIP)
    {
      /* GZIP: Keep original file open, then open as gzipped file as well */
#ifdef HAVE_ZLIB_H
      if (!gz_lib)
        fatal("Files compressed with gzip are not supported");
      if (! (h->fp_gz = (*gzdopen_p)(fileno(h->fp), "rb")))
        fatal("Unable to open gzip compressed fasta file (%s)", filename);
#else
      fatal("Files compressed with gzip are not supported");
#endif
    }

  if (h->format == FORMAT_BZIP)
    {
      /* BZIP2: Keep original file open, then open as bzipped file as well */
#ifdef HAVE_ZLIB_H
      if (!bz2_lib)
        fatal("Files compressed with bzip2 are not supported");
      int bzError;
      if (! (h->fp_bz = (*BZ2_bzReadOpen_p)(& bzError, h->fp,
                                            BZ_VERBOSE_0, BZ_MORE_MEM,
                                            NULL, 0)))
        fatal("Unable to open bzip2 compressed fasta file (%s)", filename);
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

  h->lineno = 1;
  h->seqno = -1;

  return h;
}

void fasta_close(fasta_handle h)
{
  /* Warn about stripped chars */

  if (h->stripped_all)
    {
      fprintf(stderr, "WARNING: invalid characters stripped from fasta file:");
      for (int i=0; i<256;i++)
        if (h->stripped[i])
          fprintf(stderr, " %c(%lu)", i, h->stripped[i]);
      fprintf(stderr, "\n");

      if (opt_log)
        {
          fprintf(fp_log, "WARNING: invalid characters stripped from fasta file:");
          for (int i=0; i<256;i++)
            if (h->stripped[i])
              fprintf(fp_log, " %c(%lu)", i, h->stripped[i]);
          fprintf(fp_log, "\n");
        }
    }

#ifdef HAVE_BZLIB_H
  int bz_error;
#endif
  
  switch(h->format)
    {
    case FORMAT_PLAIN:
      fclose(h->fp);
      h->fp = 0;
      break;

    case FORMAT_GZIP:
#ifdef HAVE_ZLIB_H
      (*gzclose_p)(h->fp_gz);
      h->fp_gz = 0;
      h->fp = 0;
      break;
#endif
      
    case FORMAT_BZIP:
#ifdef HAVE_BZLIB_H
      (*BZ2_bzReadClose_p)(&bz_error, h->fp_bz);
      h->fp_bz = 0;
      h->fp = 0;
      break;
#endif

    default:
      fatal("Internal error");
    }
  
  buffer_free(& h->file_buffer);
  buffer_free(& h->header_buffer);
  buffer_free(& h->sequence_buffer);

  h->file_size = 0;
  h->file_position = 0;

  h->lineno = 0;
  h->seqno = -1;

  free(h);
}


unsigned long fasta_file_fill_buffer(fasta_handle h)
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
            fatal("Error reading gzip compressed fasta file");
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
            fatal("Error reading bzip2 compressed fasta file");
          break;
#endif
          
        default:
          fatal("Internal error");
        }
      
      h->file_buffer.length += bytes_read;
      return bytes_read;
    }
}

void buffer_extend(struct fasta_buffer_s * buffer, char * buf, unsigned long len)
{
  if (buffer->length + len + 1 > buffer->alloc)
    {
      /* alloc space for len more characters + terminating zero,
         but round up to nearest block size */
      buffer->alloc = 
        FASTA_BUFFER_ALLOC * 
        (((buffer->length + len) / FASTA_BUFFER_ALLOC) + 1);
      buffer->data = (char*) xrealloc(buffer->data, buffer->alloc);
    }

  /* copy string */
  memcpy(buffer->data + buffer->length, buf, len);
  buffer->length += len;

  /* add terminator */
  buffer->data[buffer->length] = 0;
}

void fasta_truncate_header(fasta_handle h, bool truncateatspace)
{
  /* Truncate the zero-terminated header string by inserting a new
     terminator (zero byte) at the first space (if truncateatspace)
     or first linefeed. */
  
  if (truncateatspace)
    h->header_buffer.length = strcspn(h->header_buffer.data, " \n");
  else
    h->header_buffer.length = strcspn(h->header_buffer.data, "\n");
  
  h->header_buffer.data[h->header_buffer.length] = 0;
}

void fasta_filter_sequence(fasta_handle h,
                           unsigned int * char_action,
                           char * char_mapping)
{
  /* Strip unwanted characters from the sequence and raise warnings or
     errors on certain characters. */

  char * p = h->sequence_buffer.data;
  char * q = p;
  char c;
  char msg[200];

  while ((c = *p++))
    {
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
                     "illegal character '%c' on line %lu in fasta file",
                     c,
                     h->lineno);
          else
            snprintf(msg,
                     200,
                     "illegal unprintable character %#.2x (hexadecimal) on line %lu in fasta file",
                     c,
                     h->lineno);
          fatal(msg);
          break;
              
        case 3:
          /* silently stripped chars (whitespace) */
          break;

        case 4:
          /* newline (silently stripped) */
          h->lineno++;
          break;
        }
    }

  /* add zero after sequence */
  *q = 0;
  h->sequence_buffer.length = q - h->sequence_buffer.data;
}

bool fasta_next(fasta_handle h,
                bool truncateatspace,
                char * char_mapping)
{
  h->header_buffer.length = 0;
  h->sequence_buffer.length = 0;

  unsigned long rest = fasta_file_fill_buffer(h);

  if (rest == 0)
    return 0;

  /* read header */

  /* check initial > character */
  
  if (h->file_buffer.data[h->file_buffer.position] != '>')
    fatal("Invalid FASTA - header must start with > character");
  h->file_buffer.position++;
  rest--;

  char * lf = 0;
  while (lf == 0)
    {
      /* get more data if buffer empty*/
      rest = fasta_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTA - header must be terminated with newline");
      
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

  /* read one or more sequence lines */

  while (1)
    {
      /* get more data, if necessary */
      rest = fasta_file_fill_buffer(h);

      /* end if no more data */
      if (rest == 0)
        break;

      /* end if new sequence starts */
      if (lf && (h->file_buffer.data[h->file_buffer.position] == '>'))
        break;

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);

      unsigned long len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
        }
      buffer_extend(& h->sequence_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  h->seqno++;

  fasta_truncate_header(h, truncateatspace);
  fasta_filter_sequence(h, char_fasta_action, char_mapping);

#ifdef HAVE_ZLIB_H
  if (h->format == FORMAT_GZIP)
    {
      /* Circumvent the missing gzoffset function in zlib 1.2.3 and earlier */
      int fd = dup(fileno(h->fp));
      h->file_position = lseek(fd, 0, SEEK_CUR);
      close(fd);
    }
  else
#endif
    h->file_position = ftell(h->fp);
  
  return 1;
}

long fasta_get_abundance(fasta_handle h)
{
  return abundance_get(global_abundance, h->header_buffer.data);
}

unsigned long fasta_get_position(fasta_handle h)
{
  return h->file_position;
}

unsigned long fasta_get_size(fasta_handle h)
{
  return h->file_size;
}

unsigned long fasta_get_lineno(fasta_handle h)
{
  return h->lineno;
}

unsigned long fasta_get_seqno(fasta_handle h)
{
  return h->seqno;
}

unsigned long fasta_get_header_length(fasta_handle h)
{
  return h->header_buffer.length;
}

unsigned long fasta_get_sequence_length(fasta_handle h)
{
  return h->sequence_buffer.length;
}

char * fasta_get_header(fasta_handle h)
{
  return h->header_buffer.data;
}

char * fasta_get_sequence(fasta_handle h)
{
  return h->sequence_buffer.data;
}


/* fasta output */

void fasta_print_header(FILE * fp, const char * hdr)
{
  fprintf(fp, ">%s\n", hdr);
}

void fasta_print_sequence(FILE * fp, char * seq, unsigned long len, int width)
{
  /*
     The actual length of the sequence may be longer than "len", but only
     "len" characters are printed.

     Specify width of lines - zero (or <1)  means linearize (all on one line).
  */

  if (width < 1)
    fprintf(fp, "%.*s\n", (int)(len), seq);
  else
    {
      long rest = len;
      for(unsigned long i=0; i<len; i += width)
        {
          fprintf(fp, "%.*s\n", (int)(MIN(rest,width)), seq+i);
          rest -= width;
        }
    }
}

void fasta_print(FILE * fp, const char * hdr,
                 char * seq, unsigned long len)
{
  fasta_print_header(fp, hdr);
  fasta_print_sequence(fp, seq, len, opt_fasta_width);
}


void fasta_print_relabel(FILE * fp,
                         char * seq,
                         int len,
                         char * header,
                         int header_len,
                         int abundance,
                         int ordinal)
{
  if (opt_relabel || opt_relabel_sha1 || opt_relabel_md5)
    {
      fprintf(fp, ">");

      if (opt_relabel_sha1)
        fprint_seq_digest_sha1(fp, seq, len);
      else if (opt_relabel_md5)
        fprint_seq_digest_md5(fp, seq, len);
      else
        fprintf(fp, "%s%d", opt_relabel, ordinal);

      if (opt_sizeout)
        fprintf(fp, ";size=%u;", abundance);

      if (opt_relabel_keep)
        fprintf(fp, " %s", header);

      fprintf(fp, "\n");
    }
  else if (opt_sizeout)
    {
      abundance_fprint_header_with_size(global_abundance,
                                        fp,
                                        header,
                                        header_len,
                                        abundance);
    }
  else if (opt_xsize)
    {
      abundance_fprint_header_strip_size(global_abundance,
                                         fp,
                                         header,
                                         header_len);
    }
  else
    {
      fasta_print_header(fp, header);
    }

  fasta_print_sequence(fp, seq, len, opt_fasta_width);
}

void fasta_print_db_relabel(FILE * fp,
                            unsigned long seqno,
                            int ordinal)
{
  fasta_print_relabel(fp,
                      db_getsequence(seqno),
                      db_getsequencelen(seqno),
                      db_getheader(seqno),
                      db_getheaderlen(seqno),
                      db_getabundance(seqno),
                      ordinal);
}

void fasta_print_db_sequence(FILE * fp, unsigned long seqno)
{
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

void fasta_print_db(FILE * fp, unsigned long seqno)
{
  char * hdr = db_getheader(seqno);
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);

  fasta_print_header(fp, hdr);
  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

void fasta_print_db_size(FILE * fp, unsigned long seqno, unsigned long size)
{
  char * hdr = db_getheader(seqno);
  int hdrlen = db_getheaderlen(seqno);
  
  abundance_fprint_header_with_size(global_abundance,
                                    fp,
                                    hdr,
                                    hdrlen,
                                    size);

  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);

  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

void fasta_print_db_strip_size(FILE * fp, unsigned long seqno)
{
  /* write FASTA but remove abundance information, as with --xsize option */

  char * hdr = db_getheader(seqno);
  int hdrlen = db_getheaderlen(seqno);

  abundance_fprint_header_strip_size(global_abundance,
                                     fp,
                                     hdr,
                                     hdrlen);

  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);

  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

