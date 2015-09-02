/*
    Copyright (C) 2014-2015 Torbjorn Rognes & Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
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


size_t fastq_file_fill_buffer(fastq_handle h)
{
  /* read more data if necessary */
  size_t rest = h->file_buffer.length - h->file_buffer.position;
  
  if (rest > 0)
    return rest;
  else
    {
      size_t space = h->file_buffer.alloc - h->file_buffer.length;

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

void buffer_extend(struct fastq_buffer_s * buffer, char * buf, size_t len)
{
  if (buffer->length + len + 1 > buffer->alloc)
    {
      /* alloc space for len more characters + terminating zero,
         but round up to nearest block size */
      buffer->alloc = 
        (FASTQ_BUFFER_ALLOC * 
         ((buffer->length + len) / FASTQ_BUFFER_ALLOC) + 1);
      buffer->data = (char*) xrealloc(buffer->data, buffer->alloc);
    }

  /* copy string */
  memcpy(buffer->data + buffer->length, buf, len);
  buffer->length += len;

  /* add terminator */
  buffer->data[buffer->length] = 0;
}

void fastq_truncate_header(fastq_handle h, bool truncateatspace)
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

void fastq_filter_string(fastq_handle h,
                         struct fastq_buffer_s * buffer,
                         unsigned int * char_action,
                         char * char_mapping)
{
  /* Strip unwanted characters from the string and raise warnings or
     errors on certain characters. */

  char * p = buffer->data;
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
          if (c>=32)
            snprintf(msg,
                     200,
                     "illegal character '%c' on line %lu in fastq file",
                     c,
                     h->lineno);
          else
            snprintf(msg,
                     200,
                     "illegal character %#.2x (hexadecimal) on line %lu in fastq file",
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
  buffer->length = q - buffer->data;
}

bool fastq_next(fastq_handle h,
                bool truncateatspace,
                unsigned int * char_action,
                char * char_mapping)
{
  h->header_buffer.length = 0;
  h->sequence_buffer.length = 0;
  h->quality_buffer.length = 0;

  size_t rest = fastq_file_fill_buffer(h);

  /* check end of file */

  if (rest == 0)
    return 0;

  /* read header */

  /* check initial > character */
  
  if (h->file_buffer.data[h->file_buffer.position] != '@')
    fatal("Invalid FASTQ - header must start with '@' character");
  h->file_buffer.position++;
  rest--;

  char * lf = 0;
  while (lf == 0)
    {
      /* get more data if buffer empty */
      rest = fastq_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTQ - header must be terminated with newline");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);

      /* copy to header buffer */
      size_t len = rest;
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

  /* read sequence line */

  lf = 0;
  while (lf == 0)
    {
      /* get more data, if necessary */
      rest = fastq_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTQ file: sequence line must be terminated with newline");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);
      
      /* copy to sequence buffer */
      size_t len = rest;
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

  /* read + line */
  
  rest = fastq_file_fill_buffer(h);
  if (rest == 0)
    fatal("Invalid FASTQ - lacking third line in entry (+ line)");
  
  if (h->file_buffer.data[h->file_buffer.position] != '+')
    fatal("Invalid FASTQ - third line in entry must start with '+' character");
  h->file_buffer.position++;
  rest--;
  
  lf = 0;
  while (lf == 0)
    {
      /* get more data if buffer empty */
      rest = fastq_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTQ - third entry line must be terminated with newline");
      
      /* ignore rest of line */
      /* find LF */
      size_t len = rest;
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);
      if (lf)
        {
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }
      h->file_buffer.position += len;
      rest -= len;
    }

  /* read quality line */

  lf = 0;
  while (lf == 0)
    {
      /* get more data, if necessary */
      rest = fastq_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTQ file: quality line must be terminated with newline");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);
      
      /* copy to quality buffer */
      size_t len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
        }
      buffer_extend(& h->quality_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  if (h->sequence_buffer.length != h->quality_buffer.length)
    fatal("Invalid FASTQ - Different length of sequence and quality lines");

  h->seqno++;

  fastq_truncate_header(h, truncateatspace);
  fastq_filter_string(h, & h->sequence_buffer, char_action, char_mapping);
  fastq_filter_string(h, & h->quality_buffer, char_qual_action, map_identity);

#ifdef HAVE_ZLIB
  if (h->format == FORMAT_GZIP)
    h->file_position = gzoffset(h->fp_gz);
  else
#endif
    h->file_position = ftell(h->fp);
  
  return 1;
}

char * fastq_get_quality(fastq_handle h)
{
  return h->quality_buffer.data;
}

size_t fastq_get_quality_length(fastq_handle h)
{
  return h->quality_buffer.length;
}

size_t fastq_get_position(fastq_handle h)
{
  return h->file_position;
}

size_t fastq_get_size(fastq_handle h)
{
  return h->file_size;
}

size_t fastq_get_lineno(fastq_handle h)
{
  return h->lineno;
}

size_t fastq_get_seqno(fastq_handle h)
{
  return h->seqno;
}

size_t fastq_get_header_length(fastq_handle h)
{
  return h->header_buffer.length;
}

size_t fastq_get_sequence_length(fastq_handle h)
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
