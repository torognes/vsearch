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

int fastx_detect(const char * filename)
{
#ifdef HAVE_ZLIB
  gzFile fp_gz = 0;
#endif

#ifdef HAVE_BZLIB
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

  rewind(fp);

  if (format == FORMAT_GZIP)
    {
      /* GZIP: Close ordinary file and open again as gzipped file */
#ifdef HAVE_ZLIB
      fclose(fp);
      if (! (fp_gz = gzopen(filename, "rb")))
        fatal("Unable to open gzip compressed file (%s)", filename);
#else
      fatal("Files compressed with gzip are not supported");
#endif
    }

  if (format == FORMAT_BZIP)
    {
      /* BZIP2: Keep original file open, then open as bzipped file as well */
#ifdef HAVE_ZLIB
      int bzError;
      if (! (fp_bz = BZ2_bzReadOpen(& bzError, fp,
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
  
#ifdef HAVE_BZLIB
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
#ifdef HAVE_ZLIB
      bytes_read = gzread(fp_gz,
                          buffer,
                          BUFFERLEN);
      if (bytes_read < 0)
        fatal("Error reading gzip compressed file (%s)", filename);
      break;
#endif
      
    case FORMAT_BZIP:
#ifdef HAVE_BZLIB
      bytes_read = BZ2_bzRead(& bzError,
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

#ifdef HAVE_BZLIB
  int bz_error;
#endif
  
  switch(format)
    {
    case FORMAT_PLAIN:
      fclose(fp);
      fp = 0;
      break;

    case FORMAT_GZIP:
#ifdef HAVE_ZLIB
      gzclose(fp_gz);
      fp_gz = 0;
      break;
#endif
      
    case FORMAT_BZIP:
#ifdef HAVE_BZLIB
      BZ2_bzReadClose(&bz_error, fp_bz);
      fp_bz = 0;
      break;
#endif

    default:
      fatal("Internal error");
    }

  return filetype;
}
