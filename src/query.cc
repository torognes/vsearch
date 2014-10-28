/*
    Copyright (C) 2014 Torbjorn Rognes & Tomas Flouri

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

/* please note that these functions will return a pointer to a buffer
   allocated here for the query header and sequence. This buffers will
   be overwritten on the next call of query_getnext. */

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#define MEMCHUNK 4096
#define LINEALLOC LINE_MAX

extern unsigned int chrstatus[256];

static FILE * query_fp;
static char query_line[LINEALLOC];

static int query_no = -1;

static char * query_head = 0;
static char * query_seq = 0;

static long query_head_len = 0;
static long query_seq_len = 0;

static long query_head_alloc = 0;
static long query_seq_alloc = 0;

static long query_filesize = 0;
static int query_format = FORMAT_PLAIN;

static int query_lineno;

static int query_stripped_count;
static int query_stripped[256];

#ifdef HAVE_BZLIB
static BZFILE * bz_query_fp;
static int bz_error;
static char bz_buffer[LINEALLOC];
static long bz_buffer_len = 0;
#endif

#ifdef HAVE_ZLIB
gzFile gz_query_fp;
#endif

regex_t q_regexp;

long query_getfilesize()
{
  return query_filesize;
}

long query_getfilepos()
{
  return ftell(query_fp);
}

static char * FGETS(char * query_line, int size)
{
  switch (query_format)
   {
     case FORMAT_PLAIN:
       fgets(query_line, size, query_fp);
       break;
     case FORMAT_BZIP:
#ifdef HAVE_BZLIB
       bz_fgets(query_line, size, bz_query_fp, size,
                &bz_error, bz_buffer, &bz_buffer_len);
#else
       fatal("Error: Query file seems to be BZIPx compressed, but %s was not compiled with BZLIB support", PROG_NAME);
#endif
       break;
     case FORMAT_GZIP:
#ifdef HAVE_ZLIB
       gzgets(gz_query_fp, query_line, size);
       break;
#else
       fatal("Error: Database file seems to be GZIP compressed, but %s was not "
             "compiled with ZLIB support", PROG_NAME);
#endif
     default:
       fatal("Error: Unknown compression type detected");
   }
  return query_line;
}


void query_open(const char * filename)
{
  if (regcomp(&q_regexp, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
    fatal("Regular expression compilation failed");
  
  //unsigned long query_line_len;
  /* allocate space */

  query_head = NULL;
  query_seq = NULL;

  query_head_len = 0;
  query_seq_len = 0;

  query_head_alloc = MEMCHUNK;
  query_seq_alloc = MEMCHUNK;

  query_head = (char *) xmalloc((size_t)query_head_alloc);
  query_seq = (char *) xmalloc((size_t)query_seq_alloc);

  query_no = -1;

  /* detect compression type (if any) */
  query_format = detect_compress_format(filename);
  if (!query_format)
    fatal("Error: Unable to read from query file (%s)", filename);

  /* open query file */
  query_fp = NULL;
  query_fp = fopen(filename, "r");
  if (!query_fp)
    fatal("Error: Unable to open query file (%s)", filename);

  if (fseek(query_fp, 0, SEEK_END))
    fatal("Error: Unable to seek in query file (%s)", filename);

  query_filesize = ftell(query_fp);
  
  rewind(query_fp);

#ifdef HAVE_BZLIB
  /* open appropriate data steam if input file was compressed with bzip */
  if (query_format == FORMAT_BZIP)
   {
     bz_query_fp = BZ2_bzReadOpen(&bz_error, query_fp, 
                                  BZ_VERBOSE_0, BZ_MORE_MEM, NULL, 0);
     if (!bz_query_fp)
       fatal("Error: Unable to open query file (%s)", filename);
   }
#endif
#ifdef HAVE_ZLIB
  if (query_format == FORMAT_GZIP)
   {
     fclose(query_fp);
     gz_query_fp = gzopen(filename, "r");
     if (!gz_query_fp)
       fatal("Error: Unable to open query file (%s)", filename);
   }
#endif
  
  query_line[0] = 0;
  FGETS(query_line, LINEALLOC);
  query_lineno = 1;

  query_stripped_count = 0;
  for(int i=0; i<256; i++)
    query_stripped[i] = 0;
}

void query_close()
{
  /* Warn about stripped chars */

  if (query_stripped_count)
    {
      fprintf(stderr, "Warning: invalid characters stripped from query:");
      for (int i=0; i<256;i++)
        if (query_stripped[i])
          fprintf(stderr, " %c(%d)", i, query_stripped[i]);
      fprintf(stderr, "\n");
    }
#ifdef HAVE_BZLIB
  if (query_format == FORMAT_BZIP)
    BZ2_bzReadClose(&bz_error, bz_query_fp);
#endif
#ifdef HAVE_ZLIB
  if (query_format == FORMAT_GZIP)
    gzclose(gz_query_fp);
#endif

  if (query_format != FORMAT_GZIP)
    fclose(query_fp);
  
  if (query_seq)
    free(query_seq);
  if (query_head)
    free(query_head);

  regfree(&q_regexp);

  query_head = 0;
  query_seq = 0;
}

int query_getnext(char ** head, int * head_len,
                  char ** seq, int * seq_len, int * qno,
                  int * qsize, int upcase)
{
  while (query_line[0])
    {
      /* read header */

      if (query_line[0] != '>')
        fatal("Illegal header line in query fasta file");
      
      long headerlen = xstrchrnul(query_line+1, '\n') - (query_line+1);
      query_head_len = headerlen;

      if (headerlen + 1 > query_head_alloc)
        {
          query_head_alloc = headerlen + 1;
          query_head = (char *) xrealloc(query_head, (size_t)query_head_alloc);
        }

      memcpy(query_head, query_line + 1, (size_t)headerlen);
      query_head[headerlen] = 0;

      /* read size/abundance annotation */

      regmatch_t pmatch[4];

      if (!regexec(&q_regexp, query_head, 4, pmatch, 0))
        {
          unsigned long size = atol(query_head + pmatch[2].rm_so);
          if (size > 0)
            * qsize = size;
          else
            fatal("Size annotation zero in query sequence");
        }
      else
        *qsize = 1;

      /* get next line */

      query_line[0] = 0;
      FGETS(query_line, LINEALLOC);
      query_lineno++;

      /* read sequence */

      query_seq_len = 0;

      while (query_line[0] && (query_line[0] != '>'))
        {
          char c;
          char m;
          char * p = query_line;

          while((c = *p++))
            {
              m = chrstatus[(int)c];
              switch(m)
                {
                case 0:
                  /* character to be stripped */
                  query_stripped_count++;
                  query_stripped[(int)c]++;
                  break;

                case 1:
                  /* legal character */
                  if (query_seq_len + 1 > query_seq_alloc)
                    {
                      query_seq_alloc += MEMCHUNK;
                      query_seq = (char *) xrealloc(query_seq, (size_t)query_seq_alloc);
                    }
                  if (upcase)
                    c &= 0xdf;
                  *(query_seq + query_seq_len) = c;
                  query_seq_len++;

                  break;

                case 2:
                  /* fatal character */
                  char msg[200];
                  if (c>=32)
                    snprintf(msg, 200, "illegal character '%c' on line %d in the query file", c, query_lineno);
                  else
                    snprintf(msg, 200, "illegal unprintable character %#.2x (hexadecimal) on line %d in the query file", c, query_lineno);
                  fatal(msg);
                  break;

                case 3:
                  /* silently stripped chars */
                  break;

                }
            }

          query_line[0] = 0;
          FGETS(query_line, LINEALLOC);
          query_lineno++;
        }

      /* add zero after sequence */

      if (query_seq_len + 1 > query_seq_alloc)
        {
          query_seq_alloc += MEMCHUNK;
          query_seq = (char *) xrealloc(query_seq, (size_t)query_seq_alloc);
        }
      *(query_seq + query_seq_len) = 0;




      query_no++;
      *head = query_head;
      *seq = query_seq;
      *head_len = query_head_len;
      *seq_len = query_seq_len;
      *qno = query_no;

      return 1;
    }
  
  return 0;
}


