/*
    Copyright (C) 2014 Torbjorn Rognes

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


/*

char map_nt[256] =
  {
    // N = A

    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,  0, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,  0, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

*/

#define MEMCHUNK 4096
#define LINEALLOC LINE_MAX

extern char map_nt[256];

FILE * query_fp;
char query_line[LINEALLOC];

long query_no = -1;

char * query_head = 0;
char * query_seq = 0;

unsigned long query_head_len = 0;
unsigned long query_seq_len = 0;

unsigned long query_head_alloc = 0;
unsigned long query_seq_alloc = 0;

long query_filesize = 0;

long query_getfilesize()
{
  return query_filesize;
};

long query_getfilepos()
{
  return ftell(query_fp);
};

void query_open(const char * filename)
{
  /* allocate space */

  query_head = NULL;
  query_seq = NULL;

  query_head_len = 0;
  query_seq_len = 0;

  query_head_alloc = MEMCHUNK;
  query_seq_alloc = MEMCHUNK;

  query_head = (char *) xmalloc(query_head_alloc);
  query_seq = (char *) xmalloc(query_seq_alloc);

  query_no = -1;

  /* open fasta file */

  query_fp = NULL;
  query_fp = fopen(filename, "r");
  if (!query_fp)
    fatal("Error: Unable to open query file (%s)", filename);

  if (fseek(query_fp, 0, SEEK_END))
    fatal("Error: Unable to seek in query file (%s)", filename);

  query_filesize = ftell(query_fp);
  
  if (fseek(query_fp, 0, SEEK_SET))
    fatal("Error: Unable to seek in query file (%s)", filename);

  query_line[0] = 0;
  fgets(query_line, LINEALLOC, query_fp);
}

void query_close()
{
  fclose(query_fp);
  
  if (query_seq)
    free(query_seq);
  if (query_head)
    free(query_head);

  query_head = 0;
  query_seq = 0;
}

int query_getnext(char ** head, long * head_len,
                  char ** seq, long * seq_len, long * qno)
{
  if(query_line[0])
    {
      /* read header */

      if (query_line[0] != '>')
        fatal("Illegal header line in query fasta file");
      
      unsigned long headerlen = xstrchrnul(query_line+1, '\n') - (query_line+1);
      query_head_len = headerlen;

      if (headerlen + 1 > query_head_alloc)
        {
          query_head_alloc = headerlen + 1;
          query_head = (char *) xrealloc(query_head, query_head_alloc);
        }

      memcpy(query_head, query_line + 1, headerlen);
      query_head[headerlen] = 0;


      /* get next line */

      query_line[0] = 0;
      fgets(query_line, LINEALLOC, query_fp);


      /* read sequence */

      query_seq_len = 0;

      while (query_line[0] && (query_line[0] != '>'))
        {
          char c;
          char m;

          char * p = query_line;

          while((c = *p++))
            if ((m = map_nt[(int)c]) >= 0)
            {
              if (query_seq_len + 1 > query_seq_alloc)
              {
                query_seq_alloc += MEMCHUNK;
                query_seq = (char *) xrealloc(query_seq, query_seq_alloc);
              }

              *(query_seq + query_seq_len) = m;
              query_seq_len++;
            }
            else if (c != '\n')
              fatal("Illegal character in sequence.");

          query_line[0] = 0;
          fgets(query_line, LINEALLOC, query_fp);
        }

      query_no++;

      *head = query_head;
      *seq = query_seq;
      *head_len = query_head_len;
      *seq_len = query_seq_len;
      *qno = query_no;

      return 1;
    }
  else
    return 0;
}


