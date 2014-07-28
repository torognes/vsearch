/*
    Copyright (C) 2014 Torbjorn Rognes and Tomas Flouri

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

#define MEMCHUNK 4096
#define LINEALLOC LINE_MAX

extern char map_nt[256];

extern FILE * query_fp;

static BZFILE * query_bz_fp;

extern char query_line[LINEALLOC];

extern long query_no;

extern char * query_head;
extern char * query_seq;

extern unsigned long query_head_len;
extern unsigned long query_seq_len;

extern unsigned long query_head_alloc;
extern unsigned long query_seq_alloc;

static unsigned long query_line_len;
static int bzerror;

void query_bz_open(const char * filename)
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

  /* open bz2 compressed fasta file */

  query_fp = NULL;
  query_fp = fopen(filename, "rb");
  if (!query_fp)
    fatal("Error: Unable to open query file (%s)", filename);

  query_bz_fp = BZ2_bzReadOpen(&bzerror, query_fp, BZ_VERBOSE_0, BZ_MORE_MEM, NULL, 0);
  if (!query_bz_fp)
    fatal("Error: Unable to open query file (%s)", filename);

  query_line_len = BZ2_bzRead(&bzerror, query_bz_fp, query_line, LINEALLOC - 1);
  query_line[query_line_len] = 0;

}

void query_bz_close()
{
  BZ2_bzReadClose(&bzerror, query_bz_fp);
  fclose(query_fp);
  
  if (query_seq)
    free(query_seq);
  if (query_head)
    free(query_head);

  query_head = 0;
  query_seq = 0;
}

int query_bz_getnext(char ** head, long * head_len,
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

      /* read the remaining characters from the current (fasta header) line */
      while (!strchr(query_line,'\n') && (query_line_len == LINEALLOC - 1))
         query_line_len = BZ2_bzRead(&bzerror, query_bz_fp, query_line, LINEALLOC - 1);

      /* move the contents of the next line (fasta data) to the beginning of 
         query_line  and fill the remaining empty buffer with new data */
      headerlen = xstrchrnul(query_line, '\n') - (query_line);
      unsigned long remaining = query_line_len - headerlen - 1;
      memmove(query_line, query_line + headerlen + 1, remaining);
      if (query_line_len == LINEALLOC - 1)
       {
         query_line_len = BZ2_bzRead(&bzerror, query_bz_fp, query_line + remaining, LINEALLOC - remaining - 1);
         query_line_len += remaining;
       }
      else
        query_line_len = remaining;
      query_line[query_line_len] = 0;

      
      /* read sequence */

      query_seq_len = 0;

      char c;
      char m;


      char prev = 0;
      while (1)
       {
         char * p = query_line;
         while ((c = *p++))
          {
            if ((m = map_nt[(int)c]) >= 0)
             {
               if (query_seq_len + 1 > query_seq_alloc)
                {
                  query_seq_alloc += MEMCHUNK;
                  query_seq = (char *) xrealloc(query_seq, query_seq_alloc);
                }

               *(query_seq + query_seq_len) = c; //m;
               query_seq_len++;
             }
            else 
             {
               if (c != '\n' && c != '>')
                 fatal("Illegal character in sequence.");
               if (c == '>' && prev != '\n')
                 fatal("Illegal character in sequence.");
             }
            prev = c;
            if (c == '>') break;
          }

         /* two cases here - either we read more data from the bz2 file, or we
            don't in case we reached the end-of-file already */
         if (query_line_len == LINEALLOC - 1)
          {
            if (!c)
             {
               query_line_len = BZ2_bzRead(&bzerror, query_bz_fp, query_line, LINEALLOC - 1);
               query_line[query_line_len] = 0;
             }
            else
             {
               unsigned long datalen = (p-1) - query_line;
               memmove (query_line, query_line + datalen, query_line_len - datalen);
               query_line_len -= datalen;
               query_line_len += BZ2_bzRead(&bzerror, query_bz_fp, query_line + query_line_len, LINEALLOC - query_line_len - 1);
               break;
             }
          }
         else
          {
            unsigned long datalen = (p-1) - query_line;
            memmove (query_line, query_line + datalen, query_line_len - datalen);
            query_line_len -= datalen;
            break;
          }
       }
      if (!query_line_len) query_line[0] = 0;
      query_seq[query_seq_len] = 0;
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
