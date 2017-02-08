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

static int64_t line_pos;

static char * q_seq;
static char * d_seq;

static int64_t q_start;
static int64_t d_start;

static int64_t q_pos;
static int64_t d_pos;

static int64_t q_strand;

static int64_t alignlen;

static char * q_line;
static char * a_line;
static char * d_line;

static FILE * out;

static int poswidth = 3;
static int headwidth = 5;

static const char * q_name;
static const char * d_name;

static int64_t q_len;
static int64_t d_len;

inline int nt_identical(char a, char b)
{
  if (chrmap_4bit[(int)a] == chrmap_4bit[(int)b])
    return 1;
  else
    return 0;
}

inline void putop(char c, int64_t len)
{
  int64_t delta = q_strand ? -1 : +1;

  int64_t count = len;
  while(count)
    {
      if (line_pos == 0)
        {
          q_start = q_pos;
          d_start = d_pos;
        }

      char qs;
      char ds;

      switch(c)
        {
        case 'M':
          qs = q_strand ? chrmap_complement[(int)(q_seq[q_pos])] : q_seq[q_pos];
          ds = d_seq[d_pos];
          q_pos += delta;
          d_pos += 1;
          q_line[line_pos] = qs;
          a_line[line_pos] = nt_identical(qs, ds) ? '|' : ' ';
          d_line[line_pos] = ds;
          line_pos++;
          break;

        case 'D':
          qs = q_strand ? chrmap_complement[(int)(q_seq[q_pos])] : q_seq[q_pos];
          q_pos += delta;
          q_line[line_pos] = qs;
          a_line[line_pos] = ' ';
          d_line[line_pos] = '-';
          line_pos++;
          break;

        case 'I':
          ds = d_seq[d_pos];
          d_pos += 1;
          q_line[line_pos] = '-';
          a_line[line_pos] = ' ';
          d_line[line_pos] = ds;
          line_pos++;
          break;
        }

      if ((line_pos == alignlen) || ((c == 0) && (line_pos > 0)))
        {
          q_line[line_pos] = 0;
          a_line[line_pos] = 0;
          d_line[line_pos] = 0;

          int64_t q1 = q_start + 1;
          if (q1 > q_len)
            q1 = q_len;

          int64_t q2 = q_strand ? q_pos +2 : q_pos;

          int64_t d1 = d_start + 1;
          if (d1 > d_len)
            d1 = d_len;

          int64_t d2 = d_pos;

          fprintf(out, "\n");
          fprintf(out, "%*s %*" PRId64 " %c %s %" PRId64 "\n", headwidth, q_name, poswidth,
                  q1, q_strand ? '-' : '+', q_line, q2);
          fprintf(out, "%*s %*s   %s\n",      headwidth, "",     poswidth,
                  "", a_line);
          fprintf(out, "%*s %*" PRId64 " %c %s %" PRId64 "\n", headwidth, d_name, poswidth,
                  d1, '+', d_line, d2);

          line_pos = 0;
        }
      count--;
    }
}

void align_show(FILE * f,
                char * seq1,
                int64_t seq1len,
                int64_t seq1off,
                const char * seq1name,
                char * seq2,
                int64_t seq2len,
                int64_t seq2off,
                const char * seq2name,
                char * cigar,
                int64_t cigarlen,
                int numwidth,
                int namewidth,
                int alignwidth,
                int strand)
{
  out = f;

  q_seq = seq1;
  q_len = seq1len;
  q_name = seq1name;
  q_strand = strand;
  
  d_seq = seq2;
  d_len = seq2len;
  d_name = seq2name;
  
  char * p = cigar;
  char * e = p + cigarlen;
  
  poswidth = numwidth;
  headwidth = namewidth;
  alignlen = alignwidth;

  q_line = (char*) xmalloc(alignwidth+1);
  a_line = (char*) xmalloc(alignwidth+1);
  d_line = (char*) xmalloc(alignwidth+1);

  q_pos = strand ? seq1len - 1 - seq1off : seq1off;
  d_pos = seq2off;
  
  line_pos = 0;

  while(p < e)
    {
      int64_t len;
      int n;
      if (!sscanf(p, "%" PRId64 "%n", & len, & n))
        {
          n = 0;
          len = 1;
        }
      p += n;
      char op = *p++;
      putop(op, len);
    }
  
  putop(0, 1);

  xfree(q_line);
  xfree(a_line);
  xfree(d_line);
}
               
char * align_getrow(char * seq, char * cigar, int alen, int origin)
{
  char * row = (char*) xmalloc(alen+1);
  char * r = row;
  char * p = cigar;
  char * s = seq;

  while(*p)
    {
      int64_t len;
      int n;
      if (!sscanf(p, "%" PRId64 "%n", & len, & n))
        {
          n = 0;
          len = 1;
        }
      p += n;
      char op = *p++;
      
      if ((op == 'M') || 
          ((op == 'D') && (origin == 0)) ||
          ((op == 'I') && (origin == 1)))
        {
          strncpy(r, s, len);
          r += len;
          s += len;
        }
      else
        {
          /* insert len gap symbols */
          for(int64_t i = 0; i < len; i++)
            *r++ = '-';
        }
    }

  *r = 0;
  return row;
}

void align_fprint_uncompressed_alignment(FILE * f, char * cigar)
{
  char * p = cigar;
  while(*p)
    {
      if (*p > '9')
        fprintf(f, "%c", *p++);
      else
        {
          int n = 0;
          char c = 0;
          int x = 0;
          if (sscanf(p, "%d%c%n", &n, &c, &x) == 2)
            {
              for(int i = 0; i<n; i++)
                fprintf(f, "%c", c);
              p += x;
            }
          else
            fatal("bad alignment string");
        }
    }
}
