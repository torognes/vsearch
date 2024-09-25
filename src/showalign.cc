/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "maps.h"
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t
#include <cstdio>  // FILE
#include <cstring>  // std::strncpy


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

static std::FILE * out;

constexpr int poswidth_default {3};
static int poswidth = poswidth_default;
constexpr int headwidth_default {5};
static int headwidth = headwidth_default;

static const char * q_name;
static const char * d_name;

static int64_t q_len;
static int64_t d_len;

inline auto putop(char c, int64_t len) -> void
{
  const int64_t delta = q_strand != 0 ? -1 : +1;

  int64_t count = len;
  while (count != 0)
    {
      if (line_pos == 0)
        {
          q_start = q_pos;
          d_start = d_pos;
        }

      char qs = '\0';
      char ds = '\0';
      unsigned int qs4 = 0;
      unsigned int ds4 = 0;

      switch(c)
        {
        case 'M':
          qs = q_strand != 0 ? chrmap_complement[static_cast<int>(q_seq[q_pos])] : q_seq[q_pos];
          ds = d_seq[d_pos];
          q_pos += delta;
          d_pos += 1;
          q_line[line_pos] = qs;

          qs4 = chrmap_4bit[static_cast<int>(qs)];
          ds4 = chrmap_4bit[static_cast<int>(ds)];
          if ((qs4 == ds4) and (ambiguous_4bit[qs4] == 0U))
            {
              a_line[line_pos] = '|';
            }
          else if ((qs4 & ds4) != 0U)
            {
              a_line[line_pos] = '+';
            }
          else
            {
              a_line[line_pos] = ' ';
            }

          d_line[line_pos] = ds;
          ++line_pos;
          break;

        case 'D':
          qs = q_strand != 0 ? chrmap_complement[static_cast<int>(q_seq[q_pos])] : q_seq[q_pos];
          q_pos += delta;
          q_line[line_pos] = qs;
          a_line[line_pos] = ' ';
          d_line[line_pos] = '-';
          ++line_pos;
          break;

        case 'I':
          ds = d_seq[d_pos];
          d_pos += 1;
          q_line[line_pos] = '-';
          a_line[line_pos] = ' ';
          d_line[line_pos] = ds;
          ++line_pos;
          break;
        }

      if ((line_pos == alignlen) or ((c == 0) and (line_pos > 0)))
        {
          q_line[line_pos] = 0;
          a_line[line_pos] = 0;
          d_line[line_pos] = 0;

          const int64_t q1 = q_start + 1 > q_len ? q_len : q_start + 1;
          const int64_t q2 = q_strand != 0 ? q_pos + 2 : q_pos;
          const int64_t d1 = d_start + 1 > d_len ? d_len : d_start + 1;
          const int64_t d2 = d_pos;

          fprintf(out, "\n");
          fprintf(out, "%*s %*" PRId64 " %c %s %" PRId64 "\n", headwidth, q_name, poswidth,
                  q1, q_strand != 0 ? '-' : '+', q_line, q2);
          fprintf(out, "%*s %*s   %s\n",      headwidth, "",     poswidth,
                  "", a_line);
          fprintf(out, "%*s %*" PRId64 " %c %s %" PRId64 "\n", headwidth, d_name, poswidth,
                  d1, '+', d_line, d2);

          line_pos = 0;
        }
      --count;
    }
}

auto align_show(std::FILE * output_handle,
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
                int strand) -> void
{
  out = output_handle;

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

  q_line = (char *) xmalloc(alignwidth + 1);
  a_line = (char *) xmalloc(alignwidth + 1);
  d_line = (char *) xmalloc(alignwidth + 1);

  q_pos = strand != 0 ? seq1len - 1 - seq1off : seq1off;
  d_pos = seq2off;

  line_pos = 0;

  while (p < e)
    {
      int64_t len = 0;
      int n = 0;
      if (sscanf(p, "%" PRId64 "%n", & len, & n) == 0)
        {
          n = 0;
          len = 1;
        }
      p += n;
      const char op = *p++;
      putop(op, len);
    }

  putop(0, 1);

  xfree(q_line);
  xfree(a_line);
  xfree(d_line);
}

auto align_getrow(char * seq, char * cigar, int alignlen, int origin) -> char *
{
  char * row = (char *) xmalloc(alignlen + 1);
  char * r = row;
  char * p = cigar;
  char * s = seq;

  while (*p != 0)
    {
      int64_t len = 0;
      int n = 0;
      if (sscanf(p, "%" PRId64 "%n", & len, & n) == 0)
        {
          n = 0;
          len = 1;
        }
      p += n;
      const char op = *p++;

      if ((op == 'M') or
          ((op == 'D') and (origin == 0)) or
          ((op == 'I') and (origin == 1)))
        {
          strncpy(r, s, len);
          r += len;
          s += len;
        }
      else
        {
          /* insert len gap symbols */
          for (int64_t i = 0; i < len; i++)
            {
              *r++ = '-';
            }
        }
    }

  *r = 0;
  return row;
}

auto align_fprint_uncompressed_alignment(std::FILE * output_handle, char * cigar) -> void
{
  char * p = cigar;
  while (*p != 0)
    {
      if (*p > '9')
        {
          fprintf(output_handle, "%c", *p++);
        }
      else
        {
          int n = 0;
          char c = 0;
          int x = 0;
          if (sscanf(p, "%d%c%n", &n, &c, &x) == 2)
            {
              for (int i = 0; i < n; i++)
                {
                  fprintf(output_handle, "%c", c);
                }
              p += x;
            }
          else
            {
              fatal("bad alignment string");
            }
        }
    }
}
