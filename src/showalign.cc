/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "utils/fatal.hpp"
#include <algorithm>  // std::copy, std::fill_n, std::min
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t
#include <cstdio>  // FILE
#include <iterator>  // std::next
#include <vector>


static int64_t line_pos;

static int64_t q_start;
static int64_t d_start;

static int64_t q_pos;
static int64_t d_pos;

static int64_t alignlen;

std::vector<char> q_line;
std::vector<char> a_line;
std::vector<char> d_line;

static std::FILE * out;

constexpr auto poswidth_default = 3;
static auto poswidth = poswidth_default;
constexpr auto headwidth_default = 5;
static auto headwidth = headwidth_default;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {


  struct Sequence {
    char const * sequence = nullptr;
    int64_t length = 0;
    int64_t offset = 0;
    char const * name = nullptr;
  };


  struct Alignment {
    Sequence query;
    Sequence target;
    int64_t strand = 0;
  };


  inline auto putop(Alignment const & alignment, char const c, int64_t const len) -> void
  {
    int64_t const delta = alignment.strand != 0 ? -1 : +1;

    auto count = len;
    while (count != 0)
      {
        if (line_pos == 0)
          {
            q_start = q_pos;
            d_start = d_pos;
          }

        auto qs = '\0';
        auto ds = '\0';
        auto qs4 = 0U;
        auto ds4 = 0U;

        switch (c)
          {
          case 'M':
            qs = alignment.strand != 0 ? chrmap_complement[static_cast<int>(alignment.query.sequence[q_pos])] : alignment.query.sequence[q_pos];
            ds = alignment.target.sequence[d_pos];
            q_pos += delta;
            d_pos += 1;
            q_line[line_pos] = qs;

            qs4 = chrmap_4bit[static_cast<int>(qs)];
            ds4 = chrmap_4bit[static_cast<int>(ds)];
            if (opt_n_mismatch and ((qs4 == 15) or (ds4 == 15)))
              {
                a_line[line_pos] = ' ';
              }
            else if ((qs4 == ds4) and (ambiguous_4bit[qs4] == 0U))
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
            qs = alignment.strand != 0 ? chrmap_complement[static_cast<int>(alignment.query.sequence[q_pos])] : alignment.query.sequence[q_pos];
            q_pos += delta;
            q_line[line_pos] = qs;
            a_line[line_pos] = ' ';
            d_line[line_pos] = '-';
            ++line_pos;
            break;

          case 'I':
            ds = alignment.target.sequence[d_pos];
            d_pos += 1;
            q_line[line_pos] = '-';
            a_line[line_pos] = ' ';
            d_line[line_pos] = ds;
            ++line_pos;
            break;
          }

        if ((line_pos == alignlen) or ((c == '\0') and (line_pos > 0)))
          {
            q_line[line_pos] = '\0';
            a_line[line_pos] = '\0';
            d_line[line_pos] = '\0';

            int64_t const q1 = std::min(q_start + 1, alignment.query.length);
            int64_t const q2 = alignment.strand != 0 ? q_pos + 2 : q_pos;
            int64_t const d1 = std::min(d_start + 1, alignment.target.length);
            int64_t const d2 = d_pos;

            fprintf(out, "\n");
            fprintf(out, "%*s %*" PRId64 " %c %s %" PRId64 "\n", headwidth, alignment.query.name, poswidth,
                    q1, alignment.strand != 0 ? '-' : '+', q_line.data(), q2);
            fprintf(out, "%*s %*s   %s\n",      headwidth, "",     poswidth,
                    "", a_line.data());
            fprintf(out, "%*s %*" PRId64 " %c %s %" PRId64 "\n", headwidth, alignment.target.name, poswidth,
                    d1, '+', d_line.data(), d2);

            line_pos = 0;
          }
        --count;
      }
  }

}  // end of anonymous namespace


auto align_show(std::FILE * output_handle,
                char const * seq1,
                int64_t const seq1len,
                int64_t const seq1off,
                char const * seq1name,
                char const * seq2,
                int64_t const seq2len,
                int64_t const seq2off,
                char const * seq2name,
                char const * cigar,
                int64_t const cigarlen,
                int const numwidth,
                int const namewidth,
                int const alignwidth,
                int const strand) -> void
{
  out = output_handle;

  Alignment alignment;
  alignment.query.sequence = seq1;
  alignment.query.length = seq1len;
  alignment.query.offset = seq1off;
  alignment.query.name = seq1name;
  alignment.target.sequence = seq2;
  alignment.target.length = seq2len;
  alignment.target.offset = seq2off;
  alignment.target.name = seq2name;
  alignment.strand = strand;

  auto * pos = cigar;
  auto * end = pos + cigarlen;

  poswidth = numwidth;
  headwidth = namewidth;
  alignlen = alignwidth;

  q_line.resize(alignwidth + 1);
  a_line.resize(alignwidth + 1);
  d_line.resize(alignwidth + 1);

  q_pos = strand != 0 ? seq1len - 1 - seq1off : seq1off;
  d_pos = seq2off;

  line_pos = 0;

  while (pos < end)
    {
      int64_t len = 0;
      auto n = 0;
      if (sscanf(pos, "%" PRId64 "%n", &len, &n) == 0)
        {
          n = 0;
          len = 1;
        }
      pos += n;
      auto const op = *pos;
      ++pos;
      putop(alignment, op, len);
    }

  putop(alignment, '\0', 1);

  q_line.clear();
  a_line.clear();
  d_line.clear();
}


auto align_getrow(char const * seq, char const * cigar, int const alignlen, int const origin) -> std::vector<char>
{
  std::vector<char> row(alignlen + 1);
  auto * row_cursor = row.data();
  auto const * cigar_cursor = cigar;
  auto const * seq_cursor = seq;

  while (*cigar_cursor != 0)
    {
      int64_t len = 0;
      auto n = 0;
      if (sscanf(cigar_cursor, "%" PRId64 "%n", &len, &n) == 0)
        {
          n = 0;
          len = 1;
        }
      cigar_cursor += n;
      auto const op = *cigar_cursor;
      ++cigar_cursor;

      if ((op == 'M') or
          ((op == 'D') and (origin == 0)) or
          ((op == 'I') and (origin == 1)))
        {
          std::copy(seq_cursor, std::next(seq_cursor, len), row_cursor);
          row_cursor += len;
          seq_cursor += len;
        }
      else
        {
          /* insert len gap symbols */
          std::fill_n(row_cursor, len, '-');
        }
    }

  // assert(*row_cursor == '\0');  // is not always true
  *row_cursor = '\0';  // needed, not already initialized to null
  return row;
}


auto align_fprint_uncompressed_alignment(std::FILE * output_handle, char const * cigar) -> void
{
  auto * pos = cigar;
  while (*pos != '\0')
    {
      if (*pos > '9')
        {
          fprintf(output_handle, "%c", *pos);
          ++pos;
        }
      else
        {
          auto n = 0;
          auto c = '\0';
          auto x = 0;
          if (sscanf(pos, "%d%c%n", &n, &c, &x) == 2)
            {
              for (auto i = 0; i < n; ++i)
                {
                  fprintf(output_handle, "%c", c);
                }
              pos += x;
            }
          else
            {
              fatal("bad alignment string");
            }
        }
    }
}
