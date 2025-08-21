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
#include "utils/cigar.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/span.hpp"
#include <algorithm>  // std::copy, std::fill_n, std::min
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE
#include <cstring>  // std::strlen
#include <vector>


std::vector<char> q_line;
std::vector<char> a_line;
std::vector<char> d_line;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {


  struct Position {
    int64_t line = 0;
    int64_t query = 0;
    int64_t target = 0;
    int64_t query_start = 0;
    int64_t target_start = 0;
  };


  struct Sequence {
    char const * sequence = nullptr;
    int64_t length = 0;
    int64_t offset = 0;
    char const * name = nullptr;
  };


  struct Alignment {
    static constexpr auto poswidth_default = 3;
    static constexpr auto headwidth_default = 5;
    std::FILE * output_handle = nullptr;
    Sequence query;
    Sequence target;
    int poswidth = poswidth_default;
    int headwidth = headwidth_default;
    int64_t width = 0;
    int64_t strand = 0;
  };


  inline auto putop(Alignment const & alignment, Position & position, char const operation, int64_t const len) -> void
  {
    static constexpr auto is_N = 15U;
    int64_t const delta = alignment.strand != 0 ? -1 : +1;

    auto count = len;
    while (count != 0)
      {
        if (position.line == 0)
          {
            position.query_start = position.query;
            position.target_start = position.target;
          }

        auto qs = '\0';
        auto ds = '\0';
        auto qs4 = 0U;
        auto ds4 = 0U;

        switch (operation)
          {
          case 'M':
            qs = alignment.strand != 0 ? map_complement(alignment.query.sequence[position.query]) : alignment.query.sequence[position.query];
            ds = alignment.target.sequence[position.target];
            position.query += delta;
            position.target += 1;
            q_line[position.line] = qs;

            qs4 = map_4bit(qs);
            ds4 = map_4bit(ds);
            if (opt_n_mismatch and ((qs4 == is_N) or (ds4 == is_N)))
              {
                a_line[position.line] = ' ';
              }
            else if ((qs4 == ds4) and not is_ambiguous_4bit[qs4])
              {
                a_line[position.line] = '|';
              }
            else if ((qs4 & ds4) != 0U)
              {
                a_line[position.line] = '+';
              }
            else
              {
                a_line[position.line] = ' ';
              }

            d_line[position.line] = ds;
            ++position.line;
            break;

          case 'D':
            qs = alignment.strand != 0 ? map_complement(alignment.query.sequence[position.query]) : alignment.query.sequence[position.query];
            position.query += delta;
            q_line[position.line] = qs;
            a_line[position.line] = ' ';
            d_line[position.line] = '-';
            ++position.line;
            break;

          case 'I':
            ds = alignment.target.sequence[position.target];
            position.target += 1;
            q_line[position.line] = '-';
            a_line[position.line] = ' ';
            d_line[position.line] = ds;
            ++position.line;
            break;
          }

        // refactor: extract to a putop_final()?
        if ((position.line == alignment.width) or ((operation == '\0') and (position.line > 0)))
          {
            q_line[position.line] = '\0';
            a_line[position.line] = '\0';
            d_line[position.line] = '\0';

            int64_t const q1 = std::min(position.query_start + 1, alignment.query.length);
            int64_t const q2 = alignment.strand != 0 ? position.query + 2 : position.query;
            int64_t const d1 = std::min(position.target_start + 1, alignment.target.length);
            int64_t const d2 = position.target;

            fprintf(alignment.output_handle, "\n");
            fprintf(alignment.output_handle, "%*s %*" PRId64 " %c %s %" PRId64 "\n", alignment.headwidth, alignment.query.name, alignment.poswidth,
                    q1, alignment.strand != 0 ? '-' : '+', q_line.data(), q2);
            fprintf(alignment.output_handle, "%*s %*s   %s\n",      alignment.headwidth, "",     alignment.poswidth,
                    "", a_line.data());
            fprintf(alignment.output_handle, "%*s %*" PRId64 " %c %s %" PRId64 "\n", alignment.headwidth, alignment.target.name, alignment.poswidth,
                    d1, '+', d_line.data(), d2);

            position.line = 0;
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

  Alignment alignment;
  alignment.output_handle = output_handle;
  alignment.query.sequence = seq1;
  alignment.query.length = seq1len;
  alignment.query.offset = seq1off;
  alignment.query.name = seq1name;
  alignment.target.sequence = seq2;
  alignment.target.length = seq2len;
  alignment.target.offset = seq2off;
  alignment.target.name = seq2name;
  alignment.poswidth = numwidth;
  alignment.headwidth = namewidth;
  alignment.width = alignwidth;
  alignment.strand = strand;

  // C++14 refactoring: aggregate initialization of a struct with
  // default member initializers
  // Alignment const alignment = {
  //   output_handle,
  //   {seq1, seq1len, seq1off, seq1name},
  //   {seq2, seq2len, seq2off, seq2name},
  //   numwidth,
  //   namewidth,
  //   alignwidth,
  //   strand
  // };

  Position position;
  position.query = alignment.strand != 0 ? alignment.query.length - 1 - alignment.query.offset : alignment.query.offset;
  position.target = alignment.target.offset;
  position.query_start = position.query;
  position.target_start = position.target;

  q_line.resize(alignment.width + 1);
  a_line.resize(alignment.width + 1);
  d_line.resize(alignment.width + 1);

  // cigar string can be trimmed (left and right): cigarlen maybe != std::strlen(cigar)
  auto const cigar_pairs = parse_cigar_string_char(Span<char>{cigar, static_cast<size_t>(cigarlen)});
  for (auto const & a_pair: cigar_pairs) {
    auto const operation = a_pair.first;
    auto const runlength = a_pair.second;
    putop(alignment, position, operation, runlength);
  }

  putop(alignment, position, '\0', 1);

  q_line.clear();
  a_line.clear();
  d_line.clear();
}


auto align_getrow(char const * seq, char const * cigar, int const alignlen, bool const is_target) -> std::vector<char> {
  std::vector<char> row(alignlen + 1);
  auto cursor = size_t{0};
  auto const cigar_pairs = parse_cigar_string(Span<char>{cigar, std::strlen(cigar)});

  for (auto const & a_pair: cigar_pairs) {
    auto const operation = a_pair.first;
    auto const runlength = a_pair.second;
    assert(static_cast<size_t>(runlength) < row.size() - cursor);
    auto const is_query = not is_target;
    auto const is_not_a_gap = (operation == Operation::match) // a match, all good
      or ((operation == Operation::deletion) and is_query)    // seq = query, insertion in seq
      or ((operation == Operation::insertion) and is_target); // seq = target, insertion in seq
    if (is_not_a_gap) {
      std::copy(&seq[cursor], &seq[cursor + runlength], &row[cursor]);
    } else {
      /* deletion in sequence: insert gap symbols */
      std::fill_n(&row[cursor], runlength, '-');
    }
    cursor += runlength;
  }

  assert(row[cursor] == '\0');
  return row;
}
