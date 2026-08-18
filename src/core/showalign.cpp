/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "vsearch.hpp"
#include "utils/cigar.hpp"
#include "utils/maps.hpp"
#include "utils/print_record.hpp"  // OutputRecord, fprint
#include "utils/print_view.hpp"  // fprint
#include "utils/view.hpp"
#include <algorithm>  // std::copy, std::fill_n, std::max, std::min
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE
#include <cstring>  // std::strlen
#include <iterator>  // std::next
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // The three rows of a printed alignment block, each sized to the alignment
  // width (+1 for the NUL terminator). Owned by align_show() and threaded
  // through putop()/putop_final()/print_alignment_block() (formerly the
  // file-scope q_line/a_line/d_line buffers).
  struct AlignmentRows {
    std::vector<char> query;    // query
    std::vector<char> symbols;  // alignment symbols (|)
    std::vector<char> target;   // target
  };


  enum struct Viewpoint : char {
    target,
    query,
  };


  struct Position {
    int64_t line = 0;
    int64_t query = 0;
    int64_t target = 0;
    int64_t query_start = 0;
    int64_t target_start = 0;
  };


  struct Sequence {
    View<char> sequence;
    int64_t offset = 0;
    // The row label, as a view rather than a pointer, because the label is
    // written into a field of alignment.headwidth characters and the padding
    // needs its length. align_show() measures the caller's char const * once
    // per alignment, not once per printed block.
    View<char> name;
  };


  struct Alignment {
    static constexpr auto poswidth_default = 3;
    static constexpr auto headwidth_default = 5;
    std::FILE * output_handle = nullptr;
    Sequence query;
    Sequence target;
    int64_t width = 0;
    int poswidth = poswidth_default;
    int headwidth = headwidth_default;
    bool is_reverse_strand = false;
    bool n_mismatch = false;  // treat alignment against N as a mismatch (opt_n_mismatch)
  };


  auto get_aligment_symbol(char const query_nuc, char const target_nuc, bool const n_mismatch) -> char {
    static constexpr auto is_N = 15U;
    auto const query_coded = map_4bit(query_nuc);
    auto const target_coded = map_4bit(target_nuc);

    if (n_mismatch and ((query_coded == is_N) or (target_coded == is_N))) {
      return ' ';  // N are mismatches
    }
    if ((query_coded == target_coded) and not is_ambiguous_4bit(query_coded)) {
      return '|';  // a perfect match
    }
    if ((query_coded & target_coded) != 0U) {
     return '+';  // an equivalence (ambiguous nucleotides)
    }
    return ' ';
  }


  constexpr
  auto adapt_to_viewpoint(Viewpoint const viewpoint) noexcept -> Operation {
    // cigar operations are relative to target (see issue #618)
    //  - D is a deletion in target, an insertion in query
    //  - I is an insertion in target, a deletion in query
    //
    // Which operation is an insertion depends on viewpoint:
    return (viewpoint == Viewpoint::target) ? Operation::insertion : Operation::deletion;
  }


  // 'insertion_equivalent' is a template parameter rather than a function
  // one: the two wrappers below are its only callers, each passing the enum
  // constant adapt_to_viewpoint() already folds at compile time, so the
  // runtime parameter only ever carried a compile-time fact. GCC 13 at -O2
  // did not exploit that on its own -- it emitted one generic body that
  // compared against the spilled parameter on every cigar pair
  // (cmp 0xf(%rsp),%dl) -- whereas each specialisation compares against
  // immediates, and the target viewpoint collapses to one AND plus one
  // compare because 'M' and 'I' differ only in bit 2. The trade is
  // duplication: two 681-byte bodies inlined into the wrappers where one
  // 705-byte body and two 67-byte wrappers stood before, +672 bytes of
  // .text in this TU. Cheap, because the cigar walk runs once per printed
  // alignment (--alnout), not per byte.
  template <Operation insertion_equivalent>
  auto get_alignment_row(View<char> const seq_view, View<char> const cigar_view,
                         int const alignlen) -> std::vector<char> {
    std::vector<char> row(static_cast<size_t>(alignlen) + 1);
    auto cursor_src = size_t{0};
    auto cursor_dest = size_t{0};

    for (auto const & a_pair: parse_cigar_string(cigar_view)) {
      auto const operation = a_pair.first;
      auto const runlength = a_pair.second;
      assert(static_cast<size_t>(runlength) < row.size() - cursor_dest);
      if ((operation == Operation::match) or
          (operation == insertion_equivalent)) {
        auto const subsequence = seq_view.subspan(cursor_src, static_cast<size_t>(runlength));
        std::copy(subsequence.cbegin(), subsequence.cend(), &row[cursor_dest]);
        cursor_src += static_cast<size_t>(runlength);
      } else {
        // viewpoint_deletion = fill-in with gap symbols
        std::fill_n(&row[cursor_dest], runlength, '-');
      }
      cursor_dest += static_cast<size_t>(runlength);
    }

    assert(row[cursor_dest] == '\0');
    return row;
  }


  auto get_query_nucleotide(Alignment const & alignment, Position const & position) -> char {
    auto const nucleotide = alignment.query.sequence[static_cast<std::size_t>(position.query)];
    if (alignment.is_reverse_strand) {
      return map_complement(nucleotide);
    }
    return nucleotide;
  }


  auto get_target_nucleotide(Alignment const & alignment, Position const & position) -> char {
    return alignment.target.sequence[static_cast<std::size_t>(position.target)];
  }


  // What a "%*s" conversion did: pad on the left to 'width' characters, and
  // never truncate a label wider than the field.
  auto print_padded(OutputRecord & record, View<char> const text,
                    int const width) -> void {
    auto const field = static_cast<std::size_t>(std::max(width, 0));
    if (text.size() < field) { fprint_spaces(record, field - text.size()); }
    fprint(record, text);
  }


  auto print_alignment_block(Alignment const & alignment, Position const & position,
                             AlignmentRows const & rows) -> void {
    // current query and target starting and ending positions
    auto const query_start = std::min(position.query_start + 1, static_cast<int64_t>(alignment.query.sequence.size()));
    auto const query_end = alignment.is_reverse_strand ? position.query + 2 : position.query;
    auto const target_start = std::min(position.target_start + 1, static_cast<int64_t>(alignment.target.sequence.size()));
    auto const target_end = position.target;

    // The three rows are filled up to position.line and NUL-terminated there
    // by the caller; that length is what the "%s" conversions used to walk.
    // The three rows are filled up to position.line and NUL-terminated there
    // by the caller; that length is what the "%s" conversions used to walk.
    auto const row_length = static_cast<std::size_t>(position.line);

    // One record for the whole block. A block is three rows of alignwidth
    // characters, so it does not fit in the buffer and flushes a few times --
    // still far fewer writes than the ~20 the block needs field by field, and
    // the padding and the positions are what dominate the call count.
    OutputRecord record {alignment.output_handle};

    fprint(record, '\n');
    print_padded(record, alignment.query.name, alignment.headwidth);
    fprint(record, ' ');
    fprint_integer(record, query_start, static_cast<std::size_t>(alignment.poswidth));
    fprint(record, ' ');
    fprint(record, alignment.is_reverse_strand ? '-' : '+');
    fprint(record, ' ');
    fprint(record, make_view(rows.query).first(row_length));
    fprint(record, ' ');
    fprint_integer(record, query_end);
    fprint(record, '\n');

    fprint_spaces(record, static_cast<std::size_t>(std::max(alignment.headwidth, 0)));
    fprint(record, ' ');
    fprint_spaces(record, static_cast<std::size_t>(std::max(alignment.poswidth, 0)));
    fprint(record, "   ");
    fprint(record, make_view(rows.symbols).first(row_length));
    fprint(record, '\n');

    print_padded(record, alignment.target.name, alignment.headwidth);
    fprint(record, ' ');
    fprint_integer(record, target_start, static_cast<std::size_t>(alignment.poswidth));
    fprint(record, " + ");
    fprint(record, make_view(rows.target).first(row_length));
    fprint(record, ' ');
    fprint_integer(record, target_end);
    fprint(record, '\n');
  }


  inline auto putop(Alignment const & alignment, Position & position, AlignmentRows & rows,
                    Operation const operation, int64_t const runlength) -> void {
    int64_t const delta = alignment.is_reverse_strand ? -1 : +1;

    for (auto count = runlength; count != 0; --count) {

      if (position.line == 0) {
        position.query_start = position.query;
        position.target_start = position.target;
      }

      auto const query_nuc = get_query_nucleotide(alignment, position);
      auto const target_nuc = get_target_nucleotide(alignment, position);

      auto const line_index = static_cast<size_t>(position.line);

      switch (operation) {
      case Operation::match:
        position.query += delta;
        position.target += 1;
        rows.query[line_index] = query_nuc;
        rows.symbols[line_index] = get_aligment_symbol(query_nuc, target_nuc, alignment.n_mismatch);
        rows.target[line_index] = target_nuc;
        ++position.line;
        break;

      case Operation::deletion:  // gap in target (insertion in query)
        position.query += delta;
        rows.query[line_index] = query_nuc;
        rows.symbols[line_index] = ' ';
        rows.target[line_index] = '-';
        ++position.line;
        break;

      case Operation::insertion:  // insertion in target (gap in query)
        position.target += 1;
        rows.query[line_index] = '-';
        rows.symbols[line_index] = ' ';
        rows.target[line_index] = target_nuc;
        ++position.line;
        break;
      }

      if (position.line == alignment.width) {
        // maximal alignment width is reached, print alignment block
        auto const terminator_index = static_cast<size_t>(position.line);
        rows.query[terminator_index] = '\0';
        rows.symbols[terminator_index] = '\0';
        rows.target[terminator_index] = '\0';
        print_alignment_block(alignment, position, rows);
        position.line = 0;  // needed to avoid out-of-bounds
      }
    }
  }


  auto putop_final(Alignment const & alignment, Position const & position,
                   AlignmentRows & rows) -> void {
    if (position.line == 0) { return; }  // final block already printed
    auto const terminator_index = static_cast<size_t>(position.line);
    rows.query[terminator_index] = '\0';
    rows.symbols[terminator_index] = '\0';
    rows.target[terminator_index] = '\0';
    print_alignment_block(alignment, position, rows);
  }

}  // end of anonymous namespace


auto align_show(std::FILE * output_handle,
                View<char> const seq1,
                int64_t const seq1off,
                char const * seq1name,
                View<char> const seq2,
                int64_t const seq2off,
                char const * seq2name,
                View<char> const cigar,
                int const numwidth,
                int const namewidth,
                int64_t const alignwidth,
                int const strand,
                struct Parameters const & parameters) -> void
{

  Alignment alignment;
  alignment.output_handle = output_handle;
  alignment.query.sequence = seq1;
  alignment.query.offset = seq1off;
  alignment.query.name = View<char>{seq1name, std::strlen(seq1name)};
  alignment.target.sequence = seq2;
  alignment.target.offset = seq2off;
  alignment.target.name = View<char>{seq2name, std::strlen(seq2name)};
  alignment.width = alignwidth;
  alignment.poswidth = numwidth;
  alignment.headwidth = namewidth;
  alignment.is_reverse_strand = strand != 0;
  alignment.n_mismatch = parameters.opt_n_mismatch;

  // C++14 refactoring: aggregate initialization of a struct with
  // default member initializers
  // Alignment const alignment = {
  //   output_handle,
  //   {seq1, seq1off, seq1name},
  //   {seq2, seq2off, seq2name},
  //   numwidth,
  //   namewidth,
  //   alignwidth,
  //   strand
  // };

  Position position;
  position.query = alignment.is_reverse_strand
    ? static_cast<int64_t>(alignment.query.sequence.size()) - 1 - alignment.query.offset
    : alignment.query.offset;
  position.target = alignment.target.offset;
  position.query_start = position.query;
  position.target_start = position.target;

  AlignmentRows rows;
  rows.query.resize(static_cast<size_t>(alignment.width) + 1);
  rows.symbols.resize(static_cast<size_t>(alignment.width) + 1);
  rows.target.resize(static_cast<size_t>(alignment.width) + 1);

  // cigar string can be trimmed (left and right): cigar.size() maybe != std::strlen(cigar.data())
  auto const cigar_pairs = parse_cigar_string(cigar);
  for (auto const & a_pair: cigar_pairs) {
    auto const operation = a_pair.first;
    auto const runlength = a_pair.second;
    putop(alignment, position, rows, operation, runlength);
  }

  putop_final(alignment, position, rows);
}


auto get_alignment_qrow(View<char> const seq_view, View<char> const cigar_view,
                        int const alignlen) -> std::vector<char> {
  return get_alignment_row<adapt_to_viewpoint(Viewpoint::query)>(seq_view, cigar_view, alignlen);
}


auto get_alignment_trow(View<char> const seq_view, View<char> const cigar_view,
                        int const alignlen) -> std::vector<char> {
  return get_alignment_row<adapt_to_viewpoint(Viewpoint::target)>(seq_view, cigar_view, alignlen);
}
