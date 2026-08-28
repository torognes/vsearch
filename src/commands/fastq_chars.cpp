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
#include "core/fastq.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/maps.hpp"
#include "utils/quality_encoding.hpp"  // sanger_ascii_offset, solexa_ascii_offset
#include "utils/view.hpp"
#include <algorithm>  // std::find_if
#include <cassert>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <iterator>  // std::distance
#include <vector>


#ifndef NDEBUG
#include <limits>
constexpr long int char_max = std::numeric_limits<char>::max();
#endif

constexpr unsigned int n_characters = 256;

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct statistics {
    std::vector<uint64_t> sequence_chars;
    std::vector<uint64_t> quality_chars;
    std::vector<uint64_t> tail_chars;
    std::vector<int> maxrun;
    uint64_t total_chars = 0;
    uint64_t seq_count = 0;
    unsigned char qmin_n = 255;
    unsigned char qmax_n = 0;
    char qmin = '\0';
    char qmax = '\0';
    char fastq_ascii = '\0';
    char fastq_qmin = '\0';
    char fastq_qmax = '\0';
    FastqEncoding encoding = FastqEncoding::sanger;  // set by guess_quality_offset()
  };


  auto guess_quality_offset(struct statistics & stats) -> void {
    stats.encoding = classify_encoding(stats.qmin, stats.qmax);
    stats.fastq_ascii = static_cast<char>(offset_of(stats.encoding));
    stats.fastq_qmax = static_cast<char>(stats.qmax - stats.fastq_ascii);
    stats.fastq_qmin = static_cast<char>(stats.qmin - stats.fastq_ascii);
  }


  auto find_lowest_quality_symbol(struct statistics & stats) -> void {
    auto lowest = std::find_if(stats.quality_chars.cbegin(),
                               stats.quality_chars.cend(),
                               [](uint64_t const counter) -> bool {
                                 return counter != 0;
                               });
    if (lowest == stats.quality_chars.cend()) {
      return;
    }
    auto const index = std::distance(stats.quality_chars.cbegin(), lowest);
    assert(index >= 0);
    assert(index <= char_max);
    stats.qmin = static_cast<char>(index);
  }


  auto find_highest_quality_symbol(struct statistics & stats) -> void {
    // note: searching using reverse iterators
    auto highest = std::find_if(stats.quality_chars.rbegin(),
                                stats.quality_chars.rend(),
                                [](uint64_t const counter) -> bool {
                                  return counter != 0; }
                                );
    if (highest == stats.quality_chars.rend()) {
      return;
    }
    auto const index = std::distance(highest, stats.quality_chars.rend()) - 1;
    assert(index >= 0);
    assert(index <= char_max);
    stats.qmax = static_cast<char>(index);
  }


  auto search_trailing_homopolymers(View<char> const symbols, int64_t const tail_length_signed) -> char {
    // search for trailing homopolymers of length >= 'tail_length'
    assert(tail_length_signed >= 0);
    auto const tail_length = static_cast<std::size_t>(tail_length_signed);
    if (symbols.size() < tail_length) {
      return '\0';
    }
    auto const last_symbol = symbols.back();
    auto const tail = symbols.last(tail_length);
    if (std::all_of(
            tail.begin(), tail.end(),
            [last_symbol](char const symbol) -> bool { return symbol == last_symbol; })
        ) {
      return last_symbol;
    }
    return '\0';
  }


  auto stats_message(std::FILE * output_stream,
                     struct statistics const & stats) -> void {
    assert(stats.sequence_chars['n'] == 0);  // sequences are uppercased, no results for lowercase symbols
    fprint(output_stream, "Read ");
    fprint_integer(output_stream, stats.seq_count);
    fprint(output_stream, " sequences.\n");

    if (stats.seq_count == 0) {
      return;
    }

    fprint(output_stream, "Qmin ");
    fprint_integer(output_stream, static_cast<int>(stats.qmin));
    fprint(output_stream, ", Qmax ");
    fprint_integer(output_stream, static_cast<int>(stats.qmax));
    fprint(output_stream, ", Range ");
    fprint_integer(output_stream, stats.qmax - stats.qmin + 1);
    fprint(output_stream, '\n');

    fprint(output_stream, "Guess: -fastq_qmin ");
    fprint_integer(output_stream, static_cast<int>(stats.fastq_qmin));
    fprint(output_stream, " -fastq_qmax ");
    fprint_integer(output_stream, static_cast<int>(stats.fastq_qmax));
    fprint(output_stream, " -fastq_ascii ");
    fprint_integer(output_stream, static_cast<int>(stats.fastq_ascii));
    fprint(output_stream, '\n');

    // exhaustive switch, deliberately without a default case: -Wswitch then
    // flags a newly added encoding that forgot its label here
    switch (stats.encoding)
      {
      case FastqEncoding::sanger:
        // Sanger Phred+33, quality values ranging from 0 to 40 (ascii: 33 to 73)
        fprint(output_stream, "Guess: Original Sanger format (phred+33)\n");
        break;
      case FastqEncoding::illumina_1_8:
        fprint(output_stream, "Guess: Illumina 1.8+ format (phred+33)\n");
        break;
      case FastqEncoding::solexa:
        /* solexa+64, not phred+64: the offset is shared with Illumina 1.3+,
           the score definition is not. A Solexa score is -10 log10(p / (1 - p))
           and vsearch has only the Phred formula, so the options suggested
           above will read the file but overstate the error probability at low
           scores (see the fastq(5) manual page) */
        fprint(output_stream, "Guess: Solexa format (solexa+64, not supported)\n");
        break;
      case FastqEncoding::illumina_1_3:
        fprint(output_stream, "Guess: Illumina 1.3+ format (phred+64)\n");
        break;
      case FastqEncoding::illumina_1_5:
        // Illumina 1.5+ Phred+64, quality values ranging from 3 to 41 (ascii: 67 to 105)
        // Q2 (ascii 66, 'B') is the Read Segment Quality Control Indicator
        fprint(output_stream, "Guess: Illumina 1.5+ format (phred+64)\n");
        break;
      }

    fprint(output_stream, '\n');
    fprint(output_stream, "Letter          N   Freq MaxRun\n");
    fprint(output_stream, "------ ---------- ------ ------\n");

    double const percentage_factor = 100.0 / static_cast<double>(stats.total_chars);
    unsigned char index = 0;
    for (auto const counter: stats.sequence_chars)
      {
        if (counter == 0) { ++index ; continue; }
        fprint(output_stream, "     ");
        fprint(output_stream, static_cast<char>(index));
        fprint(output_stream, ' ');
        fprint_integer(output_stream, counter, 10);
        fprint(output_stream, ' ');
        std::fprintf(output_stream, "%5.1f", static_cast<double>(counter) * percentage_factor);
        fprint(output_stream, "% ");
        fprint_integer(output_stream, stats.maxrun[index], 6);
        if (index == 'N')
          {
            if (stats.qmin_n < stats.qmax_n)
              {
                fprint(output_stream, "  Q=");
                fprint(output_stream, static_cast<char>(stats.qmin_n));
                fprint(output_stream, "..");
                fprint(output_stream, static_cast<char>(stats.qmax_n));
              }
            else
              {
                fprint(output_stream, "  Q=");
                fprint(output_stream, static_cast<char>(stats.qmin_n));
              }
          }
        fprint(output_stream, '\n');
        ++index;
      }

    fprint(output_stream, '\n');
    fprint(output_stream, "Char  ASCII    Freq       Tails\n");
    fprint(output_stream, "----  -----  ------  ----------\n");

    for (char i = stats.qmin; i <= stats.qmax; ++i)
      {
        auto const quality_index = static_cast<unsigned char>(i);
        if (stats.quality_chars[quality_index] == 0) {
          continue;
        }
        fprint(output_stream, " '");
        fprint(output_stream, i);
        fprint(output_stream, "'  ");
        fprint_integer(output_stream, static_cast<int>(i), 5);
        fprint(output_stream, "  ");
        std::fprintf(output_stream, "%5.1f", static_cast<double>(stats.quality_chars[quality_index]) * percentage_factor);
        fprint(output_stream, "%  ");
        fprint_integer(output_stream, stats.tail_chars[quality_index], 10);
        fprint(output_stream, '\n');
      }
  }
}  // end of anonymous namespace


auto fastq_chars(struct Parameters const & parameters) -> void
{
  struct statistics stats;
  stats.sequence_chars.resize(n_characters);
  stats.quality_chars.resize(n_characters);
  stats.tail_chars.resize(n_characters);
  stats.maxrun.resize(n_characters);

  auto fastq_handle = fastq_open(parameters.opt_fastq_chars, parameters);
  /* this command *is* the encoding diagnostic and prints its own guess below,
     so the reader must not say the same thing first */
  fastq_handle->silence_offset_warning();

  auto const filesize = fastq_handle->get_size();

  {
    Progress progress("Reading FASTQ file", filesize, parameters);

    while (fastq_handle->next(false, chrmap_upcase()))
      {
        auto const seq_length = fastq_handle->get_sequence_length();
        auto const * seq_ptr = fastq_handle->get_sequence();
        auto const * qual_ptr = fastq_handle->get_quality();

        ++stats.seq_count;
        stats.total_chars += seq_length;

        auto run_char = -1;
        auto run = 0;

        for (auto i = 0ULL ; i < seq_length ; ++i)
          {
            auto const seq_symbol = static_cast<unsigned char>(*seq_ptr);
            std::advance(seq_ptr, 1);
            auto const qual_symbol = static_cast<unsigned char>(*qual_ptr);
            std::advance(qual_ptr, 1);
            ++stats.sequence_chars[seq_symbol];
            ++stats.quality_chars[qual_symbol];

            if (seq_symbol == 'N')
              {
                stats.qmin_n = std::min(qual_symbol, stats.qmin_n);
                stats.qmax_n = std::max(qual_symbol, stats.qmax_n);
              }

            if (seq_symbol == run_char)
              {
                ++run;
              }
            else
              {
                // a run just ended: its length was maximal, record it once
                if (run != 0)
                  {
                    auto const run_index = static_cast<unsigned char>(run_char);
                    stats.maxrun[run_index] = std::max(run, stats.maxrun[run_index]);
                  }
                run_char = seq_symbol;
                run = 0;
              }
          }

        if (run != 0)
          {
            auto const run_index = static_cast<unsigned char>(run_char);
            stats.maxrun[run_index] = std::max(run, stats.maxrun[run_index]);
          }

        // search for trailing homopolymers in quality strings
        auto const tail_char =
          search_trailing_homopolymers(View<char>{fastq_handle->get_quality(), seq_length},
                                       parameters.opt_fastq_tail);
        if (tail_char != '\0') {
          ++stats.tail_chars[static_cast<unsigned char>(tail_char)];
        }

        progress.update(fastq_handle->get_position());
      }
  }

  fastq_handle->report_stripped_warning(parameters);

  find_lowest_quality_symbol(stats);
  find_highest_quality_symbol(stats);
  guess_quality_offset(stats);

  if (not parameters.opt_quiet) {
    stats_message(stderr, stats);
  }
  if (parameters.fp_log != nullptr) {
    stats_message(parameters.fp_log, stats);
  }
}
