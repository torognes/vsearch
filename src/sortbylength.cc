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
#include <algorithm>  // std::sort, std::min
#include <cassert>
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <cstdlib>  // std::ldiv
#include <cstring>  // std::strcmp
#include <vector>

#ifndef NDEBUG
#include <limits>
#endif


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct sortinfo_length_s
  {
    unsigned int length = 0;
    unsigned int size = 0;
    unsigned int seqno = 0;
  };


  auto open_output_file(struct Parameters const & parameters) -> std::FILE * {
    if (parameters.opt_output == nullptr) {
      fatal("FASTA output file for sortbylength must be specified with --output");
    }
    auto * output_handle = fopen_output(parameters.opt_output);
    if (output_handle == nullptr) {
      fatal("Unable to open sortbylength output file for writing");
    }
    return output_handle;
  }


  auto create_deck() -> std::vector<struct sortinfo_length_s> {
    auto const dbsequencecount = db_getsequencecount();
    assert(dbsequencecount < std::numeric_limits<std::size_t>::max());
    std::vector<struct sortinfo_length_s> deck(dbsequencecount);
    progress_init("Getting lengths", deck.size());
    auto counter = std::size_t{0};
    for (auto & sequence: deck) {
      sequence.seqno = counter;
      sequence.length = db_getsequencelen(counter);
      sequence.size = db_getabundance(counter);
      progress_update(counter);
      ++counter;
    }
    progress_done();
    return deck;
  }


  auto sort_deck(std::vector<sortinfo_length_s> & deck) -> void {
    auto compare_sequences = [](struct sortinfo_length_s const & lhs,
                                struct sortinfo_length_s const & rhs) -> bool {
      // longest first...
      if (lhs.length < rhs.length) {
        return false;
      }
      if (lhs.length > rhs.length) {
        return true;
      }
      // ... then ties are sorted by decreasing abundance values...
      if (lhs.size < rhs.size) {
        return false;
      }
      if (lhs.size > rhs.size) {
        return true;
      }
      // ...then ties are sorted by sequence labels (alpha-numerical ordering),
      // preserve input order
      auto const result = std::strcmp(db_getheader(lhs.seqno), db_getheader(rhs.seqno));
      return result < 0;
    };

    static constexpr auto one_hundred_percent = 100ULL;
    progress_init("Sorting", one_hundred_percent);
    std::stable_sort(deck.begin(), deck.end(), compare_sequences);
    progress_done();
  }


  // refactoring C++17 [[nodiscard]]
  auto find_median_length(std::vector<sortinfo_length_s> const &deck) -> double {
    // function returns a round value or a value with a remainder of 0.5
    static constexpr double half = 0.5;

    if (deck.empty()) {
      return 0.0;
    }

    // refactoring C++11: use const& std::vector.size()
    auto const midarray = std::ldiv(static_cast<long>(deck.size()), 2L);

    // odd number of valid amplicons
    if (deck.size() % 2 != 0)  {
      return deck[midarray.quot].length * 1.0;  // a round value
    }

    // even number of valid amplicons
    // (average of two ints is either round or has a remainder of .5)
    // avoid risk of silent overflow for large abundance values:
    // a >= b ; (a + b) / 2 == b + (a - b) / 2
    return deck[midarray.quot].length +
      ((deck[midarray.quot - 1].length - deck[midarray.quot].length) * half);
  }


  auto output_median_length(std::vector<struct sortinfo_length_s> const & deck,
                            struct Parameters const & parameters) -> void {
    // Banker's rounding (round half to even)
    auto const median = find_median_length(deck);
    if (not parameters.opt_quiet)
      {
        std::fprintf(stderr, "Median length: %.0f\n", median);
      }
    if (parameters.opt_log != nullptr)
      {
        std::fprintf(fp_log, "Median length: %.0f\n", median);
      }
  }


  // refactoring: extract as a template
  auto truncate_deck(std::vector<struct sortinfo_length_s> &deck,
                     long int const n_first_sequences) -> void {
    if (deck.size() > static_cast<unsigned long>(n_first_sequences)) {
      deck.resize(n_first_sequences);
    }
  }


  // refactoring: extract as a template
  auto output_sorted_fasta(std::vector<struct sortinfo_length_s> const & deck,
                           std::FILE * output_file) -> void {
    progress_init("Writing output", deck.size());
    auto counter = std::size_t{0};
    for (auto const & sequence: deck) {
      fasta_print_db_relabel(output_file, sequence.seqno, counter + 1);
      progress_update(counter);
      ++counter;
    }
    progress_done();
  }

}  // end of anonymous namespace


auto sortbylength(struct Parameters const & parameters) -> void {
  auto * output_handle = open_output_file(parameters);
  db_read(parameters.opt_sortbylength, 0);
  show_rusage();

  auto deck = create_deck();
  show_rusage();

  sort_deck(deck);

  output_median_length(deck, parameters);
  show_rusage();

  truncate_deck(deck, parameters.opt_topn);
  output_sorted_fasta(deck, output_handle);
  show_rusage();

  db_free();
  if (output_handle != nullptr) {
    static_cast<void>(std::fclose(output_handle));
  }
}
