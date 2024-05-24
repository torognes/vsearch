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
#include <algorithm>  // std::min, std::sort
#include <cassert>
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <cstdlib>   // std::ldiv
#include <cstring>  // std::strcmp
#include <vector>


struct sortinfo_size_s
{
  unsigned int size = 0;
  unsigned int seqno = 0;
};


namespace {
  // anonymous namespace to avoid linker error (multiple definitions
  // of function with identical names and parameters)
  auto create_deck() -> std::vector<struct sortinfo_size_s> {
    auto const dbsequencecount = db_getsequencecount();
    assert(dbsequencecount < std::numeric_limits<std::size_t>::max());
    std::vector<struct sortinfo_size_s> deck(dbsequencecount);
    progress_init("Getting sizes", deck.size());
    auto counter = std::size_t{0};
    for (auto seqno = 0U; seqno < dbsequencecount; ++seqno) {
      auto const size = static_cast<int64_t>(db_getabundance(seqno));
      if ((size < opt_minsize) or (size > opt_maxsize)) {
        continue;
      }
      deck[counter].seqno = seqno;
      deck[counter].size = static_cast<unsigned int>(size);
      progress_update(seqno);
      ++counter;
    }
    progress_done();
    deck.resize(counter);
    return deck;
  }
}


auto sort_deck(std::vector<sortinfo_size_s> & deck) -> void {
  auto compare_sequences = [](struct sortinfo_size_s const & lhs,
                              struct sortinfo_size_s const & rhs) -> bool {
    // highest abundance first...
    if (lhs.size < rhs.size) {
      return false;
    }
    if (lhs.size > rhs.size) {
      return true;
    }
    // ...then ties are sorted by sequence labels (alpha-numerical ordering)...
    auto const result = std::strcmp(db_getheader(lhs.seqno), db_getheader(rhs.seqno));
    if (result > 0) {
      return false;
    }
    if (result < 0) {
      return true;
    }
    // ... then ties are sorted by input order (seqno)
    if (lhs.seqno < rhs.seqno) {
      return true;
    }
    return false;
  };

  static constexpr auto one_hundred_percent = 100ULL;
  progress_init("Sorting", one_hundred_percent);
  std::sort(deck.begin(), deck.end(), compare_sequences);
  progress_done();
}


// refactoring C++17 [[nodiscard]]
auto find_median_abundance(std::vector<sortinfo_size_s> const & deck) -> double
{
  // function returns a round value or a value with a remainder of 0.5
  static constexpr double half = 0.5;

  if (deck.empty()) {
    return 0.0;
  }

  // refactoring C++11: use const& std::vector.size()
  auto const midarray = std::ldiv(static_cast<long>(deck.size()), 2L);

  // odd number of valid amplicons
  if (deck.size() % 2 != 0)  {
    return deck[midarray.quot].size * 1.0;  // a round value
  }

  // even number of valid amplicons
  // (average of two ints is either round or has a remainder of .5)
  // avoid risk of silent overflow for large abundance values:
  // a >= b ; (a + b) / 2 == b + (a - b) / 2
  return deck[midarray.quot].size +
    ((deck[midarray.quot - 1].size - deck[midarray.quot].size) * half);
}


auto output_median_abundance(std::vector<sortinfo_size_s> const & deck) -> void {
  // Banker's rounding (round half to even)
  auto const median = find_median_abundance(deck);
  if (not opt_quiet) {
      static_cast<void>(fprintf(stderr, "Median abundance: %.0f\n", median));
  }
  if (opt_log != nullptr) {
    static_cast<void>(fprintf(fp_log, "Median abundance: %.0f\n", median));
  }
}


// refactoring: extract as a template
auto output_sorted_fasta(std::vector<struct sortinfo_size_s> & deck,
                           long int const n_first_sequences,
                           std::FILE * output_file) -> void {
  auto const final_size = std::min(deck.size(),
                                   static_cast<unsigned long>(n_first_sequences));
  deck.resize(final_size);
  progress_init("Writing output", deck.size());
  auto counter = std::size_t{0};
  for (auto const & sequence: deck) {
    fasta_print_db_relabel(output_file, sequence.seqno, counter + 1);
    progress_update(counter);
    ++counter;
  }
  progress_done();
}


// refactoring: trim misize and maxsize with a free function
// https://stackoverflow.com/questions/26719144/how-to-erase-a-value-efficiently-from-a-sorted-vector
// auto erase_high_abundances(std::vector<int> & vec, int value) -> void
// {
//     auto lb = std::lower_bound(std::begin(vec), std::end(vec), value);
//     if (lb != std::end(vec) and *lb == value) {
//         auto ub = std::upper_bound(lb, std::end(vec), value);
//         vec.erase(lb, ub);
//     }
// }


// refactoring:
// - create vector (no branch)
// - stable_sort vector (by increasing size, then label)
// - find lower_bound(comp(opt_minsize)),
// - deck.resize()
// - find upper_bound(comp(opt_maxsize)),
// - std::vector<S> subdeck = {deck.begin() + upper_bound, deck.end()};  // view?
// - opt_minsize = 0 by default
// - opt_maxsize = LONG_MAX by default
// - top_n = LONG_MAX by default
// - mediane, etc...
// - std::min(subdeck.size(), topn);


auto sortbysize() -> void
{
  if (opt_output == nullptr) {
    fatal("FASTA output file for sortbysize must be specified with --output");
  }

  auto * fp_output = fopen_output(opt_output);
  if (fp_output == nullptr) {
    fatal("Unable to open sortbysize output file for writing");
  }

  db_read(opt_sortbysize, 0);
  show_rusage();

  auto deck = create_deck();
  show_rusage();

  sort_deck(deck);

  output_median_abundance(deck);
  show_rusage();

  output_sorted_fasta(deck, opt_topn, fp_output);
  show_rusage();  // refactoring: why three calls to show_rusage()?

  db_free();

  if (fp_output != nullptr) {
    static_cast<void>(fclose(fp_output));
  }
}
