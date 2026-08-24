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
#include "core/db.hpp"
#include "core/fasta.hpp"
#include "utils/median.hpp"
#include "utils/open_file.hpp"
#include "utils/progress.hpp"
#include "utils/view.hpp"
#include <algorithm>  // std::min, std::sort
#include <cassert>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <vector>

#ifndef NDEBUG
#include <limits>
#endif


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct sortinfo_size_s
  {
    /* abundance, as wide as the ;size= annotation it comes from
       (Database::getabundance() returns uint64_t). A 32-bit field truncated
       any annotation above 4294967295 silently, which sorted the most
       abundant sequences to the bottom and skewed the reported median. */
    uint64_t size = 0;
    unsigned int seqno = 0;
  };


  auto create_deck(Database const & db, struct Parameters const & parameters) -> std::vector<struct sortinfo_size_s> {
    auto const dbsequencecount = db.getsequencecount();
    assert(dbsequencecount < std::numeric_limits<std::size_t>::max());
    std::vector<struct sortinfo_size_s> deck(dbsequencecount);
    Progress progress("Getting sizes", deck.size(), parameters);
    auto counter = std::size_t{0};
    for (auto seqno = 0U; seqno < dbsequencecount; ++seqno) {
      auto const abundance = db.getabundance(seqno);
      /* compared as int64_t, the type of the two bounds */
      auto const size = static_cast<int64_t>(abundance);
      if ((size < parameters.opt_minsize) or (size > parameters.opt_maxsize)) {
        continue;
      }
      deck[counter].seqno = seqno;
      deck[counter].size = abundance;
      progress.update(seqno);
      ++counter;
    }
    deck.resize(counter);
    return deck;
  }


  auto sort_deck(std::vector<sortinfo_size_s> & deck, Database const & db, struct Parameters const & parameters) -> void {
    auto compare_sequences = [&db](struct sortinfo_size_s const & lhs,
                                struct sortinfo_size_s const & rhs) noexcept -> bool {
      // highest abundance first...
      if (lhs.size < rhs.size) {
        return false;
      }
      if (lhs.size > rhs.size) {
        return true;
      }
      // ...then ties are sorted by sequence labels (alpha-numerical ordering),
      // preserve input order
      auto const order = db.header_view(lhs.seqno).compare(db.header_view(rhs.seqno));
      if (order != 0) {
        return order < 0;
      }
      // seqno is the input order, and makes the ordering total, so std::sort
      // returns exactly what std::stable_sort returned without the label
      // tie-break: no merge buffer, and fewer comparison passes
      return lhs.seqno < rhs.seqno;
    };

    static constexpr auto one_hundred_percent = 100ULL;
    Progress const progress("Sorting", one_hundred_percent, parameters);
    std::sort(deck.begin(), deck.end(), compare_sequences);
  }


  auto output_median_abundance(std::vector<sortinfo_size_s> const & deck,
                               struct Parameters const & parameters) -> void {
    // Banker's rounding (round half to even)
    auto const median = median_of_descending(
        make_view(deck),
        [](sortinfo_size_s const & entry) { return entry.size; });
    if (not parameters.opt_quiet) {
      static_cast<void>(std::fprintf(stderr, "Median abundance: %.0f\n", median));
    }
    if (parameters.fp_log != nullptr) {
      static_cast<void>(std::fprintf(parameters.fp_log, "Median abundance: %.0f\n", median));
    }
  }


  // auto trim_deck(std::vector<struct sortinfo_size_s> & deck)
  //     -> std::vector<struct sortinfo_size_s> {
  //   // assume deck is sorted by decreasing abundance
  //   // - opt_minsize = 0 by default
  //   // - opt_maxsize = LONG_MAX by default
  //   // - size is unsigned int
  //   auto begin = std::upper_bound(deck.begin(), deck.end(), opt_maxsize,
  //                                 [](int64_t maxsize, struct sortinfo_size_s & seq) -> bool {
  //                                   return seq.size > maxsize;
  //                                 });
  //   auto end = std::lower_bound(deck.begin(), deck.end(), opt_minsize,
  //                               [](int64_t minsize, struct sortinfo_size_s & seq) -> bool {
  //                                 return seq.size <= minsize;
  //                               });
  //   return std::vector<struct sortinfo_size_s>{begin, end};
  // }


  auto truncate_deck(std::vector<struct sortinfo_size_s> & deck,
                     long int const n_first_sequences) -> void {
    if (deck.size() > static_cast<unsigned long>(n_first_sequences)) {
      deck.resize(static_cast<std::size_t>(n_first_sequences));
    }
  }


  // refactoring: extract as a template
  auto output_sorted_fasta(std::vector<struct sortinfo_size_s> const & deck,
                           std::FILE * output_file,
                           Database const & db,
                           struct Parameters const & parameters) -> void {
    Progress progress("Writing output", deck.size(), parameters);
    auto counter = std::size_t{0};
    for (auto const & sequence: deck) {
      fasta_print_db_relabel(output_file, sequence.seqno, counter + 1, db, parameters);
      progress.update(counter);
      ++counter;
    }
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

}  // end of anonymous namespace


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


auto sortbysize(struct Parameters const & parameters) -> void
{
  auto const output_handle = open_mandatory_output_file(parameters.opt_output, OutputOption{"--output"});
  Database db;
  db.read(parameters.opt_sortbysize, 0, parameters);
  // memory-intensive: the entire database is now held in memory

  auto deck = create_deck(db, parameters);

  sort_deck(deck, db, parameters);

  output_median_abundance(deck, parameters);

  truncate_deck(deck, parameters.opt_topn);
  output_sorted_fasta(deck, output_handle.get(), db, parameters);

  db.clear();
}
