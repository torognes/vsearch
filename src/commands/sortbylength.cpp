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
#include "utils/open_file.hpp"
#include "utils/progress.hpp"
#include "utils/view.hpp"
#include <algorithm>  // std::sort, std::min
#include <cassert>
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <vector>

#ifndef NDEBUG
#include <limits>
#endif


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct sortinfo_length_s
  {
    /* abundance, as wide as the ;size= annotation it comes from
       (Database::getabundance() returns uint64_t). A 32-bit field truncated
       any annotation above 4294967295 silently, which broke the documented
       tie-break by decreasing abundance (sortbysize carries the same fix). */
    uint64_t size = 0;
    uint64_t label_prefix = 0;
    unsigned int length = 0;
    unsigned int seqno = 0;
  };


  /* The first eight header bytes packed big-endian into a uint64_t, a shorter
     header zero-padded. Comparing two prefixes as integers matches comparing
     the same bytes lexicographically (as unsigned char, like View's ordering):
     the padding can tie with a further-differing or equally short header, but
     never inverts, so a prefix tie falls back to the full header comparison.
     Cached in the deck because the sort's label tie-break is hot: lengths and
     sizes are heavily tied (amplicons cluster on a few lengths, and size=1
     dominates a typical set), and without the cache each of the O(n log n)
     comparisons chases two random header pointers into the database. */
  auto label_prefix_of(View<char> const header) -> uint64_t {
    static constexpr auto prefix_bytes = std::size_t{8};
    static constexpr auto bits_per_byte = 8U;
    auto const length = std::min(header.size(), prefix_bytes);
    auto prefix = uint64_t{0};
    for (auto position = std::size_t{0}; position < length; ++position) {
      prefix = (prefix << bits_per_byte)
        | static_cast<unsigned char>(header[position]);
    }
    prefix <<= bits_per_byte * (prefix_bytes - length);
    return prefix;
  }


  auto create_deck(Database const & db, struct Parameters const & parameters) -> std::vector<struct sortinfo_length_s> {
    auto const dbsequencecount = db.getsequencecount();
    assert(dbsequencecount < std::numeric_limits<std::size_t>::max());
    std::vector<struct sortinfo_length_s> deck(dbsequencecount);
    Progress progress("Getting lengths", deck.size(), parameters);
    auto counter = std::size_t{0};
    for (auto & sequence: deck) {
      sequence.seqno = static_cast<unsigned int>(counter);
      sequence.length = static_cast<unsigned int>(db.getsequencelen(counter));
      sequence.size = db.getabundance(counter);
      sequence.label_prefix = label_prefix_of(db.header_view(counter));
      progress.update(counter);
      ++counter;
    }
    return deck;
  }


  auto sort_deck(std::vector<sortinfo_length_s> & deck,
                 Database const & db,
                 struct Parameters const & parameters) -> void {
    auto compare_sequences = [&db](struct sortinfo_length_s const & lhs,
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
      if (lhs.label_prefix != rhs.label_prefix) {
        return lhs.label_prefix < rhs.label_prefix;
      }
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


  // refactoring C++17 [[nodiscard]]
  auto find_median_length(std::vector<sortinfo_length_s> const & deck) -> double {
    // function returns a round value or a value with a remainder of 0.5
    static constexpr double half = 0.5;

    if (deck.empty()) {
      return 0.0;
    }

    // refactoring C++11: use const& std::vector.size()
    // plain division on std::size_t: ldiv would have needed a narrowing cast
    // to long (32-bit on the Windows target), and its remainder is recomputed
    // with % just below anyway
    auto const mid = deck.size() / 2;

    // odd number of valid amplicons
    if (deck.size() % 2 != 0)  {
      return deck[mid].length * 1.0;  // a round value
    }

    // even number of valid amplicons
    // (average of two ints is either round or has a remainder of .5)
    // avoid risk of silent overflow for large abundance values:
    // a >= b ; (a + b) / 2 == b + (a - b) / 2
    return deck[mid].length +
      ((deck[mid - 1].length - deck[mid].length) * half);
  }


  auto output_median_length(std::vector<struct sortinfo_length_s> const & deck,
                            struct Parameters const & parameters) -> void {
    // Banker's rounding (round half to even)
    auto const median = find_median_length(deck);
    if (not parameters.opt_quiet)
      {
        std::fprintf(stderr, "Median length: %.0f\n", median);
      }
    if (parameters.fp_log != nullptr)
      {
        std::fprintf(parameters.fp_log, "Median length: %.0f\n", median);
      }
  }


  // refactoring: extract as a template
  auto truncate_deck(std::vector<struct sortinfo_length_s> & deck,
                     long int const n_first_sequences) -> void {
    if (deck.size() > static_cast<unsigned long>(n_first_sequences)) {
      deck.resize(static_cast<std::size_t>(n_first_sequences));
    }
  }


  // refactoring: extract as a template
  auto output_sorted_fasta(std::vector<struct sortinfo_length_s> const & deck,
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

}  // end of anonymous namespace


auto sortbylength(struct Parameters const & parameters) -> void {
  auto const output_handle = open_mandatory_output_file(parameters.opt_output, OutputOption{"--output"});
  Database db;
  db.read(parameters.opt_sortbylength, 0, parameters);
  // memory-intensive: the entire database is now held in memory

  auto deck = create_deck(db, parameters);

  sort_deck(deck, db, parameters);

  output_median_length(deck, parameters);

  truncate_deck(deck, parameters.opt_topn);
  output_sorted_fasta(deck, output_handle.get(), db, parameters);

  db.clear();
}
