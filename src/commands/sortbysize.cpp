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
#include "commands/deck.hpp"
#include "core/db.hpp"
#include "utils/median.hpp"
#include "utils/open_file.hpp"
#include "utils/progress.hpp"
#include "utils/view.hpp"
#include <algorithm>  // std::sort
#include <cassert>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::size_t
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
    uint64_t label_prefix = 0;  // see label_prefix_of() in deck.hpp
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
      deck[counter].label_prefix = label_prefix_of(db.header_view(seqno));
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


  auto output_median_abundance(std::vector<sortinfo_size_s> const & deck,
                               struct Parameters const & parameters) -> void {
    // Banker's rounding (round half to even)
    auto const median = median_of_descending(
        make_view(deck),
        [](sortinfo_size_s const & entry) { return entry.size; });
    report_median(median, "Median abundance", parameters);
  }


}  // end of anonymous namespace


auto sortbysize(struct Parameters const & parameters) -> void
{
  auto const output_handle = open_mandatory_output_file(parameters.opt_output, OutputOption{"--output"});
  Database db;
  db.read(parameters.input_filename, 0, parameters);
  // memory-intensive: the entire database is now held in memory

  auto deck = create_deck(db, parameters);

  sort_deck(deck, db, parameters);

  output_median_abundance(deck, parameters);

  truncate_deck(deck, parameters.opt_topn);
  output_deck_fasta(deck, output_handle.get(), db, parameters);

  db.clear();
}
