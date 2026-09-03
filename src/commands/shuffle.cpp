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
#include "utils/progress.hpp"
#include "utils/open_file.hpp"
#include "utils/random.hpp"
#include "utils/span.hpp"  // make_span
#include <numeric>  // std::iota
#include <random>  // std::mt19937_64
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* the deck holds bare sequence numbers, in database order to start with */
  auto create_deck(Database const & db) -> std::vector<unsigned int> {
    auto const dbsequencecount = db.getsequencecount();
    std::vector<unsigned int> deck(dbsequencecount);
    std::iota(deck.begin(), deck.end(), 0U);
    return deck;
  }


  auto shuffle_deck(std::vector<unsigned int> & deck, struct Parameters const & parameters) -> void {
    static constexpr auto one_hundred_percent = 100ULL;
    Progress const progress("Shuffling", one_hundred_percent, parameters);
    /* the RandomSeed carries the full 64-bit --randseed (or an OS value
       when 0); random_shuffle() is a portable Fisher-Yates so the order is
       identical across platforms for a given seed (see util.h) */
    RandomSeed const seed(parameters);
    std::mt19937_64 uniform_generator(seed.value());
    random_shuffle(make_span(deck), uniform_generator);
  }

}  // end of anonymous namespace


auto shuffle(struct Parameters const & parameters) -> void {
  auto const output_handle = open_mandatory_output_file(parameters.opt_output, OutputOption{"--output"});
  Database db;
  db.read(parameters.input_filename, 0, parameters);
  // memory-intensive: the entire database is now held in memory

  auto deck = create_deck(db);
  shuffle_deck(deck, parameters);

  truncate_deck(deck, parameters.opt_topn);
  output_deck_fasta(deck, output_handle.get(), db, parameters);

  db.clear();
}
