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

#pragma once

/* Helpers shared by the sortby* commands (sortbysize.cpp, sortbylength.cpp):
   each builds a small per-record "deck", sorts it, and writes the database
   records back out in deck order. The deck record types differ (one carries a
   length, the other does not), so the helpers that walk a deck are templates
   over the record type; each command keeps its own create_deck and sort_deck. */

#include "vsearch.hpp"  // struct Parameters
#include "core/fasta.hpp"  // fasta_print_db_relabel
#include "utils/progress.hpp"
#include "utils/view.hpp"
#include <algorithm>  // std::min
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <vector>

struct Database;


/* The first eight header bytes packed big-endian into a uint64_t, a shorter
   header zero-padded. Comparing two prefixes as integers matches comparing
   the same bytes lexicographically (as unsigned char, like View's ordering):
   the padding can tie with a further-differing or equally short header, but
   never inverts, so a prefix tie falls back to the full header comparison.
   Cached in the deck because the sort's label tie-break is hot: the primary
   keys are heavily tied (size=1 dominates a typical amplicon set, and lengths
   cluster on a few values), and without the cache each of the O(n log n)
   comparisons chases two random header pointers into the database. */
inline auto label_prefix_of(View<char> const header) -> uint64_t {
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


/* Report the deck's median (of abundances or lengths), already computed by
   median_of_descending; the printf rounding is banker's (round half to even) */
inline auto report_median(double const median,
                          char const * const quantity,
                          struct Parameters const & parameters) -> void {
  if (not parameters.opt_quiet) {
    static_cast<void>(std::fprintf(stderr, "%s: %.0f\n", quantity, median));
  }
  if (parameters.fp_log != nullptr) {
    static_cast<void>(std::fprintf(parameters.fp_log, "%s: %.0f\n", quantity, median));
  }
}


/* Keep only the first n records (--topn), n being at least 1 (CLI-checked) */
template <typename Record>
auto truncate_deck(std::vector<Record> & deck,
                   long int const n_first_sequences) -> void {
  if (deck.size() > static_cast<unsigned long>(n_first_sequences)) {
    deck.resize(static_cast<std::size_t>(n_first_sequences));
  }
}


/* Write the database records back out in deck order, relabelled if requested;
   ordinals are 1-based */
template <typename Record>
auto output_sorted_fasta(std::vector<Record> const & deck,
                         std::FILE * output_file,
                         Database const & db,
                         struct Parameters const & parameters) -> void {
  Progress progress("Writing output", deck.size(), parameters);
  auto counter = std::size_t{0};
  for (auto const & record: deck) {
    fasta_print_db_relabel(output_file, record.seqno, counter + 1, db, parameters);
    progress.update(counter);
    ++counter;
  }
}
