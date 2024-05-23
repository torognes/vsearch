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
#include <algorithm>  // std::min, std::shuffle
#include <cstdio>  // std::FILE
#include <numeric>  // std::iota
#include <random>
#include <vector>


auto create_deck() -> std::vector<int> {
  auto const dbsequencecount = db_getsequencecount();
  std::vector<int> deck(dbsequencecount);
  std::iota(deck.begin(), deck.end(), 0);
  return deck;
}


auto generate_seed(long int const user_seed) -> unsigned int {
  if (user_seed != 0) {
    return static_cast<unsigned int>(user_seed);
  }
  std::random_device number_generator;
  return number_generator();
}


auto shuffle_deck(std::vector<int> & deck) -> void {
  static constexpr auto one_hundred_percent = 100ULL;
  progress_init("Shuffling", one_hundred_percent);
  auto const seed = generate_seed(opt_randseed);
  std::mt19937_64 uniform_generator(seed);
  std::shuffle(deck.begin(), deck.end(), uniform_generator);
  progress_done();
}


auto output_shuffled_fasta(std::vector<int> &deck,
                           long int const n_first_sequences,
                           std::FILE * output_file) -> void {
  auto const final_size = std::min(deck.size(),
                                   static_cast<unsigned long>(n_first_sequences));
  deck.resize(final_size);
  progress_init("Writing output", deck.size());
  auto counter = 0;
  for(auto const sequence_id: deck) {
    fasta_print_db_relabel(output_file, sequence_id, counter + 1);
    progress_update(counter);
    ++counter;
  }
  progress_done();
}


auto shuffle() -> void
{
  if (opt_output == nullptr) {
    fatal("Output file for shuffling must be specified with --output");
  }

  auto * fp_output = fopen_output(opt_output);
  if (fp_output == nullptr) {
    fatal("Unable to open shuffle output file for writing");
  }

  db_read(opt_shuffle, 0);
  show_rusage();

  auto deck = create_deck();
  shuffle_deck(deck);
  show_rusage();

  output_shuffled_fasta(deck, opt_topn, fp_output);
  show_rusage();

  db_free();
  if (fp_output != nullptr) {
    static_cast<void>(fclose(fp_output));
  }
}
