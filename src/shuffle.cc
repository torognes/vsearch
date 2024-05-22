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
#include <cstdio>  // std::FILE
#include <vector>


// refactoring:
// - std::vector<int> deck(dbsequencecount);
// - std::iota(deck.begin(), deck.end(), 0);
// - std::random_shuffle(deck.begin(), deck.end());  // uniform by default
// STOP: how to produce predictable results? (--randseed: seed fixing?)
// - auto const new_size = std::min(deck.size(), top_n);
// - deck.resize(new_size)
// - range for-loop

auto shuffle() -> void
{
  if (opt_output == nullptr) {
    fatal("Output file for shuffling must be specified with --output");
    return;
  }

  auto * fp_output = fopen_output(opt_output);
  if (fp_output == nullptr) {
    fatal("Unable to open shuffle output file for writing");
    return;
  }

  db_read(opt_shuffle, 0);
  show_rusage();

  int dbsequencecount = db_getsequencecount();
  std::vector<int> deck_v(dbsequencecount);
  auto * deck = deck_v.data();

  for(int i = 0; i < dbsequencecount; i++)
    {
      deck[i] = i;
    }

  int passed = 0;
  progress_init("Shuffling", dbsequencecount - 1);
  for(int i = dbsequencecount - 1; i > 0; i--)
    {
      /* generate a random number j in the range 0 to i, inclusive */
      int j = random_int(i + 1);

      /* exchange elements i and j */
      int t = deck[i];
      deck[i] = deck[j];
      deck[j] = t;

      ++passed;
      progress_update(passed);
    }
  progress_done();
  show_rusage();

  passed = MIN(dbsequencecount, opt_topn);

  progress_init("Writing output", passed);
  for(int i = 0; i < passed; i++)
    {
      fasta_print_db_relabel(fp_output, deck[i], i + 1);
      progress_update(i);
    }
  progress_done();

  show_rusage();

  db_free();
  fclose(fp_output);
}
