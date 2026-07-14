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

#include "vsearch.h"
#include "utils/progress.hpp"
#include "core/mask.hpp"
#include "core/db.hpp"
#include "utils/maps.hpp"
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include <array>
#include <cctype>  // std::toupper
#include <cstdint>  // int64_t, uint64_t
#include <cstring>  // std::strcpy
#include <mutex>  // std::mutex, std::unique_lock
// #include <string>
#include <vector>


constexpr int dust_window = 64;


auto wo(int const len, const char *s, int *beg, int *end) -> int
{
  static constexpr auto dust_word = 3;
  static constexpr auto word_count = 1U << (2U * dust_word);  // 64
  static constexpr auto bitmask = word_count - 1;
  const auto l1 = len - dust_word + 1 - 5; /* smallest possible region is 8 */
  if (l1 < 0)
    {
      return 0;
    }

  auto bestv = 0;
  auto besti = 0;
  auto bestj = 0;
  std::array<int, word_count> counts {{}};
  std::array<int, dust_window> words {{}};
  auto word = 0U;

  for (auto j = 0; j < len; j++)
    {
      word <<= 2U;
      word |= map_2bit(s[j]);
      words[static_cast<std::size_t>(j)] = static_cast<int>(word & bitmask);
    }

  for (auto i = 0; i < l1; i++)
    {
      counts.fill(0);  // reset counts to zero

      auto sum = 0;

      for (auto j = dust_word - 1; j < len - i; j++)
        {
          word = static_cast<unsigned int>(words[static_cast<std::size_t>(i + j)]);
          const auto c = counts[word];
          if (c != 0)
            {
              sum += c;
              const auto v = 10 * sum / j;

              if (v > bestv)
                {
                  bestv = v;
                  besti = i;
                  bestj = j;
                }
            }
          ++counts[word];
        }
    }

  *beg = besti;
  *end = besti + bestj;

  return bestv;
}


/* Core DUST implementation with explicit hardmask parameter.
   Thread-safe: does not read any globals. */
static auto dust_core(char * seq, int const len, bool const use_hardmask) -> void
{
  static constexpr auto dust_level = 20;
  static constexpr auto half_dust_window = dust_window / 2;
  auto a = 0;
  auto b = 0;

  /* make a local copy of the original sequence */
  std::vector<char> local_seq(static_cast<std::size_t>(len) + 1);
  std::strcpy(local_seq.data(), seq);

  if (!use_hardmask)
    {
      /* convert sequence to upper case unless hardmask in effect */
      for (auto i = 0; i < len; i++)
        {
          seq[i] = static_cast<char>(std::toupper(seq[i]));
        }
      seq[len] = 0;
    }

  for (auto i = 0; i < len; i += half_dust_window)
    {
      const auto l = (len > i + dust_window) ? dust_window : len - i;
      const auto v = wo(l, &local_seq[static_cast<std::size_t>(i)], &a, &b);

      if (v > dust_level)
        {
          if (use_hardmask)
            {
              for (auto j = a + i; j <= b + i; j++)
                {
                  seq[j] = 'N';
                }
            }
          else
            {
              for (auto j = a + i; j <= b + i; j++)
                {
                  seq[j] = local_seq[static_cast<std::size_t>(j)] | 32U;  // check_5th_bit (0x20)
                }
            }

          if (b < half_dust_window)
            {
              i += half_dust_window - b;
            }
        }
    }
}


auto dust(char * seq, int const len, struct Parameters const & parameters) -> void
{
  dust_core(seq, len, parameters.opt_hardmask);
}


/* Per-invocation work-distribution state for dust_all(). This was three
   file-static globals (mutex / nextseq / seqcount); folding them into a local
   struct passed to the workers makes dust_all() reentrant and removes the
   shared mutable state, so a library caller can mask across sessions (or a
   future caller concurrently) without the counters bleeding between runs (E4). */
struct dust_state_s
{
  std::mutex mutex;
  uint64_t nextseq = 0;
  uint64_t seqcount = 0;
  Progress * progress = nullptr;  /* owner progress bar; worker updates it under state.mutex */
  Parameters const * parameters = nullptr;  /* set by dust_all(); read by dust() via the worker */
};


static auto dust_all_worker(struct dust_state_s & state, struct Database & db) -> void
{
  uint64_t seqno = 0;

  auto const has_work_to_claim = [&]() -> bool {
    if (state.nextseq >= state.seqcount) { return false; }
    seqno = state.nextseq++;
    state.progress->update(seqno);
    return true;
  };

  auto const process_sequence = [&]() {
    dust(db.mutatesequence(seqno),
         static_cast<int>(db.getsequencelen(seqno)),
         *state.parameters);
  };

  run_worker_loop(state.mutex, has_work_to_claim, process_sequence);
}


auto dust_all(struct Database & db, struct Parameters const & parameters) -> void
{
  struct dust_state_s state;
  state.seqcount = db.getsequencecount();
  state.parameters = &parameters;
  Progress progress("Masking", state.seqcount, parameters);
  state.progress = &progress;

  ThreadRunner threadrunner(static_cast<std::size_t>(parameters.opt_threads),
                            [&state, &db](uint64_t /*nth_thread*/)
                            { dust_all_worker(state, db); });
  threadrunner.run();
}


auto hardmask(char * seq, int const len) -> void
{
  /* convert all lower case letters in seq to N */
  // auto const * const end = std::next(seq, len);
  // refactoring: std::transform(seq, end, seq, [](unsigned char nuc){ if (std::islower(nuc) != 0) { return hardmask_char; } return nuc;});
  static constexpr auto check_5th_bit = 32U; // 0x20
  static constexpr auto hardmask_char = 'N';
  for (auto i = 0; i < len; i++)
    {
      if ((static_cast<unsigned int>(static_cast<unsigned char>(seq[i])) & check_5th_bit) != 0U)
        {
          seq[i] = hardmask_char;
        }
    }
}


auto hardmask_all(struct Database & db) -> void
{
  for (uint64_t i = 0; i < db.getsequencecount(); i++)
    {
      hardmask(db.mutatesequence(i), static_cast<int>(db.getsequencelen(i)));
    }
}


auto dust_single(char * seq, int const len, bool const use_hardmask) -> void
{
  dust_core(seq, len, use_hardmask);
}
