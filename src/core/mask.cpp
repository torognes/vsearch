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

#include "utils/span.hpp"
#include "utils/view.hpp"  // View, make_view
#include "vsearch.hpp"
#include "utils/progress.hpp"
#include "core/mask.hpp"
#include "core/db.hpp"
#include "utils/ascii_case.hpp"  // to_upper
#include "utils/maps.hpp"
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include <algorithm>  // std::copy_n, std::fill, std::transform
#include <cassert>
#include <iterator>  // std::next
#include <array>
#include <cstddef>
#include <cstdint>  // int64_t, uint64_t
#include <limits>  // std::numeric_limits
#include <mutex>  // std::mutex, std::unique_lock
// #include <string>
#include <vector>


constexpr int dust_window = 64;


namespace {
/* The lowest-complexity region wo() found in the window it was given, as an
   offset pair relative to that window. A zero score means no region: either
   the window is too short to hold one, or nothing in it repeated. The two
   int out-parameters this replaces were only ever read when the score cleared
   the DUST level, so the caller no longer keeps them alive between windows.

   A plain aggregate, deliberately: in C++11 a member initializer would make it
   a non-aggregate, and the returns below would no longer brace-initialize.
   Every DustRegion here is built with braces, so every member is set. */
struct DustRegion {
  int score;
  int begin;
  int end;
};


/* Ceiling for the overflow contract wo() asserts in its inner loop. At file
   scope so that a release build, where the assert compiles away, does not see
   it as an unused local. */
constexpr auto max_sum = dust_window * dust_window / 2;  // 2048

auto wo(View<char> const window) -> DustRegion
{
  static constexpr auto dust_word = 3;
  static constexpr auto word_count = 1U << (2U * dust_word);  // 64
  static constexpr auto bitmask = word_count - 1;
  /* words[] is indexed by j < len below, so a longer window would run off the
     array; dust_core() passes at most dust_window by construction */
  assert(window.size() <= static_cast<std::size_t>(dust_window));
  auto const len = static_cast<int>(window.size());
  const auto l1 = len - dust_word + 1 - 5; /* smallest possible region is 8 */
  if (l1 < 0)
    {
      return DustRegion{};
    }

  auto bestv = 0;
  auto besti = 0;
  auto bestj = 0;
  /* both hold 6-bit quantities -- words[] is masked to bitmask, and counts[]
     rises by at most one per inner iteration, so it peaks at len - i - 2 <= 62.

     unsigned char rather than int is worth 1.09x on --fastx_mask, and the
     reason is the reset below, not cache footprint: 256 bytes would fit L1
     several times over. At 64 bytes GCC clears counts[] with four movaps,
     while at 128 or 256 it emits a rep stos whose microcoded startup costs
     about 23 cycles here -- paid on every one of the 45 M outer iterations a
     40 MB input runs. Widening counts[] back to unsigned short (128 bytes,
     still half of int) hands the entire gain back, measured: it is a
     threshold, not a gradient, so do not assume a smaller type is
     proportionally better. The remaining quarter of the gain is words[],
     which the inner loop streams through 1.4 G times. */
  std::array<unsigned char, word_count> counts {{}};
  std::array<unsigned char, dust_window> words {{}};
  auto word = 0U;

  for (auto j = 0; j < len; j++)
    {
      word <<= 2U;
      word |= map_2bit(window[static_cast<std::size_t>(j)]);
      words[static_cast<std::size_t>(j)] = static_cast<unsigned char>(word & bitmask);
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
              /* 10 * sum is the one product in this loop; sum counts pairs
                 among at most dust_window window positions, so it stays four
                 orders of magnitude below INT_MAX. The assert states that
                 bound rather than leaving it to be re-derived. */
              assert(sum >= 0 and sum <= max_sum);
              const auto v = 10 * sum / j;

              if (v > bestv)
                {
                  bestv = v;
                  besti = i;
                  bestj = j;
                }
            }
          /* c is counts[word] read above, and nothing has touched the array
             since, so the increment below cannot wrap the byte */
          assert(c < std::numeric_limits<unsigned char>::max());
          ++counts[word];
        }
    }

  return DustRegion{bestv, besti, besti + bestj};
}
}  // anonymous namespace


/* Core DUST implementation with explicit hardmask parameter.
   Thread-safe: does not read any globals. */
static auto dust_core(Span<char> const sequence, bool const use_hardmask) -> void
{
  static constexpr auto dust_level = 20;
  static constexpr auto half_dust_window = dust_window / 2;

  auto const len = static_cast<int>(sequence.size());

  /* make a local copy of the original sequence, including the terminator that
     sits just past the span (see the write below) */
  std::vector<char> local_seq(static_cast<std::size_t>(len) + 1);
  std::copy_n(sequence.data(), static_cast<std::size_t>(len) + 1, local_seq.data());

  if (!use_hardmask)
    {
      /* convert sequence to upper case unless hardmask in effect */
      std::transform(sequence.begin(), sequence.end(), sequence.begin(), to_upper);
      sequence.data()[len] = 0;  /* the terminator, which sits just past the span */
    }

  /* indexed, and the index is mutated in the body: a masked region short
     enough to end inside this window rewinds the next window's start
     (i += half_dust_window - b), so this is not a traversal */
  for (auto i = 0; i < len; i += half_dust_window)
    {
      const auto l = (len > i + dust_window) ? dust_window : len - i;
      auto const worst = wo(make_view(local_seq).subspan(static_cast<std::size_t>(i),
                                                         static_cast<std::size_t>(l)));

      if (worst.score > dust_level)
        {
          /* the low-complexity region wo() found, in sequence coordinates.
             Each bound is widened before the arithmetic rather than after it,
             and the subtraction is safe unsigned because wo() returns an end
             at or past its begin (it is begin + bestj, with bestj >= 0). */
          assert(worst.end >= worst.begin);
          auto const region = sequence.subspan(static_cast<std::size_t>(worst.begin)
                                               + static_cast<std::size_t>(i),
                                               static_cast<std::size_t>(worst.end)
                                               - static_cast<std::size_t>(worst.begin) + 1);
          if (use_hardmask)
            {
              std::fill(region.begin(), region.end(), 'N');
            }
          else
            {
              auto const * const original = std::next(local_seq.data(), worst.begin + i);
              std::transform(original, std::next(original, worst.end - worst.begin + 1),
                             region.begin(),
                             [](char const nucleotide) -> char {
                               return static_cast<char>(nucleotide | 32U);  // check_5th_bit (0x20)
                             });
            }

          if (worst.end < half_dust_window)
            {
              i += half_dust_window - worst.end;
            }
        }
    }
}


auto dust(Span<char> const seq, struct Parameters const & parameters) -> void
{
  dust_core(seq, parameters.opt_hardmask);
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

  auto const process_sequence = [&]() -> void {
    dust(db.mutable_sequence(seqno), *state.parameters);
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
                            [&state, &db](uint64_t /*nth_thread*/) -> void
                            { dust_all_worker(state, db); });
  threadrunner.run();
}


auto hardmask(Span<char> const seq) -> void
{
  /* convert all lower case letters in seq to N */
  static constexpr auto check_5th_bit = 32U; // 0x20
  static constexpr auto hardmask_char = 'N';
  std::transform(seq.begin(), seq.end(), seq.begin(),
                 [](char const nucleotide) -> char {
                   if ((static_cast<unsigned int>(static_cast<unsigned char>(nucleotide)) & check_5th_bit) != 0U)
                     {
                       return hardmask_char;
                     }
                   return nucleotide;
                 });
}


auto hardmask_all(struct Database & db) -> void
{
  for (uint64_t i = 0; i < db.getsequencecount(); i++)
    {
      hardmask(db.mutable_sequence(i));
    }
}


auto dust_single(char * seq, int const len, bool const use_hardmask) -> void
{
  dust_core(Span<char>{seq, static_cast<std::size_t>(len)}, use_hardmask);
}
