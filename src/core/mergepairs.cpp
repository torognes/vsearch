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
#include "core/mergepairs.hpp"
#include "core/mergepairs_internal.hpp"
#include "core/kmerhash.hpp"
#include "os/system.hpp"  // xmalloc, xfree
#include "utils/fatal.hpp"
#include "utils/kmer_hash_struct.hpp"
#include "utils/maps.hpp"
#include <algorithm>  // std::min, std::max, std::copy
#include <array>
#include <atomic>  // std::atomic, std::memory_order
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow, std::sqrt, std::round, std::log10, std::log2
#include <cstdint>  // int64_t, uint64_t
#include <cstring>  // std::strlen, std::strcpy, std::memcpy
#include <vector>


/* scores in bits */

constexpr auto k                          = 5;
constexpr auto merge_dropmax         = 16.0;
constexpr auto merge_mismatchmax     = -4.0;

/* static variables */

constexpr auto n_quality_symbols = 128U;
static std::array<std::array<char, n_quality_symbols>, n_quality_symbols> merge_qual_same {{}};
static std::array<std::array<char, n_quality_symbols>, n_quality_symbols> merge_qual_diff {{}};
static std::array<std::array<double, n_quality_symbols>, n_quality_symbols> match_score {{}};
static std::array<std::array<double, n_quality_symbols>, n_quality_symbols> mism_score {{}};
static std::array<double, n_quality_symbols> q2p {{}};


/* A worker must never call std::exit() (e.g. via fatal()) while sibling
   workers are still running: std::exit() flushes and closes the shared
   output streams and runs static destructors concurrently with threads
   that are still writing to those streams, which is a data race that
   intermittently corrupts libc state and crashes (observed as SIGILL on
   FreeBSD). Instead, an out-of-range FASTQ quality value records the
   error here and requests a cooperative abort; every worker then unwinds
   its loop and pair_all() reports the error and exits from the main
   thread, after all workers have joined. The error details are written
   once (first worker to claim wins) and read by pair_all() after the
   join, which establishes the needed happens-before. */

static std::atomic<bool> merge_abort {false};
static std::atomic<bool> merge_error_claimed {false};
static MergeAbortReason merge_error_reason = MergeAbortReason::quality_below_qmin;
static int merge_error_value = 0;

/* mutex_chunks and cond_chunks are owned as locals in pair_all(), not at
   file scope; see the comment there. */


/* Request a cooperative abort from a worker thread. Records the first
   error seen and signals every worker to stop; the actual message and
   std::exit() happen in pair_all() on the main thread after all workers
   have joined (see merge_abort above). */
auto request_merge_abort(MergeAbortReason const reason, int const value) -> void
{
  if (not merge_error_claimed.exchange(true))
    {
      merge_error_reason = reason;
      merge_error_value = value;
    }
  merge_abort.store(true, std::memory_order_release);
}


/* Poll the cooperative-abort flag. The default relaxed order matches the
   worker-loop hint checks; callers that need to synchronise with the recorded
   error (e.g. the post-join check in pair_all) pass a stronger order. */
auto merge_aborted(std::memory_order const order) -> bool
{
  return merge_abort.load(order);
}


/* Clear the cooperative-abort state before a run. */
auto merge_abort_reset() -> void
{
  merge_abort.store(false);
  merge_error_claimed.store(false);
}


/* Report the recorded worker error and terminate. Must be called from the
   main thread only, after all workers have joined. */
auto report_merge_abort(struct Parameters const & parameters) -> void
{
  switch (merge_error_reason)
    {
    case MergeAbortReason::quality_below_qmin:
      std::fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
              PRId64 ")\n",
              merge_error_value, parameters.opt_fastq_qmin);
      if (parameters.fp_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
                  PRId64 ")\n",
                  merge_error_value, parameters.opt_fastq_qmin);
        }
      break;

    case MergeAbortReason::quality_above_qmax:
      std::fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
              PRId64 ")\n",
              merge_error_value, parameters.opt_fastq_qmax);
      std::fprintf(stderr,
              "By default, quality values range from 0 to 41.\n"
              "To allow higher quality values, "
              "please use the option --fastq_qmax %d\n", merge_error_value);
      if (parameters.fp_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
                  PRId64 ")\n",
                  merge_error_value, parameters.opt_fastq_qmax);
          std::fprintf(parameters.fp_log,
                  "By default, quality values range from 0 to 41.\n"
                  "To allow higher quality values, "
                  "please use the option --fastq_qmax %d\n", merge_error_value);
        }
      break;

    case MergeAbortReason::more_fwd_than_rev:
      fatal("More forward reads than reverse reads");
      break;
    }
  std::exit(EXIT_FAILURE);
}


inline auto get_qual(char const quality_symbol, struct Parameters const & parameters) -> int
{
  assert(quality_symbol >= 33);
  assert(quality_symbol <= 126);

  auto const quality_value = static_cast<int>(quality_symbol - parameters.opt_fastq_ascii);

  if (quality_value < parameters.opt_fastq_qmin)
    {
      request_merge_abort(MergeAbortReason::quality_below_qmin, quality_value);
    }
  else if (quality_value > parameters.opt_fastq_qmax)
    {
      request_merge_abort(MergeAbortReason::quality_above_qmax, quality_value);
    }
  return quality_value;
}


inline auto q_to_p(int const quality_symbol, struct Parameters const & parameters) -> double
{
  static constexpr auto low_quality_threshold = 2;
  static constexpr auto max_probability = 0.75;
  static constexpr auto quality_divider = 10.0;
  static constexpr auto power_base = 10.0;

  assert(quality_symbol >= 33);
  assert(quality_symbol <= 126);

  auto const quality_value = static_cast<int>(quality_symbol - parameters.opt_fastq_ascii);

  // refactor: extract branch to a separate operation
  if (quality_value < low_quality_threshold) {
    return max_probability;
  }
  // probability = 10^-(quality / 10)
  return std::pow(power_base, -quality_value / quality_divider);
}


auto precompute_qual(struct Parameters const & parameters) -> void
{
  /* Precompute tables of scores etc */
  auto const qmaxout = static_cast<double>(parameters.opt_fastq_qmaxout);
  auto const qminout = static_cast<double>(parameters.opt_fastq_qminout);

  for (auto x = 33U; x <= 126U; x++)
    {
      auto const px = q_to_p(static_cast<int>(x), parameters);
      q2p[x] = px;

      for (auto y = 33U; y <= 126U; y++)
        {
          auto const py = q_to_p(static_cast<int>(y), parameters);

          auto p = 0.0;
          auto q = 0.0;

          /* Quality score equations from Edgar & Flyvbjerg (2015) */

          /* Match */
          p = px * py / 3.0 / (1.0 - px - py + (4.0 * px * py / 3.0));
          q = std::round(-10.0 * std::log10(p));
          q = std::min(q, qmaxout);
          q = std::max(q, qminout);
          merge_qual_same[x][y] = static_cast<char>(static_cast<double>(parameters.opt_fastq_ascii) + q);

          /* Mismatch, x is highest quality */
          p = px * (1.0 - (py / 3.0)) / (px + py - (4.0 * px * py / 3.0));
          q = std::round(-10.0 * std::log10(p));
          q = std::min(q, qmaxout);
          q = std::max(q, qminout);
          merge_qual_diff[x][y] = static_cast<char>(static_cast<double>(parameters.opt_fastq_ascii) + q);

          /*
            observed match,
            p = probability that they truly are identical,
            given error probabilites of px and py, resp.
          */

          // Given two initially identical aligned bases, and
          // the error probabilities px and py,
          // what is the probability of observing a match (or a mismatch)?

          p = 1.0 - px - py + (px * py * 4.0 / 3.0);
          match_score[x][y] = std::log2(p / 0.25);

          // Use a minimum mismatch penalty

          mism_score[x][y] = std::min(std::log2((1.0 - p) / 0.75), merge_mismatchmax);
        }
    }
}


auto merge_sym(char * sym,       char * qual,
               char const fwd_sym,     char const rev_sym,
               char const fwd_qual,    char const rev_qual) -> void
{
  if (rev_sym == 'N')
    {
      * sym = fwd_sym;
      * qual = fwd_qual;
    }
  else if (fwd_sym == 'N')
    {
      * sym = rev_sym;
      * qual = rev_qual;
    }
  else if (fwd_sym == rev_sym)
    {
      /* agreement */
      * sym = fwd_sym;
      * qual = merge_qual_same[static_cast<std::size_t>(fwd_qual)][static_cast<std::size_t>(rev_qual)];
    }
  else
    {
      /* disagreement */
      if (fwd_qual > rev_qual)
        {
          * sym = fwd_sym;
          * qual = merge_qual_diff[static_cast<std::size_t>(fwd_qual)][static_cast<std::size_t>(rev_qual)];
        }
      else
        {
          * sym = rev_sym;
          * qual = merge_qual_diff[static_cast<std::size_t>(rev_qual)][static_cast<std::size_t>(fwd_qual)];
        }
    }
}


auto merge(merge_data_t & a_read_pair, struct Parameters const & parameters) -> void
{
  /* length of 5' overhang of the forward sequence not merged
     with the reverse sequence */

  auto const fwd_5prime_overhang = (a_read_pair.fwd_trunc > a_read_pair.offset) ?
    a_read_pair.fwd_trunc - a_read_pair.offset : 0;

  // reset struct members
  a_read_pair.ee_merged = 0.0;
  a_read_pair.ee_fwd = 0.0;
  a_read_pair.ee_rev = 0.0;
  a_read_pair.fwd_errors = 0;
  a_read_pair.rev_errors = 0;
  auto sym = '\0';
  auto qual = '\0';
  int64_t fwd_pos = 0;
  int64_t rev_pos = 0;
  int64_t merged_pos = 0;
  auto ee = 0.0;

  merged_pos = 0;

  // 5' overhang in forward sequence

  fwd_pos = 0;

  while (fwd_pos < fwd_5prime_overhang)
    {
      sym = a_read_pair.fwd_sequence[static_cast<std::size_t>(fwd_pos)];
      qual = a_read_pair.fwd_quality[static_cast<std::size_t>(fwd_pos)];

      a_read_pair.merged_sequence[static_cast<std::size_t>(merged_pos)] = sym;
      a_read_pair.merged_quality_v[static_cast<std::size_t>(merged_pos)] = qual;

      ee = q2p[static_cast<std::size_t>(qual)];
      a_read_pair.ee_merged += ee;
      a_read_pair.ee_fwd += ee;

      ++fwd_pos;
      ++merged_pos;
    }

  // Merged region

  auto const rev_3prime_overhang = a_read_pair.offset > a_read_pair.fwd_trunc ?
    a_read_pair.offset - a_read_pair.fwd_trunc : 0;

  rev_pos = a_read_pair.rev_trunc - 1 - rev_3prime_overhang;

  while ((fwd_pos < a_read_pair.fwd_trunc) and (rev_pos >= 0))
    {
      auto fwd_sym = a_read_pair.fwd_sequence[static_cast<std::size_t>(fwd_pos)];
      auto rev_sym = map_complement(a_read_pair.rev_sequence[static_cast<std::size_t>(rev_pos)]);
      auto fwd_qual = a_read_pair.fwd_quality[static_cast<std::size_t>(fwd_pos)];
      auto rev_qual = a_read_pair.rev_quality[static_cast<std::size_t>(rev_pos)];

      merge_sym(& sym,
                & qual,
                fwd_qual < 2 ? 'N' : fwd_sym,
                rev_qual < 2 ? 'N' : rev_sym,
                fwd_qual,
                rev_qual);

      if (sym != fwd_sym)
        {
          ++a_read_pair.fwd_errors;
        }
      if (sym != rev_sym)
        {
          ++a_read_pair.rev_errors;
        }

      a_read_pair.merged_sequence[static_cast<std::size_t>(merged_pos)] = sym;
      a_read_pair.merged_quality_v[static_cast<std::size_t>(merged_pos)] = qual;
      a_read_pair.ee_merged += q2p[static_cast<std::size_t>(qual)];
      a_read_pair.ee_fwd += q2p[static_cast<std::size_t>(fwd_qual)];
      a_read_pair.ee_rev += q2p[static_cast<std::size_t>(rev_qual)];

      ++fwd_pos;
      --rev_pos;
      ++merged_pos;
    }

  // 5' overhang in reverse sequence

  while (rev_pos >= 0)
    {
      sym = map_complement(a_read_pair.rev_sequence[static_cast<std::size_t>(rev_pos)]);
      qual = a_read_pair.rev_quality[static_cast<std::size_t>(rev_pos)];

      a_read_pair.merged_sequence[static_cast<std::size_t>(merged_pos)] = sym;
      a_read_pair.merged_quality_v[static_cast<std::size_t>(merged_pos)] = qual;
      ++merged_pos;

      ee = q2p[static_cast<std::size_t>(qual)];
      a_read_pair.ee_merged += ee;
      a_read_pair.ee_rev += ee;

      --rev_pos;
    }

  auto const mergelen = merged_pos;
  a_read_pair.merged_length = mergelen;

  a_read_pair.merged_sequence[static_cast<std::size_t>(mergelen)] = 0;
  a_read_pair.merged_quality_v[static_cast<std::size_t>(mergelen)] = 0;

  if (a_read_pair.ee_merged <= parameters.opt_fastq_maxee)
    {
      a_read_pair.reason = Reason::ok;
      a_read_pair.merged = true;
    }
  else
    {
      a_read_pair.reason = Reason::maxee;
    }
}


auto optimize(merge_data_t & a_read_pair,
              struct kh_handle_s & kmerhash,
              struct Parameters const & parameters) -> int64_t
{
  /* Merge-acceptance thresholds, relaxed for short overlaps. Derived here from
     the run's opt_fastq_minovlen (threaded via parameters; clamped to >= 5 by
     both the CLI and library entry points) instead of from file-static globals
     the caller sets. Values match the former merge_mindiagcount/merge_minscore
     globals: defaults 4 / 16.0, relaxed when opt_fastq_minovlen < 9. */
  int const merge_mindiagcount = (parameters.opt_fastq_minovlen < 9)
      ? static_cast<int>(parameters.opt_fastq_minovlen - 4) : 4;
  double const merge_minscore = (parameters.opt_fastq_minovlen < 9)
      ? 1.6 * static_cast<double>(parameters.opt_fastq_minovlen) : 16.0;

  /* ungapped alignment in each diagonal */

  int64_t const i1 = 1;
  auto const i2 = a_read_pair.fwd_trunc + a_read_pair.rev_trunc - 1;

  auto best_score = 0.0;
  int64_t best_i = 0;
  int64_t best_diffs = 0;

  auto hits = 0;

  auto kmers = 0;

  std::vector<int> diags(static_cast<std::size_t>(a_read_pair.fwd_trunc + a_read_pair.rev_trunc), 0);

  kh_insert_kmers(kmerhash, k, a_read_pair.fwd_sequence.data(), static_cast<int>(a_read_pair.fwd_trunc));
  kh_find_diagonals(kmerhash, k, a_read_pair.rev_sequence.data(), static_cast<int>(a_read_pair.rev_trunc), diags);

  for (int64_t i = i1; i <= i2; i++)
    {
      int const diag = static_cast<int>(a_read_pair.rev_trunc + a_read_pair.fwd_trunc - i);
      auto const diagcount = diags[static_cast<std::size_t>(diag)];

      if (diagcount >= merge_mindiagcount)
        {
          kmers = 1;

          /* for each interesting diagonal */

          auto const fwd_3prime_overhang
            = i > a_read_pair.rev_trunc ? i - a_read_pair.rev_trunc : 0;
          auto const rev_3prime_overhang
            = i > a_read_pair.fwd_trunc ? i - a_read_pair.fwd_trunc : 0;
          auto const overlap
            = i - fwd_3prime_overhang - rev_3prime_overhang;
          auto const fwd_pos_start
            = a_read_pair.fwd_trunc - fwd_3prime_overhang - 1;
          auto const rev_pos_start
            = a_read_pair.rev_trunc - rev_3prime_overhang - overlap;

          auto fwd_pos = fwd_pos_start;
          auto rev_pos = rev_pos_start;
          auto score = 0.0;

          int64_t diffs = 0;
          auto score_high = 0.0;
          auto dropmax = 0.0;

          for (int64_t j = 0; j < overlap; j++)
            {
              /* for each pair of bases in the overlap */

              auto const fwd_sym = a_read_pair.fwd_sequence[static_cast<std::size_t>(fwd_pos)];
              auto const rev_sym = map_complement(a_read_pair.rev_sequence[static_cast<std::size_t>(rev_pos)]);

              unsigned int const fwd_qual = static_cast<unsigned int>(a_read_pair.fwd_quality[static_cast<std::size_t>(fwd_pos)]);
              unsigned int const rev_qual = static_cast<unsigned int>(a_read_pair.rev_quality[static_cast<std::size_t>(rev_pos)]);

              --fwd_pos;
              ++rev_pos;

              if (fwd_sym == rev_sym)
                {
                  score += match_score[fwd_qual][rev_qual];
                  score_high = std::max(score, score_high);
                }
              else
                {
                  score += mism_score[fwd_qual][rev_qual];
                  ++diffs;
                  if (score < score_high - dropmax)
                    {
                      dropmax = score_high - score;
                    }
                }
            }

          if (dropmax >= merge_dropmax)
            {
              score = 0.0;
            }

          if (score >= merge_minscore)
            {
              ++hits;
            }

          if (score > best_score)
            {
              best_score = score;
              best_i = i;
              best_diffs = diffs;
            }
        }
    }

  if (hits > 1)
    {
      a_read_pair.reason = Reason::repeat;
      return 0;
    }

  if ((not parameters.opt_fastq_allowmergestagger) and (best_i > a_read_pair.fwd_trunc))
    {
      a_read_pair.reason = Reason::staggered;
      return 0;
    }

  if (best_diffs > parameters.opt_fastq_maxdiffs)
    {
      a_read_pair.reason = Reason::maxdiffs;
      return 0;
    }

  if ((100.0 * static_cast<double>(best_diffs) / static_cast<double>(best_i)) > parameters.opt_fastq_maxdiffpct)
    {
      a_read_pair.reason = Reason::maxdiffpct;
      return 0;
    }

  if (kmers == 0)
    {
      a_read_pair.reason = Reason::nokmers;
      return 0;
    }

  if (best_score < merge_minscore)
    {
      a_read_pair.reason = Reason::minscore;
      return 0;
    }

  /* the effective minimum overlap is at least 5: the CLI path requires it and
     the library entry (mergepairs_single) threads a Parameters copy clamped to
     >= 5, so this reads the effective value from parameters (E1). */
  if (best_i < parameters.opt_fastq_minovlen)
    {
      a_read_pair.reason = Reason::minovlen;
      return 0;
    }

  int const mergelen = static_cast<int>(a_read_pair.fwd_trunc + a_read_pair.rev_trunc - best_i);

  if (mergelen < parameters.opt_fastq_minmergelen)
    {
      a_read_pair.reason = Reason::minmergelen;
      return 0;
    }

  if (mergelen > parameters.opt_fastq_maxmergelen)
    {
      a_read_pair.reason = Reason::maxmergelen;
      return 0;
    }

  return best_i;
}


auto process(merge_data_t & a_read_pair,
             struct kh_handle_s & kmerhash,
             struct Parameters const & parameters) -> void
{
  a_read_pair.merged = false;

  /* another worker may have hit an out-of-range quality value and
     requested a cooperative abort; stop doing work in that case */
  if (merge_abort.load(std::memory_order_acquire))
    {
      return;
    }

  auto skip = false;

  /* check length */

  if ((a_read_pair.fwd_length < parameters.opt_fastq_minlen) or
      (a_read_pair.rev_length < parameters.opt_fastq_minlen))
    {
      a_read_pair.reason = Reason::minlen;
      skip = true;
    }

  if ((a_read_pair.fwd_length > parameters.opt_fastq_maxlen) or
      (a_read_pair.rev_length > parameters.opt_fastq_maxlen))
    {
      a_read_pair.reason = Reason::maxlen;
      skip = true;
    }

  /* truncate sequences by quality */

  int64_t fwd_trunc = a_read_pair.fwd_length;

  if (not skip)
    {
      for (int64_t i = 0; i < a_read_pair.fwd_length; i++)
        {
          auto const quality_value = get_qual(a_read_pair.fwd_quality[static_cast<std::size_t>(i)], parameters);
          if (merge_abort.load(std::memory_order_relaxed))
            {
              return;
            }
          if (quality_value <= parameters.opt_fastq_truncqual)
            {
              fwd_trunc = i;
              break;
            }
        }
      if (fwd_trunc < parameters.opt_fastq_minlen)
        {
          a_read_pair.reason = Reason::minlen;
          skip = true;
        }
    }

  a_read_pair.fwd_trunc = fwd_trunc;

  auto rev_trunc = a_read_pair.rev_length;

  if (not skip)
    {
      for (int64_t i = 0; i < a_read_pair.rev_length; i++)
        {
          auto const quality_value = get_qual(a_read_pair.rev_quality[static_cast<std::size_t>(i)], parameters);
          if (merge_abort.load(std::memory_order_relaxed))
            {
              return;
            }
          if (quality_value <= parameters.opt_fastq_truncqual)
            {
              rev_trunc = i;
              break;
            }
        }
      if (rev_trunc < parameters.opt_fastq_minlen)
        {
          a_read_pair.reason = Reason::minlen;
          skip = true;
        }
    }

  a_read_pair.rev_trunc = rev_trunc;

  /* count n's */

  /* replace quality of N's by zero */

  if (not skip)
    {
      int64_t fwd_ncount = 0;
      for (int64_t i = 0; i < fwd_trunc; i++)
        {
          if (a_read_pair.fwd_sequence[static_cast<std::size_t>(i)] == 'N')
            {
              a_read_pair.fwd_quality[static_cast<std::size_t>(i)] = static_cast<char>(parameters.opt_fastq_ascii);
              ++fwd_ncount;
            }
        }
      if (fwd_ncount > parameters.opt_fastq_maxns)
        {
          a_read_pair.reason = Reason::maxns;
          skip = true;
        }
    }

  if (not skip)
    {
      int64_t rev_ncount = 0;
      for (int64_t i = 0; i < rev_trunc; i++)
        {
          if (a_read_pair.rev_sequence[static_cast<std::size_t>(i)] == 'N')
            {
              a_read_pair.rev_quality[static_cast<std::size_t>(i)] = static_cast<char>(parameters.opt_fastq_ascii);
              ++rev_ncount;
            }
        }
      if (rev_ncount > parameters.opt_fastq_maxns)
        {
          a_read_pair.reason = Reason::maxns;
          skip = true;
        }
    }

  a_read_pair.offset = 0;

  if (not skip)
    {
      a_read_pair.offset = optimize(a_read_pair, kmerhash, parameters);
    }

  if (a_read_pair.offset > 0)
    {
      merge(a_read_pair, parameters);
    }

  a_read_pair.state = State::processed;
}


/* === Library API for embedding paired-end merging === */


auto mergepairs_init(struct Parameters const & parameters) -> void
{
  /* The short-overlap relaxation of the merge-acceptance thresholds is now
     derived inside optimize() from opt_fastq_minovlen. mergepairs_single()
     threads a Parameters copy clamped to the >= 5 minimum the merge core
     requires, so optimize() sees the relaxed thresholds for short overlaps
     with no tunables to set here (matching the CLI path). */
  precompute_qual(parameters);
}


auto mergepairs_single(struct Parameters const & parameters,
                        const char * fwd_seq,
                        const char * fwd_qual,
                        int const fwd_len,
                        const char * rev_seq,
                        const char * rev_qual,
                        int const rev_len,
                        const char * fwd_header,
                        const char * rev_header,
                        struct merge_result_s * result) -> int
{
  /* Populate merge_data_t from caller's buffers */
  merge_data_t md {};

  md.fwd_length = fwd_len;
  md.rev_length = rev_len;
  md.fwd_trunc = fwd_len;
  md.rev_trunc = rev_len;

  /* Ensure buffers are large enough */
  int64_t max_len = std::max(md.fwd_length, md.rev_length);
  md.fwd_header.resize(std::strlen(fwd_header) + 1);
  md.rev_header.resize(std::strlen(rev_header) + 1);
  md.fwd_sequence.resize(static_cast<std::size_t>(max_len + 1));
  md.rev_sequence.resize(static_cast<std::size_t>(max_len + 1));
  md.fwd_quality.resize(static_cast<std::size_t>(max_len + 1));
  md.rev_quality.resize(static_cast<std::size_t>(max_len + 1));
  md.merged_sequence.resize(static_cast<std::size_t>(fwd_len + rev_len + 1));
  md.merged_quality_v.resize(static_cast<std::size_t>(fwd_len + rev_len + 1));

  std::strcpy(md.fwd_header.data(), fwd_header);
  std::strcpy(md.rev_header.data(), rev_header);
  std::memcpy(md.fwd_sequence.data(), fwd_seq, static_cast<std::size_t>(fwd_len));
  md.fwd_sequence[static_cast<std::size_t>(fwd_len)] = '\0';
  std::memcpy(md.rev_sequence.data(), rev_seq, static_cast<std::size_t>(rev_len));
  md.rev_sequence[static_cast<std::size_t>(rev_len)] = '\0';
  std::memcpy(md.fwd_quality.data(), fwd_qual, static_cast<std::size_t>(fwd_len));
  md.fwd_quality[static_cast<std::size_t>(fwd_len)] = '\0';
  std::memcpy(md.rev_quality.data(), rev_qual, static_cast<std::size_t>(rev_len));
  md.rev_quality[static_cast<std::size_t>(rev_len)] = '\0';

  /* Run the merge pipeline. The merge core enforces a minimum overlap of 5
     (see mergepairs_init); a library caller may pass a smaller value, so thread
     a clamped local copy rather than mutating the shared config global (E1). */
  struct Parameters clamped = parameters;
  if (clamped.opt_fastq_minovlen < 5)
    {
      clamped.opt_fastq_minovlen = 5;
    }
  struct kh_handle_s kmerhash;
  process(md, kmerhash, clamped);

  /* Populate result. Zero all fields including the pointers so that a
     failed merge leaves nullptr pointers for the caller. On success
     the buffers are xmalloc'd to the exact merged length; the caller
     owns them and must release via merge_result_free(). */
  *result = {};
  result->merged = md.merged;

  if (md.merged)
    {
      int const len = static_cast<int>(md.merged_length);
      result->merged_length = len;
      result->merged_sequence = static_cast<char *>(xmalloc(static_cast<std::size_t>(len) + 1));
      result->merged_quality = static_cast<char *>(xmalloc(static_cast<std::size_t>(len) + 1));
      std::memcpy(result->merged_sequence, md.merged_sequence.data(), static_cast<std::size_t>(len));
      result->merged_sequence[len] = '\0';
      std::memcpy(result->merged_quality, md.merged_quality_v.data(), static_cast<std::size_t>(len));
      result->merged_quality[len] = '\0';
      result->ee_merged = md.ee_merged;
      result->ee_fwd = md.ee_fwd;
      result->ee_rev = md.ee_rev;
      result->fwd_errors = static_cast<int>(md.fwd_errors);
      result->rev_errors = static_cast<int>(md.rev_errors);
      result->overlap_length = static_cast<int>(md.fwd_trunc + md.rev_trunc - md.merged_length);
      return 0;
    }

  return -1;  /* merge failed */
}


auto merge_result_free(struct merge_result_s * result) -> void
{
  if (result == nullptr)
    {
      return;
    }
  if (result->merged_sequence != nullptr)
    {
      xfree(result->merged_sequence);
      result->merged_sequence = nullptr;
    }
  if (result->merged_quality != nullptr)
    {
      xfree(result->merged_quality);
      result->merged_quality = nullptr;
    }
}
