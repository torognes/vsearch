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
#include "fastq_mergepairs.h"
#include "kmerhash.h"
#include "utils/fatal.hpp"
#include "utils/kmer_hash_struct.hpp"
#include "utils/maps.hpp"
#include "utils/span.hpp"
#include "utils/threads.hpp"
#include <algorithm>  // std::copy, std::min, std::max
#include <array>
#include <atomic>  // std::atomic
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow, std::sqrt, std::round, std::log10, std::log2
#include <condition_variable>  // std::condition_variable
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::strlen
#include <mutex>  // std::mutex, std::unique_lock
#include <vector>


/* chunk constants */

constexpr auto chunk_size = 500; /* read pairs per chunk */
constexpr auto chunk_factor = 2; /* chunks per thread */

/* scores in bits */

constexpr auto k                          = 5;
static int merge_mindiagcount             = 4;
static double merge_minscore              = 16.0;
constexpr auto merge_dropmax         = 16.0;
constexpr auto merge_mismatchmax     = -4.0;

/* static variables */

constexpr auto n_quality_symbols = 128U;
std::array<std::array<char, n_quality_symbols>, n_quality_symbols> merge_qual_same {{}};
std::array<std::array<char, n_quality_symbols>, n_quality_symbols> merge_qual_diff {{}};
std::array<std::array<double, n_quality_symbols>, n_quality_symbols> match_score {{}};
std::array<std::array<double, n_quality_symbols>, n_quality_symbols> mism_score {{}};
std::array<double, n_quality_symbols> q2p {{}};

/* reasons for not merging:
   - undefined
   - ok
   - input seq too short (after truncation)
   - input seq too long
   - too many Ns in input
   - overlap too short
   - too many differences (maxdiffs)
   - too high percentage of differences (maxdiffpct)
   - staggered
   - indels in overlap region
   - potential repeats in overlap region / multiple overlaps
   - merged sequence too short
   - merged sequence too long
   - expected error too high
   - alignment score too low, insignificant, potential indel
   - too few kmers on same diag found
*/

enum struct Reason : char {
  undefined,
  ok,
  minlen,
  maxlen,
  maxns,
  minovlen,
  maxdiffs,
  maxdiffpct,
  staggered,
  indel,
  repeat,
  minmergelen,
  maxmergelen,
  maxee,
  minscore,
  nokmers
};

enum struct State: char {
  empty,
  filled,
  inprogress,
  processed
};

struct merge_data_s
{
  std::vector<char> fwd_header;
  std::vector<char> rev_header;
  std::vector<char> fwd_sequence;
  std::vector<char> rev_sequence;
  std::vector<char> fwd_quality;
  std::vector<char> rev_quality;
  int64_t header_alloc = 0;
  int64_t seq_alloc = 0;
  int64_t fwd_length = 0;
  int64_t rev_length = 0;
  int64_t fwd_trunc = 0;
  int64_t rev_trunc = 0;
  int64_t fwd_abundance = 1;
  int64_t rev_abundance = 1;
  int64_t pair_no = 0;
  std::vector<char> merged_sequence;
  std::vector<char> merged_quality_v;
  int64_t merged_length = 0;
  int64_t merged_seq_alloc = 0;
  double ee_merged = 0;
  double ee_fwd = 0;
  double ee_rev = 0;
  int64_t fwd_errors = 0;
  int64_t rev_errors = 0;
  int64_t offset = 0;
  bool merged = false;
  Reason reason = Reason::undefined;
  State state = State::empty;
};

using merge_data_t = struct merge_data_s;

struct chunk_s
{
  int size = 0; /* size of merge_data = number of pairs of reads */
  State state = State::empty; /* state of chunk: empty, read, processed */
  std::vector<struct merge_data_s> merge_data = std::vector<struct merge_data_s>(chunk_size);
};


/* Per-invocation state for a fastq_mergepairs run — previously the file-static
   output/input handles, the statistics counters, and the worker-pool chunk
   coordination block. Folding them into a struct that fastq_mergepairs() owns
   and threads through the output helpers (keep/discard), the reader
   (read_pair), the chunk drivers and the worker pool makes the command
   reentrant and removes the shared mutable state (E4). The library API
   (mergepairs_single) runs a single pair through the shared merge core
   (process) and writes into a caller-owned merge_result_s, so it uses none of
   this struct; only the per-run config (merge_mindiagcount/merge_minscore),
   the quality lookup tables and the cooperative-abort atomics remain shared,
   as they are read by that shared core. */
struct mergepairs_cli_state_s
{
  /* the run configuration, threaded through the CLI-path helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter. The
     shared merge core (get_qual/q_to_p/precompute_qual/merge/optimize/process)
     now takes a Parameters const & directly (E1 shared-infra phase), so both the
     CLI path and the library entry mergepairs_single() feed it the same config;
     only opt_fastq_minovlen stays on the global there (clamped in place by
     mergepairs_init). */
  struct Parameters const & parameters;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastqout_notmerged_fwd = nullptr;
  std::FILE * fp_fastqout_notmerged_rev = nullptr;
  std::FILE * fp_fastaout_notmerged_fwd = nullptr;
  std::FILE * fp_fastaout_notmerged_rev = nullptr;
  std::FILE * fp_eetabbedout = nullptr;

  fastx_handle fastq_fwd = nullptr;
  fastx_handle fastq_rev = nullptr;

  int64_t merged = 0;
  int64_t notmerged = 0;
  int64_t total = 0;

  double sum_read_length = 0.0;
  double sum_squared_fragment_length = 0.0;
  double sum_fragment_length = 0.0;

  double sum_ee_fwd = 0.0;
  double sum_ee_rev = 0.0;
  double sum_ee_merged = 0.0;
  uint64_t sum_errors_fwd = 0;
  uint64_t sum_errors_rev = 0;

  uint64_t failed_undefined = 0;
  uint64_t failed_minlen = 0;
  uint64_t failed_maxlen = 0;
  uint64_t failed_maxns = 0;
  uint64_t failed_minovlen = 0;
  uint64_t failed_maxdiffs = 0;
  uint64_t failed_maxdiffpct = 0;
  uint64_t failed_staggered = 0;
  uint64_t failed_indel = 0;
  uint64_t failed_repeat = 0;
  uint64_t failed_minmergelen = 0;
  uint64_t failed_maxmergelen = 0;
  uint64_t failed_maxee = 0;
  uint64_t failed_minscore = 0;
  uint64_t failed_nokmers = 0;

  std::vector<struct chunk_s> chunks;
  int chunk_count = 0;
  int chunk_read_next = 0;
  int chunk_process_next = 0;
  int chunk_write_next = 0;
  bool finished_reading = false;
  bool finished_all = false;
  int pairs_read = 0;
  int pairs_written = 0;

  explicit mergepairs_cli_state_s(struct Parameters const & params) : parameters(params) {}
};

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

enum class MergeAbortReason { quality_below_qmin, quality_above_qmax, more_fwd_than_rev };
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
inline auto request_merge_abort(MergeAbortReason const reason, int const value) -> void
{
  if (not merge_error_claimed.exchange(true))
    {
      merge_error_reason = reason;
      merge_error_value = value;
    }
  merge_abort.store(true, std::memory_order_release);
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
      if (fp_log != nullptr)
        {
          std::fprintf(fp_log,
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
      if (fp_log != nullptr)
        {
          std::fprintf(fp_log,
                  "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
                  PRId64 ")\n",
                  merge_error_value, parameters.opt_fastq_qmax);
          std::fprintf(fp_log,
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
               char fwd_sym,     char rev_sym,
               char fwd_qual,    char rev_qual) -> void
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


auto fprintf_ee_value(std::FILE * output_handle, double const expected_error) -> void
{
  /* mirror the variable-precision format used in fasta/fastq output
     (see fasta_print_general) so eetabbedout preserves small EE values */
  if (expected_error < 0.000000001) {
    std::fprintf(output_handle, "%.13lf", expected_error);
  } else if (expected_error < 0.00000001) {
    std::fprintf(output_handle, "%.12lf", expected_error);
  } else if (expected_error < 0.0000001) {
    std::fprintf(output_handle, "%.11lf", expected_error);
  } else if (expected_error < 0.000001) {
    std::fprintf(output_handle, "%.10lf", expected_error);
  } else if (expected_error < 0.00001) {
    std::fprintf(output_handle, "%.9lf", expected_error);
  } else if (expected_error < 0.0001) {
    std::fprintf(output_handle, "%.8lf", expected_error);
  } else if (expected_error < 0.001) {
    std::fprintf(output_handle, "%.7lf", expected_error);
  } else if (expected_error < 0.01) {
    std::fprintf(output_handle, "%.6lf", expected_error);
  } else if (expected_error < 0.1) {
    std::fprintf(output_handle, "%.5lf", expected_error);
  } else {
    std::fprintf(output_handle, "%.4lf", expected_error);
  }
}


auto keep(struct mergepairs_cli_state_s & state, merge_data_t const & a_read_pair) -> void
{
  ++state.merged;

  state.sum_fragment_length += static_cast<double>(a_read_pair.merged_length);
  state.sum_squared_fragment_length += static_cast<double>(a_read_pair.merged_length * a_read_pair.merged_length);

  state.sum_ee_merged += a_read_pair.ee_merged;
  state.sum_ee_fwd += a_read_pair.ee_fwd;
  state.sum_ee_rev += a_read_pair.ee_rev;
  state.sum_errors_fwd += static_cast<uint64_t>(a_read_pair.fwd_errors);
  state.sum_errors_rev += static_cast<uint64_t>(a_read_pair.rev_errors);

  if (state.parameters.opt_fastqout != nullptr)
    {
      fastq_print_general(state.fp_fastqout,
                          a_read_pair.merged_sequence.data(),
                          static_cast<int>(a_read_pair.merged_length),
                          a_read_pair.fwd_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.fwd_header.data())),
                          a_read_pair.merged_quality_v.data(),
                          static_cast<uint64_t>(static_cast<int>(a_read_pair.fwd_abundance)),
                          state.merged,
                          a_read_pair.ee_merged);
    }

  if (state.parameters.opt_fastaout != nullptr)
    {
      fasta_print_general(state.fp_fastaout,
                          nullptr,
                          a_read_pair.merged_sequence.data(),
                          static_cast<int>(a_read_pair.merged_length),
                          a_read_pair.fwd_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.fwd_header.data())),
                          static_cast<uint64_t>(static_cast<int>(a_read_pair.fwd_abundance)),
                          state.merged,
                          a_read_pair.ee_merged,
                          -1,
                          -1,
                          nullptr,
                          0.0,
                          0);
    }

  if (state.parameters.opt_eetabbedout != nullptr)
    {
      fprintf_ee_value(state.fp_eetabbedout, a_read_pair.ee_fwd);
      std::fprintf(state.fp_eetabbedout, "\t");
      fprintf_ee_value(state.fp_eetabbedout, a_read_pair.ee_rev);
      std::fprintf(state.fp_eetabbedout, "\t%" PRId64 "\t%" PRId64 "\n",
                   a_read_pair.fwd_errors, a_read_pair.rev_errors);
    }
}


auto discard(struct mergepairs_cli_state_s & state, merge_data_t const & a_read_pair) -> void
{
  switch (a_read_pair.reason)
    {
    case Reason::undefined:
      ++state.failed_undefined;
      break;

    case Reason::ok:
      break;

    case Reason::minlen:
      ++state.failed_minlen;
      break;

    case Reason::maxlen:
      ++state.failed_maxlen;
      break;

    case Reason::maxns:
      ++state.failed_maxns;
      break;

    case Reason::minovlen:
      ++state.failed_minovlen;
      break;

    case Reason::maxdiffs:
      ++state.failed_maxdiffs;
      break;

    case Reason::maxdiffpct:
      ++state.failed_maxdiffpct;
      break;

    case Reason::staggered:
      ++state.failed_staggered;
      break;

    case Reason::indel:
      ++state.failed_indel;
      break;

    case Reason::repeat:
      ++state.failed_repeat;
      break;

    case Reason::minmergelen:
      ++state.failed_minmergelen;
      break;

    case Reason::maxmergelen:
      ++state.failed_maxmergelen;
      break;

    case Reason::maxee:
      ++state.failed_maxee;
      break;

    case Reason::minscore:
      ++state.failed_minscore;
      break;

    case Reason::nokmers:
      ++state.failed_nokmers;
      break;
    }

  ++state.notmerged;

  if (state.parameters.opt_fastqout_notmerged_fwd != nullptr)
    {
      fastq_print_general(state.fp_fastqout_notmerged_fwd,
                          a_read_pair.fwd_sequence.data(),
                          static_cast<int>(a_read_pair.fwd_length),
                          a_read_pair.fwd_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.fwd_header.data())),
                          a_read_pair.fwd_quality.data(),
                          static_cast<uint64_t>(static_cast<int>(a_read_pair.fwd_abundance)),
                          state.notmerged,
                          -1.0);
    }

  if (state.parameters.opt_fastqout_notmerged_rev != nullptr)
    {
      fastq_print_general(state.fp_fastqout_notmerged_rev,
                          a_read_pair.rev_sequence.data(),
                          static_cast<int>(a_read_pair.rev_length),
                          a_read_pair.rev_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.rev_header.data())),
                          a_read_pair.rev_quality.data(),
                          static_cast<uint64_t>(static_cast<int>(a_read_pair.rev_abundance)),
                          state.notmerged,
                          -1.0);
    }

  if (state.parameters.opt_fastaout_notmerged_fwd != nullptr)
    {
      fasta_print_general(state.fp_fastaout_notmerged_fwd,
                          nullptr,
                          a_read_pair.fwd_sequence.data(),
                          static_cast<int>(a_read_pair.fwd_length),
                          a_read_pair.fwd_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.fwd_header.data())),
                          static_cast<uint64_t>(static_cast<int>(a_read_pair.fwd_abundance)),
                          state.notmerged,
                          -1.0,
                          -1, -1,
                          nullptr, 0.0,
                          0);
    }

  if (state.parameters.opt_fastaout_notmerged_rev != nullptr)
    {
      fasta_print_general(state.fp_fastaout_notmerged_rev,
                          nullptr,
                          a_read_pair.rev_sequence.data(),
                          static_cast<int>(a_read_pair.rev_length),
                          a_read_pair.rev_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.rev_header.data())),
                          static_cast<uint64_t>(static_cast<int>(a_read_pair.rev_abundance)),
                          state.notmerged,
                          -1.0,
                          -1, -1,
                          nullptr, 0.0,
                          0);
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

  /* opt_fastq_minovlen is clamped at run time by mergepairs_init()
     (library path) which writes the global, so this comparison stays on
     the global to see the clamped value; do not migrate to parameters. */
  if (best_i < opt_fastq_minovlen)
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


auto read_pair(struct mergepairs_cli_state_s & state, merge_data_t & a_read_pair) -> bool
{
  auto const fastq_fwd = state.fastq_fwd;
  auto const fastq_rev = state.fastq_rev;

  if (fastq_next(fastq_fwd, false, chrmap_upcase_vector.data()))
    {
      if (not fastq_next(fastq_rev, false, chrmap_upcase_vector.data()))
        {
          /* runs in a worker thread with the chunk lock released; request
             a cooperative abort instead of exiting here, and stop reading
             (pair_all() reports it from the main thread after join) */
          request_merge_abort(MergeAbortReason::more_fwd_than_rev, 0);
          return false;
        }

      /* allocate more memory if necessary */

      int64_t const fwd_header_len = static_cast<int64_t>(fastq_get_header_length(fastq_fwd));
      int64_t const rev_header_len = static_cast<int64_t>(fastq_get_header_length(fastq_rev));
      int64_t const header_needed = std::max(fwd_header_len, rev_header_len) + 1;

      if (header_needed > a_read_pair.header_alloc)
        {
          a_read_pair.header_alloc = header_needed;
          a_read_pair.fwd_header.resize(static_cast<std::size_t>(header_needed));
          a_read_pair.rev_header.resize(static_cast<std::size_t>(header_needed));
        }

      a_read_pair.fwd_length = static_cast<int64_t>(fastq_get_sequence_length(fastq_fwd));
      a_read_pair.rev_length = static_cast<int64_t>(fastq_get_sequence_length(fastq_rev));
      int64_t const seq_needed = std::max(a_read_pair.fwd_length, a_read_pair.rev_length) + 1;

      state.sum_read_length += static_cast<double>(a_read_pair.fwd_length + a_read_pair.rev_length);

      if (seq_needed > a_read_pair.seq_alloc)
        {
          a_read_pair.seq_alloc = seq_needed;
          a_read_pair.fwd_sequence.resize(static_cast<std::size_t>(seq_needed));
          a_read_pair.rev_sequence.resize(static_cast<std::size_t>(seq_needed));
          a_read_pair.fwd_quality.resize(static_cast<std::size_t>(seq_needed));
          a_read_pair.rev_quality.resize(static_cast<std::size_t>(seq_needed));
        }


      int64_t const merged_seq_needed = a_read_pair.fwd_length + a_read_pair.rev_length + 1;

      if (merged_seq_needed > a_read_pair.merged_seq_alloc)
        {
          a_read_pair.merged_seq_alloc = merged_seq_needed;
          a_read_pair.merged_sequence.resize(static_cast<std::size_t>(merged_seq_needed));
          a_read_pair.merged_quality_v.resize(static_cast<std::size_t>(merged_seq_needed));
        }

      /* make local copies of the seq, header and qual */

      auto const fwd_header_view = Span<char> {
        fastq_get_header(fastq_fwd),
        fastq_get_header_length(fastq_fwd)};
      std::copy(fwd_header_view.cbegin(), fwd_header_view.cend(), a_read_pair.fwd_header.begin());
      a_read_pair.fwd_header[fwd_header_view.size()] = '\0';  // fix issue when reusing allocated mem

      auto const rev_header_view = Span<char> {
        fastq_get_header(fastq_rev),
        fastq_get_header_length(fastq_rev)};
      std::copy(rev_header_view.cbegin(), rev_header_view.cend(), a_read_pair.rev_header.begin());
      a_read_pair.rev_header[rev_header_view.size()] = '\0';  // fix issue when reusing allocated mem

      auto const fwd_sequence_view = Span<char> {
        fastq_get_sequence(fastq_fwd),
        fastq_get_sequence_length(fastq_fwd)};
      std::copy(fwd_sequence_view.cbegin(), fwd_sequence_view.cend(), a_read_pair.fwd_sequence.begin());

      auto const rev_sequence_view = Span<char> {
        fastq_get_sequence(fastq_rev),
        fastq_get_sequence_length(fastq_rev)};
      std::copy(rev_sequence_view.cbegin(), rev_sequence_view.cend(), a_read_pair.rev_sequence.begin());

      auto const fwd_quality_view = Span<char> {
        fastq_get_quality(fastq_fwd),
        fastq_get_quality_length(fastq_fwd)};
      std::copy(fwd_quality_view.cbegin(), fwd_quality_view.cend(), a_read_pair.fwd_quality.begin());

      auto const rev_quality_view = Span<char> {
        fastq_get_quality(fastq_rev),
        fastq_get_quality_length(fastq_rev)};
      std::copy(rev_quality_view.cbegin(), rev_quality_view.cend(), a_read_pair.rev_quality.begin());

      a_read_pair.fwd_abundance = fastq_get_abundance(fastq_fwd);
      a_read_pair.rev_abundance = fastq_get_abundance(fastq_rev);

      a_read_pair.merged_sequence[0] = 0;
      a_read_pair.merged_quality_v[0] = 0;
      a_read_pair.merged = false;
      a_read_pair.pair_no = state.total++;

      return true;
    }
  return false;
}


auto keep_or_discard(struct mergepairs_cli_state_s & state, merge_data_t const & a_read_pair) -> void
{
  if (a_read_pair.merged)
    {
      keep(state, a_read_pair);
    }
  else
    {
      discard(state, a_read_pair);
    }
}


inline auto chunk_perform_read(struct mergepairs_cli_state_s & state,
                               std::unique_lock<std::mutex> & lock,
                               std::condition_variable & cond_chunks) -> void
{
  while ((not state.finished_reading) and (state.chunks[static_cast<std::size_t>(state.chunk_read_next)].state == State::empty))
    {
      lock.unlock();
      progress_update(fastq_get_position(state.fastq_fwd));
      auto r = 0;
      while ((r < chunk_size) and
             read_pair(state, state.chunks[static_cast<std::size_t>(state.chunk_read_next)].merge_data[static_cast<std::size_t>(r)]))
        {
          ++r;
        }
      state.chunks[static_cast<std::size_t>(state.chunk_read_next)].size = r;
      lock.lock();
      state.pairs_read += r;
      if (r > 0)
        {
          state.chunks[static_cast<std::size_t>(state.chunk_read_next)].state = State::filled;
          state.chunk_read_next = (state.chunk_read_next + 1) % state.chunk_count;
        }
      if (r < chunk_size)
        {
          state.finished_reading = true;
          if (state.pairs_written >= state.pairs_read)
            {
              state.finished_all = true;
            }
        }
      cond_chunks.notify_all();
    }
}


inline auto chunk_perform_write(struct mergepairs_cli_state_s & state,
                                std::unique_lock<std::mutex> & lock,
                                std::condition_variable & cond_chunks) -> void
{
  while (state.chunks[static_cast<std::size_t>(state.chunk_write_next)].state == State::processed)
    {
      lock.unlock();
      for (auto i = 0; i < state.chunks[static_cast<std::size_t>(state.chunk_write_next)].size; i++)
        {
          keep_or_discard(state, state.chunks[static_cast<std::size_t>(state.chunk_write_next)].merge_data[static_cast<std::size_t>(i)]);
        }
      lock.lock();
      state.pairs_written += state.chunks[static_cast<std::size_t>(state.chunk_write_next)].size;
      state.chunks[static_cast<std::size_t>(state.chunk_write_next)].state = State::empty;
      if (state.finished_reading and (state.pairs_written >= state.pairs_read))
        {
          state.finished_all = true;
        }
      state.chunk_write_next = (state.chunk_write_next + 1) % state.chunk_count;
      cond_chunks.notify_all();
    }
}


inline auto chunk_perform_process(struct mergepairs_cli_state_s & state,
                                  struct kh_handle_s & kmerhash,
                                  std::unique_lock<std::mutex> & lock,
                                  std::condition_variable & cond_chunks) -> void
{
  auto const chunk_current = state.chunk_process_next;
  if (state.chunks[static_cast<std::size_t>(chunk_current)].state == State::filled)
    {
      state.chunks[static_cast<std::size_t>(chunk_current)].state = State::inprogress;
      state.chunk_process_next = (chunk_current + 1) % state.chunk_count;
      cond_chunks.notify_all();
      lock.unlock();
      for (auto i = 0; i < state.chunks[static_cast<std::size_t>(chunk_current)].size; i++)
        {
          if (merge_abort.load(std::memory_order_relaxed))
            {
              break;
            }
          process(state.chunks[static_cast<std::size_t>(chunk_current)].merge_data[static_cast<std::size_t>(i)], kmerhash, state.parameters);
        }
      lock.lock();
      state.chunks[static_cast<std::size_t>(chunk_current)].state = State::processed;
      cond_chunks.notify_all();
    }
}


auto pair_worker(struct mergepairs_cli_state_s & state,
                 uint64_t t,
                 std::mutex & mutex_chunks,
                 std::condition_variable & cond_chunks) -> void
{
  /* new */

  struct kh_handle_s kmerhash;

  std::unique_lock<std::mutex> lock(mutex_chunks);

  while (not state.finished_all)
    {
      /* a worker hit an out-of-range quality value: stop the whole pool.
         finished_all is set under the lock so the wait predicates below
         (which test it) release, and notify_all wakes any sleepers. The
         error is reported from the main thread in pair_all() after join. */
      if (merge_abort.load(std::memory_order_relaxed))
        {
          state.finished_all = true;
          cond_chunks.notify_all();
          break;
        }

      if (state.parameters.opt_threads == 1)
        {
          /* One thread does it all */
          chunk_perform_read(state, lock, cond_chunks);
          chunk_perform_process(state, kmerhash, lock, cond_chunks);
          chunk_perform_write(state, lock, cond_chunks);
        }
      else if (state.parameters.opt_threads == 2)
        {
          if (t == 0)
            {
              /* first thread reads and processes */
              while (not
                     (
                      state.finished_all
                      or
                      (state.chunks[static_cast<std::size_t>(state.chunk_process_next)].state == State::filled)
                      or
                      ((not state.finished_reading) and
                       state.chunks[static_cast<std::size_t>(state.chunk_read_next)].state == State::empty)))
                {
                  cond_chunks.wait(lock);
                }

              chunk_perform_read(state, lock, cond_chunks);
              chunk_perform_process(state, kmerhash, lock, cond_chunks);
            }
          else /* t == 1 */
            {
              /* second thread writes and processes */
              while (not
                     (
                      state.finished_all
                      or
                      (state.chunks[static_cast<std::size_t>(state.chunk_process_next)].state == State::filled)
                      or
                      (state.chunks[static_cast<std::size_t>(state.chunk_write_next)].state == State::processed)
                      )
                     )
                {
                  cond_chunks.wait(lock);
                }

              chunk_perform_write(state, lock, cond_chunks);
              chunk_perform_process(state, kmerhash, lock, cond_chunks);
            }
        }
      else
        {
          if (t == 0)
            {
              /* first thread reads and processes */
              while (not
                     (
                      state.finished_all
                      or
                      ((not state.finished_reading) and
                       (state.chunks[static_cast<std::size_t>(state.chunk_read_next)].state == State::empty))
                      or
                      (state.chunks[static_cast<std::size_t>(state.chunk_process_next)].state == State::filled)
                      )
                     )
                {
                  cond_chunks.wait(lock);
                }

              chunk_perform_read(state, lock, cond_chunks);
              chunk_perform_process(state, kmerhash, lock, cond_chunks);
            }
          else if (t == static_cast<uint64_t>(state.parameters.opt_threads) - 1)
            {
              /* last thread writes and processes */
              while (not
                     (
                      state.finished_all
                      or
                      (state.chunks[static_cast<std::size_t>(state.chunk_write_next)].state == State::processed)
                      or
                      (state.chunks[static_cast<std::size_t>(state.chunk_process_next)].state == State::filled)
                      )
                     )
                {
                  cond_chunks.wait(lock);
                }

              chunk_perform_write(state, lock, cond_chunks);
              chunk_perform_process(state, kmerhash, lock, cond_chunks);
            }
          else
            {
              /* the other threads are only processing */
              while (not
                     (
                      state.finished_all
                      or
                      (state.chunks[static_cast<std::size_t>(state.chunk_process_next)].state == State::filled)
                      )
                     )
                {
                  cond_chunks.wait(lock);
                }

              chunk_perform_process(state, kmerhash, lock, cond_chunks);
            }
        }
    }

  /* mutex_chunks released by RAII on return */
}


auto pair_all(struct mergepairs_cli_state_s & state) -> void
{
  /* prepare chunks */

  state.chunk_count = static_cast<int>(chunk_factor * state.parameters.opt_threads);
  state.chunk_read_next = 0;
  state.chunk_process_next = 0;
  state.chunk_write_next = 0;

  /* reset the cooperative-abort state (file statics persist across
     library-API sessions) */
  merge_abort.store(false);
  merge_error_claimed.store(false);

  state.chunks.resize(static_cast<std::size_t>(state.chunk_count));

  /* The chunk mutex and condition variable are locals (not file scope) so
     their lifetime is scoped to the worker pool. Combined with the
     cooperative abort (see merge_abort), no worker ever calls std::exit():
     the only exit happens in report_merge_abort() on the main thread after
     ThreadRunner has joined every worker, so the condition variable is
     never destroyed (or left) with waiters present. */
  std::mutex mutex_chunks;
  std::condition_variable cond_chunks;

  /* run the worker pool; the workers coordinate through mutex_chunks and
     cond_chunks until all chunks have been read, processed and written */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(state.parameters.opt_threads),
                              [&state, &mutex_chunks, &cond_chunks](uint64_t nth_thread) {
                                pair_worker(state, nth_thread, mutex_chunks, cond_chunks);
                              });
    threadrunner.run();
  }

  /* all workers have joined; if one hit an out-of-range quality value,
     report it and exit now, single-threaded, so the message reliably
     reaches stderr and the --log file and no stdio teardown races a live
     worker thread */
  if (merge_abort.load())
    {
      report_merge_abort(state.parameters);
    }
}


auto print_stats(struct mergepairs_cli_state_s const & state, std::FILE * output_handle) -> void
{
  /* read-only aliases for the counters used repeatedly below; the
     single-use counters are referenced as state.<member> directly */
  auto const & total = state.total;
  auto const & merged = state.merged;
  auto const & notmerged = state.notmerged;
  auto const & sum_fragment_length = state.sum_fragment_length;
  auto const & sum_errors_fwd = state.sum_errors_fwd;
  auto const & sum_errors_rev = state.sum_errors_rev;

  std::fprintf(output_handle,
          "%10" PRId64 "  Pairs\n",
          total);

  std::fprintf(output_handle,
          "%10" PRId64 "  Merged",
          merged);
  if (total > 0)
    {
      std::fprintf(output_handle,
              " (%.1lf%%)",
              100.0 * static_cast<double>(merged) / static_cast<double>(total));
    }
  std::fprintf(output_handle, "\n");

  std::fprintf(output_handle,
          "%10" PRId64 "  Not merged",
          notmerged);
  if (total > 0)
    {
      std::fprintf(output_handle,
              " (%.1lf%%)",
              100.0 * static_cast<double>(notmerged) / static_cast<double>(total));
    }
  std::fprintf(output_handle, "\n");

  if (notmerged > 0)
    {
      std::fprintf(output_handle, "\nPairs that failed merging due to various reasons:\n");
    }

  if (state.failed_undefined != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  undefined reason\n",
              state.failed_undefined);
    }

  if (state.failed_minlen != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  reads too short (after truncation)\n",
              state.failed_minlen);
    }

  if (state.failed_maxlen != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  reads too long (after truncation)\n",
              state.failed_maxlen);
    }

  if (state.failed_maxns != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  too many N's\n",
              state.failed_maxns);
    }

  if (state.failed_nokmers != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  too few kmers found on same diagonal\n",
              state.failed_nokmers);
    }

  if (state.failed_repeat != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  multiple potential alignments\n",
              state.failed_repeat);
    }

  if (state.failed_maxdiffs != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  too many differences\n",
              state.failed_maxdiffs);
    }

  if (state.failed_maxdiffpct != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  too high percentage of differences\n",
              state.failed_maxdiffpct);
    }

  if (state.failed_minscore != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  alignment score too low, or score drop too high\n",
              state.failed_minscore);
    }

  if (state.failed_minovlen != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  overlap too short\n",
              state.failed_minovlen);
    }

  if (state.failed_maxee != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  expected error too high\n",
              state.failed_maxee);
    }

  if (state.failed_minmergelen != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  merged fragment too short\n",
              state.failed_minmergelen);
    }

  if (state.failed_maxmergelen != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  merged fragment too long\n",
              state.failed_maxmergelen);
    }

  if (state.failed_staggered != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  staggered read pairs\n",
              state.failed_staggered);
    }

  if (state.failed_indel != 0U)
    {
      std::fprintf(output_handle,
              "%10" PRIu64 "  indel errors\n",
              state.failed_indel);
    }

  std::fprintf(output_handle, "\n");

  if (total > 0)
    {
      std::fprintf(output_handle, "Statistics of all reads:\n");

      auto const mean_read_length = state.sum_read_length / (2.0 * state.pairs_read);

      std::fprintf(output_handle,
              "%10.2f  Mean read length\n",
              mean_read_length);
    }

  if (merged > 0)
    {
      std::fprintf(output_handle, "\n");

      std::fprintf(output_handle, "Statistics of merged reads:\n");

      auto const mean = sum_fragment_length / static_cast<double>(merged);

      std::fprintf(output_handle,
              "%10.2f  Mean fragment length\n",
              mean);

      auto const stdev = std::sqrt((state.sum_squared_fragment_length
                               - (2.0 * mean * sum_fragment_length)
                               + (mean * mean * static_cast<double>(merged)))
                              / (static_cast<double>(merged) + 0.0));

      std::fprintf(output_handle,
              "%10.2f  Standard deviation of fragment length\n",
              stdev);

      std::fprintf(output_handle,
              "%10.2f  Mean expected error in forward sequences\n",
              state.sum_ee_fwd / static_cast<double>(merged));

      std::fprintf(output_handle,
              "%10.2f  Mean expected error in reverse sequences\n",
              state.sum_ee_rev / static_cast<double>(merged));

      std::fprintf(output_handle,
              "%10.2f  Mean expected error in merged sequences\n",
              state.sum_ee_merged / static_cast<double>(merged));

      std::fprintf(output_handle,
              "%10.2f  Mean observed errors in merged region of forward sequences\n",
              1.0 * static_cast<double>(sum_errors_fwd) / static_cast<double>(merged));

      std::fprintf(output_handle,
              "%10.2f  Mean observed errors in merged region of reverse sequences\n",
              1.0 * static_cast<double>(sum_errors_rev) / static_cast<double>(merged));

      std::fprintf(output_handle,
              "%10.2f  Mean observed errors in merged region\n",
              1.0 * static_cast<double>(sum_errors_fwd + sum_errors_rev) / static_cast<double>(merged));
    }
}


auto fastq_mergepairs(struct Parameters const & parameters) -> void
{
  /* Per-invocation state, owned here and threaded through the worker pool and
     the output helpers (E4). Aliased by reference so the body below reads
     unchanged; the workers receive `state`, not file-static globals. */
  struct mergepairs_cli_state_s state(parameters);
  auto & fastq_fwd = state.fastq_fwd;
  auto & fastq_rev = state.fastq_rev;
  auto & fp_fastqout = state.fp_fastqout;
  auto & fp_fastaout = state.fp_fastaout;
  auto & fp_fastqout_notmerged_fwd = state.fp_fastqout_notmerged_fwd;
  auto & fp_fastqout_notmerged_rev = state.fp_fastqout_notmerged_rev;
  auto & fp_fastaout_notmerged_fwd = state.fp_fastaout_notmerged_fwd;
  auto & fp_fastaout_notmerged_rev = state.fp_fastaout_notmerged_rev;
  auto & fp_eetabbedout = state.fp_eetabbedout;

  /* fatal error if specified overlap is too small */

  if (opt_fastq_minovlen < 5)
    {
      fatal("Overlap specified with --fastq_minovlen must be at least 5");
    }

  /* relax default parameters in case of short overlaps */

  if (opt_fastq_minovlen < 9)
    {
      merge_mindiagcount = static_cast<int>(opt_fastq_minovlen - 4);
      merge_minscore = 1.6 * static_cast<double>(opt_fastq_minovlen);
    }

  /* open input files */

  fastq_fwd = fastq_open(parameters.opt_fastq_mergepairs);
  fastq_rev = fastq_open(parameters.opt_reverse);

  /* open output files */

  fp_fastqout = open_optional_output(parameters.opt_fastqout, "fastqout");
  fp_fastaout = open_optional_output(parameters.opt_fastaout, "fastaout");
  fp_fastqout_notmerged_fwd = open_optional_output(parameters.opt_fastqout_notmerged_fwd, "fastqout_notmerged_fwd");
  fp_fastqout_notmerged_rev = open_optional_output(parameters.opt_fastqout_notmerged_rev, "fastqout_notmerged_rev");
  fp_fastaout_notmerged_fwd = open_optional_output(parameters.opt_fastaout_notmerged_fwd, "fastaout_notmerged_fwd");
  fp_fastaout_notmerged_rev = open_optional_output(parameters.opt_fastaout_notmerged_rev, "fastaout_notmerged_rev");
  fp_eetabbedout = open_optional_output(parameters.opt_eetabbedout, "eetabbedout");

  /* precompute merged quality values */

  precompute_qual(parameters);

  /* main */

  uint64_t const filesize = fastq_get_size(fastq_fwd);
  progress_init("Merging reads", filesize);

  if (not fastq_fwd->is_empty)
    {
      pair_all(state);
    }

  progress_done();

  if (fastq_next(fastq_rev, true, chrmap_upcase_vector.data()))
    {
      fatal("More reverse reads than forward reads");
    }

  if (fp_log != nullptr) {
    print_stats(state, fp_log);
  }
  else {
    print_stats(state, stderr);
  }

  /* clean up */

  /* fclose_output() is a no-op on a null handle, so unopened outputs need
     no guard. */
  fclose_output(fp_eetabbedout);
  fclose_output(fp_fastaout_notmerged_rev);
  fclose_output(fp_fastaout_notmerged_fwd);
  fclose_output(fp_fastqout_notmerged_rev);
  fclose_output(fp_fastqout_notmerged_fwd);
  fclose_output(fp_fastaout);
  fclose_output(fp_fastqout);

  fastq_close(fastq_rev);
  fastq_rev = nullptr;
  fastq_close(fastq_fwd);
  fastq_fwd = nullptr;
}


/* === Library API for embedding paired-end merging === */


auto mergepairs_init(struct Parameters const & parameters) -> void
{
  /* Relax default parameters for short overlaps, matching the CLI path
     (fastq_mergepairs lines 1560-1571). Without this, short-overlap
     merges that work via the CLI will silently fail via the API.
     opt_fastq_minovlen is clamped in place on the global here: optimize()
     reads the clamped value from the global, so it stays global (do not
     migrate to parameters until this writer is refactored). */
  if (opt_fastq_minovlen < 5)
    {
      opt_fastq_minovlen = 5;
    }
  if (opt_fastq_minovlen < 9)
    {
      merge_mindiagcount = static_cast<int>(opt_fastq_minovlen - 4);
      merge_minscore = 1.6 * static_cast<double>(opt_fastq_minovlen);
    }

  precompute_qual(parameters);
}


auto mergepairs_single(struct Parameters const & parameters,
                        const char * fwd_seq,
                        const char * fwd_qual,
                        int fwd_len,
                        const char * rev_seq,
                        const char * rev_qual,
                        int rev_len,
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

  /* Run the merge pipeline */
  struct kh_handle_s kmerhash;
  process(md, kmerhash, parameters);

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
