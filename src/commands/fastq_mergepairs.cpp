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
#include "commands/fastq_mergepairs.hpp"
#include "core/mergepairs_internal.hpp"
#include "core/kmerhash.hpp"
#include "utils/fatal.hpp"
#include "utils/kmer_hash_struct.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
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
   this struct. The merge-acceptance thresholds are derived inside optimize()
   from the threaded Parameters; only the quality lookup tables and the
   cooperative-abort atomics remain shared file-static state read by that core. */
struct mergepairs_cli_state_s
{
  /* the run configuration, threaded through the CLI-path helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter. The
     shared merge core (get_qual/q_to_p/precompute_qual/merge/optimize/process)
     now takes a Parameters const & directly (E1 shared-infra phase), so both the
     CLI path and the library entry mergepairs_single() feed it the same config.
     The <5 minimum-overlap clamp is applied to a local Parameters copy that is
     threaded to the merge core, so opt_fastq_minovlen no longer needs the global. */
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

  Progress * progress = nullptr;  /* owner progress bar; the chunk reader updates it */

  explicit mergepairs_cli_state_s(struct Parameters const & params) : parameters(params) {}
};


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
                          static_cast<uint64_t>(a_read_pair.fwd_abundance),
                          state.merged,
                          a_read_pair.ee_merged,
                          state.parameters);
    }

  if (state.parameters.opt_fastaout != nullptr)
    {
      fasta_print_general(state.fp_fastaout,
                          nullptr,
                          a_read_pair.merged_sequence.data(),
                          static_cast<int>(a_read_pair.merged_length),
                          a_read_pair.fwd_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.fwd_header.data())),
                          static_cast<uint64_t>(a_read_pair.fwd_abundance),
                          state.merged,
                          a_read_pair.ee_merged,
                          -1,
                          -1,
                          nullptr,
                          0.0,
                          0,
                          state.parameters);
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
                          static_cast<uint64_t>(a_read_pair.fwd_abundance),
                          state.notmerged,
                          -1.0,
                          state.parameters);
    }

  if (state.parameters.opt_fastqout_notmerged_rev != nullptr)
    {
      fastq_print_general(state.fp_fastqout_notmerged_rev,
                          a_read_pair.rev_sequence.data(),
                          static_cast<int>(a_read_pair.rev_length),
                          a_read_pair.rev_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.rev_header.data())),
                          a_read_pair.rev_quality.data(),
                          static_cast<uint64_t>(a_read_pair.rev_abundance),
                          state.notmerged,
                          -1.0,
                          state.parameters);
    }

  if (state.parameters.opt_fastaout_notmerged_fwd != nullptr)
    {
      fasta_print_general(state.fp_fastaout_notmerged_fwd,
                          nullptr,
                          a_read_pair.fwd_sequence.data(),
                          static_cast<int>(a_read_pair.fwd_length),
                          a_read_pair.fwd_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.fwd_header.data())),
                          static_cast<uint64_t>(a_read_pair.fwd_abundance),
                          state.notmerged,
                          -1.0,
                          -1, -1,
                          nullptr, 0.0,
                          0,
                          state.parameters);
    }

  if (state.parameters.opt_fastaout_notmerged_rev != nullptr)
    {
      fasta_print_general(state.fp_fastaout_notmerged_rev,
                          nullptr,
                          a_read_pair.rev_sequence.data(),
                          static_cast<int>(a_read_pair.rev_length),
                          a_read_pair.rev_header.data(),
                          static_cast<int>(std::strlen(a_read_pair.rev_header.data())),
                          static_cast<uint64_t>(a_read_pair.rev_abundance),
                          state.notmerged,
                          -1.0,
                          -1, -1,
                          nullptr, 0.0,
                          0,
                          state.parameters);
    }
}


auto read_pair(struct mergepairs_cli_state_s & state, merge_data_t & a_read_pair) -> bool
{
  auto const fastq_fwd = state.fastq_fwd;
  auto const fastq_rev = state.fastq_rev;

  if (fastq_next(fastq_fwd, false, chrmap_upcase()))
    {
      if (not fastq_next(fastq_rev, false, chrmap_upcase()))
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
      state.progress->update(fastq_get_position(state.fastq_fwd));
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
          if (merge_aborted())
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
      if (merge_aborted())
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
  merge_abort_reset();

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
  if (merge_aborted(std::memory_order_seq_cst))
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

  if (parameters.opt_fastq_minovlen < 5)
    {
      fatal("Overlap specified with --fastq_minovlen must be at least 5");
    }

  /* The short-overlap relaxation of the merge-acceptance thresholds is derived
     inside optimize() from parameters.opt_fastq_minovlen (threaded through the
     worker pool), so nothing to set here. */

  /* open input files */

  fastq_fwd = fastq_open(parameters.opt_fastq_mergepairs, parameters);
  fastq_rev = fastq_open(parameters.opt_reverse, parameters);

  /* open output files */

  OutputFileHandle fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
  fp_fastqout = fastqout_handle.get();
  OutputFileHandle fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  fp_fastaout = fastaout_handle.get();
  OutputFileHandle fastqout_notmerged_fwd_handle = open_optional_output_file(parameters.opt_fastqout_notmerged_fwd, OutputOption{"--fastqout_notmerged_fwd"});
  fp_fastqout_notmerged_fwd = fastqout_notmerged_fwd_handle.get();
  OutputFileHandle fastqout_notmerged_rev_handle = open_optional_output_file(parameters.opt_fastqout_notmerged_rev, OutputOption{"--fastqout_notmerged_rev"});
  fp_fastqout_notmerged_rev = fastqout_notmerged_rev_handle.get();
  OutputFileHandle fastaout_notmerged_fwd_handle = open_optional_output_file(parameters.opt_fastaout_notmerged_fwd, OutputOption{"--fastaout_notmerged_fwd"});
  fp_fastaout_notmerged_fwd = fastaout_notmerged_fwd_handle.get();
  OutputFileHandle fastaout_notmerged_rev_handle = open_optional_output_file(parameters.opt_fastaout_notmerged_rev, OutputOption{"--fastaout_notmerged_rev"});
  fp_fastaout_notmerged_rev = fastaout_notmerged_rev_handle.get();
  OutputFileHandle eetabbedout_handle = open_optional_output_file(parameters.opt_eetabbedout, OutputOption{"--eetabbedout"});
  fp_eetabbedout = eetabbedout_handle.get();

  /* precompute merged quality values */

  precompute_qual(parameters);

  /* main */

  uint64_t const filesize = fastq_get_size(fastq_fwd);
  {
    Progress progress("Merging reads", filesize, parameters);
    state.progress = &progress;

    if (not fastq_fwd->is_empty)
      {
        pair_all(state);
      }
    state.progress = nullptr;  // clear before the Progress it points to is destroyed
  }

  if (fastq_next(fastq_rev, true, chrmap_upcase()))
    {
      fatal("More reverse reads than forward reads");
    }

  if (parameters.fp_log != nullptr) {
    print_stats(state, parameters.fp_log);
  }
  else {
    print_stats(state, stderr);
  }

  /* clean up */

  /* reset() is a no-op on an empty handle, so unopened outputs need
     no guard. */
  eetabbedout_handle.reset();
  fastaout_notmerged_rev_handle.reset();
  fastaout_notmerged_fwd_handle.reset();
  fastqout_notmerged_rev_handle.reset();
  fastqout_notmerged_fwd_handle.reset();
  fastaout_handle.reset();
  fastqout_handle.reset();

  fastq_close(fastq_rev, parameters);
  fastq_rev = nullptr;
  fastq_close(fastq_fwd, parameters);
  fastq_fwd = nullptr;
}
