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

#include "utils/view.hpp"
#include "utils/span.hpp"
#include "vsearch.hpp"
#include "core/buffer_headroom.hpp"
#include "core/db.hpp"  // Database
#include "core/search.hpp"
#include "core/search_internal.hpp"
#include "core/searchcore.hpp"
#include "core/align_simd.hpp"
#include "core/dbindex.hpp"
#include "core/mask.hpp"
#include "core/minheap.hpp"
#include "core/unique.hpp"
#include "utils/make_unique.hpp"
#include "utils/number_of_strands.hpp"
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include "utils/reverse_complement.hpp"
#include <algorithm>  // std::copy_n
#include <cstdint>  // uint64_t, int64_t
#include <cstring>  // std::strlen
#include <memory>  // std::unique_ptr
#include <mutex>  // std::mutex
#include <vector>


auto populate_si(struct searchinfo_s * si,
                 const char * head,
                 int const head_len,
                 const char * seq,
                 int const seq_len,
                 int const query_no,
                 int64_t const qsize,
                 int const strand) -> void
{
  si->query_no = query_no;
  si->qsize = qsize;
  si->strand = strand;

  /* allocate more memory for the sequence, if necessary */

  if (seq_len + 1 > si->seq_alloc)
    {
      si->seq_alloc = seq_len + buffer_headroom;
      si->qsequence_v.resize(static_cast<size_t>(si->seq_alloc));
    }

  /* copy the header into owned storage, then point the read-only view at it */
  si->query_head_v.resize(static_cast<std::size_t>(head_len) + 1);
  std::copy_n(head, static_cast<std::size_t>(head_len) + 1, si->query_head_v.data());
  si->query_head = View<char>{si->query_head_v.data(), static_cast<std::size_t>(head_len)};

  /* copy or reverse-complement sequence into the owned buffer, then point the
     span at it (length == seq_len; the NUL at [seq_len] sits just past the span) */
  if (strand == 0)
    {
      std::copy_n(seq, static_cast<std::size_t>(seq_len) + 1, si->qsequence_v.data());
    }
  else
    {
      reverse_complement(Span<char>{si->qsequence_v.data(), static_cast<std::size_t>(seq_len) + 1}, View<char>{seq, static_cast<std::size_t>(seq_len)});
    }
  si->qsequence = Span<char>{si->qsequence_v.data(), static_cast<std::size_t>(seq_len)};
}


auto search_thread_init(struct searchinfo_s * si, int const seqcount, int const tophits,
                        struct Parameters const & parameters,
                        struct Dbindex const & dbindex,
                        struct Database const & db) -> void
{
  /* thread specific initialiation */
  si->parameters = &parameters;  /* searchcore reads config through the si (E1) */
  si->dbindex = &dbindex;  /* searchcore reads the k-mer index through the si */
  si->db = &db;  /* searchcore reads the sequences through the si */
  /* si->uh (a Uniquer value member) is ready to use as default-constructed */
  /* kmers/hits/qsequence are backed by the searchinfo_s vectors (RAII), so a
     fatal() unwinding out of a partial init or a query frees them; the raw
     pointers below are views into that owned storage. */
  static constexpr auto overflow_padding = 16U;  // 16 * sizeof(count_t) = 32 bytes headroom
  si->kmers_v.reserve(static_cast<size_t>(seqcount) + overflow_padding);
  si->kmers_v.resize(static_cast<size_t>(seqcount));
  si->kmers = si->kmers_v.data();
  si->m = Minheap(tophits);
  si->hits_v.resize(static_cast<size_t>(tophits) * static_cast<size_t>(number_of_strands(parameters.opt_strand)));
  si->hits = si->hits_v.data();
  si->qsize = 1;
  si->query_head = View<char>{nullptr, 0};
  si->seq_alloc = 0;
  si->qsequence = Span<char>{};
  si->s.reset(search16_init(parameters.opt_match,
                        parameters.opt_mismatch,
                        parameters.opt_gap_open_query_left,
                        parameters.opt_gap_open_target_left,
                        parameters.opt_gap_open_query_interior,
                        parameters.opt_gap_open_target_interior,
                        parameters.opt_gap_open_query_right,
                        parameters.opt_gap_open_target_right,
                        parameters.opt_gap_extension_query_left,
                        parameters.opt_gap_extension_target_left,
                        parameters.opt_gap_extension_query_interior,
                        parameters.opt_gap_extension_target_interior,
                        parameters.opt_gap_extension_query_right,
                        parameters.opt_gap_extension_target_right,
                        parameters.opt_n_mismatch));
}


auto search_thread_exit(struct searchinfo_s * si) -> void
{
  /* thread specific clean up. The handles are also freed by ~searchinfo_s if an
     exception unwinds before this runs. */
  si->s.reset();
  si->uh = Uniquer();
  si->m = Minheap();
  /* kmers/hits/qsequence and query_head are views into the searchinfo_s vectors
     (kmers_v/hits_v/qsequence_v/query_head_v), which free their own storage. */
}


/* === Session-based search API (supports both-strand search) === */


struct search_session_s {
  std::unique_ptr<struct searchinfo_s> si_plus;
  std::unique_ptr<struct searchinfo_s> si_minus;  /* non-null when searching both strands */
  int seqcount = 0;  /* number of database sequences */
  int tophits = 0;   /* the maximum number of hits to keep */
  /* run configuration, set at search_session_init and read by the per-query
     path instead of the opt_* globals (E1 shared-infra phase); a pointer so the
     session stays default-constructible. */
  struct Parameters const * parameters = nullptr;
  /* the k-mer index this session searches, supplied by the caller at
     search_session_init and installed on the si's; must outlive the session. */
  struct Dbindex const * dbindex = nullptr;
  /* the sequence database this session searches, supplied by the caller at
     search_session_init and installed on the si's; must outlive the session. */
  struct Database const * db = nullptr;
};


auto search_session_alloc() -> struct search_session_s *
{
  return new search_session_s {};
}


auto search_session_free(struct search_session_s * ss) -> void
{
  if (ss != nullptr)
    {
      search_session_cleanup(ss);
      delete ss;
    }
}


auto search_session_init(struct search_session_s * ss, struct Parameters const & parameters,
                         struct Dbindex const & dbindex,
                         struct Database const & db) -> void
{
  /* Initialize search session for library use.
     Mirrors cluster_session_init: stores seqcount/tophits in the session,
     allocates si_plus (always) and si_minus (when searching both strands). */
  ss->parameters = &parameters;
  ss->dbindex = &dbindex;
  ss->db = &db;
  ss->seqcount = static_cast<int>(db.getsequencecount());
  /* The library path does not clamp to the database size (only the CLI
     search_prep does), so the sizing uses the configured values from
     parameters; si->parameters (set in search_thread_init) carries the same. */
  ss->tophits = static_cast<int>(parameters.opt_maxaccepts + parameters.opt_maxrejects + MAXDELAYED);
  ss->tophits = std::min(ss->tophits, ss->seqcount);

  ss->si_plus = make_unique<searchinfo_s>();
  search_thread_init(ss->si_plus.get(), ss->seqcount, ss->tophits, parameters, *ss->dbindex, *ss->db);
  ss->si_plus->strand = 0;

  if (parameters.opt_strand)
    {
      ss->si_minus = make_unique<searchinfo_s>();
      search_thread_init(ss->si_minus.get(), ss->seqcount, ss->tophits, parameters, *ss->dbindex, *ss->db);
      ss->si_minus->strand = 1;
    }
}


auto search_session_single(struct search_session_s * ss,
                           const char * query_seq,
                           const char * query_head,
                           int const query_len,
                           int64_t const query_size,
                           struct search_result_s * results,
                           int const max_results,
                           int * result_count) -> void
{
  int const head_len = static_cast<int>(std::strlen(query_head));
  struct searchinfo_s const * si = ss->si_plus.get();
  struct Parameters const & parameters = *ss->parameters;

  populate_si(ss->si_plus.get(),
              query_head,
              head_len,
              query_seq,
              query_len,
              0,
              query_size,
              0);

  if (parameters.opt_strand)
    {
      populate_si(ss->si_minus.get(),
                  query_head,
                  head_len,
                  query_seq,
                  query_len,
                  0,
                  query_size,
                  1);
    }

  /* Mask and search each strand independently */
  for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
    {
      struct searchinfo_s * strand_si =
        (s != 0) ? ss->si_minus.get() : ss->si_plus.get();

      if (parameters.opt_qmask == Masking::dust)
        {
          dust(strand_si->qsequence, parameters);
        }
      else if ((parameters.opt_qmask == Masking::soft) && parameters.opt_hardmask)
        {
          hardmask(strand_si->qsequence);
        }

      search_onequery(strand_si, parameters.opt_qmask);
    }

  /* Merge hits from both strands */
  std::vector<struct hit> hits;
  search_joinhits(ss->si_plus.get(),
                  parameters.opt_strand ? ss->si_minus.get() : nullptr,
                  hits);

  /* Populate results (search_joinhits returns only accepted/weak hits) */
  int count = 0;
  for (auto const & h : hits)
    {
      if (count >= max_results)
        {
          break;
        }
      auto & r = results[count];
      r.target = h.target;
      r.id = h.id;
      r.matches = h.matches;
      r.mismatches = h.mismatches;
      r.gaps = h.nwgaps;
      r.alignment_length = h.nwalignmentlength;
      r.query_length = static_cast<int>(si->qsequence.size());
      r.target_length = static_cast<int>(si->db->getsequencelen(static_cast<uint64_t>(h.target)));
      r.accepted = h.accepted;
      r.strand = h.strand;
      ++count;
    }
  *result_count = count;

  /* Free alignment strings directly from si->hits (not the joinhits copy)
     to avoid dangling pointers. Follows cluster_assign_single pattern. */
  for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
    {
      struct searchinfo_s * strand_si =
        (s != 0) ? ss->si_minus.get() : ss->si_plus.get();
      for (int i = 0; i < strand_si->hit_count; ++i)
        {
          if (strand_si->hits[i].aligned)
            {
              strand_si->hits[i].nwalignment.clear();  // std::string; free after use
            }
        }
    }
}


auto search_session_cleanup(struct search_session_s * ss) -> void
{
  if (ss->si_plus)
    {
      search_thread_exit(ss->si_plus.get());
      ss->si_plus.reset();
    }
  if (ss->si_minus)
    {
      search_thread_exit(ss->si_minus.get());
      ss->si_minus.reset();
    }
}


/* === Batch search API === */


/* Shared state for batch search worker threads */
struct search_batch_context_s {
  const char ** query_seqs;
  const char ** query_heads;
  const int * query_lens;
  const int64_t * query_sizes;
  int query_count;
  struct search_result_s * results;
  int max_results_per_query;
  int * result_counts;

  /* per-thread search state arrays (sized to opt_threads). Owned vectors so a
     fatal() during per-thread init unwinds them, running each searchinfo_s
     destructor (freeing partially-initialised handles) and releasing the array. */
  std::vector<struct searchinfo_s> batch_si_plus;
  std::vector<struct searchinfo_s> batch_si_minus;  /* empty when searching the plus strand only */

  /* run configuration, set in search_batch and read by the workers instead of
     the opt_* globals (E1 shared-infra phase). */
  struct Parameters const * parameters;

  /* work-stealing counter */
  std::mutex mutex;
  int next_query;
};


static auto search_batch_worker_fn(struct search_batch_context_s & ctx,
                                   uint64_t tid) -> void
{
  struct searchinfo_s * my_si_plus = &ctx.batch_si_plus[tid];
  struct searchinfo_s * my_si_minus =
    (not ctx.batch_si_minus.empty()) ? &ctx.batch_si_minus[tid] : nullptr;
  struct Parameters const & parameters = *ctx.parameters;

  /* grab next query */
  int qi {0};

  auto const has_work_to_claim = [&]() -> bool {
    qi = ctx.next_query++;
    return qi < ctx.query_count;
  };

  auto const process_query = [&]() -> void {
    char const * qseq = ctx.query_seqs[qi];
    char const * qhead = ctx.query_heads[qi];
    int const qlen = ctx.query_lens[qi];
    int64_t const qsize = ctx.query_sizes[qi];
    int const head_len = static_cast<int>(std::strlen(qhead));

    populate_si(my_si_plus,
                qhead,
                head_len,
                qseq,
                qlen,
                qi,
                qsize,
                0);

    if (my_si_minus != nullptr)
      {
        populate_si(my_si_minus,
                    qhead,
                    head_len,
                    qseq,
                    qlen,
                    qi,
                    qsize,
                    1);
      }

    /* Mask and search each strand independently */
    for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
      {
        struct searchinfo_s * strand_si =
          (s != 0) ? my_si_minus : my_si_plus;

        if (parameters.opt_qmask == Masking::dust)
          {
            dust(strand_si->qsequence, parameters);
          }
        else if ((parameters.opt_qmask == Masking::soft) && parameters.opt_hardmask)
          {
            hardmask(strand_si->qsequence);
          }

        search_onequery(strand_si, parameters.opt_qmask);
      }

    /* Merge hits from both strands */
    std::vector<struct hit> hits;
    search_joinhits(my_si_plus,
                    parameters.opt_strand ? my_si_minus : nullptr,
                    hits);

    /* Populate results for this query */
    struct search_result_s * qresults =
      ctx.results + (static_cast<std::ptrdiff_t>(qi) * ctx.max_results_per_query);
    int count = 0;
    for (auto const & h : hits)
      {
        if (count >= ctx.max_results_per_query)
          {
            break;
          }
        auto & r = qresults[count];
        r.target = h.target;
        r.id = h.id;
        r.matches = h.matches;
        r.mismatches = h.mismatches;
        r.gaps = h.nwgaps;
        r.alignment_length = h.nwalignmentlength;
        r.query_length = qlen;
        r.target_length = static_cast<int>(my_si_plus->db->getsequencelen(static_cast<uint64_t>(h.target)));
        r.accepted = h.accepted;
        r.strand = h.strand;
        ++count;
      }
    ctx.result_counts[qi] = count;

    /* Free alignment strings from si->hits directly */
    for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
      {
        struct searchinfo_s * strand_si =
          (s != 0) ? my_si_minus : my_si_plus;
        for (int i = 0; i < strand_si->hit_count; ++i)
          {
            if (strand_si->hits[i].aligned)
              {
                strand_si->hits[i].nwalignment.clear();  // std::string; free after use
              }
          }
      }
  };

  run_worker_loop(ctx.mutex, has_work_to_claim, process_query);
}


auto search_batch(struct Parameters const & parameters,
                  struct Dbindex const & dbindex,
                  struct Database const & db,
                  const char ** query_seqs,
                  const char ** query_heads,
                  const int * query_lens,
                  const int64_t * query_sizes,
                  int const query_count,
                  struct search_result_s * results,
                  int const max_results_per_query,
                  int * result_counts) -> void
{
  /* per-thread buffer sizes for search_thread_init (formerly file-statics).
     The library path does not clamp to the database size (only the CLI
     search_prep does), so the sizing uses the configured values from
     parameters; si->parameters (set in search_thread_init) carries the same. */
  int const seqcount = static_cast<int>(db.getsequencecount());
  int tophits = static_cast<int>(parameters.opt_maxaccepts + parameters.opt_maxrejects + MAXDELAYED);
  tophits = std::min(tophits, seqcount);

  int const nthreads = static_cast<int>(parameters.opt_threads);

  /* Allocate per-thread search state */
  struct search_batch_context_s ctx;
  ctx.query_seqs = query_seqs;
  ctx.query_heads = query_heads;
  ctx.query_lens = query_lens;
  ctx.query_sizes = query_sizes;
  ctx.query_count = query_count;
  ctx.results = results;
  ctx.max_results_per_query = max_results_per_query;
  ctx.result_counts = result_counts;
  ctx.parameters = &parameters;
  ctx.next_query = 0;

  ctx.batch_si_plus.resize(static_cast<std::size_t>(nthreads));
  if (parameters.opt_strand)
    {
      ctx.batch_si_minus.resize(static_cast<std::size_t>(nthreads));
    }

  /* Init per-thread search state before the workers start */
  for (int t = 0; t < nthreads; t++)
    {
      search_thread_init(&ctx.batch_si_plus[static_cast<std::size_t>(t)], seqcount, tophits, parameters, dbindex, db);
      if (not ctx.batch_si_minus.empty())
        {
          search_thread_init(&ctx.batch_si_minus[static_cast<std::size_t>(t)], seqcount, tophits, parameters, dbindex, db);
        }
    }

  /* run all queries through the worker pool (work-stealing on next_query) */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(nthreads),
                              [&ctx](uint64_t tid) -> void {
                                search_batch_worker_fn(ctx, tid);
                              });
    threadrunner.run();
  }

  /* clean up per-thread search state (the vectors also free themselves, and
     would run these searchinfo_s destructors on an exception unwind). */
  for (int t = 0; t < nthreads; t++)
    {
      search_thread_exit(&ctx.batch_si_plus[static_cast<std::size_t>(t)]);
      if (not ctx.batch_si_minus.empty())
        {
          search_thread_exit(&ctx.batch_si_minus[static_cast<std::size_t>(t)]);
        }
    }
}
