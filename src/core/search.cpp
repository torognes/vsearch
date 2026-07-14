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
#include <cstdint>  // uint64_t, int64_t
#include <cstring>  // std::strlen, std::strcpy
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
  si->query_head_len = head_len;
  si->qseqlen = seq_len;
  si->query_no = query_no;
  si->qsize = qsize;
  si->strand = strand;

  /* allocate more memory for the sequence, if necessary */

  if (si->qseqlen + 1 > si->seq_alloc)
    {
      si->seq_alloc = si->qseqlen + buffer_headroom;
      si->qsequence = static_cast<char *>(
        xrealloc(si->qsequence, static_cast<size_t>(si->seq_alloc)));
    }

  /* copy the header into owned storage, then point the read-only view at it */
  si->query_head_v.resize(static_cast<std::size_t>(si->query_head_len) + 1);
  std::strcpy(si->query_head_v.data(), head);
  si->query_head = si->query_head_v.data();

  /* copy or reverse-complement sequence */
  if (strand == 0)
    {
      std::strcpy(si->qsequence, seq);
    }
  else
    {
      reverse_complement(si->qsequence, seq, seq_len);
    }
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
  si->uh = unique_init();
  si->kmers = static_cast<count_t *>(xmalloc((static_cast<size_t>(seqcount) * sizeof(count_t)) + 32));
  si->m = minheap_init(tophits);
  si->hits = static_cast<struct hit *>(xmalloc
    (sizeof(struct hit) * static_cast<size_t>(tophits) * static_cast<size_t>(number_of_strands(parameters.opt_strand))));
  si->qsize = 1;
  si->query_head = nullptr;
  si->seq_alloc = 0;
  si->qsequence = nullptr;
  si->s = search16_init(parameters.opt_match,
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
                        parameters.opt_n_mismatch);
}


auto search_thread_exit(struct searchinfo_s * si) -> void
{
  /* thread specific clean up */
  search16_exit(si->s);
  unique_exit(si->uh);
  xfree(si->hits);
  minheap_exit(si->m);
  xfree(si->kmers);
  /* query_head is a view (into query_head_v, the db, or a caller buffer); the
     owned storage query_head_v frees itself. */
  if (si->qsequence != nullptr)
    {
      xfree(si->qsequence);
    }
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
  if (ss->tophits > ss->seqcount)
    {
      ss->tophits = ss->seqcount;
    }

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
  struct searchinfo_s * si = ss->si_plus.get();
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
          dust(strand_si->qsequence, strand_si->qseqlen, parameters);
        }
      else if ((parameters.opt_qmask == Masking::soft) && (parameters.opt_hardmask))
        {
          hardmask(strand_si->qsequence, strand_si->qseqlen);
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
      r.query_length = si->qseqlen;
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
          if (strand_si->hits[i].aligned &&
              strand_si->hits[i].nwalignment != nullptr)
            {
              xfree(strand_si->hits[i].nwalignment);
              strand_si->hits[i].nwalignment = nullptr;
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

  /* per-thread search state arrays (sized to opt_threads) */
  struct searchinfo_s * batch_si_plus;
  struct searchinfo_s * batch_si_minus;  /* nullptr when searching the plus strand only */

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
  struct searchinfo_s * my_si_plus = ctx.batch_si_plus + tid;
  struct searchinfo_s * my_si_minus =
    (ctx.batch_si_minus != nullptr) ? ctx.batch_si_minus + tid : nullptr;
  struct Parameters const & parameters = *ctx.parameters;

  /* grab next query */
  int qi {0};

  auto const has_work_to_claim = [&]() -> bool {
    qi = ctx.next_query++;
    return qi < ctx.query_count;
  };

  auto const process_query = [&]() {
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
            dust(strand_si->qsequence, strand_si->qseqlen, parameters);
          }
        else if ((parameters.opt_qmask == Masking::soft) && (parameters.opt_hardmask))
          {
            hardmask(strand_si->qsequence, strand_si->qseqlen);
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
      ctx.results + qi * ctx.max_results_per_query;
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
            if (strand_si->hits[i].aligned &&
                strand_si->hits[i].nwalignment != nullptr)
              {
                xfree(strand_si->hits[i].nwalignment);
                strand_si->hits[i].nwalignment = nullptr;
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
  if (tophits > seqcount)
    {
      tophits = seqcount;
    }

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

  ctx.batch_si_plus = new searchinfo_s[nthreads]{};
  if (parameters.opt_strand)
    {
      ctx.batch_si_minus = new searchinfo_s[nthreads]{};
    }
  else
    {
      ctx.batch_si_minus = nullptr;
    }

  /* Init per-thread search state before the workers start */
  for (int t = 0; t < nthreads; t++)
    {
      search_thread_init(ctx.batch_si_plus + t, seqcount, tophits, parameters, dbindex, db);
      if (ctx.batch_si_minus != nullptr)
        {
          search_thread_init(ctx.batch_si_minus + t, seqcount, tophits, parameters, dbindex, db);
        }
    }

  /* run all queries through the worker pool (work-stealing on next_query) */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(nthreads),
                              [&ctx](uint64_t tid) {
                                search_batch_worker_fn(ctx, tid);
                              });
    threadrunner.run();
  }

  /* clean up per-thread search state */
  for (int t = 0; t < nthreads; t++)
    {
      search_thread_exit(ctx.batch_si_plus + t);
      if (ctx.batch_si_minus != nullptr)
        {
          search_thread_exit(ctx.batch_si_minus + t);
        }
    }

  delete [] ctx.batch_si_plus;
  if (ctx.batch_si_minus != nullptr)
    {
      delete [] ctx.batch_si_minus;
    }
}
