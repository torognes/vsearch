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
#include "utils/view.hpp"
#include "vsearch.hpp"
#include <memory>  // std::unique_ptr
#include "core/db.hpp"
#include "core/dbhash.hpp"
#include "core/fasta.hpp"
#include "core/fastx.hpp"
#include "core/results.hpp"
#include "core/searchcore.hpp"
#include "core/buffer_headroom.hpp"
#include "utils/progress.hpp"
#include "core/mask.hpp"
#include "core/otutable.hpp"
#include "utils/fatal.hpp"
#include "utils/fatal_allocator.hpp"  // FatalAllocator
#include "utils/maps.hpp"
#include "utils/number_of_strands.hpp"
#include "utils/open_file.hpp"
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/string_normalize.hpp"
#include <algorithm>  // std::min
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t
#include <cstring>  // std::strlen
#include <mutex>  // std::mutex, std::lock_guard, std::unique_lock
#include <string>  // std::string, std::to_string
#include <vector>


/* Per-invocation state for the CLI-only --search_exact command. Folds what
   used to be file-static module state into one struct owned as a local in
   search_exact() and threaded through prep / the worker pool / done, so the
   command holds no mutable module-scope state (Band 7 E4). */
struct search_exact_state_s
{
  /* the run configuration, threaded through the helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter */
  struct Parameters const & parameters;
  struct Database db;  /* the sequence database this run owns (RAII); si->db points here */
  struct Dbhash dbhash;  /* the exact-match hash index this run owns (RAII) */

  struct searchinfo_s * si_plus = nullptr;
  struct searchinfo_s * si_minus = nullptr;

  /* set once before the worker pool runs, then read-only: no synchronization needed */
  int tophits = 0; /* the maximum number of hits to keep */
  int seqcount = 0; /* number of database sequences */
  fastx_handle query_fastx_h = nullptr;

  /* accessed by the worker threads; access serialized by the mutexes */
  std::mutex mutex_input;
  std::mutex mutex_output;
  int qmatches = 0;
  uint64_t qmatches_abundance = 0;
  int queries = 0;
  uint64_t queries_abundance = 0;
  std::vector<uint64_t, FatalAllocator<uint64_t>> dbmatched;
  FILE * fp_samout = nullptr;
  FILE * fp_alnout = nullptr;
  FILE * fp_userout = nullptr;
  FILE * fp_blast6out = nullptr;
  FILE * fp_uc = nullptr;
  FILE * fp_fastapairs = nullptr;
  FILE * fp_matched = nullptr;
  FILE * fp_notmatched = nullptr;
  FILE * fp_dbmatched = nullptr;
  FILE * fp_dbnotmatched = nullptr;
  FILE * fp_otutabout = nullptr;
  FILE * fp_mothur_shared_out = nullptr;
  FILE * fp_biomout = nullptr;
  /* accumulates the OTU table for the outputs above; the workers mutate it via
     add() under mutex_output, replacing the former file-static singleton */
  OtuTable otutable;
  FILE * fp_qsegout = nullptr;
  FILE * fp_tsegout = nullptr;

  int count_matched = 0;
  int count_notmatched = 0;

  Progress * progress = nullptr;  /* the owner's progress bar; worker updates it under mutex_output */

  explicit search_exact_state_s(struct Parameters const & params) : parameters(params) {}
};


namespace {
auto add_hit(struct searchinfo_s * si, uint64_t const seqno) -> void
{
  if (search_acceptable_unaligned(*si, static_cast<int>(seqno)))
    {
      struct hit * hp = si->hits + si->hit_count;
      si->hit_count++;

      hp->target = static_cast<int>(seqno);
      hp->strand = si->strand;

      hp->count = 0;

      auto const qseqlen = static_cast<int>(si->qsequence.size());
      hp->nwscore = static_cast<int>(static_cast<int64_t>(qseqlen) * si->parameters->opt_match);
      hp->nwdiff = 0;
      hp->nwgaps = 0;
      hp->nwindels = 0;
      hp->nwalignmentlength = qseqlen;
      hp->nwid = 100.0;
      hp->matches = qseqlen;
      hp->mismatches = 0;

      hp->nwalignment = std::to_string(qseqlen) + "M";

      hp->internal_alignmentlength = qseqlen;
      hp->internal_gaps = 0;
      hp->internal_indels = 0;
      hp->trim_q_left = 0;
      hp->trim_q_right = 0;
      hp->trim_t_left = 0;
      hp->trim_t_right = 0;
      hp->trim_aln_left = 0;
      hp->trim_aln_right = 0;

      hp->id = 100.0;
      hp->id0 = 100.0;
      hp->id1 = 100.0;
      hp->id2 = 100.0;
      hp->id3 = 100.0;
      hp->id4 = 100.0;

      hp->shortest = qseqlen;
      hp->longest = qseqlen;

      hp->aligned = true;

      hp->accepted = false;
      hp->rejected = false;
      hp->weak = false;
      (void) search_acceptable_aligned(*si, hp);
    }
}

auto search_exact_onequery(struct searchinfo_s * si, struct Dbhash const & dbhash) -> void
{
  dbhash_search_info_s info;

  char const * seq = si->qsequence.data();
  uint64_t const seqlen = si->qsequence.size();
  std::vector<char> normalized(seqlen + 1);
  string_normalize(Span<char>{normalized.data(), static_cast<std::size_t>(seqlen) + 1}, View<char>{seq, static_cast<std::size_t>(seqlen)});

  si->hit_count = 0;

  int64_t ret = dbhash.search_first(normalized.data(), seqlen, & info, *si->db);
  while (ret >= 0)
    {
      add_hit(si, static_cast<uint64_t>(ret));
      ret = dbhash.search_next(&info, *si->db);
    }
}

auto search_exact_output_results(struct search_exact_state_s & state,
                                 std::vector<struct hit> const & hits,
                                 View<char> const query_head_view,
                                 int const qseqlen,
                                 char const * qsequence,
                                 char const * qsequence_rc,
                                 int64_t const qsize) -> void
{
  struct Parameters const & parameters = state.parameters;
  std::lock_guard<std::mutex> const lock(state.mutex_output);
  auto const * const query_head = query_head_view.data();

  /* show results */
  auto const n_results_to_report = std::min(parameters.opt_maxhits, static_cast<int64_t>(hits.size()));

  if (state.fp_alnout != nullptr)
    {
      results_show_alnout(state.fp_alnout,
                          hits.data(),
                          static_cast<int>(n_results_to_report),
                          query_head,
                          qsequence,
                          qseqlen,
                          state.db,
                          parameters);
    }

  if (state.fp_samout != nullptr)
    {
      results_show_samout(state.fp_samout,
                          hits.data(),
                          static_cast<int>(n_results_to_report),
                          query_head,
                          qsequence,
                          qsequence_rc,
                          state.db,
                          parameters);
    }

  if (n_results_to_report != 0)
    {
      double const top_hit_id = hits[0].id;

      if ((parameters.opt_otutabout != nullptr) || (parameters.opt_mothur_shared_out != nullptr) || (parameters.opt_biomout != nullptr))
        {
          state.otutable.add(query_head_view,
                       state.db.header_view(static_cast<uint64_t>(hits[0].target)),
                       qsize);
        }

      for (int64_t t = 0; t < n_results_to_report; t++)
        {
          auto const & hit = hits[static_cast<std::size_t>(t)];

          if ((parameters.opt_top_hits_only != 0) && (hit.id < top_hit_id))
            {
              break;
            }

          if (state.fp_fastapairs != nullptr)
            {
              results_show_fastapairs_one(state.fp_fastapairs,
                                          &hit,
                                          query_head,
                                          qsequence,
                                          qsequence_rc,
                                          state.db,
                                          parameters);
            }

          if (state.fp_qsegout != nullptr)
            {
              results_show_qsegout_one(state.fp_qsegout,
                                       &hit,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       parameters);
            }

          if (state.fp_tsegout != nullptr)
            {
              results_show_tsegout_one(state.fp_tsegout,
                                       &hit,
                                       state.db,
                                       parameters);
            }

          if (state.fp_uc != nullptr)
            {
              if ((t == 0) || (parameters.opt_uc_allhits))
                {
                  results_show_uc_one(state.fp_uc,
                                      &hit,
                                      query_head,
                                      qseqlen,
                                      hit.target,
                                      state.db,
                                      parameters);
                }
            }

          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout,
                                       &hit,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.db,
                                       parameters);
            }

          if (state.fp_blast6out != nullptr)
            {
              results_show_blast6out_one(state.fp_blast6out,
                                         &hit,
                                         query_head,
                                         qseqlen,
                                         state.db);
            }
        }
    }
  else
    {
      if ((parameters.opt_otutabout != nullptr) || (parameters.opt_mothur_shared_out != nullptr) || (parameters.opt_biomout != nullptr))
        {
          state.otutable.add(query_head_view,
                       View<char>{},
                       qsize);
        }

      if (state.fp_uc != nullptr)
        {
          results_show_uc_one(state.fp_uc,
                              nullptr,
                              query_head,
                              qseqlen,
                              0,
                              state.db,
                              parameters);
        }

      if (parameters.opt_output_no_hits != 0)
        {
          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout,
                                       nullptr,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.db,
                                       parameters);
            }

          if (state.fp_blast6out != nullptr)
            {
              results_show_blast6out_one(state.fp_blast6out,
                                         nullptr,
                                         query_head,
                                         qseqlen,
                                         state.db);
            }
        }
    }

  if (not hits.empty())
    {
      ++state.count_matched;
      if (parameters.opt_matched != nullptr)
        {
          fasta_print_general(state.fp_matched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              static_cast<int>(std::strlen(query_head)),
                              static_cast<uint64_t>(qsize),
                              state.count_matched,
                              -1.0,
                              -1, -1, nullptr, 0.0,
                              0,
                              parameters);
        }
    }
  else
    {
      ++state.count_notmatched;
      if (parameters.opt_notmatched != nullptr)
        {
          fasta_print_general(state.fp_notmatched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              static_cast<int>(std::strlen(query_head)),
                              static_cast<uint64_t>(qsize),
                              state.count_notmatched,
                              -1.0,
                              -1, -1, nullptr, 0.0,
                              0,
                              parameters);
        }
    }

  /* update matching db sequences */
    for (auto const & hit : hits) {
      if (hit.accepted) {
        state.dbmatched[static_cast<std::size_t>(hit.target)] += static_cast<uint64_t>(parameters.opt_sizein ? qsize : 1);
      }
    }
}

auto search_exact_query(uint64_t const t, struct search_exact_state_s & state) -> int
{
  struct Parameters const & parameters = state.parameters;
  for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
    {
      struct searchinfo_s * si = (s != 0) ? state.si_minus + t : state.si_plus + t;

      /* mask query */
      if (parameters.opt_qmask == Masking::dust)
        {
          dust(si->qsequence, parameters);
        }
      else if ((parameters.opt_qmask == Masking::soft) && (parameters.opt_hardmask))
        {
          hardmask(si->qsequence);
        }

      /* perform search */
      search_exact_onequery(si, state.dbhash);
    }

  std::vector<struct hit> hits;

  search_joinhits(state.si_plus + t,
                  parameters.opt_strand ? state.si_minus + t : nullptr,
                  hits);

  search_exact_output_results(state,
                              hits,
                              state.si_plus[t].query_head,
                              static_cast<int>(state.si_plus[t].qsequence.size()),
                              state.si_plus[t].qsequence.data(),
                              parameters.opt_strand ? state.si_minus[t].qsequence.data() : nullptr,
                              state.si_plus[t].qsize);

  /* alignment strings (hit.nwalignment) are std::string and free themselves */

  return static_cast<int>(hits.size());
}

auto search_exact_thread_run(uint64_t const t, struct search_exact_state_s & state) -> void
{
  struct Parameters const & parameters = state.parameters;
  int64_t qsize = 0;
  uint64_t progress = 0;

  auto const has_work_to_claim = [&]() -> bool {
    if (not state.query_fastx_h->next((not parameters.opt_notrunclabels), chrmap_no_change()))
      {
        return false;
      }

    char const * qhead = state.query_fastx_h->get_header();
    int const query_head_len = static_cast<int>(state.query_fastx_h->get_header_length());
    char const * qseq = state.query_fastx_h->get_sequence();
    int const qseqlen = static_cast<int>(state.query_fastx_h->get_sequence_length());
    int const query_no = static_cast<int>(state.query_fastx_h->get_seqno());
    qsize = state.query_fastx_h->get_abundance();

    for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
      {
        struct searchinfo_s * si = (s != 0) ? state.si_minus + t : state.si_plus + t;

        si->query_no = query_no;
        si->qsize = qsize;
        si->strand = s;

        /* allocate more memory for the sequence, if necessary */

        if (qseqlen + 1 > si->seq_alloc)
          {
            si->seq_alloc = qseqlen + buffer_headroom;
            si->qsequence_v.resize(static_cast<size_t>(si->seq_alloc));
          }
      }

    /* plus strand: copy header and sequence into owned storage, spans point at them */
    state.si_plus[t].query_head_v.resize(static_cast<std::size_t>(query_head_len) + 1);
    std::copy_n(qhead, static_cast<std::size_t>(query_head_len) + 1, state.si_plus[t].query_head_v.data());
    state.si_plus[t].query_head = View<char>{state.si_plus[t].query_head_v.data(), static_cast<std::size_t>(query_head_len)};
    std::copy_n(qseq, static_cast<std::size_t>(qseqlen) + 1, state.si_plus[t].qsequence_v.data());
    state.si_plus[t].qsequence = Span<char>{state.si_plus[t].qsequence_v.data(), static_cast<std::size_t>(qseqlen)};

    /* get progress as amount of input file read */
    progress = state.query_fastx_h->get_position();
    return true;
  };

  auto const process_query = [&]() {
    /* minus strand: copy header and reverse complementary sequence */
    if (parameters.opt_strand)
      {
        state.si_minus[t].query_head_v = state.si_plus[t].query_head_v;
        state.si_minus[t].query_head = View<char>{state.si_minus[t].query_head_v.data(), state.si_plus[t].query_head.size()};
        reverse_complement(Span<char>{state.si_minus[t].qsequence_v.data(), state.si_plus[t].qsequence.size() + 1},
                           View<char>{state.si_plus[t].qsequence.data(), state.si_plus[t].qsequence.size()});
        state.si_minus[t].qsequence = Span<char>{state.si_minus[t].qsequence_v.data(), state.si_plus[t].qsequence.size()};
      }

    int const match = search_exact_query(t, state);

    /* lock mutex for update of global data and output */
    std::lock_guard<std::mutex> const output_lock(state.mutex_output);

    /* update stats */
    state.queries++;
    state.queries_abundance += static_cast<uint64_t>(qsize);

    if (match != 0)
      {
        state.qmatches++;
        state.qmatches_abundance += static_cast<uint64_t>(qsize);
      }

    /* show progress */
    state.progress->update(progress);
  };

  run_worker_loop(state.mutex_input, has_work_to_claim, process_query);
}

auto search_exact_thread_init(struct searchinfo_s * si, struct Parameters const & parameters, int const tophits) -> void
{
  /* thread specific initialiation */
  si->parameters = &parameters;  /* searchcore reads config through the si (E1) */
  si->uh = Uniquer();
  si->kmers = nullptr;
  si->m = Minheap();
  si->hits_v.resize(static_cast<std::size_t>(tophits * number_of_strands(parameters.opt_strand)));
  si->hits = si->hits_v.data();
  si->qsize = 1;
  si->query_head = View<char>{nullptr, 0};
  si->seq_alloc = 0;
  si->qsequence = Span<char>{};
  si->nw = nullptr;
  si->s = nullptr;
}

auto search_exact_thread_exit(struct searchinfo_s * /* si */) -> void
{
  /* thread specific clean up: query_head/qsequence are views and the hit/kmer
     buffers are RAII vectors that free themselves, so nothing to release here.
     Kept for lifecycle symmetry with search_exact_thread_init. */
}

auto search_exact_thread_worker_run(struct search_exact_state_s & state) -> void
{
  struct Parameters const & parameters = state.parameters;

  /* search_exact forces 100% identity: thread a copy with opt_id = 1.0 through
     si->parameters so the shared search_acceptable_aligned accepts only exact
     matches. The forcing lives here, not in the command dispatcher, so the
     "100% identity" semantics stay local to the command that needs them (E1).
     effective outlives the worker pool below, so si->parameters stays valid for
     the duration of the run. */
  struct Parameters effective = parameters;
  effective.opt_id = 1.0;

  /* init per-thread search state before the workers start */
  for (int t = 0; t < parameters.opt_threads; t++)
    {
      search_exact_thread_init(state.si_plus + t, effective, state.tophits);
      (state.si_plus + t)->db = &state.db;  /* searchcore reads the sequences through the si */
      if (state.si_minus != nullptr)
        {
          search_exact_thread_init(state.si_minus + t, effective, state.tophits);
          (state.si_minus + t)->db = &state.db;
        }
    }

  /* run the worker pool over the input file */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(parameters.opt_threads),
                              [&state](uint64_t const t)
                              { search_exact_thread_run(t, state); });
    threadrunner.run();
  }

  /* clean up per-thread search state */
  for (int t = 0; t < parameters.opt_threads; t++)
    {
      search_exact_thread_exit(state.si_plus + t);
      if (state.si_minus != nullptr)
        {
          search_exact_thread_exit(state.si_minus + t);
        }
    }
}

auto search_exact_prep(struct search_exact_state_s & state) -> void
{
  struct Parameters const & parameters = state.parameters;

  /* the output files were opened by the caller (search_exact), which owns
     the RAII handles for the duration of the run */

  state.db.read(parameters.opt_db, 0, parameters);

  results_show_samheader(state.fp_samout, parameters.opt_db, state.db, parameters);

  if (parameters.opt_dbmask == Masking::dust)
    {
      dust_all(state.db, parameters);
    }
  else if ((parameters.opt_dbmask == Masking::soft) && (parameters.opt_hardmask))
    {
      hardmask_all(state.db);
    }

  // memory-intensive: the entire database is now held in memory

  state.seqcount = static_cast<int>(state.db.getsequencecount());

  /* tophits = the maximum number of hits we need to store */
  state.tophits = state.seqcount;

  state.dbmatched.assign(static_cast<size_t>(state.seqcount), 0);

  state.dbhash.open(static_cast<uint64_t>(state.seqcount));
  state.dbhash.add_all(state.db, parameters);
}

auto search_exact_done(struct search_exact_state_s & state) -> void
{
  /* clean up; the output files are closed by search_exact() through the RAII
     handles it owns */
  state.dbhash.clear();

  state.db.clear();
  state.dbmatched.clear();
  state.dbmatched.shrink_to_fit();
}
}  // anonymous namespace


auto search_exact(struct Parameters const & parameters) -> void
{
  search_exact_state_s state(parameters);

  /* open output files; the handles are owned here so they outlive the worker
     pool, which reads the non-owning state.fp_* under the output lock */
  OutputFileHandle alnout_handle = open_optional_output_file(parameters.opt_alnout, OutputOption{"--alnout"});
  state.fp_alnout = alnout_handle.get();
  if (state.fp_alnout != nullptr)
    {
      std::fprintf(state.fp_alnout, "%s\n", parameters.command_line.c_str());
      std::fprintf(state.fp_alnout, "%s\n", parameters.prog_header.c_str());
    }
  OutputFileHandle samout_handle = open_optional_output_file(parameters.opt_samout, OutputOption{"--samout"});
  state.fp_samout = samout_handle.get();
  OutputFileHandle userout_handle = open_optional_output_file(parameters.opt_userout, OutputOption{"--userout"});
  state.fp_userout = userout_handle.get();
  OutputFileHandle blast6out_handle = open_optional_output_file(parameters.opt_blast6out, OutputOption{"--blast6out"});
  state.fp_blast6out = blast6out_handle.get();
  OutputFileHandle uc_handle = open_optional_output_file(parameters.opt_uc, OutputOption{"--uc"});
  state.fp_uc = uc_handle.get();
  OutputFileHandle fastapairs_handle = open_optional_output_file(parameters.opt_fastapairs, OutputOption{"--fastapairs"});
  state.fp_fastapairs = fastapairs_handle.get();
  OutputFileHandle qsegout_handle = open_optional_output_file(parameters.opt_qsegout, OutputOption{"--qsegout"});
  state.fp_qsegout = qsegout_handle.get();
  OutputFileHandle tsegout_handle = open_optional_output_file(parameters.opt_tsegout, OutputOption{"--tsegout"});
  state.fp_tsegout = tsegout_handle.get();
  OutputFileHandle matched_handle = open_optional_output_file(parameters.opt_matched, OutputOption{"--matched"});
  state.fp_matched = matched_handle.get();
  OutputFileHandle notmatched_handle = open_optional_output_file(parameters.opt_notmatched, OutputOption{"--notmatched"});
  state.fp_notmatched = notmatched_handle.get();
  OutputFileHandle dbmatched_handle = open_optional_output_file(parameters.opt_dbmatched, OutputOption{"--dbmatched"});
  state.fp_dbmatched = dbmatched_handle.get();
  OutputFileHandle dbnotmatched_handle = open_optional_output_file(parameters.opt_dbnotmatched, OutputOption{"--dbnotmatched"});
  state.fp_dbnotmatched = dbnotmatched_handle.get();
  OutputFileHandle otutabout_handle = open_optional_output_file(parameters.opt_otutabout, OutputOption{"--otutabout"});
  state.fp_otutabout = otutabout_handle.get();
  OutputFileHandle mothur_shared_out_handle = open_optional_output_file(parameters.opt_mothur_shared_out, OutputOption{"--mothur_shared_out"});
  state.fp_mothur_shared_out = mothur_shared_out_handle.get();
  OutputFileHandle biomout_handle = open_optional_output_file(parameters.opt_biomout, OutputOption{"--biomout"});
  state.fp_biomout = biomout_handle.get();

  search_exact_prep(state);

  /* prepare reading of queries */
  state.qmatches = 0;
  state.qmatches_abundance = 0;
  state.queries = 0;
  state.queries_abundance = 0;
  auto const query_owner = fastx_open(parameters.opt_search_exact, parameters);
  state.query_fastx_h = query_owner.get();  // workers borrow the raw handle

  /* The query file is parsed inside the worker threads
     (search_exact_thread_run). Defer parse errors so a malformed query
     stops the pool cooperatively instead of calling fatal()/std::exit()
     from a worker while siblings are writing output (CC3); reported below
     from the main thread after join. */
  state.query_fastx_h->enable_deferred_errors();

  /* allocate memory for thread info */
  std::vector<struct searchinfo_s> si_plus_v(static_cast<std::size_t>(parameters.opt_threads));
  state.si_plus = si_plus_v.data();
  std::vector<struct searchinfo_s> si_minus_v;
  if (parameters.opt_strand)
    {
      si_minus_v.resize(static_cast<std::size_t>(parameters.opt_threads));
      state.si_minus = si_minus_v.data();
    }

  {
    Progress progress_bar("Searching", state.query_fastx_h->get_size(), parameters);
    state.progress = &progress_bar;
    search_exact_thread_worker_run(state);
  }

  /* all workers joined; report a deferred query parse error (CC3) from the
     main thread so it does not race a worker's output */
  if (state.query_fastx_h->get_error())
    {
      fatal("%s", state.query_fastx_h->get_errmsg());
    }

  // si_plus not used below that point
  // si_minus not used below that point


  query_owner->report_stripped_warning(parameters);  // query_owner (unique_ptr) frees the reader on return

  if (! parameters.opt_quiet)
    {
      std::fprintf(stderr, "Matching unique query sequences: %d of %d",
              state.qmatches, state.queries);
      if (state.queries > 0)
        {
          std::fprintf(stderr, " (%.2f%%)", 100.0 * state.qmatches / state.queries);
        }
      std::fprintf(stderr, "\n");
      if (parameters.opt_sizein)
        {
          std::fprintf(stderr, "Matching total query sequences: %" PRIu64 " of %"
                  PRIu64,
                  state.qmatches_abundance, state.queries_abundance);
          if (state.queries_abundance > 0)
            {
              std::fprintf(stderr, " (%.2f%%)",
                      100.0 * static_cast<double>(state.qmatches_abundance) / static_cast<double>(state.queries_abundance));
            }
          std::fprintf(stderr, "\n");
        }
    }

  if (parameters.opt_log != nullptr)
    {
      std::fprintf(parameters.fp_log, "Matching unique query sequences: %d of %d",
              state.qmatches, state.queries);
      if (state.queries > 0)
        {
          std::fprintf(parameters.fp_log, " (%.2f%%)", 100.0 * state.qmatches / state.queries);
        }
      std::fprintf(parameters.fp_log, "\n");
      if (parameters.opt_sizein)
        {
          std::fprintf(parameters.fp_log, "Matching total query sequences: %" PRIu64 " of %"
                  PRIu64,
                  state.qmatches_abundance, state.queries_abundance);
          if (state.queries_abundance > 0)
            {
              std::fprintf(parameters.fp_log, " (%.2f%%)",
                      100.0 * static_cast<double>(state.qmatches_abundance) / static_cast<double>(state.queries_abundance));
            }
          std::fprintf(parameters.fp_log, "\n");
        }
    }

  // Add OTUs with no matches to OTU table
  if ((parameters.opt_otutabout != nullptr) || (parameters.opt_mothur_shared_out != nullptr) || (parameters.opt_biomout != nullptr)) {
    for (int64_t i = 0; i < state.seqcount; i++) {
      if (state.dbmatched[static_cast<std::size_t>(i)] == 0U) {
        state.otutable.add(View<char>{}, state.db.header_view(static_cast<uint64_t>(i)), 0);
      }
    }
  }

  if (state.fp_biomout != nullptr)
    {
      state.otutable.print_biomout(state.fp_biomout, parameters);
      biomout_handle.reset();
    }

  if (state.fp_otutabout != nullptr)
    {
      state.otutable.print_otutabout(state.fp_otutabout, parameters);
      otutabout_handle.reset();
    }

  if (state.fp_mothur_shared_out != nullptr)
    {
      state.otutable.print_mothur_shared_out(state.fp_mothur_shared_out, parameters);
      mothur_shared_out_handle.reset();
    }

  if ((parameters.opt_dbmatched != nullptr) || (parameters.opt_dbnotmatched != nullptr))
    {
      int count_dbmatched = 0;
      int count_dbnotmatched = 0;

      for (int64_t i = 0; i < state.seqcount; i++)
        {
          if (state.dbmatched[static_cast<std::size_t>(i)] != 0U)
            {
              ++count_dbmatched;
              if (parameters.opt_dbmatched != nullptr)
                {
                  fasta_print_general(state.fp_dbmatched,
                                      nullptr,
                                      state.db.record(static_cast<uint64_t>(i)),
                                      state.dbmatched[static_cast<std::size_t>(i)],
                                      count_dbmatched,
                                      -1.0,
                                      -1, -1, nullptr, 0.0,
                                      0,
                                      parameters);
                }
            }
          else
            {
              ++count_dbnotmatched;
              if (parameters.opt_dbnotmatched != nullptr)
                {
                  fasta_print_general(state.fp_dbnotmatched,
                                      nullptr,
                                      state.db.record(static_cast<uint64_t>(i)),
                                      0,
                                      count_dbnotmatched,
                                      -1.0,
                                      -1, -1, nullptr, 0.0,
                                      0,
                                      parameters);
                }
            }
        }
    }

  search_exact_done(state);

  /* reset() is a no-op on an empty handle, so unopened outputs need no guard.
     The fixed order matches the legacy fclose sequence: RAII scope-exit would
     reverse it and flip the flush order for outputs that share stdout. */
  dbmatched_handle.reset();
  dbnotmatched_handle.reset();
  matched_handle.reset();
  notmatched_handle.reset();
  fastapairs_handle.reset();
  qsegout_handle.reset();
  tsegout_handle.reset();
  uc_handle.reset();
  blast6out_handle.reset();
  userout_handle.reset();
  alnout_handle.reset();
  samout_handle.reset();
}
