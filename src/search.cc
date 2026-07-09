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
#include "search.h"
#include "searchcore.h"
#include "align_simd.h"
#include "dbindex.h"
#include "mask.h"
#include "minheap.h"
#include "otutable.h"
#include "udb.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/make_unique.hpp"
#include "utils/maps.hpp"
#include "utils/number_of_strands.hpp"
#include "utils/threads.hpp"
#include <algorithm>  // std::min
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // uint64_t, int64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t
#include <cstring>  // std::strlen, std::memset, std::strcpy
#include <mutex>  // std::mutex, std::lock_guard
#include <vector>


/* Per-invocation state for a usearch_global run — previously the file-static
   globals below: the per-thread searchinfo arrays, the query file handle, the
   two mutexes, the match/abundance counters, the per-db-sequence match tally
   and the sixteen output handles. Folding them into a struct that
   usearch_global() owns and threads through the output helper and the
   streaming worker pool makes the command reentrant and removes the shared
   mutable state (E4). The library session/batch paths own their own searchinfo
   arrays and return results rather than writing files, so they do not use this
   struct; the shared per-thread sizes (seqcount/tophits) are passed to
   search_thread_init as parameters instead. */
struct search_cli_state_s
{
  /* the run configuration, threaded through the CLI-path helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter. The
     helpers shared with the library path (search_thread_init, populate_si, the
     searchcore) keep reading the globals — the library session/batch entries
     have no Parameters. */
  struct Parameters const & parameters;
  /* a copy of parameters with opt_maxaccepts/opt_maxrejects clamped to the
     database size (search_prep); si->parameters points here so the shared
     searchcore reads the clamped values without a mutated global (E1). */
  struct Parameters effective_parameters;
  struct Dbindex dbindex;  /* the k-mer index this run owns (RAII); si->dbindex points here */
  int tophits = 0;   /* the maximum number of hits to keep */
  int seqcount = 0;  /* number of database sequences */
  struct searchinfo_s * si_plus = nullptr;
  struct searchinfo_s * si_minus = nullptr;
  fastx_handle query_fastx_h = nullptr;
  std::mutex mutex_input;   /* serializes query reads */
  std::mutex mutex_output;  /* serializes output + counter updates */
  int qmatches = 0;
  uint64_t qmatches_abundance = 0;
  int queries = 0;
  uint64_t queries_abundance = 0;
  uint64_t * dbmatched = nullptr;
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
  FILE * fp_lcaout = nullptr;
  FILE * fp_qsegout = nullptr;
  FILE * fp_tsegout = nullptr;
  int count_matched = 0;
  int count_notmatched = 0;
  Progress * progress = nullptr;  /* the owner's progress bar; worker updates it under mutex_output */

  explicit search_cli_state_s(struct Parameters const & params) : parameters(params) {}
};


static auto search_output_results(struct search_cli_state_s & state,
                           std::vector<struct hit> const & hits,
                           char const * query_head,
                           int qseqlen,
                           char const * qsequence,
                           char const * qsequence_rc,
                           int64_t qsize) -> void
{
  std::lock_guard<std::mutex> const lock(state.mutex_output);

  /* show results */
  auto const toreport = std::min(state.parameters.opt_maxhits, static_cast<int64_t>(hits.size()));

  if (state.fp_alnout != nullptr)
    {
      results_show_alnout(state.fp_alnout,
                          hits.data(),
                          static_cast<int>(toreport),
                          query_head,
                          qsequence,
                          qseqlen,
                          state.parameters);
    }

  if (state.fp_lcaout != nullptr)
    {
      results_show_lcaout(state.fp_lcaout,
                          hits.data(),
                          static_cast<int>(toreport),
                          query_head,
                          state.parameters);
    }

  if (state.fp_samout != nullptr)
    {
      results_show_samout(state.fp_samout,
                          hits.data(),
                          static_cast<int>(toreport),
                          query_head,
                          qsequence,
                          qsequence_rc,
                          state.parameters);
    }

  if (toreport != 0)  // hits.size() >=1 and <= opt_maxhits
    {
      double const top_hit_id = hits[0].id;

      if ((state.parameters.opt_otutabout != nullptr) || (state.parameters.opt_mothur_shared_out != nullptr) || (state.parameters.opt_biomout != nullptr))
        {
          otutable_add(query_head,
                       db_getheader(static_cast<uint64_t>(hits[0].target)),
                       qsize);
        }

      for (int64_t t = 0; t < toreport; t++)
        {
          auto const * hp = &hits[static_cast<std::size_t>(t)];

          if ((state.parameters.opt_top_hits_only != 0) && (hp->id < top_hit_id))
            {
              break;
            }

          if (state.fp_fastapairs != nullptr)
            {
              results_show_fastapairs_one(state.fp_fastapairs,
                                          hp,
                                          query_head,
                                          qsequence,
                                          qsequence_rc,
                                          state.parameters);
            }

          if (state.fp_qsegout != nullptr)
            {
              results_show_qsegout_one(state.fp_qsegout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.parameters);
            }

          if (state.fp_tsegout != nullptr)
            {
              results_show_tsegout_one(state.fp_tsegout,
                                       hp,
                                       state.parameters);
            }

          if (state.fp_uc != nullptr)
            {
              if ((t==0) || (state.parameters.opt_uc_allhits))
                {
                  results_show_uc_one(state.fp_uc,
                                      hp,
                                      query_head,
                                      qseqlen,
                                      hp->target,
                                      state.parameters);
                }
            }

          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.parameters);
            }

          if (state.fp_blast6out != nullptr)
            {
              results_show_blast6out_one(state.fp_blast6out,
                                         hp,
                                         query_head,
                                         qseqlen);
            }
        }
    }
  else
    {
      if ((state.parameters.opt_otutabout != nullptr) || (state.parameters.opt_mothur_shared_out != nullptr) || (state.parameters.opt_biomout != nullptr))
        {
          otutable_add(query_head,
                       nullptr,
                       qsize);
        }

      if (state.fp_uc != nullptr)
        {
          results_show_uc_one(state.fp_uc,
                              nullptr,
                              query_head,
                              qseqlen,
                              0,
                              state.parameters);
        }

      if (state.parameters.opt_output_no_hits != 0)
        {
          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout,
                                       nullptr,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.parameters);
            }

          if (state.fp_blast6out != nullptr)
            {
              results_show_blast6out_one(state.fp_blast6out,
                                         nullptr,
                                         query_head,
                                         qseqlen);
            }
        }
    }

  if (not hits.empty())
    {
      state.count_matched++;
      if (state.parameters.opt_matched != nullptr)
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
                              state.parameters);
        }
    }
  else
    {
      state.count_notmatched++;
      if (state.parameters.opt_notmatched != nullptr)
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
                              state.parameters);
        }
    }

  /* update matching db sequences */
  for (auto const & hit : hits) {
    if (hit.accepted or hit.weak) {
      state.dbmatched[hit.target] += state.parameters.opt_sizein ? static_cast<uint64_t>(qsize) : 1;
    }
  }
}


static auto search_query(struct search_cli_state_s & state, uint64_t t) -> int
{
  struct searchinfo_s * const si_plus = state.si_plus;
  struct searchinfo_s * const si_minus = state.si_minus;

  for (int s = 0; s < number_of_strands(state.parameters.opt_strand); s++)
    {
      struct searchinfo_s * si = (s != 0) ? si_minus + t : si_plus + t;

      /* mask query */
      if (state.parameters.opt_qmask == MASK_DUST)
        {
          dust(si->qsequence, si->qseqlen, state.parameters);
        }
      else if ((state.parameters.opt_qmask == MASK_SOFT) && (state.parameters.opt_hardmask))
        {
          hardmask(si->qsequence, si->qseqlen);
        }

      /* perform search */
      search_onequery(si, static_cast<int>(state.parameters.opt_qmask));
    }

  std::vector<struct hit> hits;

  search_joinhits(si_plus + t,
                  state.parameters.opt_strand ? si_minus + t : nullptr,
                  hits);

  search_output_results(state,
                        hits,
                        si_plus[t].query_head,
                        si_plus[t].qseqlen,
                        si_plus[t].qsequence,
                        state.parameters.opt_strand ? si_minus[t].qsequence : nullptr,
                        si_plus[t].qsize);

  /* free memory for alignment strings */
  for (auto const & hit : hits) {
    if (hit.aligned) {
      xfree(hit.nwalignment);
    }
  }

  return static_cast<int>(hits.size());
}


static auto populate_si(struct searchinfo_s * si,
                        const char * head,
                        int head_len,
                        const char * seq,
                        int seq_len,
                        int query_no,
                        int64_t qsize,
                        int strand) -> void
{
  si->query_head_len = head_len;
  si->qseqlen = seq_len;
  si->query_no = query_no;
  si->qsize = qsize;
  si->strand = strand;

  /* allocate more memory for header and sequence, if necessary */

  if (si->query_head_len + 1 > si->query_head_alloc)
    {
      si->query_head_alloc = si->query_head_len + buffer_headroom;
      si->query_head = static_cast<char *>(
        xrealloc(si->query_head, static_cast<size_t>(si->query_head_alloc)));
    }

  if (si->qseqlen + 1 > si->seq_alloc)
    {
      si->seq_alloc = si->qseqlen + buffer_headroom;
      si->qsequence = static_cast<char *>(
        xrealloc(si->qsequence, static_cast<size_t>(si->seq_alloc)));
    }

  /* copy header */
  std::strcpy(si->query_head, head);

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


static auto search_thread_run(struct search_cli_state_s & state, uint64_t t) -> void
{
  auto const query_fastx_h = state.query_fastx_h;
  struct searchinfo_s * const si_plus = state.si_plus;

  while (true)
    {
      std::unique_lock<std::mutex> input_lock(state.mutex_input);

      if (fastx_next(query_fastx_h,
                     (not state.parameters.opt_notrunclabels),
                     chrmap_no_change()))
        {
          char const * qhead = fastx_get_header(query_fastx_h);
          int const query_head_len = static_cast<int>(fastx_get_header_length(query_fastx_h));
          char const * qseq = fastx_get_sequence(query_fastx_h);
          int const qseqlen = static_cast<int>(fastx_get_sequence_length(query_fastx_h));
          int const query_no = static_cast<int>(fastx_get_seqno(query_fastx_h));
          int64_t const qsize = fastx_get_abundance(query_fastx_h);

          populate_si(si_plus + t,
                      qhead,
                      query_head_len,
                      qseq,
                      qseqlen,
                      query_no,
                      qsize,
                      0);

          /* get progress as amount of input file read */
          uint64_t const progress = fastx_get_position(query_fastx_h);

          /* let other threads read input */
          input_lock.unlock();

          if (state.parameters.opt_strand)
            {
              populate_si(state.si_minus + t,
                          si_plus[t].query_head,
                          query_head_len,
                          si_plus[t].qsequence,
                          qseqlen,
                          query_no,
                          qsize,
                          1);
            }

          int const match = search_query(state, t);

          /* lock mutex for update of global data and output */
          std::lock_guard<std::mutex> const output_lock(state.mutex_output);

          /* update stats */
          ++state.queries;
          state.queries_abundance += static_cast<uint64_t>(qsize);

          if (match != 0)
            {
              ++state.qmatches;
              state.qmatches_abundance += static_cast<uint64_t>(qsize);
            }

          /* show progress */
          state.progress->update(progress);
        }
      else
        {
          /* input_lock released by RAII when leaving the loop body */
          break;
        }
    }
}


static auto search_thread_init(struct searchinfo_s * si, int const seqcount, int const tophits,
                               struct Parameters const & parameters,
                               struct Dbindex const & dbindex) -> void
{
  /* thread specific initialiation */
  si->parameters = &parameters;  /* searchcore reads config through the si (E1) */
  si->dbindex = &dbindex;  /* searchcore reads the k-mer index through the si */
  si->uh = unique_init();
  si->kmers = static_cast<count_t *>(xmalloc((static_cast<size_t>(seqcount) * sizeof(count_t)) + 32));
  si->m = minheap_init(tophits);
  si->hits = static_cast<struct hit *>(xmalloc
    (sizeof(struct hit) * static_cast<size_t>(tophits) * static_cast<size_t>(number_of_strands(parameters.opt_strand))));
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = nullptr;
  si->seq_alloc = 0;
  si->qsequence = nullptr;
  si->s = search16_init(static_cast<CELL>(parameters.opt_match),
                        static_cast<CELL>(parameters.opt_mismatch),
                        static_cast<CELL>(parameters.opt_gap_open_query_left),
                        static_cast<CELL>(parameters.opt_gap_open_target_left),
                        static_cast<CELL>(parameters.opt_gap_open_query_interior),
                        static_cast<CELL>(parameters.opt_gap_open_target_interior),
                        static_cast<CELL>(parameters.opt_gap_open_query_right),
                        static_cast<CELL>(parameters.opt_gap_open_target_right),
                        static_cast<CELL>(parameters.opt_gap_extension_query_left),
                        static_cast<CELL>(parameters.opt_gap_extension_target_left),
                        static_cast<CELL>(parameters.opt_gap_extension_query_interior),
                        static_cast<CELL>(parameters.opt_gap_extension_target_interior),
                        static_cast<CELL>(parameters.opt_gap_extension_query_right),
                        static_cast<CELL>(parameters.opt_gap_extension_target_right),
                        parameters.opt_n_mismatch);
}


static auto search_thread_exit(struct searchinfo_s * si) -> void
{
  /* thread specific clean up */
  search16_exit(si->s);
  unique_exit(si->uh);
  xfree(si->hits);
  minheap_exit(si->m);
  xfree(si->kmers);
  if (si->query_head != nullptr)
    {
      xfree(si->query_head);
    }
  if (si->qsequence != nullptr)
    {
      xfree(si->qsequence);
    }
}



static auto search_thread_worker_run(struct search_cli_state_s & state) -> void
{
  struct searchinfo_s * const si_plus = state.si_plus;
  struct searchinfo_s * const si_minus = state.si_minus;
  int const seqcount = state.seqcount;
  int const tophits = state.tophits;

  /* init per-thread search state before the workers start */
  for (int t = 0; t < state.parameters.opt_threads; t++)
    {
      search_thread_init(si_plus + t, seqcount, tophits, state.effective_parameters, state.dbindex);
      if (si_minus != nullptr)
        {
          search_thread_init(si_minus + t, seqcount, tophits, state.effective_parameters, state.dbindex);
        }
    }

  /* run the worker pool over the input file */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(state.parameters.opt_threads),
                              [&state](uint64_t const t)
                              { search_thread_run(state, t); });
    threadrunner.run();
  }

  /* clean up per-thread search state */
  for (int t = 0; t < state.parameters.opt_threads; t++)
    {
      search_thread_exit(si_plus + t);
      if (si_minus != nullptr)
        {
          search_thread_exit(si_minus + t);
        }
    }
}


static auto search_prep(struct search_cli_state_s & state) -> void
{
  /* open output files */

  state.fp_alnout = open_optional_output(state.parameters.opt_alnout, "alignment");
  if (state.fp_alnout != nullptr)
    {
      std::fprintf(state.fp_alnout, "%s\n", state.parameters.command_line.c_str());
      std::fprintf(state.fp_alnout, "%s\n", state.parameters.prog_header.c_str());
    }

  state.fp_lcaout = open_optional_output(state.parameters.opt_lcaout, "lca");
  state.fp_samout = open_optional_output(state.parameters.opt_samout, "SAM");
  state.fp_userout = open_optional_output(state.parameters.opt_userout, "user-defined");
  state.fp_blast6out = open_optional_output(state.parameters.opt_blast6out, "blast6-like");
  state.fp_uc = open_optional_output(state.parameters.opt_uc, "uc");
  state.fp_fastapairs = open_optional_output(state.parameters.opt_fastapairs, "fastapairs");
  state.fp_qsegout = open_optional_output(state.parameters.opt_qsegout, "qsegout");
  state.fp_tsegout = open_optional_output(state.parameters.opt_tsegout, "tsegout");
  state.fp_matched = open_optional_output(state.parameters.opt_matched, "matched");
  state.fp_notmatched = open_optional_output(state.parameters.opt_notmatched, "notmatched");
  state.fp_otutabout = open_optional_output(state.parameters.opt_otutabout, "OTU table (text format)");
  state.fp_mothur_shared_out = open_optional_output(state.parameters.opt_mothur_shared_out, "OTU table (mothur format)");
  state.fp_biomout = open_optional_output(state.parameters.opt_biomout, "OTU table (biom 1.0 format)");

  /* check if it may be an UDB file */

  bool const is_udb = udb_detect_isudb(state.parameters.opt_db);

  if (is_udb)
    {
      udb_read(state.parameters.opt_db, true, true, state.dbindex, state.parameters);
      results_show_samheader(state.fp_samout, state.parameters.opt_db, state.parameters);
      // memory-intensive: the entire database is now held in memory
      state.seqcount = static_cast<int>(db_getsequencecount());
    }
  else
    {
      db_read(state.parameters.opt_db, 0, state.parameters);
      results_show_samheader(state.fp_samout, state.parameters.opt_db, state.parameters);
      if (state.parameters.opt_dbmask == MASK_DUST)
        {
          dust_all(state.parameters);
        }
      else if ((state.parameters.opt_dbmask == MASK_SOFT) && (state.parameters.opt_hardmask))
        {
          hardmask_all();
        }
      // memory-intensive: the entire database is now held in memory
      state.seqcount = static_cast<int>(db_getsequencecount());
      state.dbindex.prepare(1, static_cast<int>(state.parameters.opt_dbmask), state.parameters);
      state.dbindex.add_all_sequences(static_cast<int>(state.parameters.opt_dbmask), state.parameters);
    }

  /* tophits = the maximum number of hits we need to store */

  /* Clamp maxrejects/maxaccepts to the database size (0 or "> seqcount" means
     "all"). Apply the clamp to a local Parameters copy that is threaded to the
     workers via si->parameters, rather than mutating the shared config globals
     (E1 trap-writer): the shared searchcore reads the clamped values through
     si->parameters. */
  state.effective_parameters = state.parameters;
  if ((state.effective_parameters.opt_maxrejects == 0) ||
      (state.effective_parameters.opt_maxrejects > state.seqcount))
    {
      state.effective_parameters.opt_maxrejects = state.seqcount;
    }

  if ((state.effective_parameters.opt_maxaccepts == 0) ||
      (state.effective_parameters.opt_maxaccepts > state.seqcount))
    {
      state.effective_parameters.opt_maxaccepts = state.seqcount;
    }

  state.tophits = static_cast<int>(state.effective_parameters.opt_maxrejects +
                                   state.effective_parameters.opt_maxaccepts + MAXDELAYED);

  state.tophits = std::min(state.tophits, state.seqcount);
}


static auto search_done(struct search_cli_state_s & state) -> void
{
  /* clean up, global */

  state.dbindex.clear();
  db_free();

  /* fclose_output() is a no-op on a null handle, so the unopened outputs
     need no guard; only userout carries extra teardown. */
  fclose_output(state.fp_lcaout);
  fclose_output(state.fp_matched);
  fclose_output(state.fp_notmatched);
  fclose_output(state.fp_fastapairs);
  fclose_output(state.fp_qsegout);
  fclose_output(state.fp_tsegout);
  fclose_output(state.fp_uc);
  fclose_output(state.fp_blast6out);
  fclose_output(state.fp_userout);
  fclose_output(state.fp_alnout);
  fclose_output(state.fp_samout);
}


auto usearch_global(struct Parameters const & parameters) -> void
{
  /* Per-invocation state, owned here and threaded through the worker pool and
     the output helper (E4). Aliased by reference so the long body below reads
     unchanged; the workers receive `state`, not file-static globals. */
  struct search_cli_state_s state(parameters);
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;
  auto & query_fastx_h = state.query_fastx_h;
  auto & seqcount = state.seqcount;
  auto & qmatches = state.qmatches;
  auto & qmatches_abundance = state.qmatches_abundance;
  auto & queries = state.queries;
  auto & queries_abundance = state.queries_abundance;
  auto & dbmatched = state.dbmatched;
  auto & fp_dbmatched = state.fp_dbmatched;
  auto & fp_dbnotmatched = state.fp_dbnotmatched;

  search_prep(state);

  fp_dbmatched = open_optional_output(parameters.opt_dbmatched, "dbmatched");
  fp_dbnotmatched = open_optional_output(parameters.opt_dbnotmatched, "dbnotmatched");

  dbmatched = static_cast<uint64_t *>(xmalloc(static_cast<size_t>(seqcount) * sizeof(uint64_t)));
  std::memset(dbmatched, 0, static_cast<size_t>(seqcount) * sizeof(uint64_t));

  otutable_init();

  /* prepare reading of queries */
  qmatches = 0;
  qmatches_abundance = 0;
  queries = 0;
  queries_abundance = 0;
  query_fastx_h = fastx_open(parameters.opt_usearch_global, parameters);

  /* The query file is parsed inside the worker threads (search_thread_run).
     Defer parse errors so a malformed query stops the pool cooperatively
     instead of calling fatal()/std::exit() from a worker while siblings are
     writing output (CC3); reported below from the main thread after join. */
  query_fastx_h->defer_errors = true;

  /* allocate memory for thread info */
  si_plus = new searchinfo_s[parameters.opt_threads]{};
  if (parameters.opt_strand)
    {
      si_minus = new searchinfo_s[parameters.opt_threads]{};
    }
  else
    {
      si_minus = nullptr;
    }

  {
    Progress progress("Searching", fastx_get_size(query_fastx_h), parameters);
    state.progress = &progress;
    search_thread_worker_run(state);
  }

  /* all workers joined; report a deferred query parse error (CC3) from the
     main thread so it does not race a worker's output */
  if (fastx_get_error(query_fastx_h))
    {
      fatal("%s", fastx_get_errmsg(query_fastx_h));
    }

  delete [] si_plus;
  if (si_minus != nullptr)
    {
      delete [] si_minus;
    }

  fastx_close(query_fastx_h, parameters);

  if (! parameters.opt_quiet)
    {
      std::fprintf(stderr, "Matching unique query sequences: %d of %d",
              qmatches, queries);
      if (queries > 0)
        {
          std::fprintf(stderr, " (%.2f%%)", 100.0 * qmatches / queries);
        }
      std::fprintf(stderr, "\n");
      if (parameters.opt_sizein)
        {
          std::fprintf(stderr, "Matching total query sequences: %" PRIu64 " of %"
                  PRIu64,
                  qmatches_abundance, queries_abundance);
          if (queries_abundance > 0)
            {
              std::fprintf(stderr, " (%.2f%%)",
                      100.0 * static_cast<double>(qmatches_abundance) / static_cast<double>(queries_abundance));
            }
          std::fprintf(stderr, "\n");
        }
    }

  if (parameters.opt_log != nullptr)
    {
      std::fprintf(parameters.fp_log, "Matching unique query sequences: %d of %d",
              qmatches, queries);
      if (queries > 0)
        {
          std::fprintf(parameters.fp_log, " (%.2f%%)", 100.0 * qmatches / queries);
        }
      std::fprintf(parameters.fp_log, "\n");
      if (parameters.opt_sizein)
        {
          std::fprintf(parameters.fp_log, "Matching total query sequences: %" PRIu64 " of %"
                  PRIu64,
                  qmatches_abundance, queries_abundance);
          if (queries_abundance > 0)
            {
              std::fprintf(parameters.fp_log, " (%.2f%%)",
                      100.0 * static_cast<double>(qmatches_abundance) / static_cast<double>(queries_abundance));
            }
          std::fprintf(parameters.fp_log, "\n");
        }
    }


  // Add OTUs with no matches to OTU table
  if ((parameters.opt_otutabout != nullptr) || (parameters.opt_mothur_shared_out != nullptr) || (parameters.opt_biomout != nullptr)) {
    for (int64_t i = 0; i < seqcount; i++) {
      if (dbmatched[i] == 0U) {
        otutable_add(nullptr, db_getheader(static_cast<uint64_t>(i)), 0);
      }
    }
  }

  if (parameters.opt_biomout != nullptr)
    {
      otutable_print_biomout(state.fp_biomout, state.parameters);
      fclose_output(state.fp_biomout);
    }

  if (parameters.opt_otutabout != nullptr)
    {
      otutable_print_otutabout(state.fp_otutabout, state.parameters);
      fclose_output(state.fp_otutabout);
    }

  if (parameters.opt_mothur_shared_out != nullptr)
    {
      otutable_print_mothur_shared_out(state.fp_mothur_shared_out, state.parameters);
      fclose_output(state.fp_mothur_shared_out);
    }

  otutable_done();

  if ((parameters.opt_dbmatched != nullptr) || (parameters.opt_dbnotmatched != nullptr))
    {
      int count_dbmatched = 0;
      int count_dbnotmatched = 0;

      for (int64_t i = 0; i < seqcount; i++)
        {
          if (dbmatched[i] != 0U)
            {
              count_dbmatched++;
              if (parameters.opt_dbmatched != nullptr)
                {
                  fasta_print_general(fp_dbmatched,
                                      nullptr,
                                      db_getsequence(static_cast<uint64_t>(i)),
                                      static_cast<int>(db_getsequencelen(static_cast<uint64_t>(i))),
                                      db_getheader(static_cast<uint64_t>(i)),
                                      static_cast<int>(db_getheaderlen(static_cast<uint64_t>(i))),
                                      dbmatched[i],
                                      count_dbmatched,
                                      -1.0,
                                      -1, -1, nullptr, 0.0,
                                      0,
                                      parameters);
                }
            }
          else
            {
              count_dbnotmatched++;
              if (parameters.opt_dbnotmatched != nullptr)
                {
                  fasta_print_general(fp_dbnotmatched,
                                      nullptr,
                                      db_getsequence(static_cast<uint64_t>(i)),
                                      static_cast<int>(db_getsequencelen(static_cast<uint64_t>(i))),
                                      db_getheader(static_cast<uint64_t>(i)),
                                      static_cast<int>(db_getheaderlen(static_cast<uint64_t>(i))),
                                      db_getabundance(static_cast<uint64_t>(i)),
                                      count_dbnotmatched,
                                      -1.0,
                                      -1, -1, nullptr, 0.0,
                                      0,
                                      parameters);
                }
            }
        }
    }

  xfree(dbmatched);

  fclose_output(fp_dbmatched);
  fclose_output(fp_dbnotmatched);

  search_done(state);
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
                         struct Dbindex const & dbindex) -> void
{
  /* Initialize search session for library use.
     Mirrors cluster_session_init: stores seqcount/tophits in the session,
     allocates si_plus (always) and si_minus (when searching both strands). */
  ss->parameters = &parameters;
  ss->dbindex = &dbindex;
  ss->seqcount = static_cast<int>(db_getsequencecount());
  /* The library path does not clamp to the database size (only the CLI
     search_prep does), so the sizing uses the configured values from
     parameters; si->parameters (set in search_thread_init) carries the same. */
  ss->tophits = static_cast<int>(parameters.opt_maxaccepts + parameters.opt_maxrejects + MAXDELAYED);
  if (ss->tophits > ss->seqcount)
    {
      ss->tophits = ss->seqcount;
    }

  ss->si_plus = make_unique<searchinfo_s>();
  search_thread_init(ss->si_plus.get(), ss->seqcount, ss->tophits, parameters, *ss->dbindex);
  ss->si_plus->strand = 0;

  if (parameters.opt_strand)
    {
      ss->si_minus = make_unique<searchinfo_s>();
      search_thread_init(ss->si_minus.get(), ss->seqcount, ss->tophits, parameters, *ss->dbindex);
      ss->si_minus->strand = 1;
    }
}


auto search_session_single(struct search_session_s * ss,
                           const char * query_seq,
                           const char * query_head,
                           int query_len,
                           int64_t query_size,
                           struct search_result_s * results,
                           int max_results,
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

      if (parameters.opt_qmask == MASK_DUST)
        {
          dust(strand_si->qsequence, strand_si->qseqlen, parameters);
        }
      else if ((parameters.opt_qmask == MASK_SOFT) && (parameters.opt_hardmask))
        {
          hardmask(strand_si->qsequence, strand_si->qseqlen);
        }

      search_onequery(strand_si, static_cast<int>(parameters.opt_qmask));
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
      r.target_length = static_cast<int>(db_getsequencelen(static_cast<uint64_t>(h.target)));
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

  while (true)
    {
      /* grab next query */
      int qi {0};
      {
        std::lock_guard<std::mutex> const lock(ctx.mutex);
        qi = ctx.next_query++;
      }

      if (qi >= ctx.query_count)
        {
          break;
        }

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

          if (parameters.opt_qmask == MASK_DUST)
            {
              dust(strand_si->qsequence, strand_si->qseqlen, parameters);
            }
          else if ((parameters.opt_qmask == MASK_SOFT) && (parameters.opt_hardmask))
            {
              hardmask(strand_si->qsequence, strand_si->qseqlen);
            }

          search_onequery(strand_si, static_cast<int>(parameters.opt_qmask));
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
          r.target_length = static_cast<int>(db_getsequencelen(static_cast<uint64_t>(h.target)));
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
    }
}


auto search_batch(struct Parameters const & parameters,
                  struct Dbindex const & dbindex,
                  const char ** query_seqs,
                  const char ** query_heads,
                  const int * query_lens,
                  const int64_t * query_sizes,
                  int query_count,
                  struct search_result_s * results,
                  int max_results_per_query,
                  int * result_counts) -> void
{
  /* per-thread buffer sizes for search_thread_init (formerly file-statics).
     The library path does not clamp to the database size (only the CLI
     search_prep does), so the sizing uses the configured values from
     parameters; si->parameters (set in search_thread_init) carries the same. */
  int const seqcount = static_cast<int>(db_getsequencecount());
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
      search_thread_init(ctx.batch_si_plus + t, seqcount, tophits, parameters, dbindex);
      if (ctx.batch_si_minus != nullptr)
        {
          search_thread_init(ctx.batch_si_minus + t, seqcount, tophits, parameters, dbindex);
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
