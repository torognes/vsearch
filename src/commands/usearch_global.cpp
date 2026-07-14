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
#include "commands/usearch_global.hpp"
#include "core/db.hpp"
#include "core/fasta.hpp"
#include "core/fastx.hpp"
#include "core/results.hpp"
#include "os/system.hpp"
#include "core/search_internal.hpp"
#include "utils/progress.hpp"
#include "core/searchcore.hpp"
#include "core/dbindex.hpp"
#include "core/mask.hpp"
#include "core/otutable.hpp"
#include "core/udb.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/number_of_strands.hpp"
#include "utils/open_file.hpp"
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include <algorithm>  // std::min
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // uint64_t, int64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstring>  // std::memset
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
  struct Database db;  /* the sequence database this run owns (RAII); si->db points here */
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
  /* RAII output handles; the workers read the raw FILE * via .get() under
     mutex_output. Closed explicitly with reset() in a fixed order (see
     search_done and the OTU/db blocks in usearch_global) so streams sharing
     stdout flush in the legacy order rather than the reverse order a struct
     destructor would use. */
  OutputFileHandle fp_samout;
  OutputFileHandle fp_alnout;
  OutputFileHandle fp_userout;
  OutputFileHandle fp_blast6out;
  OutputFileHandle fp_uc;
  OutputFileHandle fp_fastapairs;
  OutputFileHandle fp_matched;
  OutputFileHandle fp_notmatched;
  OutputFileHandle fp_dbmatched;
  OutputFileHandle fp_dbnotmatched;
  OutputFileHandle fp_otutabout;
  OutputFileHandle fp_mothur_shared_out;
  OutputFileHandle fp_biomout;
  OutputFileHandle fp_lcaout;
  OutputFileHandle fp_qsegout;
  OutputFileHandle fp_tsegout;
  int count_matched = 0;
  int count_notmatched = 0;
  Progress * progress = nullptr;  /* the owner's progress bar; worker updates it under mutex_output */

  explicit search_cli_state_s(struct Parameters const & params) : parameters(params) {}
};


static auto search_output_results(struct search_cli_state_s & state,
                           std::vector<struct hit> const & hits,
                           char const * query_head,
                           int const qseqlen,
                           char const * qsequence,
                           char const * qsequence_rc,
                           int64_t const qsize) -> void
{
  std::lock_guard<std::mutex> const lock(state.mutex_output);

  /* show results */
  auto const toreport = std::min(state.parameters.opt_maxhits, static_cast<int64_t>(hits.size()));

  if (state.fp_alnout != nullptr)
    {
      results_show_alnout(state.fp_alnout.get(),
                          hits.data(),
                          static_cast<int>(toreport),
                          query_head,
                          qsequence,
                          qseqlen,
                          state.db,
                          state.parameters);
    }

  if (state.fp_lcaout != nullptr)
    {
      results_show_lcaout(state.fp_lcaout.get(),
                          hits.data(),
                          static_cast<int>(toreport),
                          query_head,
                          state.db,
                          state.parameters);
    }

  if (state.fp_samout != nullptr)
    {
      results_show_samout(state.fp_samout.get(),
                          hits.data(),
                          static_cast<int>(toreport),
                          query_head,
                          qsequence,
                          qsequence_rc,
                          state.db,
                          state.parameters);
    }

  if (toreport != 0)  // hits.size() >=1 and <= opt_maxhits
    {
      double const top_hit_id = hits[0].id;

      if ((state.parameters.opt_otutabout != nullptr) || (state.parameters.opt_mothur_shared_out != nullptr) || (state.parameters.opt_biomout != nullptr))
        {
          otutable_add(query_head,
                       state.db.getheader(static_cast<uint64_t>(hits[0].target)),
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
              results_show_fastapairs_one(state.fp_fastapairs.get(),
                                          hp,
                                          query_head,
                                          qsequence,
                                          qsequence_rc,
                                          state.db,
                                          state.parameters);
            }

          if (state.fp_qsegout != nullptr)
            {
              results_show_qsegout_one(state.fp_qsegout.get(),
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.parameters);
            }

          if (state.fp_tsegout != nullptr)
            {
              results_show_tsegout_one(state.fp_tsegout.get(),
                                       hp,
                                       state.db,
                                       state.parameters);
            }

          if (state.fp_uc != nullptr)
            {
              if ((t==0) || (state.parameters.opt_uc_allhits))
                {
                  results_show_uc_one(state.fp_uc.get(),
                                      hp,
                                      query_head,
                                      qseqlen,
                                      hp->target,
                                      state.db,
                                      state.parameters);
                }
            }

          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout.get(),
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.db,
                                       state.parameters);
            }

          if (state.fp_blast6out != nullptr)
            {
              results_show_blast6out_one(state.fp_blast6out.get(),
                                         hp,
                                         query_head,
                                         qseqlen,
                                         state.db);
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
          results_show_uc_one(state.fp_uc.get(),
                              nullptr,
                              query_head,
                              qseqlen,
                              0,
                              state.db,
                              state.parameters);
        }

      if (state.parameters.opt_output_no_hits != 0)
        {
          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout.get(),
                                       nullptr,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc,
                                       state.db,
                                       state.parameters);
            }

          if (state.fp_blast6out != nullptr)
            {
              results_show_blast6out_one(state.fp_blast6out.get(),
                                         nullptr,
                                         query_head,
                                         qseqlen,
                                         state.db);
            }
        }
    }

  if (not hits.empty())
    {
      state.count_matched++;
      if (state.parameters.opt_matched != nullptr)
        {
          fasta_print_general(state.fp_matched.get(),
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
          fasta_print_general(state.fp_notmatched.get(),
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


static auto search_query(struct search_cli_state_s & state, uint64_t const t) -> int
{
  struct searchinfo_s * const si_plus = state.si_plus;
  struct searchinfo_s * const si_minus = state.si_minus;

  for (int s = 0; s < number_of_strands(state.parameters.opt_strand); s++)
    {
      struct searchinfo_s * si = (s != 0) ? si_minus + t : si_plus + t;

      /* mask query */
      if (state.parameters.opt_qmask == Masking::dust)
        {
          dust(si->qsequence, si->qseqlen, state.parameters);
        }
      else if ((state.parameters.opt_qmask == Masking::soft) && (state.parameters.opt_hardmask))
        {
          hardmask(si->qsequence, si->qseqlen);
        }

      /* perform search */
      search_onequery(si, state.parameters.opt_qmask);
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


static auto search_thread_run(struct search_cli_state_s & state, uint64_t const t) -> void
{
  auto const query_fastx_h = state.query_fastx_h;
  struct searchinfo_s * const si_plus = state.si_plus;

  int query_head_len = 0;
  int qseqlen = 0;
  int query_no = 0;
  int64_t qsize = 0;
  uint64_t progress = 0;

  auto const has_work_to_claim = [&]() -> bool {
    if (not fastx_next(query_fastx_h,
                       (not state.parameters.opt_notrunclabels),
                       chrmap_no_change()))
      {
        return false;
      }

    char const * qhead = fastx_get_header(query_fastx_h);
    query_head_len = static_cast<int>(fastx_get_header_length(query_fastx_h));
    char const * qseq = fastx_get_sequence(query_fastx_h);
    qseqlen = static_cast<int>(fastx_get_sequence_length(query_fastx_h));
    query_no = static_cast<int>(fastx_get_seqno(query_fastx_h));
    qsize = fastx_get_abundance(query_fastx_h);

    populate_si(si_plus + t,
                qhead,
                query_head_len,
                qseq,
                qseqlen,
                query_no,
                qsize,
                0);

    /* get progress as amount of input file read */
    progress = fastx_get_position(query_fastx_h);
    return true;
  };

  auto const process_query = [&]() {
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
  };

  run_worker_loop(state.mutex_input, has_work_to_claim, process_query);
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
      search_thread_init(si_plus + t, seqcount, tophits, state.effective_parameters, state.dbindex, state.db);
      if (si_minus != nullptr)
        {
          search_thread_init(si_minus + t, seqcount, tophits, state.effective_parameters, state.dbindex, state.db);
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

  state.fp_alnout = open_optional_output_file(state.parameters.opt_alnout, OutputOption{"--alnout"});
  if (state.fp_alnout != nullptr)
    {
      std::fprintf(state.fp_alnout.get(), "%s\n", state.parameters.command_line.c_str());
      std::fprintf(state.fp_alnout.get(), "%s\n", state.parameters.prog_header.c_str());
    }

  state.fp_lcaout = open_optional_output_file(state.parameters.opt_lcaout, OutputOption{"--lcaout"});
  state.fp_samout = open_optional_output_file(state.parameters.opt_samout, OutputOption{"--samout"});
  state.fp_userout = open_optional_output_file(state.parameters.opt_userout, OutputOption{"--userout"});
  state.fp_blast6out = open_optional_output_file(state.parameters.opt_blast6out, OutputOption{"--blast6out"});
  state.fp_uc = open_optional_output_file(state.parameters.opt_uc, OutputOption{"--uc"});
  state.fp_fastapairs = open_optional_output_file(state.parameters.opt_fastapairs, OutputOption{"--fastapairs"});
  state.fp_qsegout = open_optional_output_file(state.parameters.opt_qsegout, OutputOption{"--qsegout"});
  state.fp_tsegout = open_optional_output_file(state.parameters.opt_tsegout, OutputOption{"--tsegout"});
  state.fp_matched = open_optional_output_file(state.parameters.opt_matched, OutputOption{"--matched"});
  state.fp_notmatched = open_optional_output_file(state.parameters.opt_notmatched, OutputOption{"--notmatched"});
  state.fp_otutabout = open_optional_output_file(state.parameters.opt_otutabout, OutputOption{"--otutabout"});
  state.fp_mothur_shared_out = open_optional_output_file(state.parameters.opt_mothur_shared_out, OutputOption{"--mothur_shared_out"});
  state.fp_biomout = open_optional_output_file(state.parameters.opt_biomout, OutputOption{"--biomout"});

  /* check if it may be an UDB file */

  bool const is_udb = udb_detect_isudb(state.parameters.opt_db);

  if (is_udb)
    {
      udb_read(state.parameters.opt_db, true, true, state.dbindex, state.db, state.parameters);
      results_show_samheader(state.fp_samout.get(), state.parameters.opt_db, state.db, state.parameters);
      // memory-intensive: the entire database is now held in memory
      state.seqcount = static_cast<int>(state.db.getsequencecount());
    }
  else
    {
      state.db.read(state.parameters.opt_db, 0, state.parameters);
      results_show_samheader(state.fp_samout.get(), state.parameters.opt_db, state.db, state.parameters);
      if (state.parameters.opt_dbmask == Masking::dust)
        {
          dust_all(state.db, state.parameters);
        }
      else if ((state.parameters.opt_dbmask == Masking::soft) && (state.parameters.opt_hardmask))
        {
          hardmask_all(state.db);
        }
      // memory-intensive: the entire database is now held in memory
      state.seqcount = static_cast<int>(state.db.getsequencecount());
      state.dbindex.prepare(1, state.parameters.opt_dbmask, state.db, state.parameters);
      state.dbindex.add_all_sequences(state.parameters.opt_dbmask, state.db, state.parameters);
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
  state.db.clear();

  /* reset() is a no-op on an empty handle, so unopened outputs need no guard.
     The fixed order matches the legacy fclose sequence: RAII scope-exit would
     reverse it and flip the flush order for outputs that share stdout. */
  state.fp_lcaout.reset();
  state.fp_matched.reset();
  state.fp_notmatched.reset();
  state.fp_fastapairs.reset();
  state.fp_qsegout.reset();
  state.fp_tsegout.reset();
  state.fp_uc.reset();
  state.fp_blast6out.reset();
  state.fp_userout.reset();
  state.fp_alnout.reset();
  state.fp_samout.reset();
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

  fp_dbmatched = open_optional_output_file(parameters.opt_dbmatched, OutputOption{"--dbmatched"});
  fp_dbnotmatched = open_optional_output_file(parameters.opt_dbnotmatched, OutputOption{"--dbnotmatched"});

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
        otutable_add(nullptr, state.db.getheader(static_cast<uint64_t>(i)), 0);
      }
    }
  }

  if (parameters.opt_biomout != nullptr)
    {
      otutable_print_biomout(state.fp_biomout.get(), state.parameters);
      state.fp_biomout.reset();
    }

  if (parameters.opt_otutabout != nullptr)
    {
      otutable_print_otutabout(state.fp_otutabout.get(), state.parameters);
      state.fp_otutabout.reset();
    }

  if (parameters.opt_mothur_shared_out != nullptr)
    {
      otutable_print_mothur_shared_out(state.fp_mothur_shared_out.get(), state.parameters);
      state.fp_mothur_shared_out.reset();
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
                  fasta_print_general(fp_dbmatched.get(),
                                      nullptr,
                                      state.db.getsequence(static_cast<uint64_t>(i)),
                                      static_cast<int>(state.db.getsequencelen(static_cast<uint64_t>(i))),
                                      state.db.getheader(static_cast<uint64_t>(i)),
                                      static_cast<int>(state.db.getheaderlen(static_cast<uint64_t>(i))),
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
                  fasta_print_general(fp_dbnotmatched.get(),
                                      nullptr,
                                      state.db.getsequence(static_cast<uint64_t>(i)),
                                      static_cast<int>(state.db.getsequencelen(static_cast<uint64_t>(i))),
                                      state.db.getheader(static_cast<uint64_t>(i)),
                                      static_cast<int>(state.db.getheaderlen(static_cast<uint64_t>(i))),
                                      state.db.getabundance(static_cast<uint64_t>(i)),
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

  fp_dbmatched.reset();
  fp_dbnotmatched.reset();

  search_done(state);
}
