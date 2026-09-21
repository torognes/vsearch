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
#include "vsearch.hpp"
#include <memory>  // std::unique_ptr
#include "commands/usearch_global.hpp"
#include "commands/global_search_internal.hpp"  // enum struct Prefilter
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/db.hpp"
#include "core/match_counts.hpp"  // vsearch::MatchCounts, vsearch::print_match_counts
#include "core/fasta.hpp"
#include "core/fastx.hpp"
#include "core/results.hpp"
#include "core/search_internal.hpp"
#include "utils/progress.hpp"
#include "core/searchcore.hpp"
#include "core/dbindex.hpp"
#include "core/mask.hpp"
#include "core/otutable.hpp"
#include "core/udb.hpp"
#include "utils/base_mapping.hpp"
#include "utils/fatal.hpp"
#include "utils/fatal_allocator.hpp"  // FatalAllocator
#include "utils/number_of_strands.hpp"
#include "utils/open_file.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include <algorithm>  // std::min
#include <array>  // std::array
#include <cstdint>  // uint64_t, int64_t
#include <cstdio>  // std::FILE, std::fprintf
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
  /* which candidate set a query is compared against: the pre-filter's pick
     (--usearch_global) or the whole database (--search_global). The only
     difference between the two commands this struct serves. */
  Prefilter prefilter;
  /* a copy of parameters with opt_maxaccepts/opt_maxrejects clamped to the
     database size (search_prep); si->parameters points here so the shared
     searchcore reads the clamped values without a mutated global (E1). For
     --search_global it also carries the exhaustive configuration, which is
     what makes the options steering the pre-filter ignored rather than
     rejected (search_prep). */
  struct Parameters effective_parameters;
  struct Database db;  /* the sequence database this run owns (RAII); si->db points here */
  struct Dbindex dbindex;  /* the k-mer index this run owns (RAII); si->dbindex points here */
  int tophits = 0;   /* the maximum number of hits to keep */
  int seqcount = 0;  /* number of database sequences */
  std::vector<searchinfo_s> si_plus;
  std::vector<searchinfo_s> si_minus;  /* empty unless --strand both */
  fastx_handle query_fastx_h = nullptr;
  std::mutex mutex_input;   /* serializes query reads */
  std::mutex mutex_output;  /* serializes output + counter updates */
  int qmatches = 0;
  uint64_t qmatches_abundance = 0;
  int queries = 0;
  uint64_t queries_abundance = 0;
  std::vector<uint64_t, FatalAllocator<uint64_t>> dbmatched;
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
  /* accumulates the OTU table for the outputs above; the workers mutate it via
     add() under mutex_output, replacing the former file-static singleton */
  OtuTable otutable;
  OutputFileHandle fp_lcaout;
  OutputFileHandle fp_qsegout;
  OutputFileHandle fp_tsegout;
  int count_matched = 0;
  int count_notmatched = 0;
  Progress * progress = nullptr;  /* the owner's progress bar; worker updates it under mutex_output */

  search_cli_state_s(struct Parameters const & params, Prefilter const how)
    : parameters(params), prefilter(how) {}
};


static auto search_output_results(struct search_cli_state_s & state,
                           std::vector<struct hit> const & hits,
                           View<char> const query_head,
                           View<char> const qsequence,
                           View<char> const qsequence_rc,
                           int64_t const qsize) -> void
{
  std::lock_guard<std::mutex> const lock(state.mutex_output);
  auto const qseqlen = static_cast<int>(qsequence.size());

  /* show results: the hits --maxhits keeps, of which the per-hit writers below
     report the --top_hits_only prefix. The whole clamp is named once here
     rather than spelled at each of the four writers that used it. */
  auto const to_report = make_view(hits)
    .first(static_cast<std::size_t>(std::min(state.parameters.opt_maxhits,
                                             static_cast<int64_t>(hits.size()))));

  if (state.fp_alnout != nullptr)
    {
      results_show_alnout(state.fp_alnout.get(),
                          to_report,
                          query_head,
                          qsequence,
                          state.db,
                          state.parameters);
    }

  if (state.fp_lcaout != nullptr)
    {
      results_show_lcaout(state.fp_lcaout.get(),
                          to_report,
                          query_head,
                          state.db,
                          state.parameters);
    }

  if (state.fp_samout != nullptr)
    {
      results_show_samout(state.fp_samout.get(),
                          to_report,
                          query_head,
                          qsequence,
                          qsequence_rc,
                          state.db,
                          state.parameters);
    }

  if (not to_report.empty())  // hits.size() >=1 and <= opt_maxhits
    {
      auto const top = top_hits(to_report, state.parameters.opt_top_hits_only != 0);

      if (needs_otu_table(state.parameters))
        {
          state.otutable.add(query_head,
                       state.db.header_view(static_cast<uint64_t>(hits[0].target)),
                       qsize);
        }

      PerHitOutputFiles per_hit_files;
      per_hit_files.fastapairs = state.fp_fastapairs.get();
      per_hit_files.qsegout = state.fp_qsegout.get();
      per_hit_files.tsegout = state.fp_tsegout.get();
      per_hit_files.uc = state.fp_uc.get();
      per_hit_files.userout = state.fp_userout.get();
      per_hit_files.blast6out = state.fp_blast6out.get();
      results_show_hits(per_hit_files, top, query_head, qsequence, qsequence_rc,
                        state.db, state.parameters);
    }
  else
    {
      if (needs_otu_table(state.parameters))
        {
          state.otutable.add(query_head,
                       View<char>{},
                       qsize);
        }

      NoHitOutputFiles no_hit_files;
      no_hit_files.uc = state.fp_uc.get();
      no_hit_files.userout = state.fp_userout.get();
      no_hit_files.blast6out = state.fp_blast6out.get();
      results_show_no_hit(no_hit_files, query_head, qsequence, qsequence_rc,
                          qseqlen, state.db, state.parameters);
    }

  auto const matched = not hits.empty();
  results_show_matched_query(matched ? state.fp_matched.get() : state.fp_notmatched.get(),
                             matched ? state.count_matched : state.count_notmatched,
                             query_head, qsequence,
                             static_cast<uint64_t>(qsize), state.parameters);

  /* update matching db sequences */
  for (auto const & hit : hits) {
    if (hit.accepted or hit.weak) {
      state.dbmatched[static_cast<std::size_t>(hit.target)] += state.parameters.opt_sizein ? static_cast<uint64_t>(qsize) : 1;
    }
  }
}


static auto search_query(struct search_cli_state_s & state, uint64_t const t) -> int
{
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;

  std::array<struct searchinfo_s *, 2> const strands
    {{&si_plus[t], si_minus.empty() ? nullptr : &si_minus[t]}};
  for (auto * const si : make_view(strands)
         .first(static_cast<std::size_t>(number_of_strands(state.parameters.opt_strand))))
    {
      /* mask query */
      apply_masking(si->qsequence, state.parameters.opt_qmask, state.parameters);

      /* perform search */
      search_onequery(si, state.parameters.opt_qmask);
    }

  std::vector<struct hit> hits;

  search_joinhits(&si_plus[t],
                  state.parameters.opt_strand ? &si_minus[t] : nullptr,
                  hits);

  auto const qsequence = View<char>{si_plus[t].qsequence};
  auto const qsequence_rc = state.parameters.opt_strand
    ? View<char>{si_minus[t].qsequence}
    : View<char>{};

  search_output_results(state,
                        hits,
                        si_plus[t].query_head,
                        qsequence,
                        qsequence_rc,
                        si_plus[t].qsize);

  /* alignment strings (hit.nwalignment) are std::string and free themselves */

  return static_cast<int>(hits.size());
}


static auto search_thread_run(struct search_cli_state_s & state, uint64_t const t) -> void
{
  auto * const query_fastx_h = state.query_fastx_h;
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;

  int query_no = 0;
  int64_t qsize = 0;
  uint64_t progress = 0;

  auto const has_work_to_claim = [&]() -> bool {
    if (not query_fastx_h->next(
                       header_truncation(state.parameters.opt_notrunclabels),
                       Mapping::none))
      {
        return false;
      }

    query_no = static_cast<int>(query_fastx_h->get_seqno());
    qsize = query_fastx_h->get_abundance();

    populate_si(si_plus[t],
                query_fastx_h->header_view(),
                query_fastx_h->sequence_view(),
                query_no,
                qsize,
                0);

    /* get progress as amount of input file read */
    progress = query_fastx_h->get_position();
    return true;
  };

  auto const process_query = [&]() -> void {
    if (state.parameters.opt_strand)
      {
        populate_si(si_minus[t],
                    si_plus[t].query_head,
                    View<char>{si_plus[t].qsequence},
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
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;
  int const seqcount = state.seqcount;
  int const tophits = state.tophits;

  /* init per-thread search state before the workers start. Both vectors are
     resize()d to opt_threads (si_minus stays empty unless --strand both), so
     walking each one whole is the same set of calls the index made -- and an
     empty si_minus needs no emptiness test of its own. */
  for (auto & si : si_plus)
    {
      search_thread_init(si, seqcount, tophits, state.effective_parameters, state.dbindex, state.db);
    }
  for (auto & si : si_minus)
    {
      search_thread_init(si, seqcount, tophits, state.effective_parameters, state.dbindex, state.db);
    }

  /* run the worker pool over the input file */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(state.parameters.opt_threads),
                              [&state](uint64_t const t) -> void
                              { search_thread_run(state, t); });
    threadrunner.run();
  }

  /* clean up per-thread search state */
  for (auto & si : si_plus)
    {
      search_thread_exit(si);
    }
  for (auto & si : si_minus)
    {
      search_thread_exit(si);
    }
}


static auto search_prep(struct search_cli_state_s & state) -> void
{
  /* The configuration the search engine reads, as opposed to the one the user
     gave: resolved here, before the database is read, because the index is
     built against it too. For --usearch_global it is the user's, plus the
     database-size clamp applied at the end of this function. */
  state.effective_parameters = state.parameters;

  if (state.prefilter == Prefilter::none)
    {
      /* --search_global compares every query against every target, so the
         three options that bound the pre-filtered search are ignored: the
         word threshold that selects candidates, and the two counters that
         stop the scan early. Overriding them here, on the engine's copy
         rather than on the user's, is what "ignored" means for this command
         -- the option stays accepted so a command line moved over from
         --usearch_global keeps working, and it simply has no effect.

         Zero is the "no bound" spelling for all three: minwordmatches 0 asks
         a candidate for no shared words at all, and maxaccepts/maxrejects 0
         are clamped up to the database size below. The index is built from
         this same copy, so dbindex.minwordmatches agrees with what
         search_topscores asserts against. */
      state.effective_parameters.opt_minwordmatches = 0;
      state.effective_parameters.opt_maxaccepts = 0;
      state.effective_parameters.opt_maxrejects = 0;
    }

  /* open output files */

  state.fp_alnout = open_optional_output_file(state.parameters.opt_alnout, OutputOption{"--alnout"});
  if (state.fp_alnout != nullptr)
    {
      fprint(state.fp_alnout.get(), make_view(state.parameters.runtime.command_line));
      fprint(state.fp_alnout.get(), '\n');
      fprint(state.fp_alnout.get(), make_view(state.parameters.runtime.prog_header));
      fprint(state.fp_alnout.get(), '\n');
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
      udb_read(state.parameters.opt_db, UdbUse::search, state.dbindex, state.db, state.effective_parameters);
      results_show_samheader(state.fp_samout.get(), state.parameters.opt_db, state.db, state.parameters);
      // memory-intensive: the entire database is now held in memory
      state.seqcount = static_cast<int>(state.db.getsequencecount());
    }
  else
    {
      state.db.read(state.parameters.opt_db, 0, state.parameters);
      results_show_samheader(state.fp_samout.get(), state.parameters.opt_db, state.db, state.parameters);
      apply_masking(state.db, state.parameters.opt_dbmask, state.parameters);
      // memory-intensive: the entire database is now held in memory
      state.seqcount = static_cast<int>(state.db.getsequencecount());
      state.dbindex.prepare(state.parameters.opt_dbmask, state.db, state.effective_parameters);
      state.dbindex.add_all_sequences(state.parameters.opt_dbmask, state.db, state.effective_parameters);
    }

  /* tophits = the maximum number of hits we need to store */

  /* Clamp maxrejects/maxaccepts to the database size (0 or "> seqcount" means
     "all"). Apply the clamp to a local Parameters copy that is threaded to the
     workers via si->parameters, rather than mutating the shared config globals
     (E1 trap-writer): the shared searchcore reads the clamped values through
     si->parameters. The copy itself was made at the top of this function,
     because the index is built from it too. */
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

}


auto run_global_search(struct Parameters const & parameters,
                       Prefilter const prefilter) -> void
{
  /* Per-invocation state, owned here and threaded through the worker pool and
     the output helper (E4). Aliased by reference so the long body below reads
     unchanged; the workers receive `state`, not file-static globals. */
  struct search_cli_state_s state(parameters, prefilter);
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;
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

  dbmatched.assign(static_cast<size_t>(seqcount), 0);

  /* prepare reading of queries */
  qmatches = 0;
  qmatches_abundance = 0;
  queries = 0;
  queries_abundance = 0;
  auto const query_fastx_h = fastx_open(parameters.input_filename, parameters);
  state.query_fastx_h = query_fastx_h.get();  // workers borrow the raw handle

  /* The query file is parsed inside the worker threads (search_thread_run).
     Defer parse errors so a malformed query stops the pool cooperatively
     instead of calling fatal()/std::exit() from a worker while siblings are
     writing output (CC3); reported below from the main thread after join. */
  query_fastx_h->enable_deferred_errors();

  /* allocate memory for thread info */
  si_plus.resize(static_cast<std::size_t>(parameters.opt_threads));
  if (parameters.opt_strand)
    {
      si_minus.resize(static_cast<std::size_t>(parameters.opt_threads));
    }

  {
    Progress progress("Searching", query_fastx_h->get_size(), parameters);
    state.progress = &progress;
    search_thread_worker_run(state);
  }

  /* all workers joined; report a deferred query parse error (CC3) from the
     main thread so it does not race a worker's output */
  if (query_fastx_h->get_error())
    {
      fatal(query_fastx_h->get_errmsg());
    }

  /* si_plus/si_minus are std::vector members of state (RAII) */

  query_fastx_h->report_stripped_warning(parameters);

  auto const match_counts = vsearch::MatchCounts{state.qmatches, state.queries,
                                                 state.qmatches_abundance,
                                                 state.queries_abundance};

  if (! parameters.opt_quiet)
    {
      vsearch::print_match_counts(stderr, match_counts, parameters.opt_sizein);
    }

  if (parameters.fp_log != nullptr)
    {
      vsearch::print_match_counts(parameters.fp_log, match_counts, parameters.opt_sizein);
    }


  // Add OTUs with no matches to OTU table
  if (needs_otu_table(parameters)) {
    for (int64_t i = 0; i < seqcount; i++) {
      if (dbmatched[static_cast<std::size_t>(i)] == 0U) {
        state.otutable.add(View<char>{}, state.db.header_view(static_cast<uint64_t>(i)), 0);
      }
    }
  }

  if (parameters.opt_biomout != nullptr)
    {
      state.otutable.print_biomout(state.fp_biomout.get(), state.parameters);
      state.fp_biomout.reset();
    }

  if (parameters.opt_otutabout != nullptr)
    {
      state.otutable.print_otutabout(state.fp_otutabout.get(), state.parameters);
      state.fp_otutabout.reset();
    }

  if (parameters.opt_mothur_shared_out != nullptr)
    {
      state.otutable.print_mothur_shared_out(state.fp_mothur_shared_out.get(), state.parameters);
      state.fp_mothur_shared_out.reset();
    }

  if ((parameters.opt_dbmatched != nullptr) || (parameters.opt_dbnotmatched != nullptr))
    {
      int count_dbmatched = 0;
      int count_dbnotmatched = 0;

      for (int64_t i = 0; i < seqcount; i++)
        {
          if (dbmatched[static_cast<std::size_t>(i)] != 0U)
            {
              count_dbmatched++;
              if (parameters.opt_dbmatched != nullptr)
                {
                  fasta_print_general(fp_dbmatched.get(),
                                      state.db.record(static_cast<uint64_t>(i)),
                                      OutputAnnotations{dbmatched[static_cast<std::size_t>(i)],
                                                        count_dbmatched},
                                      parameters);
                }
            }
          else
            {
              count_dbnotmatched++;
              if (parameters.opt_dbnotmatched != nullptr)
                {
                  fasta_print_general(fp_dbnotmatched.get(),
                                      state.db.record(static_cast<uint64_t>(i)),
                                      OutputAnnotations{state.db.getabundance(static_cast<uint64_t>(i)),
                                                        count_dbnotmatched},
                                      parameters);
                }
            }
        }
    }

  fp_dbmatched.reset();
  fp_dbnotmatched.reset();

  search_done(state);
}


/* --usearch_global: the heuristic global search. The word pre-filter picks
   and ranks the candidates, and --maxaccepts / --maxrejects stop the scan
   once enough of them have been settled. The exhaustive counterpart is
   --search_global (commands/search_global.cpp), which shares the run above
   and differs only in the candidate set. */
auto usearch_global(struct Parameters const & parameters) -> void
{
  run_global_search(parameters, Prefilter::kmer);
}
