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

#include "vsearch.hpp"
#include "utils/progress.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/align_simd.hpp"
#include "core/linmemalign.hpp"
#include "core/mask.hpp"
#include "core/db.hpp"
#include "core/fasta.hpp"
#include "core/results.hpp"
#include "core/searchcore.hpp"
#include "utils/open_file.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/threads.hpp"
#include "utils/view.hpp"  // View
#include "utils/worker_loop.hpp"
#include <algorithm>  // std::min, std::max
#include <cstddef>
#include <cstdint>  // int64_t
#include <cstdio>  // std::fprintf, std::FILE, std:fclose, std::size_t
#include <iterator>  // std::next
#include <limits>
#include <string>
#include <utility>  // std::move
#include <mutex>  // std::mutex, std::lock_guard, std::unique_lock
#include <vector>


/* Per-invocation state for an allpairs_global run — previously eighteen
   file-static globals. Folding them into a struct that allpairs_global() owns
   and threads through the output helper and the worker pool makes the command
   reentrant and removes the shared mutable state (E4). */
struct allpairs_state_s
{
  /* the run configuration, threaded through the helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter */
  struct Parameters const & parameters;
  struct Database db;  /* the sequence database this run owns (RAII); si->db points here */
  int seqcount = 0;         /* number of database sequences */
  std::mutex mutex_input;   /* serializes query reads */
  std::mutex mutex_output;  /* serializes output + counter updates */
  int qmatches = 0;
  int queries = 0;
  int64_t progress = 0;
  FILE * fp_alnout = nullptr;
  FILE * fp_samout = nullptr;
  FILE * fp_userout = nullptr;
  FILE * fp_blast6out = nullptr;
  FILE * fp_uc = nullptr;
  FILE * fp_fastapairs = nullptr;
  FILE * fp_matched = nullptr;
  FILE * fp_notmatched = nullptr;
  FILE * fp_qsegout = nullptr;
  FILE * fp_tsegout = nullptr;
  int count_matched = 0;
  int count_notmatched = 0;

  Progress * progress_bar = nullptr;  /* owner progress bar; worker updates it under output_lock (state.progress is the counter) */

  explicit allpairs_state_s(struct Parameters const & params) : parameters(params) {}
};


namespace {
inline auto allpairs_hit_compare_typed(struct hit const & lhs, struct hit const & rhs) -> int
{
  // high id, then low id
  // early target, then late target

  if (lhs.id > rhs.id)
    {
      return -1;
    }
  if (lhs.id < rhs.id)
    {
      return +1;
    }
  if (lhs.target < rhs.target)
    {
      return -1;
    }
  if (lhs.target > rhs.target)
    {
      return +1;
    }
  return 0;
}
/* The state is passed whole rather than its two int counters, which would be
   two adjacent swappable arguments. */
auto print_match_count(std::FILE * output_stream,
                       struct allpairs_state_s const & state) -> void
{
  fprint(output_stream, "Matching query sequences: ");
  fprint_integer(output_stream, state.qmatches);
  fprint(output_stream, " of ");
  fprint_integer(output_stream, state.queries);
  if (state.queries > 0)
    {
      fprint(output_stream, " (");
      std::fprintf(output_stream, "%.2f", 100.0 * state.qmatches / state.queries);
      fprint(output_stream, "%)");
    }
  fprint(output_stream, '\n');
}

}  // anonymous namespace


static auto allpairs_output_results(struct allpairs_state_s & state,
                             View<struct hit> const hits,
                             View<char> const query_head,
                             View<char> const qsequence,
                             View<char> const qsequence_rc) -> void
{
  auto const qseqlen = static_cast<int>(qsequence.size());

  /* show results: the hits --maxhits keeps, of which the per-hit writers below
     report the --top_hits_only prefix. The whole clamp is named once here
     rather than spelled at each of the three writers that used it. */
  auto const to_report = hits
    .first(static_cast<std::size_t>(std::min(state.parameters.opt_maxhits,
                                             static_cast<int64_t>(hits.size()))));

  if (state.fp_alnout != nullptr)
    {
      results_show_alnout(state.fp_alnout,
                          to_report,
                          query_head,
                          qsequence,
                          state.db,
                          state.parameters);
    }

  if (state.fp_samout != nullptr)
    {
      results_show_samout(state.fp_samout,
                          to_report,
                          query_head,
                          qsequence,
                          qsequence_rc,
                          state.db,
                          state.parameters);
    }

  if (not to_report.empty())
    {
      auto const top = top_hits(to_report, state.parameters.opt_top_hits_only != 0);

      /* indexed rather than a range-for: --uc reports the best hit only,
         unless --uc_allhits, so the position is part of the output */
      for (std::size_t t = 0; t < top.size(); ++t)
        {
          struct hit const * hp = &top[t];

          if (state.fp_fastapairs != nullptr)
            {
              results_show_fastapairs_one(state.fp_fastapairs,
                                          *hp,
                                          query_head,
                                          qsequence,
                                          qsequence_rc,
                                          state.db,
                                          state.parameters);
            }

          if (state.fp_qsegout != nullptr)
            {
              results_show_qsegout_one(state.fp_qsegout,
                                       *hp,
                                       query_head,
                                       qsequence,
                                       qsequence_rc,
                                       state.parameters);
            }

          if (state.fp_tsegout != nullptr)
            {
              results_show_tsegout_one(state.fp_tsegout,
                                       *hp,
                                       state.db,
                                       state.parameters);
            }

          if ((state.fp_uc != nullptr) and ((t == 0) or state.parameters.opt_uc_allhits))
            {
              results_show_uc_one(state.fp_uc,
                                  hp,
                                  query_head,
                                  qseqlen,
                                  hp->target,
                                  state.db,
                                  state.parameters,
                                  PerfectMatch::whole_alignment);
            }

          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qsequence_rc,
                                       state.db,
                                       state.parameters);
            }

          if (state.fp_blast6out != nullptr)
            {
              results_show_blast6out_one(state.fp_blast6out,
                                         hp,
                                         query_head,
                                         qseqlen,
                                         state.db);
            }
        }
    }
  else
    {
      if (state.fp_uc != nullptr)
        {
          results_show_uc_one(state.fp_uc,
                              nullptr,
                              query_head,
                              qseqlen,
                              0,
                              state.db,
                              state.parameters,
                              PerfectMatch::whole_alignment);
        }

      if (state.parameters.opt_output_no_hits != 0)
        {
          if (state.fp_userout != nullptr)
            {
              results_show_userout_one(state.fp_userout,
                                       nullptr,
                                       query_head,
                                       qsequence,
                                       qsequence_rc,
                                       state.db,
                                       state.parameters);
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
      if (state.parameters.opt_matched != nullptr)
        {
          fasta_print_general(state.fp_matched,
                              nullptr,
                              qsequence,
                              query_head,
                              OutputAnnotations{0, state.count_matched},
                              state.parameters);
        }
    }
  else
    {
      ++state.count_notmatched;
      if (state.parameters.opt_notmatched != nullptr)
        {
          fasta_print_general(state.fp_notmatched,
                              nullptr,
                              qsequence,
                              query_head,
                              OutputAnnotations{0, state.count_notmatched},
                              state.parameters);
        }
    }
}


static auto allpairs_thread_run(struct allpairs_state_s & state, uint64_t const t) -> void
{
  (void) t;

  struct searchinfo_s searchinfo;
  searchinfo.parameters = &state.parameters;  /* searchcore reads config through the si (E1) */
  searchinfo.db = &state.db;  /* searchcore reads the sequences through the si */

  searchinfo.hits_v.resize(static_cast<std::size_t>(state.seqcount));

  searchinfo.s.reset(search16_init(state.parameters.opt_match,
                        state.parameters.opt_mismatch,
                        state.parameters.opt_gap_open_query_left,
                        state.parameters.opt_gap_open_target_left,
                        state.parameters.opt_gap_open_query_interior,
                        state.parameters.opt_gap_open_target_interior,
                        state.parameters.opt_gap_open_query_right,
                        state.parameters.opt_gap_open_target_right,
                        state.parameters.opt_gap_extension_query_left,
                        state.parameters.opt_gap_extension_target_left,
                        state.parameters.opt_gap_extension_query_interior,
                        state.parameters.opt_gap_extension_target_interior,
                        state.parameters.opt_gap_extension_query_right,
                        state.parameters.opt_gap_extension_target_right,
                        state.parameters.opt_n_mismatch));


  struct Scoring const scoring = scoring_from_options(state.parameters);


  LinearMemoryAligner lma(scoring);


  /* allocate memory for alignment results */
  auto const maxhits = static_cast<std::size_t>(state.seqcount);
  std::vector<unsigned int> pseqnos(maxhits);
  std::vector<CELL> pscores(maxhits);
  std::vector<unsigned short> paligned(maxhits);
  std::vector<unsigned short> pmatches(maxhits);
  std::vector<unsigned short> pmismatches(maxhits);
  std::vector<unsigned short> pgaps(maxhits);
  std::vector<std::string> pcigar(maxhits);
  std::vector<struct hit> finalhits(maxhits);

  int query_no = 0;

  auto const has_work_to_claim = [&]() -> bool {
    query_no = state.queries;
    if (query_no >= state.seqcount) { return false; }
    ++state.queries;
    return true;
  };

  auto const process_query = [&]() -> void {
    /* init search info */
    auto const query_no_u = static_cast<uint64_t>(query_no);
    searchinfo.query_no = query_no;
    searchinfo.qsize = static_cast<int64_t>(state.db.getabundance(query_no_u));
    searchinfo.query_head = state.db.header_view(query_no_u);
    searchinfo.qsequence = state.db.mutable_sequence(query_no_u);
    searchinfo.rejects = 0;
    searchinfo.accepts = 0;
    searchinfo.hit_count = 0;

    for (int target = searchinfo.query_no + 1; target < state.seqcount; target++)
      {
        if ((state.parameters.opt_acceptall != 0) or search_acceptable_unaligned(searchinfo, target))
          {
            pseqnos[static_cast<std::size_t>(searchinfo.hit_count)] = static_cast<unsigned int>(target);
            ++searchinfo.hit_count;
          }
      }

    if (searchinfo.hit_count != 0)
      {
        /* perform alignments */

        search16_qprep(searchinfo.s.get(), View<char>{searchinfo.qsequence});

        /* the hits accumulated above, not the whole maxhits buffers */
        auto const found = static_cast<std::size_t>(searchinfo.hit_count);
        search16(searchinfo.s.get(),
                 make_view(pseqnos).first(found),
                 make_span(pscores).first(found),
                 make_span(paligned).first(found),
                 make_span(pmatches).first(found),
                 make_span(pmismatches).first(found),
                 make_span(pgaps).first(found),
                 make_span(pcigar).first(found),
                 state.db);

        /* convert to hit structure */
        for (std::size_t h = 0; h < static_cast<std::size_t>(searchinfo.hit_count); h++)
          {
            struct hit * hit = &searchinfo.hits_v[h];

            unsigned int const target = pseqnos[h];
            int64_t nwscore = pscores[h];

            std::string nwcigar;
            int64_t nwalignmentlength {0};
            int64_t nwmatches {0};
            int64_t nwmismatches {0};
            int64_t nwgaps {0};

            if (nwscore == std::numeric_limits<short>::max())
              {
                /* In case the SIMD aligner cannot align,
                   perform a new alignment with the
                   linear memory aligner */

                auto const tseq = state.db.sequence_view(target);
                auto const qseq = View<char>{searchinfo.qsequence};

                nwcigar = lma.align(qseq, tseq);
                auto const stats = lma.alignstats(nwcigar.c_str(), qseq, tseq);
                nwscore = stats.score;
                nwalignmentlength = stats.alignmentlength;
                nwmatches = stats.matches;
                nwmismatches = stats.mismatches;
                nwgaps = stats.gaps;
              }
            else
              {
                nwcigar = std::move(pcigar[h]);
                nwalignmentlength = paligned[h];
                nwmatches = pmatches[h];
                nwmismatches = pmismatches[h];
                nwgaps = pgaps[h];
              }

            hit->target = static_cast<int>(target);
            hit->strand = 0;
            hit->count = 0;

            hit->accepted = false;
            hit->rejected = false;
            hit->aligned = true;
            hit->weak = false;

            hit->nwscore = static_cast<int>(nwscore);
            hit->nwdiff = static_cast<int>(nwalignmentlength - nwmatches);
            hit->nwgaps = static_cast<int>(nwgaps);
            hit->nwindels = static_cast<int>(nwalignmentlength - nwmatches - nwmismatches);
            hit->nwalignmentlength = static_cast<int>(nwalignmentlength);
            hit->nwid = 100.0 * static_cast<double>(nwalignmentlength - hit->nwdiff) /
              static_cast<double>(nwalignmentlength);
            hit->nwalignment = std::move(nwcigar);  // owned cigar (empty means no alignment)
            hit->matches = static_cast<int>(nwalignmentlength - hit->nwdiff);
            hit->mismatches = hit->nwdiff - hit->nwindels;

            auto const dseqlen = static_cast<int>(state.db.getsequencelen(target));
            hit->shortest = std::min(static_cast<int>(searchinfo.qsequence.size()), dseqlen);
            hit->longest = std::max(static_cast<int>(searchinfo.qsequence.size()), dseqlen);

            /* trim alignment, compute numbers excluding terminal gaps */
            align_trim(*hit, state.parameters);

            /* test accept/reject criteria after alignment */
            if ((state.parameters.opt_acceptall != 0) or search_acceptable_aligned(searchinfo, *hit))
              {
                finalhits[static_cast<std::size_t>(searchinfo.accepts)] = *hit;
                ++searchinfo.accepts;
              }
          }

        /* sort the accepted hits */
        std::sort(finalhits.begin(),
                  std::next(finalhits.begin(), static_cast<std::ptrdiff_t>(searchinfo.accepts)),
                  [](struct hit const & lhs, struct hit const & rhs) -> bool {
                    return allpairs_hit_compare_typed(lhs, rhs) < 0;
                  });
      }

    /* lock mutex for update of global data and output */
    std::unique_lock<std::mutex> output_lock(state.mutex_output);

    /* output results */
    allpairs_output_results(state,
                            make_view(finalhits).first(static_cast<std::size_t>(searchinfo.accepts)),
                            searchinfo.query_head,
                            View<char>{searchinfo.qsequence},
                            View<char>{});

    /* update stats */
    if (searchinfo.accepts != 0)
      {
        ++state.qmatches;
      }

    /* show progress */
    state.progress += state.seqcount - query_no - 1;
    state.progress_bar->update(static_cast<uint64_t>(state.progress));

    output_lock.unlock();

    /* alignment strings (hit.nwalignment) are std::string and free themselves */
  };

  run_worker_loop(state.mutex_input, has_work_to_claim, process_query);

  searchinfo.s.reset();
}


static auto allpairs_thread_worker_run(struct allpairs_state_s & state) -> void
{
  /* run the worker pool; each worker keeps its own search state and
     processes queries until the shared counter is exhausted */
  ThreadRunner threadrunner(static_cast<std::size_t>(state.parameters.opt_threads),
                            [&state](uint64_t const t) -> void
                            { allpairs_thread_run(state, t); });
  threadrunner.run();
}


auto allpairs_global(struct Parameters const & parameters) -> void
{
  /* Per-invocation state, owned here and threaded through the worker pool (E4).
     Aliased by reference so the body below reads unchanged; the workers receive
     `state`, not file-static globals. */
  struct allpairs_state_s state(parameters);
  auto & fp_alnout = state.fp_alnout;
  auto & fp_samout = state.fp_samout;
  auto & fp_userout = state.fp_userout;
  auto & fp_blast6out = state.fp_blast6out;
  auto & fp_uc = state.fp_uc;
  auto & fp_fastapairs = state.fp_fastapairs;
  auto & fp_qsegout = state.fp_qsegout;
  auto & fp_tsegout = state.fp_tsegout;
  auto & fp_matched = state.fp_matched;
  auto & fp_notmatched = state.fp_notmatched;
  auto & seqcount = state.seqcount;
  auto & qmatches = state.qmatches;
  auto & queries = state.queries;
  auto & progress = state.progress;

  /* open output files */

  OutputFileHandle alnout_handle = open_optional_output_file(parameters.opt_alnout, OutputOption{"--alnout"});
  fp_alnout = alnout_handle.get();
  if (fp_alnout != nullptr)
    {
      fprint(fp_alnout, make_view(parameters.runtime.command_line));
      fprint(fp_alnout, '\n');
      fprint(fp_alnout, make_view(parameters.runtime.prog_header));
      fprint(fp_alnout, '\n');
    }

  OutputFileHandle samout_handle = open_optional_output_file(parameters.opt_samout, OutputOption{"--samout"});
  fp_samout = samout_handle.get();
  OutputFileHandle userout_handle = open_optional_output_file(parameters.opt_userout, OutputOption{"--userout"});
  fp_userout = userout_handle.get();
  OutputFileHandle blast6out_handle = open_optional_output_file(parameters.opt_blast6out, OutputOption{"--blast6out"});
  fp_blast6out = blast6out_handle.get();
  OutputFileHandle uc_handle = open_optional_output_file(parameters.opt_uc, OutputOption{"--uc"});
  fp_uc = uc_handle.get();
  OutputFileHandle fastapairs_handle = open_optional_output_file(parameters.opt_fastapairs, OutputOption{"--fastapairs"});
  fp_fastapairs = fastapairs_handle.get();
  OutputFileHandle qsegout_handle = open_optional_output_file(parameters.opt_qsegout, OutputOption{"--qsegout"});
  fp_qsegout = qsegout_handle.get();
  OutputFileHandle tsegout_handle = open_optional_output_file(parameters.opt_tsegout, OutputOption{"--tsegout"});
  fp_tsegout = tsegout_handle.get();
  OutputFileHandle matched_handle = open_optional_output_file(parameters.opt_matched, OutputOption{"--matched"});
  fp_matched = matched_handle.get();
  OutputFileHandle notmatched_handle = open_optional_output_file(parameters.opt_notmatched, OutputOption{"--notmatched"});
  fp_notmatched = notmatched_handle.get();

  state.db.read(parameters.input_filename, 0, parameters);

  results_show_samheader(fp_samout, parameters.input_filename, state.db, parameters);

  if (parameters.opt_qmask == Masking::dust)
    {
      dust_all(state.db, parameters);
    }
  else if ((parameters.opt_qmask == Masking::soft) and parameters.opt_hardmask)
    {
      hardmask_all(state.db);
    }

  // memory-intensive: the entire database is now held in memory

  seqcount = static_cast<int>(state.db.getsequencecount());

  /* prepare reading of queries */
  qmatches = 0;
  queries = 0;

  progress = 0;
  {
    Progress progress_bar("Aligning", static_cast<uint64_t>(std::max(int64_t{0}, (static_cast<int64_t>(seqcount)) * (static_cast<int64_t>(seqcount) - 1)) / 2), parameters);
    state.progress_bar = &progress_bar;
    allpairs_thread_worker_run(state);
  }

  if (not parameters.opt_quiet)
    {
      print_match_count(stderr, state);
    }

  if (parameters.fp_log != nullptr)
    {
      /* No trailing separator: this report is the last thing written to the
         log before LogFile's footer, which opens with its own newline
         (utils/logfile.cpp). The separator other reports add exists to part
         them from the report that follows; there is none here. */
      print_match_count(parameters.fp_log, state);
    }

  /* clean up, global */
  state.db.clear();
}
