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

/*

  Implements the Sintax algorithm as described in Robert Edgar's preprint:

  Robert Edgar (2016)
  SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences
  BioRxiv, 074161
  doi: https://doi.org/10.1101/074161

  Further details:

  https://www.drive5.com/usearch/manual/cmd_sintax.html


  Note that due to the lack of details in the description, this implementation
  in vsearch is surely somewhat different from the one in usearch.

*/

#include "vsearch.hpp"
#include "core/db.hpp"
#include "core/fastx.hpp"
#include "core/searchcore.hpp"
#include "arch/increment_counters.hpp"
#include "core/buffer_headroom.hpp"
#include "utils/progress.hpp"
#include "core/bitmap.hpp"
#include "core/dbindex.hpp"
#include "core/mask.hpp"
#include "core/minheap.hpp"
#include "core/tax.hpp"
#include "core/udb.hpp"
#include "core/unique.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/number_of_strands.hpp"
#include "utils/open_file.hpp"
#include "utils/taxonomic_fields.h"
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include "utils/random.hpp"
#include "utils/reverse_complement.hpp"
#include <algorithm>  // std::min, std::max
#include <array>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t
#include <cstring>  // std::memset, std::strncmp, std::strcpy
#include <mutex>  // std::mutex, std::lock_guard, std::unique_lock


constexpr auto subset_size = 32;
constexpr auto bootstrap_count = 100;

/* Per-invocation state for a sintax run — previously ten file-static globals.
   Folding them into a struct that sintax() owns and threads through the helper
   functions and workers makes the command reentrant and removes the shared
   mutable state (E4). subset_size / bootstrap_count stay file-scope constexpr:
   they are compile-time constants, not mutable state. */
struct sintax_state_s
{
  /* the run configuration, threaded through the helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter */
  struct Parameters const & parameters;
  struct Database db;  /* the sequence database this run owns (RAII); si->db points here */
  struct Dbindex dbindex;  /* the k-mer index this run owns (RAII); si->dbindex points here */
  struct searchinfo_s * si_plus = nullptr;
  struct searchinfo_s * si_minus = nullptr;
  int tophits = 0;   /* the maximum number of hits to keep */
  int seqcount = 0;  /* number of database sequences */
  fastx_handle query_fastx_h = nullptr;
  std::mutex mutex_input;   /* serializes query reads */
  std::mutex mutex_output;  /* serializes output + counter updates */
  std::FILE * fp_tabbedout = nullptr;
  int queries = 0;
  int classified = 0;

  Progress * progress = nullptr;  /* owner progress bar; worker updates it under the output lock */

  explicit sintax_state_s(struct Parameters const & params) : parameters(params) {}
};


static auto sintax_analyse(struct sintax_state_s & state,
                    char const * query_head,
                    int const strand,
                    int const * all_seqno,
                    int const count) -> void
{
  std::FILE * const fp_tabbedout = state.fp_tabbedout;

  std::array<int, tax_levels> level_matchcount {{}};
  std::array<int, tax_levels> level_best {{}};
  std::array<std::array<char const *, tax_levels>, bootstrap_count> cand_level_name_start {{}};
  std::array<std::array<int, tax_levels>, bootstrap_count> cand_level_name_len {{}};

  /* Check number of successful bootstraps, must be at least half */

  auto const is_enough = count >= (bootstrap_count + 1) / 2;

  if (is_enough)
    {
      /* Find the most common name at each taxonomic rank,
         but with the same names at higher ranks. */

      for (auto i = 0; i < count ; i++)
        {
          /* Split headers of all candidates by taxonomy ranks */

          auto const seqno = all_seqno[i];
          std::array<int, tax_levels> new_level_name_start {{}};
          std::array<int, tax_levels> new_level_name_len {{}};
          tax_split(seqno, new_level_name_start.data(), new_level_name_len.data(), state.db);
          cand_level_name_len[static_cast<std::size_t>(i)] = new_level_name_len;
          for (auto k = 0; k < tax_levels; k++)
            {
              cand_level_name_start[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                state.db.getheader(static_cast<uint64_t>(seqno)) + new_level_name_start[static_cast<std::size_t>(k)];
            }
        }

      std::array<bool, bootstrap_count> cand_included {{}};
      cand_included.fill(true);

      /* Count matching names among candidates */

      for (auto k = 0; k < tax_levels; k++)
        {
          auto const level = static_cast<std::size_t>(k);
          level_best[level] = -1;
          level_matchcount[level] = 0;

          std::array<int, bootstrap_count> cand_match {{}};
          cand_match.fill(-1);
          std::array<int, bootstrap_count> cand_matchcount {{}};

          for (auto i = 0; i < count ; i++) {
            auto const cand_i = static_cast<std::size_t>(i);
            if (cand_included[cand_i]) {
              for (auto j = 0; j <= i ; j++) {
                auto const cand_j = static_cast<std::size_t>(j);
                if (cand_included[cand_j])
                  {
                    /* check match at current level */
                    if ((cand_level_name_len[cand_i][level] == cand_level_name_len[cand_j][level]) &&
                        (std::strncmp(cand_level_name_start[cand_i][level],
                                           cand_level_name_start[cand_j][level],
                                           static_cast<std::size_t>(cand_level_name_len[cand_i][level])) == 0))
                      {
                        cand_match[cand_i] = j;
                        cand_matchcount[cand_j]++;
                        break; /* stop at first match */
                      }
                  }
              }
            }
          }

          for (auto i = 0; i < count ; i++) {
            auto const cand_i = static_cast<std::size_t>(i);
            if (cand_matchcount[cand_i] > level_matchcount[level])
              {
                level_best[level] = i;
                level_matchcount[level] = cand_matchcount[cand_i];
              }
          }

          for (auto i = 0; i < count; i++) {
            auto const cand_i = static_cast<std::size_t>(i);
            if (cand_match[cand_i] != level_best[level]) {
              cand_included[cand_i] = false;
            }
          }
        }
    }

  /* write to tabbedout file */
  std::lock_guard<std::mutex> const output_lock(state.mutex_output);
  std::fprintf(fp_tabbedout, "%s\t", query_head);

  state.queries++;

  if (is_enough)
    {
      state.classified++;

      auto comma = false;
      for (auto j = 0; j < tax_levels; j++)
        {
          auto const level = static_cast<std::size_t>(j);
          auto const best = static_cast<std::size_t>(level_best[level]);
          if (cand_level_name_len[best][level] > 0)
            {
              std::fprintf(fp_tabbedout,
                      "%s%c:%.*s(%.2f)",
                      (comma ? "," : ""),
                      taxonomic_fields[level],
                      cand_level_name_len[best][level],
                      cand_level_name_start[best][level],
                      1.0 * level_matchcount[level] / count);
              comma = true;
            }
        }

      std::fprintf(fp_tabbedout, "\t%c", (strand != 0) ? '-' : '+');

      if (state.parameters.opt_sintax_cutoff > 0.0)
        {
          std::fprintf(fp_tabbedout, "\t");
          auto comma_cutoff = false;
          for (auto j = 0; j < tax_levels; j++)
            {
              auto const level = static_cast<std::size_t>(j);
              auto const best = static_cast<std::size_t>(level_best[level]);
              if ((cand_level_name_len[best][level] > 0) &&
                  (1.0 * level_matchcount[level] / count >= state.parameters.opt_sintax_cutoff))
                {
                  std::fprintf(fp_tabbedout,
                          "%s%c:%.*s",
                          (comma_cutoff ? "," : ""),
                          taxonomic_fields[level],
                          cand_level_name_len[best][level],
                          cand_level_name_start[best][level]);
                  comma_cutoff = true;
                }
            }
        }
    }
  else
    {
      if (state.parameters.opt_sintax_cutoff > 0.0)
        {
          std::fprintf(fp_tabbedout, "\t\t");
        }
      else
        {
          std::fprintf(fp_tabbedout, "\t");
        }
    }

  std::fprintf(fp_tabbedout, "\n");
}


auto sintax_search_topscores(struct searchinfo_s * searchinfo,
                             SplitMix64 & rng,
                             struct Parameters const & parameters) -> void
{
  /*
    Count the number of kmer hits in each database sequence and select
    the database sequence with the highest number of matching kmers.
    If several sequences have equally many kmer matches, choose one of
    them according to the following rules: By default, choose the
    shortest. If two are equally short, choose the one that comes
    first in the database.  If the sintax_random option is in effect,
    ties will instead be chosen randomly.
  */

  /* count kmer hits in the database sequences */
  unsigned int const indexed_count = searchinfo->dbindex->getcount();

  /* zero counts */
  std::memset(searchinfo->kmers, 0, indexed_count * sizeof(count_t));

  for (auto i = 0U; i < searchinfo->kmersamplecount; i++)
    {
      unsigned int const kmer = searchinfo->kmersample[i];
      auto * bitmap = searchinfo->dbindex->getbitmap(kmer);

      if (bitmap != nullptr)
        {
#ifdef __x86_64__
          if (parameters.ssse3_present != 0)
            {
              increment_counters_from_bitmap_ssse3(searchinfo->kmers,
                                                   bitmap, indexed_count);
            }
          else
            {
              increment_counters_from_bitmap_sse2(searchinfo->kmers,
                                                  bitmap, indexed_count);
            }
#else
          increment_counters_from_bitmap(searchinfo->kmers, bitmap, indexed_count);
#endif
        }
      else
        {
          auto const * list = searchinfo->dbindex->getmatchlist(kmer);
          auto const count = searchinfo->dbindex->getmatchcount(kmer);
          for (auto j = 0U; j < count; j++)
            {
              searchinfo->kmers[list[j]]++;
            }
        }
    }

  auto tophit_count = 0U;

  elem_t best;
  best.count = 0;
  best.seqno = 0;
  best.length = 0;

  for (auto i = 0U; i < indexed_count; i++)
    {
      count_t const count = searchinfo->kmers[i];
      auto const seqno = searchinfo->dbindex->getmapping(i);
      unsigned int const length = static_cast<unsigned int>(searchinfo->db->getsequencelen(seqno));

      if (count > best.count)
        {
          best.count = count;
          best.seqno = seqno;
          best.length = length;
          tophit_count = 1;
        }
      else if (count == best.count)
        {
          if (parameters.opt_sintax_random)
            {
              tophit_count++;
              if (random_bounded(rng, tophit_count) == 0)
                {
                  best.seqno = seqno;
                  best.length = length;
                }
            }
          else
            {
              if (length < best.length)
                {
                  best.seqno = seqno;
                  best.length = length;
                }
              else if (length == best.length)
                {
                  best.seqno = std::min(seqno, best.seqno);
                }
            }
        }
    }

  searchinfo->m.clear();
  if (best.count > 1) {
    searchinfo->m.add(best);
  }
}


static auto sintax_query(struct sintax_state_s & state, uint64_t const t) -> void
{
  struct searchinfo_s * const si_plus = state.si_plus;
  struct searchinfo_s * const si_minus = state.si_minus;

  std::array<std::array<int, bootstrap_count>, 2> all_seqno {{}};
  std::array<int, 2> boot_count = {0, 0};
  std::array<unsigned int, 2> best_count = {0, 0};
  int const qseqlen = si_plus[t].qseqlen;
  char const * query_head = si_plus[t].query_head;

  /* Per-query RNG: seed from the global base seed and this query's input
     number, so the random subsampling and tie-breaking are reproducible
     regardless of how many threads process the queries or in what order. */
  SplitMix64 rng(random_substream_seed(random_base_seed(),
                                       static_cast<uint64_t>(si_plus[t].query_no)));

  Bitmap b(static_cast<unsigned int>(qseqlen));

  for (auto s = 0; s < number_of_strands(state.parameters.opt_strand); s++)
    {
      struct searchinfo_s * si = (s != 0) ? si_minus + t : si_plus + t;

      /* perform search */

      auto kmersamplecount = 0U;
      unsigned int const * kmersample = nullptr;

      /* find unique kmers at dbindex.wordlength, the effective index width (set
         by Dbindex::prepare for a FASTA db, or udb_read for a UDB db); reading
         parameters.opt_wordlength would use the wrong width against a UDB index. */
      si->uh.count(static_cast<int>(si->dbindex->wordlength),
                   si->qseqlen, si->qsequence,
                   &kmersamplecount, &kmersample, Masking::none);

      /* perform 100 bootstraps */

      if (kmersamplecount >= subset_size)
        {
          for (auto i = 0; i < bootstrap_count ; i++)
            {
              /* subsample 32 kmers */
              std::array<unsigned int, subset_size> kmersample_subset {{}};
              auto subsamples = 0;
              b.reset_all();
              for (auto j = 0; j < subset_size ; j++)
                {
                  int64_t const x = static_cast<int64_t>(random_bounded(rng, kmersamplecount));
                  if (not b.is_set(static_cast<unsigned int>(x)))
                    {
                      kmersample_subset[static_cast<std::size_t>(subsamples++)] = kmersample[x];
                      b.set(static_cast<unsigned int>(x));
                    }
                }

              si->kmersamplecount = static_cast<unsigned int>(subsamples);
              si->kmersample = kmersample_subset.data();

              sintax_search_topscores(si, rng, state.parameters);

              if (! si->m.is_empty())
                {
                  auto const e = si->m.pop_last();

                  auto const strand_idx = static_cast<std::size_t>(s);
                  auto & boot = boot_count[strand_idx];
                  all_seqno[strand_idx][static_cast<std::size_t>(boot++)] =
                    static_cast<int>(e.seqno);

                  best_count[strand_idx] = std::max(e.count, best_count[strand_idx]);
                }
            }
        }
    }

  auto best_strand = 0;

  if (not state.parameters.opt_strand)
    {
      best_strand = 0;
    }
  else
    {
      if (best_count[0] > best_count[1])
        {
          best_strand = 0;
        }
      else if (best_count[1] > best_count[0])
        {
          best_strand = 1;
        }
      else
        {
          if (boot_count[0] >= boot_count[1])
            {
              best_strand = 0;
            }
          else
            {
              best_strand = 1;
            }
        }
    }

  sintax_analyse(state,
                 query_head,
                 best_strand,
                 all_seqno[static_cast<std::size_t>(best_strand)].data(),
                 boot_count[static_cast<std::size_t>(best_strand)]);
}


static auto sintax_thread_run(struct sintax_state_s & state, uint64_t const t) -> void
{
  std::mutex & mutex_input = state.mutex_input;
  std::mutex & mutex_output = state.mutex_output;
  fastx_handle const query_fastx_h = state.query_fastx_h;
  struct searchinfo_s * const si_plus = state.si_plus;
  struct searchinfo_s * const si_minus = state.si_minus;

  uint64_t progress = 0;

  auto const has_work_to_claim = [&]() -> bool {
    if (not fastx_next(query_fastx_h,
                       not state.parameters.opt_notrunclabels,
                       chrmap_no_change()))
      {
        /* End of input, or a deferred parse error was recorded (CC3):
           fastx_next() returns false in both cases, so the worker stops
           here cooperatively. The error, if any, is reported by sintax()
           from the main thread after the pool joins. */
        return false;
      }

    auto const * qhead = fastx_get_header(query_fastx_h);
    int const query_head_len = static_cast<int>(fastx_get_header_length(query_fastx_h));
    auto const * qseq = fastx_get_sequence(query_fastx_h);
    int const qseqlen = static_cast<int>(fastx_get_sequence_length(query_fastx_h));
    int const query_no = static_cast<int>(fastx_get_seqno(query_fastx_h));
    int64_t const qsize = fastx_get_abundance(query_fastx_h);

    for (auto s = 0; s < number_of_strands(state.parameters.opt_strand); s++)
      {
        struct searchinfo_s * si = (s != 0) ? si_minus + t : si_plus + t;

        si->query_head_len = query_head_len;
        si->qseqlen = qseqlen;
        si->query_no = query_no;
        si->qsize = qsize;
        si->strand = s;

        /* allocate more memory for the sequence, if necessary */

        if (si->qseqlen + 1 > si->seq_alloc)
          {
            si->seq_alloc = si->qseqlen + buffer_headroom;
            si->qsequence = static_cast<char *>(
              xrealloc(si->qsequence, static_cast<size_t>(si->seq_alloc)));
          }
      }

    /* plus strand: copy header (into owned storage, view points at it) and sequence */
    si_plus[t].query_head_v.resize(static_cast<std::size_t>(query_head_len) + 1);
    std::strcpy(si_plus[t].query_head_v.data(), qhead);
    si_plus[t].query_head = si_plus[t].query_head_v.data();
    std::strcpy(si_plus[t].qsequence, qseq);

    /* get progress as amount of input file read */
    progress = fastx_get_position(query_fastx_h);
    return true;
  };

  auto const process_query = [&]() {
    /* minus strand: copy header and reverse complementary sequence */
    if (state.parameters.opt_strand)
      {
        si_minus[t].query_head_v = si_plus[t].query_head_v;
        si_minus[t].query_head = si_minus[t].query_head_v.data();
        reverse_complement(si_minus[t].qsequence,
                           si_plus[t].qsequence,
                           si_plus[t].qseqlen);
      }

    sintax_query(state, t);

    /* lock mutex for update of global data and output */
    std::lock_guard<std::mutex> const output_lock(mutex_output);

    /* show progress */
    state.progress->update(progress);
  };

  run_worker_loop(mutex_input, has_work_to_claim, process_query);
}


static auto sintax_thread_init(struct sintax_state_s const & state, struct searchinfo_s * si) -> void
{
  /* thread specific initialiation */
  si->parameters = &state.parameters;  /* searchcore reads config through the si (E1) */
  si->dbindex = &state.dbindex;  /* searchcore reads the k-mer index through the si */
  si->db = &state.db;  /* searchcore reads the sequences through the si */
  /* si->uh (a Uniquer value member) is ready to use as default-constructed */
  si->kmers = static_cast<count_t *>(xmalloc((static_cast<size_t>(state.seqcount) * sizeof(count_t)) + 32));
  si->m = Minheap(state.tophits);
  si->hits = nullptr;
  si->qsize = 1;
  si->query_head = nullptr;
  si->seq_alloc = 0;
  si->qsequence = nullptr;
  si->nw = nullptr;
  si->s.reset();
}


static auto sintax_thread_exit(struct searchinfo_s * searchinfo) -> void
{
  /* thread specific clean up */
  searchinfo->uh = Uniquer();
  searchinfo->m = Minheap();
  xfree(searchinfo->kmers);
  /* query_head is a view; its owned storage query_head_v frees itself */
  if (searchinfo->qsequence != nullptr)
    {
      xfree(searchinfo->qsequence);
    }
}


static auto sintax_thread_worker_run(struct sintax_state_s & state) -> void
{
  struct searchinfo_s * const si_plus = state.si_plus;
  struct searchinfo_s * const si_minus = state.si_minus;

  /* init per-thread search state before the workers start */
  for (auto t = 0; t < state.parameters.opt_threads; t++)
    {
      sintax_thread_init(state, si_plus + t);
      if (si_minus != nullptr)
        {
          sintax_thread_init(state, si_minus + t);
        }
    }

  /* run the worker pool over the input file */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(state.parameters.opt_threads),
                              [&state](uint64_t const t)
                              { sintax_thread_run(state, t); });
    threadrunner.run();
  }

  /* clean up per-thread search state */
  for (auto t = 0; t < state.parameters.opt_threads; t++)
    {
      sintax_thread_exit(si_plus + t);
      if (si_minus != nullptr)
        {
          sintax_thread_exit(si_minus + t);
        }
    }
}


auto sintax(struct Parameters const & parameters) -> void
{
  /* Per-invocation state, owned here and threaded through the workers (E4).
     Aliased by reference so the body below reads unchanged; the workers
     receive `state`, not file-static globals. */
  struct sintax_state_s state(parameters);
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;
  auto & tophits = state.tophits;
  auto & seqcount = state.seqcount;
  auto & query_fastx_h = state.query_fastx_h;
  int & queries = state.queries;
  int & classified = state.classified;

  /* tophits = the maximum number of hits we need to store */

  tophits = 1;

  /* open output files */

  if (parameters.opt_db == nullptr)
    {
      fatal("No database file specified with --db");
    }

  auto const output_handle = open_mandatory_output_file(parameters.opt_tabbedout, OutputOption{"--tabbedout"});
  state.fp_tabbedout = output_handle.get();

  /* check if db may be an UDB file */

  auto const is_udb = udb_detect_isudb(parameters.opt_db);

  if (is_udb)
    {
      udb_read(parameters.opt_db, true, true, state.dbindex, state.db, parameters);
    }
  else
    {
      state.db.read(parameters.opt_db, 0, parameters);
    }

  seqcount = static_cast<int>(state.db.getsequencecount());

  if (! is_udb)
    {
      state.dbindex.prepare(1, parameters.opt_dbmask, state.db, parameters);
      state.dbindex.add_all_sequences(parameters.opt_dbmask, state.db, parameters);
    }

  /* prepare reading of queries */

  query_fastx_h = fastx_open(parameters.opt_sintax, parameters);

  /* The query file is parsed inside the worker threads (see
     sintax_thread_run). Enable deferred error reporting so a malformed
     query records the error and stops the pool cooperatively, rather than
     calling fatal()/std::exit() from a worker thread while siblings are
     still writing output (CC3). The error is reported below, from the
     main thread, after the pool has joined. */
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

  /* run */

  {
    Progress progress("Classifying sequences", fastx_get_size(query_fastx_h), parameters);
    state.progress = &progress;
    sintax_thread_worker_run(state);
  }

  /* All workers have joined. If one hit a malformed query, report it now
     from the main thread (single-threaded) so the message is emitted and
     the process exits without racing any worker (CC3). */
  if (fastx_get_error(query_fastx_h))
    {
      fatal("%s", fastx_get_errmsg(query_fastx_h));
    }

  if (! parameters.opt_quiet)
    {
      std::fprintf(stderr, "Classified %d of %d sequences", classified, queries);
      if (queries > 0)
        {
          std::fprintf(stderr, " (%.2f%%)", 100.0 * classified / queries);
        }
      std::fprintf(stderr, "\n");
    }

  if (parameters.opt_log != nullptr)
    {
      std::fprintf(parameters.fp_log, "Classified %d of %d sequences", classified, queries);
      if (queries > 0)
        {
          std::fprintf(parameters.fp_log, " (%.2f%%)", 100.0 * classified / queries);
        }
      std::fprintf(parameters.fp_log, "\n");
    }

  /* clean up */

  delete [] si_plus;
  if (si_minus != nullptr)
    {
      delete [] si_minus;
    }

  fastx_close(query_fastx_h, parameters);

  state.dbindex.clear();
  state.db.clear();
}
