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

#include "utils/view.hpp"
#include "utils/span.hpp"
#include "vsearch.hpp"
#include <memory>  // std::unique_ptr
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
#include "core/tax.hpp"  // TaxLevel, tax_levels, tax_split
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
#include "utils/print_view.hpp"  // fprint
#include <algorithm>  // std::copy_n, std::fill_n, std::max, std::max_element, std::min, std::transform
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::distance, std::next
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <mutex>  // std::mutex, std::lock_guard, std::unique_lock
#include <vector>  // std::vector


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
  RandomSeed seed;  /* per-run base seed for the reproducible per-query RNG;
                       resolved once here, read lock-free by every worker */
  struct Database db;  /* the sequence database this run owns (RAII); si->db points here */
  struct Dbindex dbindex;  /* the k-mer index this run owns (RAII); si->dbindex points here */
  std::vector<searchinfo_s> si_plus;
  std::vector<searchinfo_s> si_minus;  /* empty unless --strand both */
  int tophits = 0;   /* the maximum number of hits to keep */
  int seqcount = 0;  /* number of database sequences */
  fastx_handle query_fastx_h = nullptr;
  std::mutex mutex_input;   /* serializes query reads */
  std::mutex mutex_output;  /* serializes output + counter updates */
  std::FILE * fp_tabbedout = nullptr;
  int queries = 0;
  int classified = 0;

  Progress * progress = nullptr;  /* owner progress bar; worker updates it under the output lock */

  explicit sintax_state_s(struct Parameters const & params) : parameters(params), seed(params) {}
};


/* `candidates` is the seqno of the best hit of each successful bootstrap: the
   filled prefix of one row of the caller's all_seqno, whose fill level used to
   be handed over beside the pointer as a separate `count`. */
static auto sintax_analyse(struct sintax_state_s & state,
                    View<char> const query_head,
                    int const strand,
                    View<int> const candidates) -> void
{
  std::FILE * const fp_tabbedout = state.fp_tabbedout;

  constexpr auto levels = static_cast<std::size_t>(tax_levels);
  std::array<int, tax_levels> level_matchcount {{}};
  std::array<int, tax_levels> level_best {{}};
  /* the taxonomy name of each candidate at each rank, as a window into that
     candidate's header; a single array of views replaces the parallel
     start-pointer and length arrays this used to keep side by side */
  std::array<std::array<View<char>, tax_levels>, bootstrap_count> cand_level_name {{}};

  /* Check number of successful bootstraps, must be at least half */

  auto const is_enough =
    candidates.size() >= static_cast<std::size_t>((bootstrap_count + 1) / 2);

  if (is_enough)
    {
      /* which makes the range every algorithm below runs over non-empty, and
         no longer than the per-candidate arrays, since a bootstrap contributes
         at most one candidate */
      assert(not candidates.empty());
      assert(candidates.size() <= cand_level_name.size());

      /* Find the most common name at each taxonomic rank,
         but with the same names at higher ranks. */

      /* Split headers of all candidates by taxonomy ranks */

      std::transform(candidates.begin(), candidates.end(), cand_level_name.begin(),
                     [&state](int const seqno) -> std::array<View<char>, tax_levels> {
                       std::array<TaxLevel, tax_levels> ranks {{}};
                       tax_split(seqno, ranks, state.db);
                       auto const header = state.db.header_view(static_cast<uint64_t>(seqno));
                       std::array<View<char>, tax_levels> names {{}};
                       std::transform(ranks.begin(), ranks.end(), names.begin(),
                                      [header](TaxLevel const & rank) -> View<char> {
                                        return header.subspan(static_cast<std::size_t>(rank.start),
                                                              static_cast<std::size_t>(rank.length));
                                      });
                       return names;
                     });

      std::array<bool, bootstrap_count> cand_included {{}};
      cand_included.fill(true);

      /* Count matching names among candidates */

      /* indexed rather than a range-for: the rank is the position shared by
         level_best, level_matchcount, taxonomic_fields and the second subscript
         of cand_level_name, so here the index is data */
      for (std::size_t level = 0; level < levels; ++level)
        {
          level_best[level] = -1;
          level_matchcount[level] = 0;

          std::array<int, bootstrap_count> cand_match {{}};
          cand_match.fill(-1);
          std::array<int, bootstrap_count> cand_matchcount {{}};

          /* the first still-included candidate, at or before `cand_i`, carrying
             the same name at this rank -- the candidate that stands for the
             group. A candidate always matches itself, so for an included cand_i
             the search cannot fail; naming it is what removes the break from
             the middle of the nested loop this used to be. */
          auto const group_of = [&](std::size_t const cand_i) -> std::size_t {
            for (std::size_t cand_j = 0; cand_j <= cand_i; ++cand_j)
              {
                if (cand_included[cand_j] and
                    (cand_level_name[cand_i][level] == cand_level_name[cand_j][level]))
                  {
                    return cand_j;
                  }
              }
            return cand_i;  /* unreachable: cand_i matches itself */
          };

          for (std::size_t cand_i = 0; cand_i < candidates.size(); ++cand_i)
            {
              if (cand_included[cand_i])
                {
                  auto const group = group_of(cand_i);
                  cand_match[cand_i] = static_cast<int>(group);
                  ++cand_matchcount[group];
                }
            }

          /* the largest group at this rank. Strictly-greater kept the first of
             a tie, which is what max_element returns, so the two agree
             candidate for candidate -- but max_element over an all-zero range
             would report candidate 0, where the loop left level_best at -1,
             hence the explicit test. */
          auto * const largest_group =
            std::max_element(cand_matchcount.begin(),
                             std::next(cand_matchcount.begin(),
                                       static_cast<std::ptrdiff_t>(candidates.size())));
          if (*largest_group > 0)
            {
              level_best[level] =
                static_cast<int>(std::distance(cand_matchcount.begin(), largest_group));
              level_matchcount[level] = *largest_group;
            }

          /* indexed for the same reason: cand_match and cand_included are two
             arrays indexed by the same candidate number */
          for (std::size_t cand_i = 0; cand_i < candidates.size(); ++cand_i)
            {
              if (cand_match[cand_i] != level_best[level])
                {
                  cand_included[cand_i] = false;
                }
            }
        }
    }

  /* write to tabbedout file */
  std::lock_guard<std::mutex> const output_lock(state.mutex_output);
  fprint(fp_tabbedout, query_head);
  fprint(fp_tabbedout, '\t');

  state.queries++;

  if (is_enough)
    {
      state.classified++;

      auto comma = false;
      for (std::size_t level = 0; level < levels; ++level)
        {
          /* every rank has a largest group once is_enough holds: candidate 0
             stays included at every rank, matching itself */
          assert(level_best[level] >= 0);
          auto const best = static_cast<std::size_t>(level_best[level]);
          auto const & level_name = cand_level_name[best][level];
          if (not level_name.empty())
            {
              std::fputs((comma ? "," : ""), fp_tabbedout);
              fprint(fp_tabbedout, static_cast<char>(taxonomic_fields[level]));
              fprint(fp_tabbedout, ':');
              fprint(fp_tabbedout, level_name);
              fprint(fp_tabbedout, '(');
              std::fprintf(fp_tabbedout, "%.2f",
                           1.0 * level_matchcount[level] / static_cast<double>(candidates.size()));
              fprint(fp_tabbedout, ')');
              comma = true;
            }
        }

      fprint(fp_tabbedout, '\t');
      fprint(fp_tabbedout, (strand != 0) ? '-' : '+');

      if (state.parameters.opt_sintax_cutoff > 0.0)
        {
          fprint(fp_tabbedout, '\t');
          auto comma_cutoff = false;
          for (std::size_t level = 0; level < levels; ++level)
            {
              assert(level_best[level] >= 0);
              auto const best = static_cast<std::size_t>(level_best[level]);
              auto const & level_name = cand_level_name[best][level];
              if ((not level_name.empty()) &&
                  (1.0 * level_matchcount[level] / static_cast<double>(candidates.size())
                   >= state.parameters.opt_sintax_cutoff))
                {
                  std::fputs((comma_cutoff ? "," : ""), fp_tabbedout);
                  fprint(fp_tabbedout, static_cast<char>(taxonomic_fields[level]));
                  fprint(fp_tabbedout, ':');
                  fprint(fp_tabbedout, level_name);
                  comma_cutoff = true;
                }
            }
        }
    }
  else
    {
      if (state.parameters.opt_sintax_cutoff > 0.0)
        {
          fprint(fp_tabbedout, "\t\t");
        }
      else
        {
          fprint(fp_tabbedout, '\t');
        }
    }

  fprint(fp_tabbedout, '\n');
}


namespace {
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
  auto const kmer_counts = make_span(searchinfo->kmers_v);
  assert(indexed_count <= kmer_counts.size());

  /* zero counts */
  std::fill_n(kmer_counts.begin(), indexed_count, count_t{0});

  for (auto const kmer : searchinfo->kmersample)
    {
      auto const * bitmap = searchinfo->dbindex->getbitmap(kmer);

      if (bitmap != nullptr)
        {
#ifdef __x86_64__
          if (parameters.ssse3_present != 0)
            {
              increment_counters_from_bitmap_ssse3(kmer_counts.data(),
                                                   bitmap, indexed_count);
            }
          else
            {
              increment_counters_from_bitmap_sse2(kmer_counts.data(),
                                                  bitmap, indexed_count);
            }
#else
          increment_counters_from_bitmap(kmer_counts.data(), bitmap, indexed_count);
#endif
        }
      else
        {
          auto const * list = searchinfo->dbindex->getmatchlist(kmer);
          auto const count = searchinfo->dbindex->getmatchcount(kmer);
          for (auto j = 0U; j < count; j++)
            {
              kmer_counts[list[j]]++;
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
      count_t const count = kmer_counts[i];
      auto const seqno = searchinfo->dbindex->getmapping(i);
      auto const length = static_cast<unsigned int>(searchinfo->db->getsequencelen(seqno));

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
/* The state is passed whole rather than its two int counters, which would be
   two adjacent swappable arguments. */
auto print_classified_count(std::FILE * output_stream,
                            struct sintax_state_s const & state) -> void
{
  fprint(output_stream, "Classified ");
  fprint_integer(output_stream, state.classified);
  fprint(output_stream, " of ");
  fprint_integer(output_stream, state.queries);
  fprint(output_stream, " sequences");
  if (state.queries > 0)
    {
      fprint(output_stream, " (");
      std::fprintf(output_stream, "%.2f", 100.0 * state.classified / state.queries);
      fprint(output_stream, "%)");
    }
  fprint(output_stream, '\n');
}

}  // anonymous namespace


static auto sintax_query(struct sintax_state_s & state, uint64_t const t) -> void
{
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;

  std::array<std::array<int, bootstrap_count>, 2> all_seqno {{}};
  std::array<int, 2> boot_count = {0, 0};
  std::array<unsigned int, 2> best_count = {0, 0};
  int const qseqlen = static_cast<int>(si_plus[t].qsequence.size());
  auto const query_head = si_plus[t].query_head;

  /* Per-query RNG: seed from the run's base seed and this query's input
     number, so the random subsampling and tie-breaking are reproducible
     regardless of how many threads process the queries or in what order. */
  SplitMix64 rng(state.seed.substream(static_cast<uint64_t>(si_plus[t].query_no)));

  Bitmap b(static_cast<unsigned int>(qseqlen));

  /* indexed rather than a range-for over the two strands: the counter is also
     the row this strand's bootstrap hits are accumulated into (all_seqno,
     boot_count and best_count below), so here the index is data */
  for (auto s = 0; s < number_of_strands(state.parameters.opt_strand); s++)
    {
      struct searchinfo_s * const si = (s != 0) ? &si_minus[t] : &si_plus[t];

      /* perform search */

      /* find unique kmers at dbindex.wordlength, the effective index width (set
         by Dbindex::prepare for a FASTA db, or udb_read for a UDB db); reading
         parameters.opt_wordlength would use the wrong width against a UDB index. */
      auto const kmersample = si->uh.count(static_cast<int>(si->dbindex->wordlength),
                                           View<char>{si->qsequence}, Masking::none);

      /* perform 100 bootstraps */

      if (kmersample.size() >= subset_size)
        {
          /* declared outside the loop, so that the view si->kmersample holds
             into it cannot outlive it; refilled from index 0 every round, and
             only the first `subsamples` entries are ever read */
          std::array<unsigned int, subset_size> kmersample_subset {{}};

          for (auto i = 0; i < bootstrap_count ; i++)
            {
              /* subsample 32 kmers */
              auto subsamples = 0;
              b.reset_all();
              for (auto j = 0; j < subset_size ; j++)
                {
                  std::size_t const x = random_bounded(rng, static_cast<unsigned int>(kmersample.size()));
                  if (not b.is_set(static_cast<unsigned int>(x)))
                    {
                      kmersample_subset[static_cast<std::size_t>(subsamples++)] = kmersample[x];
                      b.set(static_cast<unsigned int>(x));
                    }
                }

              si->kmersample = make_view(kmersample_subset).first(static_cast<std::size_t>(subsamples));

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

  auto const winner = static_cast<std::size_t>(best_strand);
  sintax_analyse(state,
                 query_head,
                 best_strand,
                 make_view(all_seqno[winner])
                   .first(static_cast<std::size_t>(boot_count[winner])));
}


static auto sintax_thread_run(struct sintax_state_s & state, uint64_t const t) -> void
{
  std::mutex & mutex_input = state.mutex_input;
  std::mutex & mutex_output = state.mutex_output;
  struct fastx_s * const query_fastx_h = state.query_fastx_h;
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;

  uint64_t progress = 0;

  auto const has_work_to_claim = [&]() -> bool {
    if (not query_fastx_h->next(
                       not state.parameters.opt_notrunclabels,
                       chrmap_no_change()))
      {
        /* End of input, or a deferred parse error was recorded (CC3):
           fastx_next() returns false in both cases, so the worker stops
           here cooperatively. The error, if any, is reported by sintax()
           from the main thread after the pool joins. */
        return false;
      }

    auto const qhead = query_fastx_h->header_view();
    auto const qseq = query_fastx_h->sequence_view();
    auto const qseqlen = static_cast<int>(qseq.size());
    int const query_no = static_cast<int>(query_fastx_h->get_seqno());
    int64_t const qsize = query_fastx_h->get_abundance();

    /* indexed rather than a range-for over the two strands: the counter is
       stored as si->strand, so here the index is data */
    for (auto s = 0; s < number_of_strands(state.parameters.opt_strand); s++)
      {
        struct searchinfo_s * const si = (s != 0) ? &si_minus[t] : &si_plus[t];

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
    si_plus[t].query_head_v.resize(qhead.size() + 1);
    std::copy(qhead.cbegin(), qhead.cend(), si_plus[t].query_head_v.begin());
    si_plus[t].query_head_v[qhead.size()] = '\0';
    si_plus[t].query_head = make_view(si_plus[t].query_head_v).first(qhead.size());
    std::copy(qseq.cbegin(), qseq.cend(), si_plus[t].qsequence_v.begin());
    si_plus[t].qsequence_v[qseq.size()] = '\0';
    si_plus[t].qsequence = make_span(si_plus[t].qsequence_v).first(qseq.size());

    /* get progress as amount of input file read */
    progress = query_fastx_h->get_position();
    return true;
  };

  auto const process_query = [&]() -> void {
    /* minus strand: copy header and reverse complementary sequence */
    if (state.parameters.opt_strand)
      {
        si_minus[t].query_head_v = si_plus[t].query_head_v;
        si_minus[t].query_head = make_view(si_minus[t].query_head_v).first(si_plus[t].query_head.size());
        reverse_complement(make_span(si_minus[t].qsequence_v).first(si_plus[t].qsequence.size() + 1),
                           View<char>{si_plus[t].qsequence});
        si_minus[t].qsequence = make_span(si_minus[t].qsequence_v).first(si_plus[t].qsequence.size());
      }

    sintax_query(state, t);

    /* lock mutex for update of global data and output */
    std::lock_guard<std::mutex> const output_lock(mutex_output);

    /* show progress */
    state.progress->update(progress);
  };

  run_worker_loop(mutex_input, has_work_to_claim, process_query);
}


static auto sintax_thread_init(struct sintax_state_s const & state, struct searchinfo_s & si) -> void
{
  /* thread specific initialiation */
  si.parameters = &state.parameters;  /* searchcore reads config through the si (E1) */
  si.dbindex = &state.dbindex;  /* searchcore reads the k-mer index through the si */
  si.db = &state.db;  /* searchcore reads the sequences through the si */
  /* si->uh (a Uniquer value member) is ready to use as default-constructed */
  /* the kmer counts live in the searchinfo_s kmers_v vector (RAII), matching
     search/cluster; the reserve headroom keeps the SIMD counter stores that
     may run past the logical end in bounds. */
  static constexpr auto overflow_padding = 16U;  // 16 * sizeof(count_t) = 32 bytes headroom
  si.kmers_v.reserve(static_cast<size_t>(state.seqcount) + overflow_padding);
  si.kmers_v.resize(static_cast<size_t>(state.seqcount));
  si.m = Minheap(state.tophits);
  si.qsize = 1;
  si.query_head = View<char>{nullptr, 0};
  si.seq_alloc = 0;
  si.qsequence = Span<char>{};
  si.nw = nullptr;
  si.s.reset();
}


static auto sintax_thread_exit(struct searchinfo_s & searchinfo) -> void
{
  /* thread specific clean up */
  searchinfo.uh = Uniquer();
  searchinfo.m = Minheap();
  /* the kmer counts, the query header and the query sequence live in the
     searchinfo_s vectors (kmers_v/query_head_v/qsequence_v), which free
     their own storage */
}


static auto sintax_thread_worker_run(struct sintax_state_s & state) -> void
{
  auto & si_plus = state.si_plus;
  auto & si_minus = state.si_minus;

  /* init per-thread search state before the workers start */
  for (auto t = 0; t < state.parameters.opt_threads; t++)
    {
      sintax_thread_init(state, si_plus[static_cast<std::size_t>(t)]);
      if (not si_minus.empty())
        {
          sintax_thread_init(state, si_minus[static_cast<std::size_t>(t)]);
        }
    }

  /* run the worker pool over the input file */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(state.parameters.opt_threads),
                              [&state](uint64_t const t) -> void
                              { sintax_thread_run(state, t); });
    threadrunner.run();
  }

  /* clean up per-thread search state */
  for (auto t = 0; t < state.parameters.opt_threads; t++)
    {
      sintax_thread_exit(si_plus[static_cast<std::size_t>(t)]);
      if (not si_minus.empty())
        {
          sintax_thread_exit(si_minus[static_cast<std::size_t>(t)]);
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
      udb_read(parameters.opt_db, UdbUse::search, state.dbindex, state.db, parameters);
    }
  else
    {
      state.db.read(parameters.opt_db, 0, parameters);
    }

  seqcount = static_cast<int>(state.db.getsequencecount());

  if (! is_udb)
    {
      state.dbindex.prepare(parameters.opt_dbmask, state.db, parameters);
      state.dbindex.add_all_sequences(parameters.opt_dbmask, state.db, parameters);
    }

  /* prepare reading of queries */

  auto const query_fastx_h = fastx_open(parameters.opt_sintax, parameters);
  state.query_fastx_h = query_fastx_h.get();  // workers borrow the raw handle

  /* The query file is parsed inside the worker threads (see
     sintax_thread_run). Enable deferred error reporting so a malformed
     query records the error and stops the pool cooperatively, rather than
     calling fatal()/std::exit() from a worker thread while siblings are
     still writing output (CC3). The error is reported below, from the
     main thread, after the pool has joined. */
  query_fastx_h->enable_deferred_errors();

  /* allocate memory for thread info */

  si_plus.resize(static_cast<std::size_t>(parameters.opt_threads));
  if (parameters.opt_strand)
    {
      si_minus.resize(static_cast<std::size_t>(parameters.opt_threads));
    }

  /* run */

  {
    Progress progress("Classifying sequences", query_fastx_h->get_size(), parameters);
    state.progress = &progress;
    sintax_thread_worker_run(state);
  }

  /* All workers have joined. If one hit a malformed query, report it now
     from the main thread (single-threaded) so the message is emitted and
     the process exits without racing any worker (CC3). */
  if (query_fastx_h->get_error())
    {
      fatal(query_fastx_h->get_errmsg());
    }

  if (! parameters.opt_quiet)
    {
      print_classified_count(stderr, state);
    }

  if (parameters.fp_log != nullptr)
    {
      print_classified_count(parameters.fp_log, state);
    }

  /* clean up */

  /* si_plus/si_minus are std::vector members of state (RAII) */

  query_fastx_h->report_stripped_warning(parameters);

  state.dbindex.clear();
  state.db.clear();
}
