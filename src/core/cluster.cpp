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
#include "core/align_simd.hpp"
#include "core/cluster.hpp"
#include "core/cluster_internal.hpp"
#include "searchcore.h"
#include "core/attributes.hpp"
#include "dbindex.h"
#include "core/linmemalign.hpp"
#include "core/mask.hpp"
#include "core/minheap.hpp"
#include "core/msa.hpp"
#include "otutable.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/make_unique.hpp"
#include "utils/open_file.hpp"
#include "utils/number_of_strands.hpp"
#include "utils/threads.hpp"
#include <algorithm>  // std::count, std::minmax_element, std::max_element, std::min
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstring>  // std::strcpy, std::strlen
#include <limits>
#include <map>
#include <memory>  // std::unique_ptr
#include <utility>  // std::get
#include <vector>



struct clusterinfo_s
{
  int seqno;
  int clusterno;
  char * cigar;
  int strand;
};

using clusterinfo_t = struct clusterinfo_s;

/* Per-invocation CLI state for the cluster*() commands — the output handles,
   the per-sequence cluster-assignment array, the cluster count and the
   matched/not-matched counters. Threaded through cluster() and its core
   drivers/output helpers so the command holds no shared mutable output state
   (E4). The search/worker state lives separately in cluster_work_pool_s, and
   the library session/batch paths return results rather than writing files, so
   they do not use this struct. */
struct cluster_cli_state_s
{
  /* the run configuration, threaded through the CLI-path helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter. The
     helpers shared with the library path (cluster_query_init, evaluate_extra_hits,
     free_hit_alignments, the searchcore) keep reading the globals — the library
     session/batch entries have no Parameters. */
  struct Parameters const & parameters;
  /* a copy of parameters with opt_maxaccepts/opt_maxrejects clamped to the
     database size (cluster()); si->parameters points here so the shared
     searchcore and evaluate_extra_hits read the clamped values without a
     mutated global (E1). */
  struct Parameters effective_parameters;
  struct Dbindex dbindex;  /* the k-mer index this run owns (RAII); si->dbindex points here */
  clusterinfo_t * clusterinfo = nullptr;
  int clusters = 0;
  int count_matched = 0;
  int count_notmatched = 0;
  std::FILE * fp_centroids = nullptr;
  std::FILE * fp_uc = nullptr;
  std::FILE * fp_alnout = nullptr;
  std::FILE * fp_samout = nullptr;
  std::FILE * fp_userout = nullptr;
  std::FILE * fp_blast6out = nullptr;
  std::FILE * fp_fastapairs = nullptr;
  std::FILE * fp_matched = nullptr;
  std::FILE * fp_notmatched = nullptr;
  std::FILE * fp_otutabout = nullptr;
  std::FILE * fp_mothur_shared_out = nullptr;
  std::FILE * fp_biomout = nullptr;
  std::FILE * fp_qsegout = nullptr;
  std::FILE * fp_tsegout = nullptr;

  explicit cluster_cli_state_s(struct Parameters const & params) : parameters(params) {}
};

/* per-thread slice of queries assigned for the current round; owned by
   cluster_work_pool_s (the threading primitives live inside its ThreadRunner) */
struct thread_work_s
{
  int query_first;
  int query_count;
};


inline auto cluster_query_core(struct searchinfo_s * si, struct Parameters const & parameters) -> void
{
  /* the main core function for clustering */

  /* get sequence etc */
  const int seqno = si->query_no;
  auto const useqno = static_cast<uint64_t>(seqno);
  si->query_head_len = static_cast<int>(db_getheaderlen(useqno));
  si->query_head = db_getheader(useqno);
  si->qsize = static_cast<int64_t>(db_getabundance(useqno));
  si->qseqlen = static_cast<int>(db_getsequencelen(useqno));
  if (si->strand != 0)
    {
      reverse_complement(si->qsequence, db_getsequence(useqno), si->qseqlen);
    }
  else
    {
      std::strcpy(si->qsequence, db_getsequence(useqno));
    }

  /* perform search */
  search_onequery(si, parameters.opt_qmask);
}


auto cluster_query_init(struct searchinfo_s * si, int const seqcount, int const tophits,
                        struct Parameters const & parameters,
                        struct Dbindex const & dbindex) -> void
{
  /* initialisation of data for one thread; run once for each thread */
  /* thread specific initialiation */

  si->parameters = &parameters;  /* searchcore reads config through the si (E1) */
  si->dbindex = &dbindex;  /* searchcore reads the k-mer index through the si */
  si->qsize = 1;
  si->nw = nullptr;
  si->hit_count = 0;

  /* allocate memory for sequence */

  si->seq_alloc = static_cast<int>(db_getlongestsequence() + 1);
  si->qsequence = static_cast<char *>(xmalloc(static_cast<std::size_t>(si->seq_alloc)));

  si->kmers = static_cast<count_t *>(xmalloc((static_cast<std::size_t>(seqcount) * sizeof(count_t)) + 32));
  si->hits = static_cast<struct hit *>(xmalloc(sizeof(struct hit) * static_cast<std::size_t>(tophits)));

  si->uh = unique_init();
  si->m = minheap_init(tophits);
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


auto cluster_query_exit(struct searchinfo_s * si) -> void
{
  /* clean up after thread execution; called once per thread */

  search16_exit(si->s);
  unique_exit(si->uh);
  minheap_exit(si->m);

  if (si->qsequence != nullptr)
    {
      xfree(si->qsequence);
    }
  if (si->hits != nullptr)
    {
      xfree(si->hits);
    }
  if (si->kmers != nullptr)
    {
      xfree(si->kmers);
    }
}


/* Self-contained per-invocation worker pool for the clustering search phase.
   It owns the per-thread searchinfo_s arrays, the per-round query slices, and
   its own ThreadRunner, so a caller drives its own pool with no shared
   file-static state (E4) — this is what lets cluster_assign_batch() stop
   borrowing the CLI path's si_plus/si_minus/thread_work/cluster_threadrunner
   via a save/restore hack.

   Similar to the Scanner class in swarm (src/utils/scanner.{h,cc}), the sister
   project's equivalent abstraction: the per-thread search state is a member
   vector, the worker reads its own slice, and the ThreadRunner is held as the
   last-created member whose lambda captures `this` — so the object must keep a
   stable address and is non-copyable/non-movable. */
struct cluster_work_pool_s
{
  struct Parameters const & parameters;  // run config, read by the workers (E1)
  struct Dbindex const & dbindex;       // k-mer index the workers search
  std::vector<searchinfo_s> si_plus;    // one entry per thread
  std::vector<searchinfo_s> si_minus;   // empty unless searching both strands
  std::vector<thread_work_s> thread_work;
  std::unique_ptr<ThreadRunner> runner;  // constructed last; lambda captures this

  cluster_work_pool_s(int const nthreads, int const seqcount,
                      int const tophits, bool const need_minus,
                      struct Parameters const & params,
                      struct Dbindex const & index)
    : parameters(params),
      dbindex(index),
      si_plus(static_cast<std::size_t>(nthreads)),
      si_minus(need_minus ? static_cast<std::size_t>(nthreads) : std::size_t{0}),
      thread_work(static_cast<std::size_t>(nthreads))
  {
    for (auto & si : si_plus)
      {
        cluster_query_init(&si, seqcount, tophits, parameters, dbindex);
        si.strand = 0;
      }
    for (auto & si : si_minus)
      {
        cluster_query_init(&si, seqcount, tophits, parameters, dbindex);
        si.strand = 1;
      }
    runner = make_unique<ThreadRunner>(static_cast<std::size_t>(nthreads),
                                       [this](uint64_t const t) { worker(t); });
  }

  ~cluster_work_pool_s()
  {
    runner.reset();  // join the workers before freeing the buffers they read
    for (auto & si : si_plus) { cluster_query_exit(&si); }
    for (auto & si : si_minus) { cluster_query_exit(&si); }
  }

  cluster_work_pool_s(cluster_work_pool_s const &) = delete;
  cluster_work_pool_s(cluster_work_pool_s &&) = delete;
  auto operator=(cluster_work_pool_s const &) -> cluster_work_pool_s & = delete;
  auto operator=(cluster_work_pool_s &&) -> cluster_work_pool_s & = delete;

  /* worker body: process this thread's assigned slice of the current round */
  auto worker(uint64_t const t) -> void
  {
    auto const & work = thread_work[t];
    for (int q = 0; q < work.query_count; q++)
      {
        cluster_query_core(si_plus.data() + work.query_first + q, parameters);
        if (not si_minus.empty())
          {
            cluster_query_core(si_minus.data() + work.query_first + q, parameters);
          }
      }
  }

  /* distribute `queries` across the threads and run one round */
  auto wakeup(int const queries) -> void
  {
    int const nthreads = static_cast<int>(thread_work.size());
    int const active = queries > nthreads ? nthreads : queries;
    int queries_rest = queries;
    int threads_rest = active;
    int query_next = 0;
    for (int t = 0; t < nthreads; t++)
      {
        auto const tdx = static_cast<std::size_t>(t);
        if (t < active)
          {
            thread_work[tdx].query_first = query_next;
            thread_work[tdx].query_count = (queries_rest + threads_rest - 1) / threads_rest;
            queries_rest -= thread_work[tdx].query_count;
            query_next += thread_work[tdx].query_count;
            --threads_rest;
          }
        else
          {
            thread_work[tdx].query_first = query_next;
            thread_work[tdx].query_count = 0;
          }
      }
    runner->run();
  }
};


auto relabel_otu(int clusterno, char const * sequence, int seqlen, struct Parameters const & parameters) -> char *
{
  char * label = nullptr;
  if (parameters.opt_relabel != nullptr)
    {
      int const size = static_cast<int>(std::strlen(parameters.opt_relabel)) + 21;
      label = static_cast<char *>(xmalloc(static_cast<std::size_t>(size)));
      std::snprintf(label, static_cast<std::size_t>(size), "%s%d", parameters.opt_relabel, clusterno + 1);
    }
  else if (parameters.opt_relabel_self)
    {
      int const size = seqlen + 1;
      label = static_cast<char *>(xmalloc(static_cast<std::size_t>(size)));
      std::snprintf(label, static_cast<std::size_t>(size), "%.*s", seqlen, sequence);
    }
  else if (parameters.opt_relabel_sha1)
    {
      label = static_cast<char *>(xmalloc(len_hex_dig_sha1));
      get_hex_seq_digest_sha1(label, sequence, seqlen);
    }
  else if (parameters.opt_relabel_md5)
    {
      label = static_cast<char *>(xmalloc(len_hex_dig_md5));
      get_hex_seq_digest_md5(label, sequence, seqlen);
    }
  return label;
}


auto cluster_core_results_hit(struct cluster_cli_state_s & state,
                              struct hit const * best,
                              int clusterno,
                              char const * query_head,
                              int qseqlen,
                              char const * qsequence,
                              char const * qsequence_rc,
                              int64_t qsize) -> void
{
  ++state.count_matched;

  if ((state.parameters.opt_otutabout != nullptr) or (state.parameters.opt_mothur_shared_out != nullptr) or (state.parameters.opt_biomout != nullptr))
    {
      if ((state.parameters.opt_relabel != nullptr) or state.parameters.opt_relabel_self or state.parameters.opt_relabel_sha1 or state.parameters.opt_relabel_md5)
        {
          char * label = relabel_otu(clusterno,
                                     db_getsequence(static_cast<uint64_t>(best->target)),
                                     static_cast<int>(db_getsequencelen(static_cast<uint64_t>(best->target))),
                                     state.parameters);
          otutable_add(query_head, label, qsize);
          xfree(label);
        }
      else
        {
          otutable_add(query_head,
                       db_getheader(static_cast<uint64_t>(best->target)),
                       qsize);
        }
    }

  if (state.fp_uc != nullptr)
    {
      results_show_uc_one(state.fp_uc,
                          best, query_head,
                          qseqlen,
                          clusterno,
                          state.parameters);
    }

  if (state.fp_alnout != nullptr)
    {
      results_show_alnout(state.fp_alnout,
                          best, 1, query_head,
                          qsequence, qseqlen,
                          state.parameters);
    }

  if (state.fp_samout != nullptr)
    {
      results_show_samout(state.fp_samout,
                          best, 1, query_head,
                          qsequence, qsequence_rc,
                          state.parameters);
    }

  if (state.fp_fastapairs != nullptr)
    {
      results_show_fastapairs_one(state.fp_fastapairs,
                                  best,
                                  query_head,
                                  qsequence,
                                  qsequence_rc,
                                  state.parameters);
    }

  if (state.fp_qsegout != nullptr)
    {
      results_show_qsegout_one(state.fp_qsegout,
                               best,
                               query_head,
                               qsequence,
                               qseqlen,
                               qsequence_rc,
                               state.parameters);
    }

  if (state.fp_tsegout != nullptr)
    {
      results_show_tsegout_one(state.fp_tsegout,
                               best,
                               state.parameters);
    }

  if (state.fp_userout != nullptr)
    {
      results_show_userout_one(state.fp_userout, best, query_head,
                               qsequence, qseqlen, qsequence_rc,
                               state.parameters);
    }

  if (state.fp_blast6out != nullptr)
    {
      results_show_blast6out_one(state.fp_blast6out, best, query_head,
                                 qseqlen);
    }

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


auto cluster_core_results_nohit(struct cluster_cli_state_s & state,
                                int clusterno,
                                char const * query_head,
                                int qseqlen,
                                char const * qsequence,
                                char const * qsequence_rc,
                                int64_t qsize) -> void
{
  ++state.count_notmatched;

  if ((state.parameters.opt_otutabout != nullptr) or (state.parameters.opt_mothur_shared_out != nullptr) or (state.parameters.opt_biomout != nullptr))
    {
      if ((state.parameters.opt_relabel != nullptr) or state.parameters.opt_relabel_self or state.parameters.opt_relabel_sha1 or state.parameters.opt_relabel_md5)
        {
          char * label = relabel_otu(clusterno, qsequence, qseqlen, state.parameters);
          otutable_add(query_head, label, qsize);
          xfree(label);
        }
      else
        {
          otutable_add(query_head, query_head, qsize);
        }
    }

  if (state.parameters.opt_uc != nullptr)
    {
      std::fprintf(state.fp_uc, "S\t%d\t%d\t*\t*\t*\t*\t*\t", state.clusters, qseqlen);
      header_fprint_strip(state.fp_uc,
                          query_head,
                          static_cast<int>(std::strlen(query_head)),
                          state.parameters.opt_xsize,
                          state.parameters.opt_xee,
                          state.parameters.opt_xlength);
      std::fprintf(state.fp_uc, "\t*\n");
    }

  if (state.parameters.opt_output_no_hits != 0)
    {
      if (state.fp_userout != nullptr)
        {
          results_show_userout_one(state.fp_userout, nullptr, query_head,
                                   qsequence, qseqlen, qsequence_rc,
                                   state.parameters);
        }

      if (state.fp_blast6out != nullptr)
        {
          results_show_blast6out_one(state.fp_blast6out, nullptr, query_head,
                                     qseqlen);
        }
    }

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


auto compare_kmersample(const void * a, const void * b) -> int
{
  unsigned int const lhs = * static_cast<unsigned int const *>(a);
  unsigned int const rhs = * static_cast<unsigned int const *>(b);

  if (lhs < rhs)
    {
      return -1;
    }
  if (lhs > rhs)
    {
      return +1;
    }
  return 0;
}

static auto evaluate_extra_hits(struct searchinfo_s * si,
                                struct searchinfo_s const * si_plus,
                                const int * extra_list,
                                int extra_count,
                                LinearMemoryAligner & lma,
                                int const tophits) -> void
{
  int added = 0;

  /* Keep at most this many hits. The list is the tophits-sized si->hits buffer,
     but tophits is clamped to seqcount; on a small dataset with large
     --maxaccepts/--maxrejects the raw bound (maxaccepts + maxrejects - 1) can
     exceed tophits, so the insertion/shift below would write past the buffer.
     Clamp the bound to the actual capacity (S10). */
  int const hit_capacity =
    static_cast<int>(std::min<int64_t>(si->parameters->opt_maxaccepts + si->parameters->opt_maxrejects - 1,
                                       tophits));

  if (extra_count != 0)
    {
      /* Check if there is a hit with one of the non-matching
         extra sequences just analysed in this round */

      for (int j = 0; j < extra_count; j++)
        {
          struct searchinfo_s const * sic = si_plus + extra_list[j];

          /* find the number of shared unique kmers */
          auto const shared
            = unique_count_shared(*si->uh,
                                  static_cast<int>(si->dbindex->wordlength),
                                  static_cast<int>(sic->kmersamplecount),
                                  sic->kmersample);

          /* check if min number of shared kmers is satisfied */
          if (search_enough_kmers(*si, shared))
            {
              unsigned int const length = static_cast<unsigned int>(sic->qseqlen);

              /* Go through the list of hits and see if the current
                 match is better than any on the list in terms of
                 more shared kmers (or shorter length if equal
                 no of kmers). Determine insertion point (x). */

              int x = si->hit_count;
              while ((x > 0) and
                     ((si->hits[x - 1].count < shared) or
                      ((si->hits[x - 1].count == shared) and
                       (db_getsequencelen(static_cast<uint64_t>(si->hits[x - 1].target))
                        > length))))
                {
                  --x;
                }

              if (x < hit_capacity)
                {
                  /* insert into list at position x */

                  /* trash bottom element if no more space */
                  if (si->hit_count >= hit_capacity)
                    {
                      if (si->hits[si->hit_count-1].aligned)
                        {
                          xfree(si->hits[si->hit_count - 1].nwalignment);
                        }
                      --si->hit_count;
                    }

                  /* move the rest down */
                  for (int z = si->hit_count; z > x; z--)
                    {
                      si->hits[z] = si->hits[z - 1];
                    }

                  /* init new hit */
                  struct hit * hit = si->hits + x;
                  ++si->hit_count;

                  hit->target = sic->query_no;
                  hit->strand = si->strand;
                  hit->count = shared;
                  hit->accepted = false;
                  hit->rejected = false;
                  hit->aligned = false;
                  hit->weak = false;
                  hit->nwalignment = nullptr;

                  ++added;
                }
            }
        }
    }

  /* now go through the hits and determine final status of each */

  if (added != 0)
    {
      si->rejects = 0;
      si->accepts = 0;

      /* set all statuses to undetermined */

      for (int t = 0; t < si->hit_count; t++)
        {
          si->hits[t].accepted = false;
          si->hits[t].rejected = false;
        }

      for (int t = 0;
           (si->accepts < si->parameters->opt_maxaccepts) and
             (si->rejects < si->parameters->opt_maxrejects) and
             (t < si->hit_count);
           ++t)
        {
          struct hit * hit = si->hits + t;

          if (not hit->aligned)
            {
              /* Test accept/reject criteria before alignment */
              unsigned int const target = static_cast<unsigned int>(hit->target);
              if (search_acceptable_unaligned(*si, static_cast<int>(target)))
                {
                  /* perform vectorized alignment */
                  /* but only using 1 sequence ! */

                  unsigned int nwtarget = target;

                  int64_t nwscore = 0;
                  int64_t nwalignmentlength = 0;
                  int64_t nwmatches = 0;
                  int64_t nwmismatches = 0;
                  int64_t nwgaps = 0;
                  char * nwcigar = nullptr;

                  /* short variants for simd aligner */
                  CELL snwscore = 0;
                  unsigned short snwalignmentlength = 0;
                  unsigned short snwmatches = 0;
                  unsigned short snwmismatches = 0;
                  unsigned short snwgaps = 0;

                  search16(si->s,
                           1,
                           & nwtarget,
                           & snwscore,
                           & snwalignmentlength,
                           & snwmatches,
                           & snwmismatches,
                           & snwgaps,
                           & nwcigar);

                  int64_t const tseqlen = static_cast<int64_t>(db_getsequencelen(target));

                  if (snwscore == std::numeric_limits<short>::max())
                    {
                      /* In case the SIMD aligner cannot align,
                         perform a new alignment with the
                         linear memory aligner */

                      char * tseq = db_getsequence(target);

                      if (nwcigar != nullptr)
                        {
                          xfree(nwcigar);
                        }

                      nwcigar = xstrdup(lma.align(si->qsequence,
                                                  tseq,
                                                  si->qseqlen,
                                                  tseqlen));

                      lma.alignstats(nwcigar,
                                     si->qsequence,
                                     tseq,
                                     & nwscore,
                                     & nwalignmentlength,
                                     & nwmatches,
                                     & nwmismatches,
                                     & nwgaps);
                    }
                  else
                    {
                      nwscore = snwscore;
                      nwalignmentlength = snwalignmentlength;
                      nwmatches = snwmatches;
                      nwmismatches = snwmismatches;
                      nwgaps = snwgaps;
                    }


                  int64_t const nwdiff = nwalignmentlength - nwmatches;
                  int64_t const nwindels = nwdiff - nwmismatches;

                  hit->aligned = true;
                  hit->nwalignment = nwcigar;
                  hit->nwscore = static_cast<int>(nwscore);
                  hit->nwdiff = static_cast<int>(nwdiff);
                  hit->nwgaps = static_cast<int>(nwgaps);
                  hit->nwindels = static_cast<int>(nwindels);
                  hit->nwalignmentlength = static_cast<int>(nwalignmentlength);
                  hit->matches = static_cast<int>(nwmatches);
                  hit->mismatches = static_cast<int>(nwmismatches);

                  hit->nwid = 100.0 *
                    static_cast<double>(nwalignmentlength - hit->nwdiff) /
                    static_cast<double>(nwalignmentlength);

                  hit->shortest = std::min(si->qseqlen, static_cast<int>(tseqlen));
                  hit->longest = std::max(si->qseqlen, static_cast<int>(tseqlen));

                  /* trim alignment and compute numbers
                     excluding terminal gaps */
                  align_trim(hit, *si->parameters);
                }
              else
                {
                  /* rejection without alignment */
                  hit->rejected = true;
                  ++si->rejects;
                }
            }

          if (not hit->rejected)
            {
              /* test accept/reject criteria after alignment */
              if (search_acceptable_aligned(*si, hit))
                {
                  ++si->accepts;
                }
              else
                {
                  ++si->rejects;
                }
            }
        }

      /* delete all undetermined hits */

      int new_hit_count = si->hit_count;
      for (int t = si->hit_count - 1; t >= 0; t--)
        {
          struct hit const * hit = si->hits + t;
          if (not hit->accepted and not hit->rejected)
            {
              new_hit_count = t;
              if (hit->aligned)
                {
                  xfree(hit->nwalignment);
                }
            }
        }
      si->hit_count = new_hit_count;
    }
}

static auto free_hit_alignments(struct searchinfo_s * si_p,
                                struct searchinfo_s * si_m,
                                struct Parameters const & parameters) -> void
{
  /* free alignments */
  for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
    {
      struct searchinfo_s const * si = (s != 0) ? si_m : si_p;
      for (int j = 0; j < si->hit_count; j++)
        {
          if ((si->hits[j].aligned) && (si->hits[j].nwalignment != nullptr))
            {
              xfree(si->hits[j].nwalignment);
              si->hits[j].nwalignment = nullptr;
            }
        }
    }
}

auto cluster_core_parallel(struct cluster_cli_state_s & state,
                           int const seqcount, int const tophits) -> void
{
  constexpr static int queries_per_thread = 1;
  const int max_queries = queries_per_thread * static_cast<int>(state.parameters.opt_threads);

  /* Own worker pool + per-thread search state (E4); see cluster_work_pool_s.
     The local si_plus/si_minus aliases let the loops below read unchanged. */
  cluster_work_pool_s pool(static_cast<int>(state.parameters.opt_threads), seqcount, tophits,
                           state.parameters.opt_strand, state.effective_parameters, state.dbindex);
  searchinfo_s * const si_plus = pool.si_plus.data();
  searchinfo_s * const si_minus = pool.si_minus.empty() ? nullptr : pool.si_minus.data();

  std::vector<int> extra_list(static_cast<std::size_t>(max_queries));

  struct Scoring scoring = scoring_from_options(state.parameters);


  LinearMemoryAligner lma(scoring);


  auto lastlength = std::numeric_limits<int>::max();

  int seqno = 0;

  int64_t sum_nucleotides = 0;

  Progress progress("Clustering", db_getnucleotidecount(), state.parameters);

  while (seqno < seqcount)
    {
      /* prepare work for the threads in sia[i] */
      /* read query sequences into the search info (si) for each thread */

      int queries = 0;

      for (int i = 0; i < max_queries; i++)
        {
          if (seqno < seqcount)
            {
              int const length = static_cast<int>(db_getsequencelen(static_cast<uint64_t>(seqno)));

              if ((state.parameters.opt_cluster_smallmem != nullptr) and (state.parameters.opt_usersort == 0) and (length > lastlength))
                {
                  fatal("Sequences not sorted by length and --usersort not specified.");
                }

              lastlength = length;

              si_plus[i].query_no = seqno;
              si_plus[i].strand = 0;

              if (state.parameters.opt_strand)
                {
                  si_minus[i].query_no = seqno;
                  si_minus[i].strand = 1;
                }

              ++queries;
              ++seqno;
            }
        }

      /* perform work in threads */
      pool.wakeup(queries);

      /* analyse results */
      int extra_count = 0;

      for (int i = 0; i < queries; i++)
        {
          struct searchinfo_s * si_p = si_plus + i;
          struct searchinfo_s * si_m = state.parameters.opt_strand ? si_minus + i : nullptr;

          for (int s = 0; s < number_of_strands(state.parameters.opt_strand); s++)
            {
              struct searchinfo_s * si = (s != 0) ? si_m : si_p;

              evaluate_extra_hits(si, si_plus, extra_list.data(), extra_count, lma, tophits);
            }

          /* find best hit */
          struct hit * best = nullptr;
          if (state.parameters.opt_sizeorder)
            {
              best = search_findbest2_bysize(si_p, si_m);
            }
          else
            {
              best = search_findbest2_byid(si_p, si_m);
            }

          int const myseqno = si_p->query_no;

          if (best != nullptr)
            {
              /* a hit was found, cluster current sequence with hit */
              int const target = best->target;

              /* output intermediate results to uc etc */
              cluster_core_results_hit(state, best,
                                       state.clusterinfo[target].clusterno,
                                       si_p->query_head,
                                       si_p->qseqlen,
                                       si_p->qsequence,
                                       (best->strand != 0) ? si_m->qsequence : nullptr,
                                       si_p->qsize);

              /* update cluster info about this sequence */
              state.clusterinfo[myseqno].seqno = myseqno;
              state.clusterinfo[myseqno].clusterno = state.clusterinfo[target].clusterno;
              state.clusterinfo[myseqno].cigar = best->nwalignment;
              state.clusterinfo[myseqno].strand = best->strand;
              best->nwalignment = nullptr;
            }
          else
            {
              /* no hit found; add it to the list of extra sequences
                 that must be considered by the coming queries in this
                 round */
              extra_list[static_cast<std::size_t>(extra_count)] = i;
              ++extra_count;

              /* update cluster info about this sequence */
              state.clusterinfo[myseqno].seqno = myseqno;
              state.clusterinfo[myseqno].clusterno = state.clusters;
              state.clusterinfo[myseqno].cigar = nullptr;
              state.clusterinfo[myseqno].strand = 0;

              /* add current sequence to database */
              state.dbindex.add_sequence(static_cast<unsigned int>(myseqno), state.parameters.opt_qmask);

              /* output intermediate results to uc etc */
              cluster_core_results_nohit(state, state.clusters,
                                         si_p->query_head,
                                         si_p->qseqlen,
                                         si_p->qsequence,
                                         nullptr,
                                         si_p->qsize);
              ++state.clusters;
            }

          free_hit_alignments(si_p, si_m, state.parameters);

          sum_nucleotides += si_p->qseqlen;
        }

      progress.update(static_cast<uint64_t>(sum_nucleotides));
    }

  /* pool's destructor joins the workers and frees the per-thread search state */
}


auto cluster_core_serial(struct cluster_cli_state_s & state,
                         int const seqcount, int const tophits) -> void
{
  std::array<struct searchinfo_s, 1> si_p {{}};  // refactoring: direct initialization?
  std::array<struct searchinfo_s, 1> si_m {{}};

  cluster_query_init(si_p.data(), seqcount, tophits, state.effective_parameters, state.dbindex);
  if (state.parameters.opt_strand)
    {
      cluster_query_init(si_m.data(), seqcount, tophits, state.effective_parameters, state.dbindex);
    }

  auto lastlength = std::numeric_limits<int>::max();

  Progress progress("Clustering", static_cast<uint64_t>(seqcount), state.parameters);
  for (int seqno = 0; seqno < seqcount; seqno++)
    {
      int const length = static_cast<int>(db_getsequencelen(static_cast<uint64_t>(seqno)));

      if ((state.parameters.opt_cluster_smallmem != nullptr) and (state.parameters.opt_usersort == 0) and (length > lastlength))
        {
          fatal("Sequences not sorted by length and --usersort not specified.");
        }

      lastlength = length;

      si_p[0].query_no = seqno;
      si_p[0].strand = 0;
      cluster_query_core(si_p.data(), state.parameters);

      if (state.parameters.opt_strand)
        {
          si_m[0].query_no = seqno;
          si_m[0].strand = 1;
          cluster_query_core(si_m.data(), state.parameters);
        }

      struct hit * best = nullptr;
      if (state.parameters.opt_sizeorder)
        {
          best = search_findbest2_bysize(si_p.data(), si_m.data());
        }
      else
        {
          best = search_findbest2_byid(si_p.data(), si_m.data());
        }

      if (best != nullptr)
        {
          int const target = best->target;
          cluster_core_results_hit(state, best,
                                   state.clusterinfo[target].clusterno,
                                   si_p[0].query_head,
                                   si_p[0].qseqlen,
                                   si_p[0].qsequence,
                                   (best->strand != 0) ? si_m[0].qsequence : nullptr,
                                   si_p[0].qsize);
          state.clusterinfo[seqno].seqno = seqno;
          state.clusterinfo[seqno].clusterno = state.clusterinfo[target].clusterno;
          state.clusterinfo[seqno].cigar = best->nwalignment;
          state.clusterinfo[seqno].strand = best->strand;
          best->nwalignment = nullptr;
        }
      else
        {
          state.clusterinfo[seqno].seqno = seqno;
          state.clusterinfo[seqno].clusterno = state.clusters;
          state.clusterinfo[seqno].cigar = nullptr;
          state.clusterinfo[seqno].strand = 0;
          state.dbindex.add_sequence(static_cast<unsigned int>(seqno), state.parameters.opt_qmask);
          cluster_core_results_nohit(state, state.clusters,
                                     si_p[0].query_head,
                                     si_p[0].qseqlen,
                                     si_p[0].qsequence,
                                     nullptr,
                                     si_p[0].qsize);
          ++state.clusters;
        }

      free_hit_alignments(si_p.data(), si_m.data(), state.parameters);

      progress.update(static_cast<uint64_t>(seqno));
    }

  cluster_query_exit(si_p.data());
  if (state.parameters.opt_strand)
    {
      cluster_query_exit(si_m.data());
    }
}


auto cluster(char const * dbname,
             struct Parameters const & parameters) -> void
{
  cluster_cli_state_s state(parameters);
  auto & clusterinfo = state.clusterinfo;
  auto & clusters = state.clusters;
  auto & fp_centroids = state.fp_centroids;
  auto & fp_uc = state.fp_uc;
  auto & fp_alnout = state.fp_alnout;
  auto & fp_samout = state.fp_samout;
  auto & fp_userout = state.fp_userout;
  auto & fp_blast6out = state.fp_blast6out;
  auto & fp_fastapairs = state.fp_fastapairs;
  auto & fp_matched = state.fp_matched;
  auto & fp_notmatched = state.fp_notmatched;
  auto & fp_otutabout = state.fp_otutabout;
  auto & fp_mothur_shared_out = state.fp_mothur_shared_out;
  auto & fp_biomout = state.fp_biomout;
  auto & fp_qsegout = state.fp_qsegout;
  auto & fp_tsegout = state.fp_tsegout;

  OutputFileHandle centroids_handle = open_optional_output_file(parameters.opt_centroids, OutputOption{"--centroids"});
  fp_centroids = centroids_handle.get();
  OutputFileHandle uc_handle = open_optional_output_file(parameters.opt_uc, OutputOption{"--uc"});
  fp_uc = uc_handle.get();

  OutputFileHandle alnout_handle = open_optional_output_file(parameters.opt_alnout, OutputOption{"--alnout"});
  fp_alnout = alnout_handle.get();
  if (fp_alnout != nullptr)
    {
      std::fprintf(fp_alnout, "%s\n", parameters.command_line.c_str());
      std::fprintf(fp_alnout, "%s\n", parameters.prog_header.c_str());
    }

  OutputFileHandle samout_handle = open_optional_output_file(parameters.opt_samout, OutputOption{"--samout"});
  fp_samout = samout_handle.get();
  OutputFileHandle userout_handle = open_optional_output_file(parameters.opt_userout, OutputOption{"--userout"});
  fp_userout = userout_handle.get();
  OutputFileHandle blast6out_handle = open_optional_output_file(parameters.opt_blast6out, OutputOption{"--blast6out"});
  fp_blast6out = blast6out_handle.get();
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
  OutputFileHandle otutabout_handle = open_optional_output_file(parameters.opt_otutabout, OutputOption{"--otutabout"});
  fp_otutabout = otutabout_handle.get();
  OutputFileHandle mothur_shared_out_handle = open_optional_output_file(parameters.opt_mothur_shared_out, OutputOption{"--mothur_shared_out"});
  fp_mothur_shared_out = mothur_shared_out_handle.get();
  OutputFileHandle biomout_handle = open_optional_output_file(parameters.opt_biomout, OutputOption{"--biomout"});
  fp_biomout = biomout_handle.get();

  db_read(dbname, 0, parameters);

  otutable_init();

  results_show_samheader(fp_samout, dbname, parameters);

  if (parameters.opt_qmask == Masking::dust)
    {
      dust_all(parameters);
    }
  else if ((parameters.opt_qmask == Masking::soft) and (parameters.opt_hardmask))
    {
      hardmask_all();
    }

  // memory-intensive: the entire database is now held in memory

  int const seqcount = static_cast<int>(db_getsequencecount());

  if (parameters.opt_cluster_fast != nullptr)
    {
      db_sortbylength(parameters);
    }
  else if ((parameters.opt_cluster_size != nullptr) or (parameters.opt_cluster_unoise != nullptr))
    {
      db_sortbyabundance(parameters);
    }

  state.dbindex.prepare(1, parameters.opt_qmask, parameters);

  /* tophits = the maximum number of hits we need to store */

  /* Clamp maxrejects/maxaccepts to the database size (0 or "> seqcount" means
     "all"). Apply the clamp to a local Parameters copy threaded to cluster_core_*
     and evaluate_extra_hits via si->parameters, rather than mutating the shared
     config globals (E1 trap-writer). */
  state.effective_parameters = state.parameters;
  if ((state.effective_parameters.opt_maxrejects == 0) or
      (state.effective_parameters.opt_maxrejects > seqcount))
    {
      state.effective_parameters.opt_maxrejects = seqcount;
    }

  if ((state.effective_parameters.opt_maxaccepts == 0) or
      (state.effective_parameters.opt_maxaccepts > seqcount))
    {
      state.effective_parameters.opt_maxaccepts = seqcount;
    }

  int tophits = static_cast<int>(state.effective_parameters.opt_maxrejects +
                                 state.effective_parameters.opt_maxaccepts + MAXDELAYED);
  tophits = std::min(tophits, seqcount);

  std::vector<clusterinfo_t> clusterinfo_v(static_cast<std::size_t>(seqcount));
  clusterinfo = clusterinfo_v.data();

  if (parameters.opt_log != nullptr)
    {
      uint64_t const slots = 1ULL << (static_cast<uint64_t>(parameters.opt_wordlength) << 1ULL);
      std::fprintf(parameters.fp_log, "\n");
      std::fprintf(parameters.fp_log, "      Alphabet  nt\n");
      std::fprintf(parameters.fp_log, "    Word width  %" PRId64 "\n", parameters.opt_wordlength);
      std::fprintf(parameters.fp_log, "     Word ones  %" PRId64 "\n", parameters.opt_wordlength);
      std::fprintf(parameters.fp_log, "        Spaced  No\n");
      std::fprintf(parameters.fp_log, "        Hashed  No\n");
      std::fprintf(parameters.fp_log, "         Coded  No\n");
      std::fprintf(parameters.fp_log, "       Stepped  No\n");
      std::fprintf(parameters.fp_log, "         Slots  %" PRIu64 " (%.1fk)\n", slots, static_cast<double>(slots)/1000.0);
      std::fprintf(parameters.fp_log, "       DBAccel  100%%\n");
      std::fprintf(parameters.fp_log, "\n");
    }

  if (parameters.opt_threads == 1)
    {
      cluster_core_serial(state, seqcount, tophits);
    }
  else
    {
      cluster_core_parallel(state, seqcount, tophits);
    }


  /* find size and abundance of each cluster and save stats */

  std::vector<int64_t> cluster_abundance_v(static_cast<std::size_t>(clusters));
  std::vector<int> cluster_size_v(static_cast<std::size_t>(clusters));

  for (int i = 0; i < seqcount; i++)
    {
      int const seqno = clusterinfo_v[static_cast<std::size_t>(i)].seqno;
      int const clusterno = clusterinfo_v[static_cast<std::size_t>(i)].clusterno;
      cluster_abundance_v[static_cast<std::size_t>(clusterno)] +=
        parameters.opt_sizein ? static_cast<int64_t>(db_getabundance(static_cast<uint64_t>(seqno))) : 1;
      ++cluster_size_v[static_cast<std::size_t>(clusterno)];
    }

  // refactoring: isolate in a function (returns struct abundance_stats)
  auto const minmax_elements = std::minmax_element(cluster_abundance_v.cbegin(),
                                                   cluster_abundance_v.cend());
  auto const abundance_min = cluster_abundance_v.empty() ? 0 : *std::get<0>(minmax_elements);
  auto const abundance_max = cluster_abundance_v.empty() ? 0 : *std::get<1>(minmax_elements);
  int const singletons = static_cast<int>(std::count(cluster_abundance_v.cbegin(),
                                    cluster_abundance_v.cend(), int64_t{1}));
  auto const max_element = std::max_element(cluster_size_v.cbegin(),
                                            cluster_size_v.cend());
  auto const size_max = cluster_size_v.empty() ? 0 : *max_element;


  /* Sort sequences in clusters by their abundance or ordinal number */
  /* Sequences in same cluster must always come right after each other. */
  /* The centroid sequence must be the first in each cluster. */

  {
    Progress const progress("Sorting clusters", static_cast<uint64_t>(clusters), parameters);
    if (parameters.opt_clusterout_sort)
      {
        /* by cluster abundance (descending), then cluster number, then seqno.
           The seqno tiebreak is a strict total order, so this matches the
           previous std::qsort result exactly. */
        std::sort(clusterinfo_v.begin(), clusterinfo_v.end(),
                  [&cluster_abundance_v](clusterinfo_t const & lhs, clusterinfo_t const & rhs) -> bool {
                    auto const lhs_ab = cluster_abundance_v[static_cast<std::size_t>(lhs.clusterno)];
                    auto const rhs_ab = cluster_abundance_v[static_cast<std::size_t>(rhs.clusterno)];
                    if (lhs_ab != rhs_ab) { return lhs_ab > rhs_ab; }
                    if (lhs.clusterno != rhs.clusterno) { return lhs.clusterno < rhs.clusterno; }
                    return lhs.seqno < rhs.seqno;
                  });
      }
    else
      {
        /* by cluster number, then seqno */
        std::sort(clusterinfo_v.begin(), clusterinfo_v.end(),
                  [](clusterinfo_t const & lhs, clusterinfo_t const & rhs) -> bool {
                    if (lhs.clusterno != rhs.clusterno) { return lhs.clusterno < rhs.clusterno; }
                    return lhs.seqno < rhs.seqno;
                  });
      }
  }


  /* allocate memory for full file name of the clusters files */
  OutputFileHandle fp_clusters;
  static constexpr auto space_for_cluster_id = 25;  // up to 25 digits
  std::vector<char> fn_clusters;
  if (parameters.opt_clusters != nullptr) {
    fn_clusters.reserve(std::strlen(parameters.opt_clusters) + space_for_cluster_id);
  }

  int lastcluster = -1;
  uint64_t ordinal = 0;

  {
    Progress progress("Writing clusters", static_cast<uint64_t>(seqcount), parameters);
    for (int i = 0; i < seqcount; i++)
      {
        int const seqno = clusterinfo_v[static_cast<std::size_t>(i)].seqno;
        int const clusterno = clusterinfo_v[static_cast<std::size_t>(i)].clusterno;

        if (clusterno != lastcluster)
          {
            /* prepare for new cluster */
            /* performed with first sequence only in each cluster */
            /* the first sequence is always the centroid */

            if (parameters.opt_centroids != nullptr)
              {
                fasta_print_general(fp_centroids,
                                    nullptr,
                                    db_getsequence(static_cast<uint64_t>(seqno)),
                                    static_cast<int>(db_getsequencelen(static_cast<uint64_t>(seqno))),
                                    db_getheader(static_cast<uint64_t>(seqno)),
                                    static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno))),
                                    static_cast<uint64_t>(cluster_abundance_v[static_cast<std::size_t>(clusterno)]),
                                    clusterno + 1,
                                    -1.0,
                                    -1,
                                    parameters.opt_clusterout_id ? clusterno : -1,
                                    nullptr, 0.0,
                                    db_getabundance(static_cast<uint64_t>(seqno)),
                                    parameters);
              }

            if (parameters.opt_uc != nullptr)
              {
                std::fprintf(fp_uc, "C\t%d\t%" PRId64 "\t*\t*\t*\t*\t*\t",
                        clusterno,
                        cluster_abundance_v[static_cast<std::size_t>(clusterno)]);
                header_fprint_strip(fp_uc,
                                    db_getheader(static_cast<uint64_t>(seqno)),
                                    static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno))),
                                    parameters.opt_xsize,
                                    parameters.opt_xee,
                                    parameters.opt_xlength);
                std::fprintf(fp_uc, "\t*\n");
              }

            if (parameters.opt_clusters != nullptr)
              {
                /* close previous (except for first time) and open new file */
                if (lastcluster != -1)
                  {
                    fp_clusters.reset();
                  }

                ordinal = 0;
                std::snprintf(fn_clusters.data(),
                         fn_clusters.capacity(),
                         "%s%d",
                         parameters.opt_clusters,
                         clusterno);
                fp_clusters = open_output_file(fn_clusters.data());
                if (not fp_clusters)
                  {
                    fatal("Unable to open clusters file for writing (%s)", fn_clusters.data());
                  }
              }

            lastcluster = clusterno;
          }

        /* performed for all sequences */

        if (parameters.opt_clusters != nullptr)
          {
            ++ordinal;
            fasta_print_db_relabel(fp_clusters.get(), static_cast<uint64_t>(seqno), ordinal, parameters);
          }

        progress.update(static_cast<uint64_t>(i));
      }

    if (lastcluster != -1)
      {
        /* performed with the last sequence */
        if (parameters.opt_clusters != nullptr)
          {
            fp_clusters.reset();
          }
      }
  }


  if (clusters < 1)
    {
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr, "Clusters: 0\n");
          std::fprintf(stderr, "Singletons: 0\n");
        }
      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log, "Clusters: 0\n");
          std::fprintf(parameters.fp_log, "Singletons: 0\n");
        }
    }
  else
    {
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr,
                  "Clusters: %d Size min %" PRId64 ", max %" PRId64 ", avg %.1f\n",
                  clusters,
                  abundance_min,
                  abundance_max,
                  1.0 * seqcount / clusters);
          std::fprintf(stderr,
                  "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
                  singletons,
                  100.0 * singletons / seqcount,
                  100.0 * singletons / clusters);
        }

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "Clusters: %d Size min %" PRId64 ", max %" PRId64 ", avg %.1f\n",
                  clusters,
                  abundance_min,
                  abundance_max,
                  1.0 * seqcount / clusters);
          std::fprintf(parameters.fp_log,
                  "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
                  singletons,
                  100.0 * singletons / seqcount,
                  100.0 * singletons / clusters);
          std::fprintf(parameters.fp_log, "\n");
        }
    }

  if ((parameters.opt_msaout != nullptr) or (parameters.opt_consout != nullptr) or (parameters.opt_profile != nullptr))
    {
      int msa_target_count = 0;
      std::vector<struct msa_target_s> msa_target_list_v(static_cast<std::size_t>(size_max));

      auto msaout_handle = open_optional_output_file(parameters.opt_msaout, OutputOption{"--msaout"});
      auto consout_handle = open_optional_output_file(parameters.opt_consout, OutputOption{"--consout"});
      auto profile_handle = open_optional_output_file(parameters.opt_profile, OutputOption{"--profile"});
      std::FILE * const fp_msaout = msaout_handle.get();
      std::FILE * const fp_consout = consout_handle.get();
      std::FILE * const fp_profile = profile_handle.get();

      lastcluster = -1;

      {
        Progress progress("Multiple alignments", static_cast<uint64_t>(seqcount), parameters);
        for (int i = 0; i < seqcount; i++)
          {
            int const clusterno = clusterinfo_v[static_cast<std::size_t>(i)].clusterno;
            int const seqno = clusterinfo_v[static_cast<std::size_t>(i)].seqno;
            char * cigar = clusterinfo_v[static_cast<std::size_t>(i)].cigar;
            int const strand = clusterinfo_v[static_cast<std::size_t>(i)].strand;

            if (clusterno != lastcluster)
              {
                if (lastcluster != -1)
                  {
                    /* compute msa & consensus */
                    msa(fp_msaout, fp_consout, fp_profile,
                        lastcluster,
                        msa_target_count, msa_target_list_v,
                        cluster_abundance_v[static_cast<std::size_t>(lastcluster)],
                        parameters);
                  }

                /* start new cluster */
                msa_target_count = 0;
                lastcluster = clusterno;
              }

            /* add current sequence to the cluster */
            msa_target_list_v[static_cast<std::size_t>(msa_target_count)].seqno = seqno;
            msa_target_list_v[static_cast<std::size_t>(msa_target_count)].cigar = cigar;
            msa_target_list_v[static_cast<std::size_t>(msa_target_count)].strand = strand;
            ++msa_target_count;

            progress.update(static_cast<uint64_t>(i));
          }

        if (lastcluster != -1)
          {
            /* compute msa & consensus */
            msa(fp_msaout, fp_consout, fp_profile,
                lastcluster,
                msa_target_count, msa_target_list_v,
                cluster_abundance_v[static_cast<std::size_t>(lastcluster)],
                parameters);
          }
      }


      profile_handle.reset();
      msaout_handle.reset();
      consout_handle.reset();
    }

  // cluster_abundance not used below that point
  // cluster_size not used below that point

  /* free cigar strings for all aligned sequences */

  for (auto const & cluster_entry : clusterinfo_v) {
    if (cluster_entry.cigar != nullptr)
      {
        xfree(cluster_entry.cigar);
      }
  }

  // clusterinfo not used after this point

  if (fp_biomout != nullptr)
    {
      otutable_print_biomout(fp_biomout, parameters);
      biomout_handle.reset();
    }

  if (fp_otutabout != nullptr)
    {
      otutable_print_otutabout(fp_otutabout, parameters);
      otutabout_handle.reset();
    }

  if (fp_mothur_shared_out != nullptr)
    {
      otutable_print_mothur_shared_out(fp_mothur_shared_out, parameters);
      mothur_shared_out_handle.reset();
    }

  otutable_done();

  /* reset() is a no-op on an empty handle, so unopened outputs need
     no guard; only userout carries extra teardown. */
  matched_handle.reset();
  notmatched_handle.reset();
  fastapairs_handle.reset();
  qsegout_handle.reset();
  tsegout_handle.reset();
  blast6out_handle.reset();
  userout_handle.reset();
  alnout_handle.reset();
  samout_handle.reset();
  uc_handle.reset();
  centroids_handle.reset();

  state.dbindex.clear();
  db_free();
}


/* === Library API for embedding clustering === */


struct cluster_session_s {
  std::unique_ptr<struct searchinfo_s> si;
  std::unique_ptr<struct searchinfo_s> si_minus;  /* non-null when searching both strands */
  int cluster_count = 0;
  int seqcount = 0;
  int tophits = 0;
  std::map<int, int> centroid_cluster_ids;  /* seqno → cluster_id for centroids */
  /* run configuration, set at cluster_session_init and read by the per-query
     path instead of the opt_* globals (E1 shared-infra phase); a pointer so the
     session stays default-constructible. */
  struct Parameters const * parameters = nullptr;
  /* the k-mer index this session clusters into, supplied by the caller at
     cluster_session_init; centroids are added to it incrementally (mutable), so
     it must outlive the session. */
  struct Dbindex * dbindex = nullptr;
};


auto cluster_session_alloc() -> struct cluster_session_s *
{
  return new cluster_session_s {};
}


auto cluster_session_free(struct cluster_session_s * cs) -> void
{
  if (cs != nullptr)
    {
      delete cs;
    }
}


auto cluster_session_init(struct cluster_session_s * cs, struct Parameters const & parameters,
                          struct Dbindex & dbindex) -> void
{
  /* Initialize clustering session for library use.
     Reads configuration from the passed parameters (same one given to
     vsearch_session_begin); assumes the database is already set up
     (sequences loaded, masked, and dbindex.prepare called with
     bitmap=1 but WITHOUT add_all_sequences — centroids are indexed
     incrementally as they are discovered, into the passed dbindex). The
     session stores references to parameters and dbindex, which must
     outlive the session.

     The database must be pre-sorted:
     - by length (descending) for cluster_fast behavior
     - by abundance (descending) for cluster_size behavior */

  cs->parameters = &parameters;
  cs->dbindex = &dbindex;

  /* Set search parameters matching the CLI cluster path.
     seqcount must be set BEFORE cluster_query_init (it sizes the kmers buffer).
     MAXDELAYED (8) is needed as safety buffer for align_delayed().
     Clamp tophits to seqcount to avoid oversized allocations.
     The library path does not clamp to the database size (only the CLI
     cluster() does), so the sizing uses the configured values from parameters;
     si->parameters (set in cluster_query_init) carries the same. */
  cs->seqcount = static_cast<int>(db_getsequencecount());
  cs->tophits = static_cast<int>(parameters.opt_maxaccepts + parameters.opt_maxrejects + MAXDELAYED);
  if (cs->tophits > cs->seqcount)
    {
      cs->tophits = cs->seqcount;
    }

  cs->si = make_unique<searchinfo_s>();
  cluster_query_init(cs->si.get(), cs->seqcount, cs->tophits, parameters, *cs->dbindex);
  cs->si->strand = 0;

  if (parameters.opt_strand)
    {
      cs->si_minus = make_unique<searchinfo_s>();
      cluster_query_init(cs->si_minus.get(), cs->seqcount, cs->tophits, parameters, *cs->dbindex);
      cs->si_minus->strand = 1;
    }

  cs->cluster_count = 0;
}


auto cluster_assign_single(struct cluster_session_s * cs,
                            int seqno,
                            struct cluster_result_s * result) -> void
{
  /* Assign a single database sequence to a cluster.
     Must be called in order (seqno 0, 1, 2, ...).
     Each call either assigns the sequence to an existing cluster
     (if a match is found above opt_id threshold) or creates a new
     cluster with this sequence as the centroid.

     New centroids are automatically indexed for subsequent queries. */

  *result = {};
  struct Parameters const & parameters = *cs->parameters;

  cs->si->query_no = seqno;
  cs->si->strand = 0;
  cluster_query_core(cs->si.get(), parameters);

  if (parameters.opt_strand)
    {
      cs->si_minus->query_no = seqno;
      cs->si_minus->strand = 1;
      cluster_query_core(cs->si_minus.get(), parameters);
    }

  struct hit const * best = nullptr;
  if (parameters.opt_sizeorder)
    {
      best = search_findbest2_bysize(cs->si.get(), cs->si_minus.get());
    }
  else
    {
      best = search_findbest2_byid(cs->si.get(), cs->si_minus.get());
    }

  if (best != nullptr)
    {
      /* Match found — assign to existing cluster */
      result->is_centroid = false;
      result->cluster_id = cs->centroid_cluster_ids.at(best->target);
      result->centroid_seqno = best->target;
      result->identity = best->id;
      std::snprintf(result->centroid_label, sizeof(result->centroid_label),
                    "%.*s",
                    static_cast<int>(db_getheaderlen(static_cast<uint64_t>(best->target))),
                    db_getheader(static_cast<uint64_t>(best->target)));
      if (best->nwalignment != nullptr)
        {
          int n = std::snprintf(result->cigar, sizeof(result->cigar), "%s",
                                best->nwalignment);
          result->cigar_truncated =
            (n >= static_cast<int>(sizeof(result->cigar)));
        }
    }
  else
    {
      /* No match — this sequence becomes a new centroid */
      result->is_centroid = true;
      result->cluster_id = cs->cluster_count;
      result->centroid_seqno = seqno;
      result->identity = 100.0;
      std::snprintf(result->centroid_label, sizeof(result->centroid_label),
                    "%.*s",
                    static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno))),
                    db_getheader(static_cast<uint64_t>(seqno)));

      cs->centroid_cluster_ids[seqno] = cs->cluster_count;
      cs->dbindex->add_sequence(static_cast<unsigned int>(seqno), parameters.opt_qmask);
      ++cs->cluster_count;
    }

  /* Free ALL hit alignments for both strands.
     search_onequery / align_delayed allocates nwalignment for each aligned hit. */
  free_hit_alignments(cs->si.get(), cs->si_minus.get(), parameters);
}


auto cluster_assign_batch(struct cluster_session_s * cs,
                          int start_seqno,
                          int count,
                          struct cluster_result_s * results) -> void
{
  if (count <= 0)
    {
      return;
    }

  /* cluster_session_init() captured the database sequence count into
     cs->seqcount and sized this session's search state from the file-static
     seqcount; the per-thread buffers allocated below are sized from the same
     count via cluster_query_init(). If the database has changed since session
     init (a reload, or an interleaved cluster path), those buffers no longer
     match the data being indexed and the k-mer buffer would be sized for the
     wrong count -> out-of-bounds indexing. Fail loudly in release builds (an
     assert would be stripped under NDEBUG) rather than corrupt memory (L2d). */
  if (cs->seqcount != static_cast<int>(db_getsequencecount()))
    {
      fatal("cluster_assign_batch: the database changed since "
            "cluster_session_init(); re-initialize the clustering session.");
    }

  struct Parameters const & parameters = *cs->parameters;

  constexpr int queries_per_thread = 1;
  int const max_queries =
    queries_per_thread * static_cast<int>(parameters.opt_threads);

  /* This batch owns its worker pool and per-thread search state via
     cluster_work_pool_s, rather than borrowing the CLI path's file-static
     si_plus/si_minus/thread_work/cluster_threadrunner through a save/restore
     hack (E4). The constructor allocates and cluster_query_init()s the
     per-thread searchinfo arrays and creates the ThreadRunner; the destructor
     joins the workers and cluster_query_exit()s them. The local si_plus/
     si_minus aliases let the loop below read unchanged. */
  cluster_work_pool_s pool(static_cast<int>(parameters.opt_threads), cs->seqcount,
                           cs->tophits, parameters.opt_strand, parameters, *cs->dbindex);
  searchinfo_s * const si_plus = pool.si_plus.data();
  searchinfo_s * const si_minus = pool.si_minus.empty() ? nullptr : pool.si_minus.data();

  /* Scoring for intra-batch fixup alignment */
  struct Scoring scoring = scoring_from_options(parameters);

  LinearMemoryAligner lma(scoring);

  std::vector<int> extra_list(static_cast<std::size_t>(max_queries));

  int seqno = start_seqno;
  int const end_seqno = start_seqno + count;

  while (seqno < end_seqno)
    {
      /* Load batch of queries */
      int queries = 0;
      for (int i = 0; i < max_queries; i++)
        {
          if (seqno < end_seqno)
            {
              si_plus[i].query_no = seqno;
              si_plus[i].strand = 0;

              if (parameters.opt_strand)
                {
                  si_minus[i].query_no = seqno;
                  si_minus[i].strand = 1;
                }

              ++queries;
              ++seqno;
            }
        }

      /* Parallel search phase */
      pool.wakeup(queries);

      /* Serial analysis with intra-batch fixup
         (adapted from cluster_core_parallel lines 725-1053) */
      int extra_count = 0;

      for (int i = 0; i < queries; i++)
        {
          struct searchinfo_s * si_p = si_plus + i;
          struct searchinfo_s * si_m =
            parameters.opt_strand ? si_minus + i : nullptr;

          for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
            {
              struct searchinfo_s * si = (s != 0) ? si_m : si_p;

              evaluate_extra_hits(si, si_plus, extra_list.data(), extra_count, lma, cs->tophits);
            }

          /* Find best hit across strands */
          struct hit const * best = nullptr;
          if (parameters.opt_sizeorder)
            {
              best = search_findbest2_bysize(si_p, si_m);
            }
          else
            {
              best = search_findbest2_byid(si_p, si_m);
            }

          int const myseqno = si_p->query_no;
          int const ri = myseqno - start_seqno;
          std::memset(&results[ri], 0, sizeof(results[ri]));

          if (best != nullptr)
            {
              /* Match found — assign to existing cluster */
              results[ri].is_centroid = false;
              results[ri].cluster_id =
                cs->centroid_cluster_ids.at(best->target);
              results[ri].centroid_seqno = best->target;
              results[ri].identity = best->id;
              std::snprintf(results[ri].centroid_label,
                            sizeof(results[ri].centroid_label),
                            "%.*s",
                            static_cast<int>(db_getheaderlen(static_cast<uint64_t>(best->target))),
                            db_getheader(static_cast<uint64_t>(best->target)));
              if (best->nwalignment != nullptr)
                {
                  int n = std::snprintf(results[ri].cigar,
                                        sizeof(results[ri].cigar),
                                        "%s", best->nwalignment);
                  results[ri].cigar_truncated =
                    (n >= static_cast<int>(sizeof(results[ri].cigar)));
                }
            }
          else
            {
              /* No match — new centroid */
              extra_list[static_cast<std::size_t>(extra_count)] = i;
              ++extra_count;

              results[ri].is_centroid = true;
              results[ri].cluster_id = cs->cluster_count;
              results[ri].centroid_seqno = myseqno;
              results[ri].identity = 100.0;
              std::snprintf(results[ri].centroid_label,
                            sizeof(results[ri].centroid_label),
                            "%.*s",
                            static_cast<int>(db_getheaderlen(static_cast<uint64_t>(myseqno))),
                            db_getheader(static_cast<uint64_t>(myseqno)));

              cs->centroid_cluster_ids[myseqno] = cs->cluster_count;
              cs->dbindex->add_sequence(static_cast<unsigned int>(myseqno), parameters.opt_qmask);
              ++cs->cluster_count;
            }

          free_hit_alignments(si_p, si_m, parameters);
        }
    }

  /* pool's destructor joins the workers and frees the per-thread search state */
}


auto cluster_session_cleanup(struct cluster_session_s * cs) -> void
{
  if (cs->si)
    {
      cluster_query_exit(cs->si.get());
      cs->si.reset();
    }
  if (cs->si_minus)
    {
      cluster_query_exit(cs->si_minus.get());
      cs->si_minus.reset();
    }
}
