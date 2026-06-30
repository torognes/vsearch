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
#include "align_simd.h"
#include "cluster.h"
#include "searchcore.h"
#include "attributes.h"
#include "dbindex.h"
#include "linmemalign.h"
#include "mask.h"
#include "minheap.h"
#include "msa.h"
#include "otutable.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/threads.hpp"
#include <algorithm>  // std::count, std::minmax_element, std::max_element, std::min
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstdlib>  // std::qsort
#include <cstring>  // std::strcpy, std::strlen
#include <limits>
#include <map>
#include <memory>  // std::unique_ptr
#include <utility>  // std::get
#include <vector>


static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */

struct clusterinfo_s
{
  int seqno;
  int clusterno;
  char * cigar;
  int strand;
};

using clusterinfo_t = struct clusterinfo_s;

static clusterinfo_t * clusterinfo = nullptr;
static int clusters = 0;

static int count_matched = 0;
static int count_notmatched = 0;

static int64_t * cluster_abundance;

static std::FILE * fp_centroids = nullptr;
static std::FILE * fp_uc = nullptr;
static std::FILE * fp_alnout = nullptr;
static std::FILE * fp_samout = nullptr;
static std::FILE * fp_userout = nullptr;
static std::FILE * fp_blast6out = nullptr;
static std::FILE * fp_fastapairs = nullptr;
static std::FILE * fp_matched = nullptr;
static std::FILE * fp_notmatched = nullptr;
static std::FILE * fp_otutabout = nullptr;
static std::FILE * fp_mothur_shared_out = nullptr;
static std::FILE * fp_biomout = nullptr;
static std::FILE * fp_qsegout = nullptr;
static std::FILE * fp_tsegout = nullptr;

static struct searchinfo_s * si_plus;
static struct searchinfo_s * si_minus;

/* per-thread slice of queries assigned for the current round; the
   threading primitives themselves live inside the ThreadRunner */
struct thread_work_s
{
  int query_first;
  int query_count;
};

static std::vector<struct thread_work_s> thread_work;
static std::unique_ptr<ThreadRunner> cluster_threadrunner;


inline auto compare_byclusterno(const void * a, const void * b) -> int
{
  auto const * lhs = static_cast<clusterinfo_t const *>(a);
  auto const * rhs = static_cast<clusterinfo_t const *>(b);

  if (lhs->clusterno < rhs->clusterno)
    {
      return -1;
    }
  if (lhs->clusterno > rhs->clusterno)
    {
      return +1;
    }

  if (lhs->seqno < rhs->seqno)
    {
      return -1;
    }
  if (lhs->seqno > rhs->seqno)
    {
      return +1;
    }
  return 0;
}


inline auto compare_byclusterabundance(const void * a, const void * b) -> int
{
  auto const * lhs = static_cast<clusterinfo_t const *>(a);
  auto const * rhs = static_cast<clusterinfo_t const *>(b);

  if (cluster_abundance[lhs->clusterno] > cluster_abundance[rhs->clusterno])
    {
      return -1;
    }
  if (cluster_abundance[lhs->clusterno] < cluster_abundance[rhs->clusterno])
    {
      return +1;
    }

  if (lhs->clusterno < rhs->clusterno)
    {
      return -1;
    }
  if (lhs->clusterno > rhs->clusterno)
    {
      return +1;
    }

  if (lhs->seqno < rhs->seqno)
    {
      return -1;
    }
  if (lhs->seqno > rhs->seqno)
    {
      return +1;
    }
  return 0;
}


inline auto cluster_query_core(struct searchinfo_s * si) -> void
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
  search_onequery(si, static_cast<int>(opt_qmask));
}


inline auto cluster_worker(uint64_t t) -> void
{
  /* wrapper for the main threaded core function for clustering */
  for (int q = 0; q < thread_work[t].query_count; q++)
    {
      cluster_query_core(si_plus + thread_work[t].query_first + q);
      if (opt_strand > 1)
        {
          cluster_query_core(si_minus + thread_work[t].query_first + q);
        }
    }
}


auto threads_wakeup(int queries) -> void
{
  int const threads = queries > opt_threads ? static_cast<int>(opt_threads) : queries;
  int queries_rest = queries;
  int threads_rest = threads;
  int query_next = 0;

  /* assign a slice of the queries to each thread; threads beyond the
     number of queries get an empty slice and do nothing this round */
  for (int t = 0; t < opt_threads; t++)
    {
      auto const tdx = static_cast<std::size_t>(t);
      if (t < threads)
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

  cluster_threadrunner->run();
}


auto threads_init() -> void
{
  thread_work.resize(static_cast<std::size_t>(opt_threads));
  cluster_threadrunner.reset(new ThreadRunner(static_cast<std::size_t>(opt_threads),
                                              cluster_worker));
}


auto threads_exit() -> void
{
  /* joins and destroys the worker threads */
  cluster_threadrunner.reset();
  thread_work.clear();
}


auto cluster_query_init(struct searchinfo_s * si) -> void
{
  /* initialisation of data for one thread; run once for each thread */
  /* thread specific initialiation */

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
  si->s = search16_init(static_cast<CELL>(opt_match),
                        static_cast<CELL>(opt_mismatch),
                        static_cast<CELL>(opt_gap_open_query_left),
                        static_cast<CELL>(opt_gap_open_target_left),
                        static_cast<CELL>(opt_gap_open_query_interior),
                        static_cast<CELL>(opt_gap_open_target_interior),
                        static_cast<CELL>(opt_gap_open_query_right),
                        static_cast<CELL>(opt_gap_open_target_right),
                        static_cast<CELL>(opt_gap_extension_query_left),
                        static_cast<CELL>(opt_gap_extension_target_left),
                        static_cast<CELL>(opt_gap_extension_query_interior),
                        static_cast<CELL>(opt_gap_extension_target_interior),
                        static_cast<CELL>(opt_gap_extension_query_right),
                        static_cast<CELL>(opt_gap_extension_target_right));
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


auto relabel_otu(int clusterno, char const * sequence, int seqlen) -> char *
{
  char * label = nullptr;
  if (opt_relabel != nullptr)
    {
      int const size = static_cast<int>(std::strlen(opt_relabel)) + 21;
      label = static_cast<char *>(xmalloc(static_cast<std::size_t>(size)));
      std::snprintf(label, static_cast<std::size_t>(size), "%s%d", opt_relabel, clusterno + 1);
    }
  else if (opt_relabel_self)
    {
      int const size = seqlen + 1;
      label = static_cast<char *>(xmalloc(static_cast<std::size_t>(size)));
      std::snprintf(label, static_cast<std::size_t>(size), "%.*s", seqlen, sequence);
    }
  else if (opt_relabel_sha1)
    {
      label = static_cast<char *>(xmalloc(len_hex_dig_sha1));
      get_hex_seq_digest_sha1(label, sequence, seqlen);
    }
  else if (opt_relabel_md5)
    {
      label = static_cast<char *>(xmalloc(len_hex_dig_md5));
      get_hex_seq_digest_md5(label, sequence, seqlen);
    }
  return label;
}


auto cluster_core_results_hit(struct hit const * best,
                              int clusterno,
                              char const * query_head,
                              int qseqlen,
                              char const * qsequence,
                              char const * qsequence_rc,
                              int64_t qsize) -> void
{
  ++count_matched;

  if ((opt_otutabout != nullptr) or (opt_mothur_shared_out != nullptr) or (opt_biomout != nullptr))
    {
      if ((opt_relabel != nullptr) or opt_relabel_self or opt_relabel_sha1 or opt_relabel_md5)
        {
          char * label = relabel_otu(clusterno,
                                     db_getsequence(static_cast<uint64_t>(best->target)),
                                     static_cast<int>(db_getsequencelen(static_cast<uint64_t>(best->target))));
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

  if (fp_uc != nullptr)
    {
      results_show_uc_one(fp_uc,
                          best, query_head,
                          qseqlen,
                          clusterno);
    }

  if (fp_alnout != nullptr)
    {
      results_show_alnout(fp_alnout,
                          best, 1, query_head,
                          qsequence, qseqlen);
    }

  if (fp_samout != nullptr)
    {
      results_show_samout(fp_samout,
                          best, 1, query_head,
                          qsequence, qsequence_rc);
    }

  if (fp_fastapairs != nullptr)
    {
      results_show_fastapairs_one(fp_fastapairs,
                                  best,
                                  query_head,
                                  qsequence,
                                  qsequence_rc);
    }

  if (fp_qsegout != nullptr)
    {
      results_show_qsegout_one(fp_qsegout,
                               best,
                               query_head,
                               qsequence,
                               qseqlen,
                               qsequence_rc);
    }

  if (fp_tsegout != nullptr)
    {
      results_show_tsegout_one(fp_tsegout,
                               best);
    }

  if (fp_userout != nullptr)
    {
      results_show_userout_one(fp_userout, best, query_head,
                               qsequence, qseqlen, qsequence_rc);
    }

  if (fp_blast6out != nullptr)
    {
      results_show_blast6out_one(fp_blast6out, best, query_head,
                                 qseqlen);
    }

  if (opt_matched != nullptr)
    {
      fasta_print_general(fp_matched,
                          nullptr,
                          qsequence,
                          qseqlen,
                          query_head,
                          static_cast<int>(std::strlen(query_head)),
                          static_cast<uint64_t>(qsize),
                          count_matched,
                          -1.0,
                          -1, -1, nullptr, 0.0,
                          0);
    }
}


auto cluster_core_results_nohit(int clusterno,
                                char const * query_head,
                                int qseqlen,
                                char const * qsequence,
                                char const * qsequence_rc,
                                int64_t qsize) -> void
{
  ++count_notmatched;

  if ((opt_otutabout != nullptr) or (opt_mothur_shared_out != nullptr) or (opt_biomout != nullptr))
    {
      if ((opt_relabel != nullptr) or opt_relabel_self or opt_relabel_sha1 or opt_relabel_md5)
        {
          char * label = relabel_otu(clusterno, qsequence, qseqlen);
          otutable_add(query_head, label, qsize);
          xfree(label);
        }
      else
        {
          otutable_add(query_head, query_head, qsize);
        }
    }

  if (opt_uc != nullptr)
    {
      std::fprintf(fp_uc, "S\t%d\t%d\t*\t*\t*\t*\t*\t", clusters, qseqlen);
      header_fprint_strip(fp_uc,
                          query_head,
                          static_cast<int>(std::strlen(query_head)),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      std::fprintf(fp_uc, "\t*\n");
    }

  if (opt_output_no_hits != 0)
    {
      if (fp_userout != nullptr)
        {
          results_show_userout_one(fp_userout, nullptr, query_head,
                                   qsequence, qseqlen, qsequence_rc);
        }

      if (fp_blast6out != nullptr)
        {
          results_show_blast6out_one(fp_blast6out, nullptr, query_head,
                                     qseqlen);
        }
    }

  if (opt_notmatched != nullptr)
    {
      fasta_print_general(fp_notmatched,
                          nullptr,
                          qsequence,
                          qseqlen,
                          query_head,
                          static_cast<int>(std::strlen(query_head)),
                          static_cast<uint64_t>(qsize),
                          count_notmatched,
                          -1.0,
                          -1, -1, nullptr, 0.0,
                          0);
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
                                const int * extra_list,
                                int extra_count,
                                LinearMemoryAligner & lma) -> void
{
  int added = 0;

  /* Keep at most this many hits. The list is the tophits-sized si->hits buffer,
     but tophits is clamped to seqcount; on a small dataset with large
     --maxaccepts/--maxrejects the raw bound (maxaccepts + maxrejects - 1) can
     exceed tophits, so the insertion/shift below would write past the buffer.
     Clamp the bound to the actual capacity (S10). */
  int const hit_capacity =
    static_cast<int>(std::min<int64_t>(opt_maxaccepts + opt_maxrejects - 1,
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
                                  static_cast<int>(opt_wordlength),
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
           (si->accepts < opt_maxaccepts) and
             (si->rejects < opt_maxrejects) and
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
                  align_trim(hit);
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
                                struct searchinfo_s * si_m) -> void
{
  /* free alignments */
  for (int s = 0; s < opt_strand; s++)
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

auto cluster_core_parallel() -> void
{
  /* create threads and set them in stand-by mode */
  threads_init();

  constexpr static int queries_per_thread = 1;
  const int max_queries = queries_per_thread * static_cast<int>(opt_threads);

  /* allocate memory for the search information for each query;
     and initialize it */
  si_plus = new searchinfo_s[max_queries]{};
  if (opt_strand > 1)
    {
      si_minus = new searchinfo_s[max_queries]{};
    }
  for (int i = 0; i < max_queries; i++)
    {
      cluster_query_init(si_plus + i);
      si_plus[i].strand = 0;
      if (opt_strand > 1)
        {
          cluster_query_init(si_minus + i);
          si_minus[i].strand = 1;
        }
    }

  std::vector<int> extra_list(static_cast<std::size_t>(max_queries));

  struct Scoring scoring;
  scoring.match = opt_match;
  scoring.mismatch = opt_mismatch;
  scoring.gap_open_query_left = opt_gap_open_query_left;
  scoring.gap_open_target_left = opt_gap_open_target_left;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_open_target_interior = opt_gap_open_target_interior;
  scoring.gap_open_query_right = opt_gap_open_query_right;
  scoring.gap_open_target_right = opt_gap_open_target_right;
  scoring.gap_extension_query_left = opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opt_gap_extension_target_left;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_extension_target_interior = opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opt_gap_extension_target_right;


  LinearMemoryAligner lma(scoring);


  auto lastlength = std::numeric_limits<int>::max();

  int seqno = 0;

  int64_t sum_nucleotides = 0;

  progress_init("Clustering", db_getnucleotidecount());

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

#if 1
              if ((opt_cluster_smallmem != nullptr) and (opt_usersort == 0) and (length > lastlength))
                {
                  fatal("Sequences not sorted by length and --usersort not specified.");
                }
#endif

              lastlength = length;

              si_plus[i].query_no = seqno;
              si_plus[i].strand = 0;

              if (opt_strand > 1)
                {
                  si_minus[i].query_no = seqno;
                  si_minus[i].strand = 1;
                }

              ++queries;
              ++seqno;
            }
        }

      /* perform work in threads */
      threads_wakeup(queries);

      /* analyse results */
      int extra_count = 0;

      for (int i = 0; i < queries; i++)
        {
          struct searchinfo_s * si_p = si_plus + i;
          struct searchinfo_s * si_m = opt_strand > 1 ? si_minus + i : nullptr;

          for (int s = 0; s < opt_strand; s++)
            {
              struct searchinfo_s * si = (s != 0) ? si_m : si_p;

              evaluate_extra_hits(si, extra_list.data(), extra_count, lma);
            }

          /* find best hit */
          struct hit * best = nullptr;
          if (opt_sizeorder)
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
              cluster_core_results_hit(best,
                                       clusterinfo[target].clusterno,
                                       si_p->query_head,
                                       si_p->qseqlen,
                                       si_p->qsequence,
                                       (best->strand != 0) ? si_m->qsequence : nullptr,
                                       si_p->qsize);

              /* update cluster info about this sequence */
              clusterinfo[myseqno].seqno = myseqno;
              clusterinfo[myseqno].clusterno = clusterinfo[target].clusterno;
              clusterinfo[myseqno].cigar = best->nwalignment;
              clusterinfo[myseqno].strand = best->strand;
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
              clusterinfo[myseqno].seqno = myseqno;
              clusterinfo[myseqno].clusterno = clusters;
              clusterinfo[myseqno].cigar = nullptr;
              clusterinfo[myseqno].strand = 0;

              /* add current sequence to database */
              dbindex_addsequence(static_cast<unsigned int>(myseqno), static_cast<int>(opt_qmask));

              /* output intermediate results to uc etc */
              cluster_core_results_nohit(clusters,
                                         si_p->query_head,
                                         si_p->qseqlen,
                                         si_p->qsequence,
                                         nullptr,
                                         si_p->qsize);
              ++clusters;
            }

          free_hit_alignments(si_p, si_m);

          sum_nucleotides += si_p->qseqlen;
        }

      progress_update(static_cast<uint64_t>(sum_nucleotides));
    }
  progress_done();

  /* clean up search info */
  for (int i = 0; i < max_queries; i++)
    {
      cluster_query_exit(si_plus + i);
      if (opt_strand > 1)
        {
          cluster_query_exit(si_minus + i);
        }
    }

  // extra_list no used after that point

  delete [] si_plus;
  if (opt_strand > 1)
    {
      delete [] si_minus;
    }

  /* terminate threads and clean up */
  threads_exit();

}


auto cluster_core_serial() -> void
{
  std::array<struct searchinfo_s, 1> si_p {{}};  // refactoring: direct initialization?
  std::array<struct searchinfo_s, 1> si_m {{}};

  cluster_query_init(si_p.data());
  if (opt_strand > 1)
    {
      cluster_query_init(si_m.data());
    }

  auto lastlength = std::numeric_limits<int>::max();

  progress_init("Clustering", static_cast<uint64_t>(seqcount));
  for (int seqno = 0; seqno < seqcount; seqno++)
    {
      int const length = static_cast<int>(db_getsequencelen(static_cast<uint64_t>(seqno)));

#if 1
      if ((opt_cluster_smallmem != nullptr) and (opt_usersort == 0) and (length > lastlength))
        {
          fatal("Sequences not sorted by length and --usersort not specified.");
        }
#endif

      lastlength = length;

      si_p[0].query_no = seqno;
      si_p[0].strand = 0;
      cluster_query_core(si_p.data());

      if (opt_strand > 1)
        {
          si_m[0].query_no = seqno;
          si_m[0].strand = 1;
          cluster_query_core(si_m.data());
        }

      struct hit * best = nullptr;
      if (opt_sizeorder)
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
          cluster_core_results_hit(best,
                                   clusterinfo[target].clusterno,
                                   si_p[0].query_head,
                                   si_p[0].qseqlen,
                                   si_p[0].qsequence,
                                   (best->strand != 0) ? si_m[0].qsequence : nullptr,
                                   si_p[0].qsize);
          clusterinfo[seqno].seqno = seqno;
          clusterinfo[seqno].clusterno = clusterinfo[target].clusterno;
          clusterinfo[seqno].cigar = best->nwalignment;
          clusterinfo[seqno].strand = best->strand;
          best->nwalignment = nullptr;
        }
      else
        {
          clusterinfo[seqno].seqno = seqno;
          clusterinfo[seqno].clusterno = clusters;
          clusterinfo[seqno].cigar = nullptr;
          clusterinfo[seqno].strand = 0;
          dbindex_addsequence(static_cast<unsigned int>(seqno), static_cast<int>(opt_qmask));
          cluster_core_results_nohit(clusters,
                                     si_p[0].query_head,
                                     si_p[0].qseqlen,
                                     si_p[0].qsequence,
                                     nullptr,
                                     si_p[0].qsize);
          ++clusters;
        }

      free_hit_alignments(si_p.data(), si_m.data());

      progress_update(static_cast<uint64_t>(seqno));
    }
  progress_done();

  cluster_query_exit(si_p.data());
  if (opt_strand > 1)
    {
      cluster_query_exit(si_m.data());
    }
}


auto cluster(char const * dbname,
             char const * cmdline,
             char const * progheader) -> void
{
  if (opt_centroids != nullptr)
    {
      fp_centroids = fopen_output(opt_centroids);
      if (fp_centroids == nullptr)
        {
          fatal("Unable to open centroids file for writing");
        }
    }

  if (opt_uc != nullptr)
    {
      fp_uc = fopen_output(opt_uc);
      if (fp_uc == nullptr)
        {
          fatal("Unable to open uc file for writing");
        }
    }

  if (opt_alnout != nullptr)
    {
      fp_alnout = fopen_output(opt_alnout);
      if (fp_alnout == nullptr)
        {
          fatal("Unable to open alignment output file for writing");
        }

      std::fprintf(fp_alnout, "%s\n", cmdline);
      std::fprintf(fp_alnout, "%s\n", progheader);
    }

  if (opt_samout != nullptr)
    {
      fp_samout = fopen_output(opt_samout);
      if (fp_samout == nullptr)
        {
          fatal("Unable to open SAM output file for writing");
        }
    }

  if (opt_userout != nullptr)
    {
      fp_userout = fopen_output(opt_userout);
      if (fp_userout == nullptr)
        {
          fatal("Unable to open user-defined output file for writing");
        }
    }

  if (opt_blast6out != nullptr)
    {
      fp_blast6out = fopen_output(opt_blast6out);
      if (fp_blast6out == nullptr)
        {
          fatal("Unable to open blast6-like output file for writing");
        }
    }

  if (opt_fastapairs != nullptr)
    {
      fp_fastapairs = fopen_output(opt_fastapairs);
      if (fp_fastapairs == nullptr)
        {
          fatal("Unable to open fastapairs output file for writing");
        }
    }

  if (opt_qsegout != nullptr)
    {
      fp_qsegout = fopen_output(opt_qsegout);
      if (fp_qsegout == nullptr)
        {
          fatal("Unable to open qsegout output file for writing");
        }
    }

  if (opt_tsegout != nullptr)
    {
      fp_tsegout = fopen_output(opt_tsegout);
      if (fp_tsegout == nullptr)
        {
          fatal("Unable to open tsegout output file for writing");
        }
    }

  if (opt_matched != nullptr)
    {
      fp_matched = fopen_output(opt_matched);
      if (fp_matched == nullptr)
        {
          fatal("Unable to open matched output file for writing");
        }
    }

  if (opt_notmatched != nullptr)
    {
      fp_notmatched = fopen_output(opt_notmatched);
      if (fp_notmatched == nullptr)
        {
          fatal("Unable to open notmatched output file for writing");
        }
    }

  if (opt_otutabout != nullptr)
    {
      fp_otutabout = fopen_output(opt_otutabout);
      if (fp_otutabout == nullptr)
        {
          fatal("Unable to open OTU table (text format) output file for writing");
        }
    }

  if (opt_mothur_shared_out != nullptr)
    {
      fp_mothur_shared_out = fopen_output(opt_mothur_shared_out);
      if (fp_mothur_shared_out == nullptr)
        {
          fatal("Unable to open OTU table (mothur format) output file for writing");
        }
    }

  if (opt_biomout != nullptr)
    {
      fp_biomout = fopen_output(opt_biomout);
      if (fp_biomout == nullptr)
        {
          fatal("Unable to open OTU table (biom 1.0 format) output file for writing");
        }
    }

  db_read(dbname, 0);

  otutable_init();

  results_show_samheader(fp_samout, cmdline, dbname);

  if (opt_qmask == MASK_DUST)
    {
      dust_all();
    }
  else if ((opt_qmask == MASK_SOFT) and (opt_hardmask != 0))
    {
      hardmask_all();
    }

  show_rusage();

  seqcount = static_cast<int>(db_getsequencecount());

  if (opt_cluster_fast != nullptr)
    {
      db_sortbylength();
    }
  else if ((opt_cluster_size != nullptr) or (opt_cluster_unoise != nullptr))
    {
      db_sortbyabundance();
    }

  dbindex_prepare(1, static_cast<int>(opt_qmask));

  /* tophits = the maximum number of hits we need to store */

  if ((opt_maxrejects == 0) or (opt_maxrejects > seqcount))
    {
      opt_maxrejects = seqcount;
    }

  if ((opt_maxaccepts == 0) or (opt_maxaccepts > seqcount))
    {
      opt_maxaccepts = seqcount;
    }

  tophits = static_cast<int>(opt_maxrejects + opt_maxaccepts + MAXDELAYED);
  tophits = std::min(tophits, seqcount);

  std::vector<clusterinfo_t> clusterinfo_v(static_cast<std::size_t>(seqcount));
  clusterinfo = clusterinfo_v.data();

  if (opt_log != nullptr)
    {
      uint64_t const slots = 1ULL << (static_cast<uint64_t>(opt_wordlength) << 1ULL);
      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "      Alphabet  nt\n");
      std::fprintf(fp_log, "    Word width  %" PRId64 "\n", opt_wordlength);
      std::fprintf(fp_log, "     Word ones  %" PRId64 "\n", opt_wordlength);
      std::fprintf(fp_log, "        Spaced  No\n");
      std::fprintf(fp_log, "        Hashed  No\n");
      std::fprintf(fp_log, "         Coded  No\n");
      std::fprintf(fp_log, "       Stepped  No\n");
      std::fprintf(fp_log, "         Slots  %" PRIu64 " (%.1fk)\n", slots, static_cast<double>(slots)/1000.0);
      std::fprintf(fp_log, "       DBAccel  100%%\n");
      std::fprintf(fp_log, "\n");
    }

  if (opt_threads == 1)
    {
      cluster_core_serial();
    }
  else
    {
      cluster_core_parallel();
    }


  /* find size and abundance of each cluster and save stats */

  std::vector<int64_t> cluster_abundance_v(static_cast<std::size_t>(clusters));
  cluster_abundance = cluster_abundance_v.data();
  std::vector<int> cluster_size_v(static_cast<std::size_t>(clusters));

  for (int i = 0; i < seqcount; i++)
    {
      int const seqno = clusterinfo_v[static_cast<std::size_t>(i)].seqno;
      int const clusterno = clusterinfo_v[static_cast<std::size_t>(i)].clusterno;
      cluster_abundance_v[static_cast<std::size_t>(clusterno)] +=
        opt_sizein ? static_cast<int64_t>(db_getabundance(static_cast<uint64_t>(seqno))) : 1;
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

  progress_init("Sorting clusters", static_cast<uint64_t>(clusters));
  if (seqcount > 0)  // qsort requires a non-null pointer even for zero elements
    {
      if (opt_clusterout_sort)
        {
          std::qsort(clusterinfo_v.data(), static_cast<std::size_t>(seqcount), sizeof(clusterinfo_t),
                compare_byclusterabundance);
        }
      else
        {
          std::qsort(clusterinfo_v.data(), static_cast<std::size_t>(seqcount), sizeof(clusterinfo_t),
                compare_byclusterno);
        }
    }
  progress_done();

  progress_init("Writing clusters", static_cast<uint64_t>(seqcount));

  /* allocate memory for full file name of the clusters files */
  std::FILE * fp_clusters = nullptr;
  static constexpr auto space_for_cluster_id = 25;  // up to 25 digits
  std::vector<char> fn_clusters;
  if (opt_clusters != nullptr) {
    fn_clusters.reserve(std::strlen(opt_clusters) + space_for_cluster_id);
  }

  int lastcluster = -1;
  uint64_t ordinal = 0;

  for (int i = 0; i < seqcount; i++)
    {
      int const seqno = clusterinfo_v[static_cast<std::size_t>(i)].seqno;
      int const clusterno = clusterinfo_v[static_cast<std::size_t>(i)].clusterno;

      if (clusterno != lastcluster)
        {
          /* prepare for new cluster */
          /* performed with first sequence only in each cluster */
          /* the first sequence is always the centroid */

          if (opt_centroids != nullptr)
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
                                  opt_clusterout_id ? clusterno : -1,
                                  nullptr, 0.0,
                                  db_getabundance(static_cast<uint64_t>(seqno)));
            }

          if (opt_uc != nullptr)
            {
              std::fprintf(fp_uc, "C\t%d\t%" PRId64 "\t*\t*\t*\t*\t*\t",
                      clusterno,
                      cluster_abundance_v[static_cast<std::size_t>(clusterno)]);
              header_fprint_strip(fp_uc,
                                  db_getheader(static_cast<uint64_t>(seqno)),
                                  static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno))),
                                  opt_xsize,
                                  opt_xee,
                                  opt_xlength);
              std::fprintf(fp_uc, "\t*\n");
            }

          if (opt_clusters != nullptr)
            {
              /* close previous (except for first time) and open new file */
              if (lastcluster != -1)
                {
                  std::fclose(fp_clusters);
                }

              ordinal = 0;
              std::snprintf(fn_clusters.data(),
                       fn_clusters.capacity(),
                       "%s%d",
                       opt_clusters,
                       clusterno);
              fp_clusters = fopen_output(fn_clusters.data());
              if (fp_clusters == nullptr)
                {
                  fatal("Unable to open clusters file for writing");
                }
            }

          lastcluster = clusterno;
        }

      /* performed for all sequences */

      if (opt_clusters != nullptr)
        {
          ++ordinal;
          fasta_print_db_relabel(fp_clusters, static_cast<uint64_t>(seqno), ordinal);
        }

      progress_update(static_cast<uint64_t>(i));
    }

  if (lastcluster != -1)
    {
      /* performed with the last sequence */
      if (opt_clusters != nullptr)
        {
          std::fclose(fp_clusters);
        }
    }

  progress_done();

  if (clusters < 1)
    {
      if (not opt_quiet)
        {
          std::fprintf(stderr, "Clusters: 0\n");
          std::fprintf(stderr, "Singletons: 0\n");
        }
      if (opt_log != nullptr)
        {
          std::fprintf(fp_log, "Clusters: 0\n");
          std::fprintf(fp_log, "Singletons: 0\n");
        }
    }
  else
    {
      if (not opt_quiet)
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

      if (opt_log != nullptr)
        {
          std::fprintf(fp_log,
                  "Clusters: %d Size min %" PRId64 ", max %" PRId64 ", avg %.1f\n",
                  clusters,
                  abundance_min,
                  abundance_max,
                  1.0 * seqcount / clusters);
          std::fprintf(fp_log,
                  "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
                  singletons,
                  100.0 * singletons / seqcount,
                  100.0 * singletons / clusters);
          std::fprintf(fp_log, "\n");
        }
    }

  if ((opt_msaout != nullptr) or (opt_consout != nullptr) or (opt_profile != nullptr))
    {
      int msa_target_count = 0;
      std::vector<struct msa_target_s> msa_target_list_v(static_cast<std::size_t>(size_max));
      progress_init("Multiple alignments", static_cast<uint64_t>(seqcount));

      std::FILE * fp_msaout = nullptr;
      std::FILE * fp_consout = nullptr;
      std::FILE * fp_profile = nullptr;

      if (opt_msaout != nullptr)
        {
          fp_msaout = fopen_output(opt_msaout);
          if (fp_msaout == nullptr)
            {
              fatal("Unable to open msaout file");
            }
        }

      if (opt_consout != nullptr)
        {
          fp_consout = fopen_output(opt_consout);
          if (fp_consout == nullptr)
            {
              fatal("Unable to open consout file");
            }
        }

      if (opt_profile != nullptr)
        {
          fp_profile = fopen_output(opt_profile);
          if (fp_profile == nullptr)
            {
              fatal("Unable to open profile file");
            }
        }

      lastcluster = -1;

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
                      cluster_abundance_v[static_cast<std::size_t>(lastcluster)]);
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

          progress_update(static_cast<uint64_t>(i));
        }

      if (lastcluster != -1)
        {
          /* compute msa & consensus */
          msa(fp_msaout, fp_consout, fp_profile,
              lastcluster,
              msa_target_count, msa_target_list_v,
              cluster_abundance_v[static_cast<std::size_t>(lastcluster)]);
        }

      progress_done();

      if (fp_profile != nullptr)
        {
          std::fclose(fp_profile);
        }

      if (fp_msaout != nullptr)
        {
          std::fclose(fp_msaout);
        }

      if (fp_consout != nullptr)
        {
          std::fclose(fp_consout);
        }
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
      otutable_print_biomout(fp_biomout);
      std::fclose(fp_biomout);
    }

  if (fp_otutabout != nullptr)
    {
      otutable_print_otutabout(fp_otutabout);
      std::fclose(fp_otutabout);
    }

  if (fp_mothur_shared_out != nullptr)
    {
      otutable_print_mothur_shared_out(fp_mothur_shared_out);
      std::fclose(fp_mothur_shared_out);
    }

  otutable_done();

  if (opt_matched != nullptr)
    {
      std::fclose(fp_matched);
    }
  if (opt_notmatched != nullptr)
    {
      std::fclose(fp_notmatched);
    }
  if (opt_fastapairs != nullptr)
    {
      std::fclose(fp_fastapairs);
    }
  if (opt_qsegout != nullptr)
    {
      std::fclose(fp_qsegout);
    }
  if (opt_tsegout != nullptr)
    {
      std::fclose(fp_tsegout);
    }
  if (fp_blast6out != nullptr)
    {
      std::fclose(fp_blast6out);
    }
  if (fp_userout != nullptr)
    {
      std::fclose(fp_userout);
      clean_up(); // free userfields allocation
    }
  if (fp_alnout != nullptr)
    {
      std::fclose(fp_alnout);
    }
  if (fp_samout != nullptr)
    {
      std::fclose(fp_samout);
    }
  if (fp_uc != nullptr)
    {
      std::fclose(fp_uc);
    }
  if (fp_centroids != nullptr)
    {
      std::fclose(fp_centroids);
    }

  dbindex_free();
  db_free();
  show_rusage();
}


auto cluster_fast(char const * cmdline, char const * progheader) -> void
{
  cluster(opt_cluster_fast, cmdline, progheader);
}


auto cluster_smallmem(char const * cmdline, char const * progheader) -> void
{
  cluster(opt_cluster_smallmem, cmdline, progheader);
}


auto cluster_size(char const * cmdline, char const * progheader) -> void
{
  cluster(opt_cluster_size, cmdline, progheader);
}


auto cluster_unoise(char const * cmdline, char const * progheader) -> void
{
  cluster(opt_cluster_unoise, cmdline, progheader);
}


/* === Library API for embedding clustering === */


struct cluster_session_s {
  struct searchinfo_s * si = nullptr;
  struct searchinfo_s * si_minus = nullptr;  /* non-null when opt_strand > 1 */
  int cluster_count = 0;
  int seqcount = 0;
  std::map<int, int> centroid_cluster_ids;  /* seqno → cluster_id for centroids */
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


auto cluster_session_init(struct cluster_session_s * cs) -> void
{
  /* Initialize clustering session for library use.
     Assumes global opt_* variables and database are already set up
     (sequences loaded, masked, and dbindex_prepare called with
     bitmap=1 but WITHOUT dbindex_addallsequences — centroids are
     indexed incrementally as they are discovered).

     The database must be pre-sorted:
     - by length (descending) for cluster_fast behavior
     - by abundance (descending) for cluster_size behavior */

  /* Set search parameters matching the CLI cluster path.
     seqcount must be set BEFORE cluster_query_init (it sizes the kmers buffer).
     MAXDELAYED (8) is needed as safety buffer for align_delayed().
     Clamp tophits to seqcount to avoid oversized allocations. */
  seqcount = static_cast<int>(db_getsequencecount());
  tophits = static_cast<int>(opt_maxaccepts + opt_maxrejects + MAXDELAYED);
  if (tophits > seqcount)
    {
      tophits = seqcount;
    }

  cs->si = new searchinfo_s {};
  cluster_query_init(cs->si);
  cs->si->strand = 0;

  if (opt_strand > 1)
    {
      cs->si_minus = new searchinfo_s {};
      cluster_query_init(cs->si_minus);
      cs->si_minus->strand = 1;
    }

  cs->cluster_count = 0;
  cs->seqcount = seqcount;
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

  cs->si->query_no = seqno;
  cs->si->strand = 0;
  cluster_query_core(cs->si);

  if (opt_strand > 1)
    {
      cs->si_minus->query_no = seqno;
      cs->si_minus->strand = 1;
      cluster_query_core(cs->si_minus);
    }

  struct hit const * best = nullptr;
  if (opt_sizeorder)
    {
      best = search_findbest2_bysize(cs->si, cs->si_minus);
    }
  else
    {
      best = search_findbest2_byid(cs->si, cs->si_minus);
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
      dbindex_addsequence(static_cast<unsigned int>(seqno), static_cast<int>(opt_qmask));
      ++cs->cluster_count;
    }

  /* Free ALL hit alignments for both strands.
     search_onequery / align_delayed allocates nwalignment for each aligned hit. */
  free_hit_alignments(cs->si, cs->si_minus);
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

  constexpr int queries_per_thread = 1;
  int const max_queries =
    queries_per_thread * static_cast<int>(opt_threads);

  /* Save file-statics used by cluster_worker / threads infrastructure */
  struct searchinfo_s * saved_si_plus = si_plus;
  struct searchinfo_s * saved_si_minus = si_minus;

  /* Allocate per-thread search state for the batch */
  si_plus = new searchinfo_s[max_queries]{};
  if (opt_strand > 1)
    {
      si_minus = new searchinfo_s[max_queries]{};
    }
  else
    {
      si_minus = nullptr;
    }

  for (int i = 0; i < max_queries; i++)
    {
      cluster_query_init(si_plus + i);
      si_plus[i].strand = 0;
      if (opt_strand > 1)
        {
          cluster_query_init(si_minus + i);
          si_minus[i].strand = 1;
        }
    }

  /* Scoring for intra-batch fixup alignment */
  struct Scoring scoring;
  scoring.match = opt_match;
  scoring.mismatch = opt_mismatch;
  scoring.gap_open_query_left = opt_gap_open_query_left;
  scoring.gap_open_target_left = opt_gap_open_target_left;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_open_target_interior = opt_gap_open_target_interior;
  scoring.gap_open_query_right = opt_gap_open_query_right;
  scoring.gap_open_target_right = opt_gap_open_target_right;
  scoring.gap_extension_query_left = opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opt_gap_extension_target_left;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_extension_target_interior = opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opt_gap_extension_target_right;

  LinearMemoryAligner lma(scoring);

  std::vector<int> extra_list(static_cast<std::size_t>(max_queries));

  /* Create thread pool */
  threads_init();

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

              if (opt_strand > 1)
                {
                  si_minus[i].query_no = seqno;
                  si_minus[i].strand = 1;
                }

              ++queries;
              ++seqno;
            }
        }

      /* Parallel search phase */
      threads_wakeup(queries);

      /* Serial analysis with intra-batch fixup
         (adapted from cluster_core_parallel lines 725-1053) */
      int extra_count = 0;

      for (int i = 0; i < queries; i++)
        {
          struct searchinfo_s * si_p = si_plus + i;
          struct searchinfo_s * si_m =
            opt_strand > 1 ? si_minus + i : nullptr;

          for (int s = 0; s < opt_strand; s++)
            {
              struct searchinfo_s * si = (s != 0) ? si_m : si_p;

              evaluate_extra_hits(si, extra_list.data(), extra_count, lma);
            }

          /* Find best hit across strands */
          struct hit const * best = nullptr;
          if (opt_sizeorder)
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
              dbindex_addsequence(static_cast<unsigned int>(myseqno), static_cast<int>(opt_qmask));
              ++cs->cluster_count;
            }

          free_hit_alignments(si_p, si_m);
        }
    }

  /* Destroy thread pool */
  threads_exit();

  /* Clean up batch search state */
  for (int i = 0; i < max_queries; i++)
    {
      cluster_query_exit(si_plus + i);
      if (opt_strand > 1)
        {
          cluster_query_exit(si_minus + i);
        }
    }
  delete [] si_plus;
  if (si_minus != nullptr)
    {
      delete [] si_minus;
    }

  /* Restore file-statics */
  si_plus = saved_si_plus;
  si_minus = saved_si_minus;
}


auto cluster_session_cleanup(struct cluster_session_s * cs) -> void
{
  if (cs->si != nullptr)
    {
      cluster_query_exit(cs->si);
      delete cs->si;
      cs->si = nullptr;
    }
  if (cs->si_minus != nullptr)
    {
      cluster_query_exit(cs->si_minus);
      delete cs->si_minus;
      cs->si_minus = nullptr;
    }
}
