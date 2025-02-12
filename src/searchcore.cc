/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "dbindex.h"
#include "maps.h"
#include "minheap.h"
#include "otutable.h"
#include "unique.h"
#include <algorithm>  // std::min, std::max
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::sscanf
#include <cstdlib>  // std::qsort
#include <cstring>  // std::strlen, std::memset, std::strcmp
#include <limits>


/* per thread data */

inline auto hit_compare_byid_typed(struct hit * lhs, struct hit * rhs) -> int
{
  /*
    Order:
    accepted, then rejected (weak)
    high id, then low id
    early target, then late target
  */

  if (lhs->rejected < rhs->rejected)
    {
      return -1;
    }
  if (lhs->rejected > rhs->rejected)
    {
      return +1;
    }
  if (lhs->aligned > rhs->aligned)
    {
      return -1;
    }
  if (lhs->aligned < rhs->aligned)
    {
      return +1;
    }
  if (lhs->aligned == 0)
    {
      return 0;
    }
  if (lhs->id > rhs->id)
    {
      return -1;
    }
  if (lhs->id < rhs->id)
    {
      return +1;
    }
  if (lhs->target < rhs->target)
    {
      return -1;
    }
  if (lhs->target > rhs->target)
    {
      return +1;
    }
  return 0;
}


inline auto hit_compare_bysize_typed(struct hit * lhs, struct hit * rhs) -> int
{
  // high abundance, then low abundance
  // high id, then low id
  // early target, then late target

  if (lhs->rejected < rhs->rejected)
    {
      return -1;
    }
  if (lhs->rejected > rhs->rejected)
    {
      return +1;
    }
  if (lhs->rejected == 1)
    {
      return 0;
    }

  if (lhs->aligned > rhs->aligned)
    {
      return -1;
    }
  if (lhs->aligned < rhs->aligned)
    {
      return +1;
    }
  if (lhs->aligned == 0)
    {
      return 0;
    }

  auto const lhs_abundance = db_getabundance(lhs->target);
  auto const rhs_abundance = db_getabundance(rhs->target);
  if (lhs_abundance > rhs_abundance)
    {
      return -1;
    }
  if (lhs_abundance < rhs_abundance)
    {
      return +1;
    }

  if (lhs->id > rhs->id)
    {
      return -1;
    }
  if (lhs->id < rhs->id)
    {
      return +1;
    }

  if (lhs->target < rhs->target)
    {
      return -1;
    }
  if (lhs->target > rhs->target)
    {
      return +1;
    }
  return 0;
}


auto hit_compare_byid(const void * lhs, const void * rhs) -> int
{
  return hit_compare_byid_typed((struct hit *) lhs, (struct hit *) rhs);
}


auto hit_compare_bysize(const void * lhs, const void * rhs) -> int
{
  return hit_compare_bysize_typed((struct hit *) lhs, (struct hit *) rhs);
}


auto search_enough_kmers(struct searchinfo_s * searchinfo,
                         unsigned int const count) -> bool
{
  return (count >= opt_minwordmatches) or (count >= searchinfo->kmersamplecount);
}


auto search_topscores(struct searchinfo_s * searchinfo) -> void
{
  /*
    Count the kmer hits in each database sequence and
    make a sorted list of a given number (th)
    of the database sequences with the highest number of matching kmers.
    These are stored in the min heap array.
  */

  /* count kmer hits in the database sequences */
  const int indexed_count = dbindex_getcount();

  /* zero counts */
  std::memset(searchinfo->kmers, 0, indexed_count * sizeof(count_t));

  minheap_empty(searchinfo->m);

  for (auto i = 0U; i < searchinfo->kmersamplecount; i++)
    {
      auto const kmer = searchinfo->kmersample[i];
      auto * bitmap = dbindex_getbitmap(kmer);

      if (bitmap != nullptr)
        {
#ifdef __x86_64__
          if (ssse3_present != 0)
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
          increment_counters_from_bitmap(si->kmers, bitmap, indexed_count);
#endif
        }
      else
        {
          auto * list = dbindex_getmatchlist(kmer);
          auto const count = dbindex_getmatchcount(kmer);
          for (auto j = 0U; j < count; j++)
            {
              searchinfo->kmers[list[j]]++;
            }
        }
    }

  auto const minmatches = std::min(static_cast<unsigned int>(opt_minwordmatches), searchinfo->kmersamplecount);

  for (auto i = 0; i < indexed_count; i++)
    {
      auto const count = searchinfo->kmers[i];
      if (count >= minmatches)
        {
          auto const seqno = dbindex_getmapping(i);
          unsigned int const length = db_getsequencelen(seqno);

          elem_t novel;
          novel.count = count;
          novel.seqno = seqno;
          novel.length = length;

          minheap_add(searchinfo->m, & novel);
        }
    }

  minheap_sort(searchinfo->m);
}


auto seqncmp(char * a, char * b, uint64_t n) -> int
{
  for (auto i = 0U; i < n; i++)
    {
      auto const x = chrmap_4bit[(int)(a[i])];
      auto const y = chrmap_4bit[(int)(b[i])];
      if (x < y)
        {
          return -1;
        }
      if (x > y)
        {
          return +1;
        }
    }
  return 0;
}


auto align_trim(struct hit * hit) -> void
{
  /* trim alignment and fill in info */
  /* assumes that the hit has been aligned */

  /* info for semi-global alignment (without gaps at ends) */

  hit->trim_aln_left = 0;
  hit->trim_q_left = 0;
  hit->trim_t_left = 0;
  hit->trim_aln_right = 0;
  hit->trim_q_right = 0;
  hit->trim_t_right = 0;

  /* left trim alignment */

  auto * p = hit->nwalignment;
  auto op = '\0';
  int64_t run = 0;
  if (*p != 0)
    {
      run = 1;
      auto scanlength = 0;
      sscanf(p, "%" PRId64 "%n", &run, &scanlength);
      op = *(p + scanlength);
      if (op != 'M')
        {
          hit->trim_aln_left = 1 + scanlength;
          if (op == 'D')
            {
              hit->trim_q_left = run;
            }
          else
            {
              hit->trim_t_left = run;
            }
        }
    }

  /* right trim alignment */

  auto * e = hit->nwalignment + strlen(hit->nwalignment);
  if (e > hit->nwalignment)
    {
      p = e - 1;
      op = *p;
      if (op != 'M')
        {
          while ((p > hit->nwalignment) && (*(p-1) <= '9'))
            {
              --p;
            }
          run = 1;
          sscanf(p, "%" PRId64, &run);
          hit->trim_aln_right = e - p;
          if (op == 'D')
            {
              hit->trim_q_right = run;
            }
          else
            {
              hit->trim_t_right = run;
            }
        }
    }

  if (hit->trim_q_left >= hit->nwalignmentlength)
    {
      hit->trim_q_right = 0;
    }

  if (hit->trim_t_left >= hit->nwalignmentlength)
    {
      hit->trim_t_right = 0;
    }

  hit->internal_alignmentlength = hit->nwalignmentlength
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;

  hit->internal_indels = hit->nwindels
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;

  hit->internal_gaps = hit->nwgaps
    - ((hit->trim_q_left  + hit->trim_t_left)  > 0 ? 1 : 0)
    - ((hit->trim_q_right + hit->trim_t_right) > 0 ? 1 : 0);

  /* CD-HIT */
  hit->id0 = hit->shortest > 0 ? 100.0 * hit->matches / hit->shortest : 0.0;
  /* all diffs */
  hit->id1 = hit->nwalignmentlength > 0 ?
    100.0 * hit->matches / hit->nwalignmentlength : 0.0;
  /* internal diffs */
  hit->id2 = hit->internal_alignmentlength > 0 ?
    100.0 * hit->matches / hit->internal_alignmentlength : 0.0;
  /* Marine Biology Lab */
  hit->id3 = std::max(0.0, 100.0 * (1.0 - (1.0 * (hit->mismatches + hit->nwgaps) /
                                      hit->longest)));
  /* BLAST */
  hit->id4 = hit->nwalignmentlength > 0 ?
    100.0 * hit->matches / hit->nwalignmentlength : 0.0;

  switch (opt_iddef)
    {
    case 0:
      hit->id = hit->id0;
      break;
    case 1:
      hit->id = hit->id1;
      break;
    case 2:
      hit->id = hit->id2;
      break;
    case 3:
      hit->id = hit->id3;
      break;
    case 4:
      hit->id = hit->id4;
      break;
    }
}


auto search_acceptable_unaligned(struct searchinfo_s * searchinfo,
                                 int const target) -> bool
{
  /* consider whether a hit satisfies accept criteria before alignment */

  auto * qseq = searchinfo->qsequence;
  auto * dlabel = db_getheader(target);
  auto * dseq = db_getsequence(target);
  const int64_t dseqlen = db_getsequencelen(target);
  const int64_t tsize = db_getabundance(target);

  if (
      /* maxqsize */
      (searchinfo->qsize <= opt_maxqsize)
      &&
      /* mintsize */
      (tsize >= opt_mintsize)
      &&
      /* minsizeratio */
      (searchinfo->qsize >= opt_minsizeratio * tsize)
      &&
      /* maxsizeratio */
      (searchinfo->qsize <= opt_maxsizeratio * tsize)
      &&
      /* minqt */
      (searchinfo->qseqlen >= opt_minqt * dseqlen)
      &&
      /* maxqt */
      (searchinfo->qseqlen <= opt_maxqt * dseqlen)
      &&
      /* minsl */
      (searchinfo->qseqlen < dseqlen ?
       searchinfo->qseqlen >= opt_minsl * dseqlen :
       dseqlen >= opt_minsl * searchinfo->qseqlen)
      &&
      /* maxsl */
      (searchinfo->qseqlen < dseqlen ?
       searchinfo->qseqlen <= opt_maxsl * dseqlen :
       dseqlen <= opt_maxsl * searchinfo->qseqlen)
      &&
      /* idprefix */
      ((searchinfo->qseqlen >= opt_idprefix) &&
       (dseqlen >= opt_idprefix) &&
       (seqncmp(qseq, dseq, opt_idprefix) == 0))
      &&
      /* idsuffix */
      ((searchinfo->qseqlen >= opt_idsuffix) &&
       (dseqlen >= opt_idsuffix) &&
       (seqncmp(qseq + searchinfo->qseqlen - opt_idsuffix,
                 dseq + dseqlen - opt_idsuffix,
                 opt_idsuffix) == 0))
      &&
      /* self */
      ((opt_self == 0) or (std::strcmp(searchinfo->query_head, dlabel) != 0))
      &&
      /* selfid */
      ((opt_selfid == 0) or
       (searchinfo->qseqlen != dseqlen) or
       (seqncmp(qseq, dseq, searchinfo->qseqlen) != 0))
      )
    {
      /* needs further consideration */
      return true;
    }

  /* reject */
  return false;
}


auto search_acceptable_aligned(struct searchinfo_s * searchinfo,
                               struct hit * hit) -> bool
{
  if (/* weak_id */
      (hit->id >= 100.0 * opt_weak_id) &&
      /* maxsubs */
      (hit->mismatches <= opt_maxsubs) &&
      /* maxgaps */
      (hit->internal_gaps <= opt_maxgaps) &&
      /* mincols */
      (hit->internal_alignmentlength >= opt_mincols) &&
      /* leftjust */
      ((opt_leftjust == 0) or (hit->trim_q_left +
                           hit->trim_t_left == 0)) &&
      /* rightjust */
      ((opt_rightjust == 0) or (hit->trim_q_right +
                            hit->trim_t_right == 0)) &&
      /* query_cov */
      (hit->matches + hit->mismatches >= opt_query_cov * searchinfo->qseqlen) &&
      /* target_cov */
      (hit->matches + hit->mismatches >=
       opt_target_cov * db_getsequencelen(hit->target)) &&
      /* maxid */
      (hit->id <= 100.0 * opt_maxid) &&
      /* mid */
      (100.0 * hit->matches / (hit->matches + hit->mismatches) >= opt_mid) &&
      /* maxdiffs */
      (hit->mismatches + hit->internal_indels <= opt_maxdiffs))
    {
      if (opt_cluster_unoise != nullptr)
        {
          const auto mismatches = hit->mismatches;
          auto const skew = 1.0 * searchinfo->qsize / db_getabundance(hit->target);
          auto const beta = 1.0 / std::pow(2, (1.0 * opt_unoise_alpha * mismatches) + 1);

          if (skew <= beta or mismatches == 0)
            {
              /* accepted */
              hit->accepted = true;
              hit->weak = false;
              return true;
            }
          /* rejected, but weak hit */
          hit->rejected = true;
          hit->weak = true;
          return false;
        }

      if (hit->id >= 100.0 * opt_id)
        {
          /* accepted */
          hit->accepted = true;
          hit->weak = false;
          return true;
        }
      /* rejected, but weak hit */
      hit->rejected = true;
      hit->weak = true;
      return false;
    }

  /* rejected */
  hit->rejected = true;
  hit->weak = false;
  return false;
}


auto align_delayed(struct searchinfo_s * searchinfo) -> void
{
  /* compute global alignment */

  std::array<unsigned int, MAXDELAYED> target_list {{}};
  std::array<CELL, MAXDELAYED> nwscore_list {{}};
  std::array<unsigned short, MAXDELAYED> nwalignmentlength_list {{}};
  std::array<unsigned short, MAXDELAYED> nwmatches_list {{}};
  std::array<unsigned short, MAXDELAYED> nwmismatches_list {{}};
  std::array<unsigned short, MAXDELAYED> nwgaps_list {{}};
  std::array<char *, MAXDELAYED> nwcigar_list {{}};

  int target_count = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++)
    {
      struct hit * hit = searchinfo->hits + x;
      if (not hit->rejected)
        {
          target_list[target_count++] = hit->target;
        }
    }

  if (target_count != 0)
    {
      search16(searchinfo->s,
               target_count,
               target_list.data(),
               nwscore_list.data(),
               nwalignmentlength_list.data(),
               nwmatches_list.data(),
               nwmismatches_list.data(),
               nwgaps_list.data(),
               nwcigar_list.data());
    }

  int i = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++)
    {
      /* maxrejects or maxaccepts reached - ignore remaining hits */
      if ((searchinfo->rejects < opt_maxrejects) && (searchinfo->accepts < opt_maxaccepts))
        {
          struct hit * hit = searchinfo->hits + x;

          if (hit->rejected)
            {
              searchinfo->rejects++;
            }
          else
            {
              int64_t const target = hit->target;
              int64_t nwscore = nwscore_list[i];

              char * nwcigar = nullptr;
              int64_t nwalignmentlength = 0;
              int64_t nwmatches = 0;
              int64_t nwmismatches = 0;
              int64_t nwgaps = 0;

              int64_t const dseqlen = db_getsequencelen(target);

              if (nwscore == std::numeric_limits<short>::max())
                {
                  /* In case the SIMD aligner cannot align,
                     perform a new alignment with the
                     linear memory aligner */

                  char * dseq = db_getsequence(target);

                  if (nwcigar_list[i] != nullptr)
                    {
                      xfree(nwcigar_list[i]);
                    }

                  nwcigar = xstrdup(searchinfo->lma->align(searchinfo->qsequence,
                                                   dseq,
                                                   searchinfo->qseqlen,
                                                   dseqlen));

                  searchinfo->lma->alignstats(nwcigar,
                                      searchinfo->qsequence,
                                      dseq,
                                      & nwscore,
                                      & nwalignmentlength,
                                      & nwmatches,
                                      & nwmismatches,
                                      & nwgaps);
                }
              else
                {
                  nwalignmentlength = nwalignmentlength_list[i];
                  nwmatches = nwmatches_list[i];
                  nwmismatches = nwmismatches_list[i];
                  nwgaps = nwgaps_list[i];
                  nwcigar = nwcigar_list[i];
                }

              hit->aligned = true;
              hit->shortest = std::min(searchinfo->qseqlen, static_cast<int>(dseqlen));
              hit->longest = std::max(searchinfo->qseqlen, static_cast<int>(dseqlen));
              hit->nwalignment = nwcigar;
              hit->nwscore = nwscore;
              hit->nwdiff = nwalignmentlength - nwmatches;
              hit->nwgaps = nwgaps;
              hit->nwindels = nwalignmentlength - nwmatches - nwmismatches;
              hit->nwalignmentlength = nwalignmentlength;
              hit->nwid = 100.0 * (nwalignmentlength - hit->nwdiff) /
                nwalignmentlength;
              hit->matches = nwalignmentlength - hit->nwdiff;
              hit->mismatches = hit->nwdiff - hit->nwindels;

              /* trim alignment and compute numbers excluding terminal gaps */
              align_trim(hit);

              /* test accept/reject criteria after alignment */
              if (search_acceptable_aligned(searchinfo, hit))
                {
                  searchinfo->accepts++;
                }
              else
                {
                  searchinfo->rejects++;
                }

              ++i;
            }
        }
    }

  /* free ignored alignments */
  while (i < target_count)
    {
      xfree(nwcigar_list[i++]);
    }

  searchinfo->finalized = searchinfo->hit_count;
}


auto search_onequery(struct searchinfo_s * searchinfo, int seqmask) -> void
{
  searchinfo->hit_count = 0;

  search16_qprep(searchinfo->s, searchinfo->qsequence, searchinfo->qseqlen);

  searchinfo->lma = new LinearMemoryAligner;

  int64_t * scorematrix = searchinfo->lma->scorematrix_create(opt_match, opt_mismatch);

  searchinfo->lma->set_parameters(scorematrix,
                          opt_gap_open_query_left,
                          opt_gap_open_target_left,
                          opt_gap_open_query_interior,
                          opt_gap_open_target_interior,
                          opt_gap_open_query_right,
                          opt_gap_open_target_right,
                          opt_gap_extension_query_left,
                          opt_gap_extension_target_left,
                          opt_gap_extension_query_interior,
                          opt_gap_extension_target_interior,
                          opt_gap_extension_query_right,
                          opt_gap_extension_target_right);

  /* extract unique kmer samples from query*/
  unique_count(searchinfo->uh, opt_wordlength,
               searchinfo->qseqlen, searchinfo->qsequence,
               &searchinfo->kmersamplecount, &searchinfo->kmersample, seqmask);

  /* find database sequences with the most kmer hits */
  search_topscores(searchinfo);

  /* analyse targets with the highest number of kmer hits */
  searchinfo->accepts = 0;
  searchinfo->rejects = 0;
  searchinfo->finalized = 0;

  int delayed = 0;

  while ((searchinfo->finalized + delayed < opt_maxaccepts + opt_maxrejects - 1) &&
         (searchinfo->rejects < opt_maxrejects) &&
         (searchinfo->accepts < opt_maxaccepts) &&
         (not minheap_isempty(searchinfo->m)))
    {
      elem_t const e = minheap_poplast(searchinfo->m);

      struct hit * hit = searchinfo->hits + searchinfo->hit_count;

      hit->target = e.seqno;
      hit->count = e.count;
      hit->strand = searchinfo->strand;
      hit->rejected = false;
      hit->accepted = false;
      hit->aligned = false;
      hit->weak = false;
      hit->nwalignment = nullptr;

      /* Test some accept/reject criteria before alignment */
      if (search_acceptable_unaligned(searchinfo, e.seqno))
        {
          ++delayed;
        }
      else
        {
          hit->rejected = true;
        }

      searchinfo->hit_count++;

      if (delayed == MAXDELAYED)
        {
          align_delayed(searchinfo);
          delayed = 0;
        }
    }
  if (delayed > 0)
    {
      align_delayed(searchinfo);
    }

  delete searchinfo->lma;
  xfree(scorematrix);
}


auto search_findbest2_byid(struct searchinfo_s * si_p,
                           struct searchinfo_s * si_m) -> struct hit *
{
  struct hit * best = nullptr;

  for (int i = 0; i < si_p->hit_count; i++)
    {
      if ((best == nullptr) or (hit_compare_byid_typed(si_p->hits + i, best) < 0))
        {
          best = si_p->hits + i;
        }
    }

  if (opt_strand>1)
    {
      for (int i = 0; i < si_m->hit_count; i++)
        {
          if ((best == nullptr) or (hit_compare_byid_typed(si_m->hits + i, best) < 0))
            {
              best = si_m->hits + i;
            }
        }
    }

  if ((best != nullptr) and not best->accepted)
    {
      best = nullptr;
    }

  return best;
}


auto search_findbest2_bysize(struct searchinfo_s * si_p,
                             struct searchinfo_s * si_m) -> struct hit *
{
  struct hit * best = nullptr;

  for (int i = 0; i < si_p->hit_count; i++)
    {
      if ((best == nullptr) or (hit_compare_bysize_typed(si_p->hits + i, best) < 0))
        {
          best = si_p->hits + i;
        }
    }

  if (opt_strand>1)
    {
      for (int i = 0; i < si_m->hit_count; i++)
        {
          if ((best == nullptr) or (hit_compare_bysize_typed(si_m->hits + i, best) < 0))
            {
              best = si_m->hits + i;
            }
        }
    }

  if ((best != nullptr) and not best->accepted)
    {
      best = nullptr;
    }

  return best;
}


auto search_joinhits(struct searchinfo_s * si_p,
                     struct searchinfo_s * si_m,
                     struct hit * * hitsp,
                     int * hit_count) -> void
{
  /* join and sort accepted and weak hits from both strands */
  /* free the remaining alignments */

  /* first, just count the number of hits to keep */
  int a = 0;
  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = (s != 0) ? si_m : si_p;
      for (int i = 0; i<si->hit_count; i++)
        {
          struct hit * h = si->hits + i;
          if (h->accepted || h->weak)
            {
              ++a;
            }
        }
    }

  /* allocate new array of hits */
  auto * hits = (struct hit *) xmalloc(a * sizeof(struct hit));

  /* copy over the hits to be kept */
  a = 0;
  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = (s != 0) ? si_m : si_p;
      for (int i = 0; i<si->hit_count; i++)
        {
          struct hit * h = si->hits + i;
          if (h->accepted || h->weak)
            {
              hits[a++] = *h;
            }
          else if (h->aligned)
            {
              xfree(h->nwalignment);
            }
        }
    }

  /* last, sort the hits */
  qsort(hits, a, sizeof(struct hit), hit_compare_byid);

  *hitsp = hits;
  *hit_count = a;
}
