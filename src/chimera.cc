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
#include "attributes.h"
#include "chimera.h"
#include "dbindex.h"
#include "maps.h"
#include "mask.h"
#include "minheap.h"
#include "udb.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/xpthread.hpp"
#include <algorithm>  // std::fill, std::fill_n, std::max, std::max_element, std::min, std::transform
#include <array>
#include <cassert>
#include <cctype>  // std::tolower
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdlib>  // std::qsort
#include <cstdio>  // std::FILE, std::fprintf, std::sscanf
#include <cstring>  // std::strlen, std::strncpy, std::strcpy
#include <iterator>  // std::next
#include <limits>
#include <numeric>  // std::accumulate
#include <pthread.h>
#include <vector>


/*
  This code implements the method described in this paper:

  Robert C. Edgar, Brian J. Haas, Jose C. Clemente, Christopher Quince
  and Rob Knight (2011)
  UCHIME improves sensitivity and speed of chimera detection
  Bioinformatics, 27, 16, 2194-2200
  https://doi.org/10.1093/bioinformatics/btr381
*/

/* global constants/data, no need for synchronization */
static auto parts = 0;
constexpr auto maxparts = 100;
constexpr auto window = 32;
constexpr auto few = 4;
constexpr auto maxcandidates = few * maxparts;
constexpr auto rejects = 16;
constexpr auto chimera_id = 0.55;
static int tophits;
static pthread_attr_t attr;
static pthread_t * pthread;
static fastx_handle query_fasta_h;

/* mutexes and global data protected by mutex */
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static unsigned int seqno = 0;
static uint64_t progress = 0;
static int chimera_count = 0;
static int nonchimera_count = 0;
static int borderline_count = 0;
static int total_count = 0;
static int64_t chimera_abundance = 0;
static int64_t nonchimera_abundance = 0;
static int64_t borderline_abundance = 0;
static int64_t total_abundance = 0;
static std::FILE * fp_chimeras = nullptr;
static std::FILE * fp_nonchimeras = nullptr;
static std::FILE * fp_uchimealns = nullptr;
static std::FILE * fp_uchimeout = nullptr;
static std::FILE * fp_borderline = nullptr;

/* information for each query sequence to be checked */
struct chimera_info_s
{
  int query_alloc = 0; /* the longest query sequence allocated memory for */
  int head_alloc = 0; /* the longest header allocated memory for */

  int query_no = 0;
  std::vector<char> query_head;
  int query_head_len = 0;
  int query_size = 0;
  std::vector<char> query_seq;
  int query_len = 0;

  std::array<struct searchinfo_s, maxparts> si {{}};

  std::array<unsigned int, maxcandidates> cand_list {{}};
  int cand_count = 0;

  struct s16info_s * s = nullptr;
  std::array<CELL, maxcandidates> snwscore {{}};
  std::array<unsigned short, maxcandidates> snwalignmentlength {{}};
  std::array<unsigned short, maxcandidates> snwmatches {{}};
  std::array<unsigned short, maxcandidates> snwmismatches {{}};
  std::array<unsigned short, maxcandidates> snwgaps {{}};
  std::array<int64_t, maxcandidates> nwscore {{}};
  std::array<int64_t, maxcandidates> nwalignmentlength {{}};
  std::array<int64_t, maxcandidates> nwmatches {{}};
  std::array<int64_t, maxcandidates> nwmismatches {{}};
  std::array<int64_t, maxcandidates> nwgaps {{}};
  std::array<char *, maxcandidates> nwcigar {{}};

  int match_size = 0;
  std::vector<int> match;
  std::vector<int> insert;
  std::vector<int> smooth;
  std::vector<int> maxsmooth;

  std::vector<double> scan_p;
  std::vector<double> scan_q;

  int parents_found = 0;
  std::array<int, maxparents> best_parents {{}};
  std::array<int, maxparents> best_start {{}};
  std::array<int, maxparents> best_len {{}};

  int best_target = 0;
  char * best_cigar = nullptr;

  std::vector<int> maxi;  // longest insertion per position
  std::vector<std::vector<char>> paln;
  std::vector<char> qaln;
  std::vector<char> diffs;
  std::vector<char> votes;
  std::vector<char> model;
  std::vector<char> ignore;

  struct hit * all_hits = nullptr;  // unused?
  double best_h = 0;
};


static struct chimera_info_s * cia;


enum struct Status {
  no_parents = 0,   // non-chimeric
  no_alignment = 1, // score < 0, non-chimeric
  low_score = 2,    // score < minh, non-chimeric
  suspicious = 3,   // score >= minh, not available with uchime2_denovo and uchime3_denovo
  chimeric = 4      // score >= minh && divdiff >= opt_mindiv && ...
};


auto realloc_arrays(struct chimera_info_s * ci) -> void
{
  if (opt_chimeras_denovo != nullptr)
    {
      if (opt_chimeras_parts == 0) {
        parts = (ci->query_len + maxparts - 1) / maxparts;  // bug fix: std::max(expr, parts);?
      }
      else {
        parts = opt_chimeras_parts;
      }
      if (parts < 2) {
        parts = 2;
      }
      else if (parts > maxparts) {
        parts = maxparts;
      }
    }
  else
    {
      /* default for uchime, uchime2, and uchime3 */
      parts = 4;
    }

  const int maxhlen = std::max(ci->query_head_len, 1);
  if (maxhlen > ci->head_alloc)
    {
      ci->head_alloc = maxhlen;
      ci->query_head.resize(maxhlen + 1);
    }

  /* realloc arrays based on query length */

  const int maxqlen = std::max(ci->query_len, 1);
  const int maxpartlen = (maxqlen + parts - 1) / parts;

  if (maxqlen > ci->query_alloc)
    {
      ci->query_alloc = maxqlen;

      ci->query_seq.resize(maxqlen + 1);

      for (auto & query_info: ci->si)
        {
          query_info.qsequence_v.resize(maxpartlen + 1);
          query_info.qsequence = query_info.qsequence_v.data();
        }

      ci->maxi.resize(maxqlen + 1);
      ci->maxsmooth.resize(maxqlen);
      ci->match.resize(maxcandidates * maxqlen);
      ci->insert.resize(maxcandidates * maxqlen);
      ci->smooth.resize(maxcandidates * maxqlen);

      ci->scan_p.resize(maxqlen + 1);
      ci->scan_q.resize(maxqlen + 1);

      const int maxalnlen = maxqlen + (2 * db_getlongestsequence());
      ci->paln.resize(maxparents);
      for (auto & a_parent_alignment : ci->paln) {
        a_parent_alignment.resize(maxalnlen + 1);
      }
      ci->qaln.resize(maxalnlen + 1);
      ci->diffs.resize(maxalnlen + 1);
      ci->votes.resize(maxalnlen + 1);
      ci->model.resize(maxalnlen + 1);
      ci->ignore.resize(maxalnlen + 1);
    }
}


auto reset_matches(struct chimera_info_s * a_chimera_info) -> void {
  // refactoring: initialization to zero? (useless), or reset to zero??
  std::fill(a_chimera_info->match.begin(), a_chimera_info->match.end(), 0);
  std::fill(a_chimera_info->insert.begin(), a_chimera_info->insert.end(), 0);
}


auto find_matches(struct chimera_info_s * ci) -> void
{
  /* find the positions with matches for each potential parent */
  /* also note the positions with inserts in front */

  char * qseq = ci->query_seq.data();

  for (int i = 0; i < ci->cand_count; ++i)
    {
      char * tseq = db_getsequence(ci->cand_list[i]);

      int qpos = 0;
      int tpos = 0;

      char * p = ci->nwcigar[i];
      char * e = p + std::strlen(p);

      while (p < e)
        {
          int run = 1;
          int scanlength = 0;
          std::sscanf(p, "%d%n", &run, &scanlength);
          p += scanlength;
          char const op = *p;
          ++p;
          switch (op)
            {
            case 'M':
              for (int k = 0; k < run; ++k)
                {
                  if ((map_4bit(qseq[qpos]) &
                       map_4bit(tseq[tpos])) != 0U)
                    {
                      ci->match[(i * ci->query_len) + qpos] = 1;
                    }
                  ++qpos;
                  ++tpos;
                }
              break;

            case 'I':
              ci->insert[(i * ci->query_len) + qpos] = run;
              tpos += run;
              break;

            case 'D':
              qpos += run;
              break;
            }
        }
    }
}


struct parents_info_s
{
  int cand = 0;
  int start = 0;
  int len = 0;
};


auto compare_positions(const void * a, const void * b) -> int
{
  const int lhs = ((const parents_info_s *) a)->start;
  const int rhs = ((const parents_info_s *) b)->start;

  if (lhs < rhs) {
    return -1;
  }
  if (lhs > rhs) {
    return +1;
  }
  return 0;
}


auto scan_matches(struct chimera_info_s * ci,
                  int * matches,
                  int len,
                  double const percentage,
                  int * best_start,
                  int * best_len) -> bool
{
  /*
    Scan matches array of zeros and ones, and find the longest subsequence
    having a match fraction above or equal to the given percentage (e.g. 2%).
    Based on an idea of finding the longest positive sum substring:
    https://stackoverflow.com/questions/28356453/longest-positive-sum-substring
    If the percentage is 2%, matches are given a score of 2 and mismatches -98.
  */

  auto const score_match = percentage;
  auto const score_mismatch = percentage - 100.0;

  auto & p = ci->scan_p;
  auto & q = ci->scan_q;

  p[0] = 0.0;
  for (auto i = 0; i < len; ++i) {
    p[i + 1] = p[i] + ((matches[i] != 0) ? score_match : score_mismatch);
  }

  q[len] = p[len];
  for (auto i = len - 1; i >= 0; --i) {
    q[i] = std::max(q[i + 1], p[i]);
  }

  auto best_i = 0;
  auto best_d = -1;
  auto best_c = -1.0;
  auto i = 1;
  auto j = 1;
  while (j <= len)
    {
      auto const c = q[j] - p[i - 1];
      if (c >= 0.0)
        {
          auto const d = j - i + 1;
          if (d > best_d)
            {
              best_i = i;
              best_d = d;
              best_c = c;
            }
          j += 1;
        }
      else
        {
          i += 1;
        }
    }

  if (best_c >= 0.0)
    {
      *best_start = best_i - 1;
      *best_len = best_d;
      return true;
    }
  return false;
}


auto find_best_parents_long(struct chimera_info_s * ci) -> int
{
  /* Find parents with longest matching regions, without indels, allowing
     a given percentage of mismatches (specified with --chimeras_diff_pct),
     and excluding regions matched by previously identified parents. */

  reset_matches(ci);
  find_matches(ci);

  std::array<struct parents_info_s, maxparents> best_parents {{}};

  for (int f = 0; f < maxparents; ++f)
    {
      best_parents[f].cand = -1;
      best_parents[f].start = -1;
    }

  std::vector<bool> position_used(ci->query_len, false);

  int pos_remaining = ci->query_len;
  int parents_found = 0;

  for (int f = 0; f < opt_chimeras_parents_max; ++f)
    {
      /* scan each candidate and find longest matching region */

      int best_start = 0;
      int best_len = 0;
      int best_cand = -1;

      for (int i = 0; i < ci->cand_count; ++i)
        {
          int start = 0;
          int len = 0;
          int j = 0;
          while (j < ci->query_len)
            {
              start = j;
              len = 0;
              while ((j < ci->query_len) &&
                     (not position_used[j]) &&
                     ((len == 0) or (ci->insert[(i * ci->query_len) + j] == 0)))
                {
                  ++len;
                  ++j;
                }
              if (len > best_len)
                {
                  int scan_best_start = 0;
                  int scan_best_len = 0;
                  if (scan_matches(ci,
                                   &ci->match[(i * ci->query_len) + start],
                                   len,
                                   opt_chimeras_diff_pct,
                                   & scan_best_start,
                                   & scan_best_len))
                    {
                      if (scan_best_len > best_len)
                        {
                          best_cand = i;
                          best_start = start + scan_best_start;
                          best_len = scan_best_len;
                        }
                    }
                }
              ++j;
            }
        }

      if (best_len >= opt_chimeras_length_min)
        {
          best_parents[f].cand = best_cand;
          best_parents[f].start = best_start;
          best_parents[f].len = best_len;
          ++parents_found;

#if 0
          if (f == 0)
            printf("\n");
          printf("Best parents long: %d %d %d %d %s %s\n",
                 f,
                 best_cand,
                 best_start,
                 best_len,
                 ci->query_head.data(),
                 db_getheader(ci->cand_list[best_cand]));
#endif

          /* mark positions used */
          for (int j = best_start; j < best_start + best_len; ++j)
            {
              position_used[j] = true;
            }
          pos_remaining -= best_len;
        }
      else {
        break;
      }
    }

  /* sort parents by position */
  std::qsort(best_parents.data(),
             parents_found,
             sizeof(struct parents_info_s),
             compare_positions);

  ci->parents_found = parents_found;

  for (int f = 0; f < parents_found; ++f)
    {
      ci->best_parents[f] = best_parents[f].cand;
      ci->best_start[f] = best_parents[f].start;
      ci->best_len[f] = best_parents[f].len;
    }

#if 0
  if (pos_remaining == 0)
    printf("Fully covered!\n");
  else
    printf("Not covered completely (%d).\n", pos_remaining);
#endif

  return static_cast<int>((parents_found > 1) and (pos_remaining == 0));
}


auto find_best_parents(struct chimera_info_s * ci) -> int
{
  reset_matches(ci);
  find_matches(ci);

  std::array<int, maxparents> best_parent_cand {{}};

  for (int f = 0; f < 2; ++f)
    {
      best_parent_cand[f] = -1;
      ci->best_parents[f] = -1;
    }

  std::vector<bool> cand_selected(ci->cand_count, false);

  for (int f = 0; f < 2; ++f)
    {
      if (f > 0)
        {
          /* for all parents except the first */

          /* wipe out matches for all candidates in positions
             covered by the previous parent */

          for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
            {
              int const z = (best_parent_cand[f - 1] * ci->query_len) + qpos;
              if (ci->smooth[z] == ci->maxsmooth[qpos])
                {
                  for (int i = qpos + 1 - window; i <= qpos; ++i)
                    {
                      for (int j = 0; j < ci->cand_count; ++j)
                        {
                          ci->match[(j * ci->query_len) + i] = 0;
                        }
                    }
                }
            }
        }


      /* Compute smoothed score in a 32bp window for each candidate. */
      /* Record max smoothed score for each position among candidates left. */

      // refactoring: reset or initialize?
      std::fill(ci->maxsmooth.begin(), ci->maxsmooth.end(), 0);

      for (int i = 0; i < ci->cand_count; ++i)
        {
          if (not cand_selected[i])
            {
              int sum = 0;
              for (int qpos = 0; qpos < ci->query_len; ++qpos)
                {
                  int const z = (i * ci->query_len) + qpos;
                  sum += ci->match[z];
                  if (qpos >= window)
                    {
                      sum -= ci->match[z - window];
                    }
                  if (qpos >= window - 1)
                    {
                      ci->smooth[z] = sum;
                      ci->maxsmooth[qpos] = std::max(ci->smooth[z], ci->maxsmooth[qpos]);
                    }
                }
            }
        }


      /* find parent with the most wins */

      std::vector<int> wins(ci->cand_count, 0);

      for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
        {
          if (ci->maxsmooth[qpos] != 0)
            {
              for (int i = 0; i < ci->cand_count; ++i)
                {
                  if (not cand_selected[i])
                    {
                      int const z = (i * ci->query_len) + qpos;
                      if (ci->smooth[z] == ci->maxsmooth[qpos])
                        {
                          ++wins[i];
                        }
                    }
                }
            }
        }

      /* select best parent based on most wins */

      int maxwins = 0;
      for (int i = 0; i < ci->cand_count; ++i)
        {
          int const w = wins[i];
          if (w > maxwins)
            {
              maxwins = w;
              best_parent_cand[f] = i;
            }
        }

      /* terminate loop if no parent found */

      if (best_parent_cand[f] < 0) {
        break;
      }

#if 0
      printf("Query %d: Best parent (%d) candidate: %d. Wins: %d\n",
             ci->query_no, f, best_parent_cand[f], maxwins);
#endif

      ci->best_parents[f] = best_parent_cand[f];
      cand_selected[best_parent_cand[f]] = true;
    }

  /* Check if at least 2 candidates selected */

  return static_cast<int>((best_parent_cand[0] >= 0) and (best_parent_cand[1] >= 0));
}


auto find_total_alignment_length(struct chimera_info_s const * chimera_info) -> int {
  // query_len, plus the sum of the longest insertion runs (I) for each position
  return std::accumulate(chimera_info->maxi.begin(),
                         chimera_info->maxi.end(),
                         chimera_info->query_len);
}


auto fill_max_alignment_length(struct chimera_info_s * ci) -> void
{
  /* find max insertions in front of each position in the query sequence */

  std::fill(ci->maxi.begin(), ci->maxi.end(), 0);

  for (int f = 0; f < ci->parents_found; ++f)
    {
      int const best_parent = ci->best_parents[f];
      char * p = ci->nwcigar[best_parent];
      char * e = p + std::strlen(p);
      int pos = 0;
      while (p < e)
        {
          int run = 1;
          int scanlength = 0;
          std::sscanf(p, "%d%n", &run, &scanlength);
          p += scanlength;
          char const op = *p;
          ++p;
          switch (op)
            {
            case 'M':
            case 'D':
              pos += run;
              break;

            case 'I':
              ci->maxi[pos] = std::max(run, ci->maxi[pos]);
              break;
            }
        }
    }
}


auto fill_alignment_parents(struct chimera_info_s * ci) -> void
{
  /* fill in alignment strings for the parents */

  for (int j = 0; j < ci->parents_found; ++j)
    {
      auto & alignment = ci->paln[j];
      int const cand = ci->best_parents[j];
      int const target_seqno = ci->cand_list[cand];
      char * target_seq = db_getsequence(target_seqno);

      auto is_inserted = false;
      int qpos = 0;
      int tpos = 0;
      int alnpos = 0;

      char * p = ci->nwcigar[cand];
      char * e = p + std::strlen(p);

      while (p < e)
        {
          int run = 1;
          int scanlength = 0;
          std::sscanf(p, "%d%n", &run, &scanlength);
          p += scanlength;
          char const op = *p;
          ++p;

          if (op == 'I')
            {
              for (int x = 0; x < ci->maxi[qpos]; ++x)
                {
                  if (x < run)
                    {
                      alignment[alnpos] = map_uppercase(target_seq[tpos]);
                      ++tpos;
                      ++alnpos;
                    }
                  else
                    {
                      alignment[alnpos] = '-';
                      ++alnpos;
                    }
                }
              is_inserted = true;
            }
          else
            {
              for (int x = 0; x < run; ++x)
                {
                  if (not is_inserted)
                    {
                      std::fill_n(&alignment[alnpos], ci->maxi[qpos], '-');
                      alnpos += ci->maxi[qpos];
                    }

                  if (op == 'M')
                    {
                      alignment[alnpos] = map_uppercase(target_seq[tpos]);
                      ++tpos;
                      ++alnpos;
                    }
                  else
                    {
                      alignment[alnpos] = '-';
                      ++alnpos;
                    }

                  ++qpos;
                  is_inserted = false;
                }
            }
        }

      /* add any gaps at the end */

      if (not is_inserted)
        {
          std::fill_n(&alignment[alnpos], ci->maxi[qpos], '-');
          alnpos += ci->maxi[qpos];
        }

      /* end of sequence string */
      alignment[alnpos] = '\0';
    }
}


auto fill_in_alignment_string_for_query(struct chimera_info_s * chimera_info) -> void {
  auto alnpos = 0;
  auto qpos = 0;
  for (auto const nucleotide: chimera_info->query_seq) {
    // add insertion (if any):
    auto const insert_length = chimera_info->maxi[qpos];
    std::fill_n(&chimera_info->qaln[alnpos], insert_length, '-');
    alnpos += insert_length;

    // add (mis-)matching position:
    chimera_info->qaln[alnpos] = map_uppercase(nucleotide);
    ++alnpos;
    ++qpos;
  }
  // add terminal gap (if any):
  auto const insert_length = chimera_info->maxi[chimera_info->query_len];
  std::fill_n(&chimera_info->qaln[alnpos], insert_length, '-');
  alnpos += insert_length;
  chimera_info->qaln[alnpos] = '\0';
}


auto fill_in_model_string_for_query(struct chimera_info_s * chimera_info) -> void {
  int nth_parent = 0;
  auto alnpos = 0;
  for (int qpos = 0; qpos < chimera_info->query_len; ++qpos)
    {
      if (qpos >= (chimera_info->best_start[nth_parent] + chimera_info->best_len[nth_parent])) {
        ++nth_parent;
      }
      // add insertion (if any):
      auto const insert_length = chimera_info->maxi[qpos];
      std::fill_n(&chimera_info->model[alnpos], insert_length, 'A' + nth_parent);
      alnpos += insert_length;

      // add (mis-)matching position:
      chimera_info->model[alnpos] = 'A' + nth_parent;
      ++alnpos;
    }
  // add terminal gap (if any):
  auto const insert_length = chimera_info->maxi[chimera_info->query_len];
  std::fill_n(&chimera_info->model[alnpos], insert_length, 'A' + nth_parent);
  alnpos += insert_length;
  chimera_info->model[alnpos] = '\0';
}


auto count_matches_with_parents(struct chimera_info_s const * chimera_info,
                                int const alignment_length) -> std::array<int, maxparents> {
  std::array<int, maxparents> matches {{}};

  for (auto i = 0; i < alignment_length; ++i)
    {
      auto const qsym = map_4bit(chimera_info->qaln[i]);

      for (auto f = 0; f < chimera_info->parents_found; ++f)
        {
          auto const psym = map_4bit(chimera_info->paln[f][i]);
          if (qsym == psym) {
            ++matches[f];
          }
        }
    }
  return matches;
}


auto compute_global_similarities_with_parents(
    std::array<int, maxparents> const & match_counts,
    int const alignment_length) -> std::array<double, maxparents> {
  std::array<double, maxparents> similarities {{}};
  auto compute_percentage = [alignment_length](int const match_count) -> double {
    return 100.0 * match_count / alignment_length;
  };
  std::transform(match_counts.begin(), match_counts.end(),
                 similarities.begin(), compute_percentage);
  return similarities;
}


auto compute_diffs(struct chimera_info_s const * ci,
                   std::vector<char> const & psym,
                   char const qsym) -> char {
  auto const all_defined = (qsym != '\0') and
    std::all_of(psym.begin(),
                psym.end(),
                [](char const symbol) { return symbol != '\0'; });

  char diff = ' ';

  if (not all_defined) { return diff; }

  auto z = 0;
  for (auto f = 0; f < ci->parents_found; ++f) {
    if (psym[f] == qsym) {
      diff = 'A' + f;
      ++z;
    }
  }
  if (z > 1) {
    diff = ' ';
  }
  return diff;
}


auto eval_parents_long(struct chimera_info_s * ci) -> Status
{
  /* always chimeric if called */
  auto const status = Status::chimeric;

  fill_max_alignment_length(ci);
  auto const alnlen = find_total_alignment_length(ci);

  fill_alignment_parents(ci);

  fill_in_alignment_string_for_query(ci);
  fill_in_model_string_for_query(ci);

  std::vector<char> psym;
  psym.reserve(maxparents);

  for (int i = 0; i < alnlen; ++i)
    {
      auto const qsym = map_4bit(ci->qaln[i]);
      for (int f = 0; f < ci->parents_found; ++f) {
        psym.emplace_back(map_4bit(ci->paln[f][i]));
      }

      /* lower case parent symbols that differ from query */

      for (int f = 0; f < ci->parents_found; ++f) {
        if ((psym[f] != '\0') and (psym[f] != qsym)) {
          ci->paln[f][i] = std::tolower(ci->paln[f][i]);
        }
      }

      /* compute diffs */
      ci->diffs[i] = compute_diffs(ci, psym, qsym);
      psym.clear();
    }

  ci->diffs[alnlen] = '\0';


  auto const match_QP = count_matches_with_parents(ci, alnlen);

  int const seqno_a = ci->cand_list[ci->best_parents[0]];
  int const seqno_b = ci->cand_list[ci->best_parents[1]];
  int const seqno_c = ci->parents_found > 2 ? ci->cand_list[ci->best_parents[2]] : -1;

  auto const QP = compute_global_similarities_with_parents(match_QP, alnlen);
  auto const QT = *std::max_element(QP.begin(), QP.end());

  double const QA = QP[0];
  double const QB = QP[1];
  double const QC = ci->parents_found > 2 ? QP[2] : 0.00;
  double const QM = 100.00;
  double const divfrac = 100.00 * (QM - QT) / QT;  // divergence of the model with the closest parent

  xpthread_mutex_lock(&mutex_output);

  if ((opt_alnout != nullptr) and (status == Status::chimeric))
    {
      std::fprintf(fp_uchimealns, "\n");
      std::fprintf(fp_uchimealns, "----------------------------------------"
                   "--------------------------------\n");
      std::fprintf(fp_uchimealns, "Query   (%5d nt) ",
                   ci->query_len);
      header_fprint_strip(fp_uchimealns,
                          ci->query_head.data(),
                          ci->query_head_len,
                          opt_xsize,
                          opt_xee,
                          opt_xlength);

      assert(ci->parents_found <= 20);  // 20 parents max ('A' to 'U')
      for (int f = 0; f < ci->parents_found; ++f)
        {
          int const seqno = ci->cand_list[ci->best_parents[f]];
          std::fprintf(fp_uchimealns, "\nParent%c (%5" PRIu64 " nt) ",
                       'A' + f,
                       db_getsequencelen(seqno));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(seqno),
                              db_getheaderlen(seqno),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
        }

      std::fprintf(fp_uchimealns, "\n\n");


      int const width = opt_alignwidth > 0 ? opt_alignwidth : alnlen;
      int qpos = 0;
      std::array<int, maxparents> ppos {{}};
      int rest = alnlen;

      for (int i = 0; i < alnlen; i += width)
        {
          /* count non-gap symbols on current line */

          int qnt = 0;
          std::array<int, maxparents> pnt {{}};

          int const w = std::min(rest, width);

          for (int j = 0; j < w; ++j)
            {
              if (ci->qaln[i + j] != '-')
                {
                  ++qnt;
                }

              for (int f = 0; f < ci->parents_found; ++f) {
                if (ci->paln[f][i + j] != '-')
                  {
                    ++pnt[f];
                  }
              }
            }

          fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                  qpos + 1, w, &ci->qaln[i], qpos + qnt);

          for (int f = 0; f < ci->parents_found; ++f)
            {
              fprintf(fp_uchimealns, "%c %5d %.*s %d\n",
                      'A' + f,
                      ppos[f] + 1, w, &ci->paln[f][i], ppos[f] + pnt[f]);
            }

          fprintf(fp_uchimealns, "Diffs   %.*s\n", w, &ci->diffs[i]);
          fprintf(fp_uchimealns, "Model   %.*s\n", w, &ci->model[i]);
          fprintf(fp_uchimealns, "\n");

          rest -= width;
          qpos += qnt;
          for (int f = 0; f < ci->parents_found; ++f) {
            ppos[f] += pnt[f];
          }
        }

      fprintf(fp_uchimealns, "Ids.  QA %.2f%%, QB %.2f%%, QC %.2f%%, "
              "QT %.2f%%, QModel %.2f%%, Div. %+.2f%%\n",
              QA, QB, QC, QT, QM, divfrac);
    }

  if (opt_tabbedout != nullptr)
    {
      fprintf(fp_uchimeout, "%.4f\t", 99.9999);

      header_fprint_strip(fp_uchimeout,
                          ci->query_head.data(),
                          ci->query_head_len,
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(fp_uchimeout, "\t");
      header_fprint_strip(fp_uchimeout,
                          db_getheader(seqno_a),
                          db_getheaderlen(seqno_a),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(fp_uchimeout, "\t");
      header_fprint_strip(fp_uchimeout,
                          db_getheader(seqno_b),
                          db_getheaderlen(seqno_b),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(fp_uchimeout, "\t");
      if (seqno_c >= 0)
        {
          header_fprint_strip(fp_uchimeout,
                              db_getheader(seqno_c),
                              db_getheaderlen(seqno_c),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
        }
      else
        {
          fprintf(fp_uchimeout, "*");
        }
      fprintf(fp_uchimeout, "\t");

      fprintf(fp_uchimeout,
              "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t"
              "%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%c\n",
              QM,
              QA,
              QB,
              QC,
              QT,
              0, /* ignore, left yes */
              0, /* ignore, left no */
              0, /* ignore, left abstain */
              0, /* ignore, right yes */
              0, /* ignore, right no */
              0, /* ignore, right abstain */
              0.00,
              status == Status::chimeric ? 'Y' : (status == Status::low_score ? 'N' : '?'));
    }

  xpthread_mutex_unlock(&mutex_output);

  return status;
}


auto eval_parents(struct chimera_info_s * ci) -> Status
{
  auto status = Status::no_alignment;
  ci->parents_found = 2;

  fill_max_alignment_length(ci);
  auto const alnlen = find_total_alignment_length(ci);

  fill_alignment_parents(ci);

  /* fill in alignment string for query */

  char * q = ci->qaln.data();
  int qpos = 0;
  for (int i = 0; i < ci->query_len; ++i)
    {
      for (int j = 0; j < ci->maxi[i]; ++j)
        {
          *q = '-';
          ++q;
        }
      *q = map_uppercase(ci->query_seq[qpos]);
      ++qpos;
      ++q;
    }
  for (int j = 0; j < ci->maxi[ci->query_len]; ++j)
    {
      *q = '-';
      ++q;
    }
  *q = 0;

  /* mark positions to ignore in voting */

  for (int i = 0; i < alnlen; ++i) {
    ci->ignore[i] = 0;
  }

  for (int i = 0; i < alnlen; ++i)
    {
      unsigned int const qsym  = chrmap_4bit[(int) (ci->qaln[i])];
      unsigned int const p1sym = chrmap_4bit[(int) (ci->paln[0][i])];
      unsigned int const p2sym = chrmap_4bit[(int) (ci->paln[1][i])];

      /* ignore gap positions and those next to the gap */
      if ((qsym == 0U) or (p1sym == 0U) or (p2sym == 0U))
        {
          ci->ignore[i] = 1;
          if (i > 0)
            {
              ci->ignore[i - 1] = 1;
            }
          if (i < alnlen - 1)
            {
              ci->ignore[i + 1] = 1;
            }
        }

      /* ignore ambiguous symbols */
      if ((ambiguous_4bit[qsym] != 0U) or
          (ambiguous_4bit[p1sym] != 0U) or
          (ambiguous_4bit[p2sym] != 0U))
        {
          ci->ignore[i] = 1;
        }

      /* lower case parent symbols that differ from query */

      if ((p1sym != 0U) and (p1sym != qsym))
        {
          ci->paln[0][i] = std::tolower(ci->paln[0][i]);
        }

      if ((p2sym != 0U) and (p2sym != qsym))
        {
          ci->paln[1][i] = std::tolower(ci->paln[1][i]);
        }

      /* compute diffs */

      char diff = '\0';

      if ((qsym != 0U) and (p1sym != 0U) and (p2sym != 0U))
        {
          if (p1sym == p2sym)
            {
              if (qsym == p1sym)
                {
                  diff = ' ';
                }
              else
                {
                  diff = 'N';
                }
            }
          else
            {
              if (qsym == p1sym)
                {
                  diff = 'A';
                }
              else if (qsym == p2sym)
                {
                  diff = 'B';
                }
              else
                {
                  diff = '?';
                }
            }
        }
      else
        {
          diff = ' ';
        }

      ci->diffs[i] = diff;
    }

  ci->diffs[alnlen] = '\0';

  /* compute score */

  int sumA = 0;
  int sumB = 0;
  int sumN = 0;

  for (int i = 0; i < alnlen; ++i)
    {
      if (ci->ignore[i] == 0)
        {
          char const diff = ci->diffs[i];

          if (diff == 'A')
            {
              ++sumA;
            }
          else if (diff == 'B')
            {
              ++sumB;
            }
          else if (diff != ' ')
            {
              ++sumN;
            }
        }

    }

  int left_n = 0;
  int left_a = 0;
  int left_y = 0;
  int right_n = sumA;
  int right_a = sumN;
  int right_y = sumB;

  double best_h = -1;
  int best_i = -1;
  auto best_is_reverse = false;

  int best_left_y = 0;
  int best_right_y = 0;
  int best_left_n = 0;
  int best_right_n = 0;
  int best_left_a = 0;
  int best_right_a = 0;

  for (int i = 0; i < alnlen; ++i)
    {
      if (ci->ignore[i] == 0)
        {
          char const diff = ci->diffs[i];
          if (diff != ' ')
            {
              if (diff == 'A')
                {
                  ++left_y;
                  --right_n;
                }
              else if (diff == 'B')
                {
                  ++left_n;
                  --right_y;
                }
              else
                {
                  ++left_a;
                  --right_a;
                }

              double left_h = 0;
              double right_h = 0;
              double h = 0;

              if ((left_y > left_n) and (right_y > right_n))
                {
                  left_h = left_y / (opt_xn * (left_n + opt_dn) + left_a);
                  right_h = right_y / (opt_xn * (right_n + opt_dn) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_is_reverse = false;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_n;
                      best_left_y = left_y;
                      best_left_a = left_a;
                      best_right_n = right_n;
                      best_right_y = right_y;
                      best_right_a = right_a;
                    }
                }
              else if ((left_n > left_y) and (right_n > right_y))
                {
                  /* swap left/right and yes/no */

                  left_h = left_n / (opt_xn * (left_y + opt_dn) + left_a);
                  right_h = right_n / (opt_xn * (right_y + opt_dn) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_is_reverse = true;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_y;
                      best_left_y = left_n;
                      best_left_a = left_a;
                      best_right_n = right_y;
                      best_right_y = right_n;
                      best_right_a = right_a;
                    }
                }
            }
        }
    }

  ci->best_h = best_h > 0 ? best_h : 0.0;

  if (best_h >= 0.0)
    {
      status = Status::low_score;

      /* flip A and B if necessary */

      if (best_is_reverse)
        {
          for (int i = 0; i < alnlen; ++i)
            {
              char const diff = ci->diffs[i];
              if (diff == 'A')
                {
                  ci->diffs[i] = 'B';
                }
              else if (diff == 'B')
                {
                  ci->diffs[i] = 'A';
                }
            }
        }

      /* fill in votes and model */

      for (int i = 0; i < alnlen; ++i)
        {
          char const m = i <= best_i ? 'A' : 'B';
          ci->model[i] = m;

          char v = ' ';
          if (ci->ignore[i] == 0)
            {
              char const d = ci->diffs[i];

              if ((d == 'A') or (d == 'B'))
                {
                  if (d == m)
                    {
                      v = '+';
                    }
                  else
                    {
                      v = '!';
                    }
                }
              else if ((d == 'N') or (d == '?'))
                {
                  v = '0';
                }
            }
          ci->votes[i] = v;

          /* lower case diffs for no votes */
          if (v == '!')
            {
              ci->diffs[i] = std::tolower(ci->diffs[i]);
            }
        }

      /* fill in crossover region */

      for (int i = best_i + 1; i < alnlen; ++i)
        {
          if ((ci->diffs[i] == ' ') or (ci->diffs[i] == 'A'))
            {
              ci->model[i] = 'x';
            }
          else
            {
              break;
            }
        }

      ci->votes[alnlen] = 0;
      ci->model[alnlen] = 0;

      /* count matches */

      int const index_a = best_is_reverse ? 1 : 0;
      int const index_b = best_is_reverse ? 0 : 1;

      int match_QA = 0;
      int match_QB = 0;
      int match_AB = 0;
      int match_QM = 0;
      int cols = 0;

      for (auto i = 0; i < alnlen; i++)
        {
          if (ci->ignore[i] == '\0')
            {
              ++cols;

              auto const qsym = map_4bit(ci->qaln[i]);
              auto const asym = map_4bit(ci->paln[index_a][i]);
              auto const bsym = map_4bit(ci->paln[index_b][i]);
              auto const msym = (i <= best_i) ? asym : bsym;

              if (qsym == asym)
                {
                  ++match_QA;
                }

              if (qsym == bsym)
                {
                  ++match_QB;
                }

              if (asym == bsym)
                {
                  ++match_AB;
                }

              if (qsym == msym)
                {
                  ++match_QM;
                }
            }
        }

      int const seqno_a = ci->cand_list[ci->best_parents[index_a]];
      int const seqno_b = ci->cand_list[ci->best_parents[index_b]];

      double const QA = 100.0 * match_QA / cols;
      double const QB = 100.0 * match_QB / cols;
      double const AB = 100.0 * match_AB / cols;
      double const QT = std::max(QA, QB);
      double const QM = 100.0 * match_QM / cols;
      double const divdiff = QM - QT;
      double const divfrac = 100.0 * divdiff / QT;

      int const sumL = best_left_n + best_left_a + best_left_y;
      int const sumR = best_right_n + best_right_a + best_right_y;

      if ((opt_uchime2_denovo != nullptr) or (opt_uchime3_denovo != nullptr))
        {
          // fix -Wfloat-equal: if match_QM == cols, then QM == 100.0
          if ((match_QM == cols) and (QT < 100.0))
            {
              status = Status::chimeric;
            }
        }
      else
        if (best_h >= opt_minh)
          {
            status = Status::suspicious;
            if ((divdiff >= opt_mindiv) and
                (sumL >= opt_mindiffs) and
                (sumR >= opt_mindiffs))
              {
                status = Status::chimeric;
              }
          }

      /* print alignment */

      xpthread_mutex_lock(&mutex_output);

      if ((opt_uchimealns != nullptr) and (status == Status::chimeric))
        {
          fprintf(fp_uchimealns, "\n");
          fprintf(fp_uchimealns, "----------------------------------------"
                  "--------------------------------\n");
          fprintf(fp_uchimealns, "Query   (%5d nt) ",
                  ci->query_len);

          header_fprint_strip(fp_uchimealns,
                              ci->query_head.data(),
                              ci->query_head_len,
                              opt_xsize,
                              opt_xee,
                              opt_xlength);

          fprintf(fp_uchimealns, "\nParentA (%5" PRIu64 " nt) ",
                  db_getsequencelen(seqno_a));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(seqno_a),
                              db_getheaderlen(seqno_a),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);

          fprintf(fp_uchimealns, "\nParentB (%5" PRIu64 " nt) ",
                  db_getsequencelen(seqno_b));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(seqno_b),
                              db_getheaderlen(seqno_b),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimealns, "\n\n");

          auto const width = opt_alignwidth > 0 ? opt_alignwidth : alnlen;
          qpos = 0;
          auto p1pos = 0;
          auto p2pos = 0;
          auto rest = alnlen;

          for (auto i = 0; i < alnlen; i += width)
            {
              /* count non-gap symbols on current line */

              auto qnt = 0;
              auto p1nt = 0;
              auto p2nt = 0;

              auto const w = std::min(rest, width);

              for (auto j = 0; j < w; ++j)
                {
                  if (ci->qaln[i + j] != '-')
                    {
                      ++qnt;
                    }
                  if (ci->paln[0][i + j] != '-')
                    {
                      ++p1nt;
                    }
                  if (ci->paln[1][i + j] != '-')
                    {
                      ++p2nt;
                    }
                }

              if (not best_is_reverse)
                {
                  fprintf(fp_uchimealns, "A %5d %.*s %d\n",
                          p1pos + 1, w, &ci->paln[0][i], p1pos + p1nt);
                  fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos + 1, w, &ci->qaln[i], qpos + qnt);
                  fprintf(fp_uchimealns, "B %5d %.*s %d\n",
                          p2pos + 1, w, &ci->paln[1][i], p2pos + p2nt);
                }
              else
                {
                  fprintf(fp_uchimealns, "A %5d %.*s %d\n",
                          p2pos + 1, w, &ci->paln[1][i], p2pos + p2nt);
                  fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos + 1, w, &ci->qaln[i], qpos + qnt);
                  fprintf(fp_uchimealns, "B %5d %.*s %d\n",
                          p1pos + 1, w, &ci->paln[0][i], p1pos + p1nt);
                }

              fprintf(fp_uchimealns, "Diffs   %.*s\n", w, &ci->diffs[i]);
              fprintf(fp_uchimealns, "Votes   %.*s\n", w, &ci->votes[i]);
              fprintf(fp_uchimealns, "Model   %.*s\n", w, &ci->model[i]);
              fprintf(fp_uchimealns, "\n");

              qpos += qnt;
              p1pos += p1nt;
              p2pos += p2nt;
              rest -= width;
            }

          fprintf(fp_uchimealns, "Ids.  QA %.1f%%, QB %.1f%%, AB %.1f%%, "
                  "QModel %.1f%%, Div. %+.1f%%\n",
                  QA, QB, AB, QM, divfrac);

          fprintf(fp_uchimealns, "Diffs Left %d: N %d, A %d, Y %d (%.1f%%); "
                  "Right %d: N %d, A %d, Y %d (%.1f%%), Score %.4f\n",
                  sumL, best_left_n, best_left_a, best_left_y,
                  100.0 * best_left_y / sumL,
                  sumR, best_right_n, best_right_a, best_right_y,
                  100.0 * best_right_y / sumR,
                  best_h);
        }

      if (opt_uchimeout != nullptr)
        {
          fprintf(fp_uchimeout, "%.4f\t", best_h);

          header_fprint_strip(fp_uchimeout,
                              ci->query_head.data(),
                              ci->query_head_len,
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimeout, "\t");
          header_fprint_strip(fp_uchimeout,
                              db_getheader(seqno_a),
                              db_getheaderlen(seqno_a),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimeout, "\t");
          header_fprint_strip(fp_uchimeout,
                              db_getheader(seqno_b),
                              db_getheaderlen(seqno_b),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimeout, "\t");

          if (opt_uchimeout5 == 0)
            {
              if (QA >= QB)
                {
                  header_fprint_strip(fp_uchimeout,
                                      db_getheader(seqno_a),
                                      db_getheaderlen(seqno_a),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                }
              else
                {
                  header_fprint_strip(fp_uchimeout,
                                      db_getheader(seqno_b),
                                      db_getheaderlen(seqno_b),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                }
              fprintf(fp_uchimeout, "\t");
            }

          fprintf(fp_uchimeout,
                  "%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t"
                  "%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%c\n",
                  QM,
                  QA,
                  QB,
                  AB,
                  QT,
                  best_left_y,
                  best_left_n,
                  best_left_a,
                  best_right_y,
                  best_right_n,
                  best_right_a,
                  divdiff,
                  status == Status::chimeric ? 'Y' : (status == Status::low_score ? 'N' : '?'));
        }
      xpthread_mutex_unlock(&mutex_output);
    }

  return status;
}


auto query_init(struct searchinfo_s * search_info) -> void
{
  static constexpr auto overflow_padding = 16U;  // 16 * sizeof(short) = 32 bytes
  search_info->hits_v.resize(tophits);
  search_info->hits = search_info->hits_v.data();
  search_info->kmers_v.reserve(db_getsequencecount() + overflow_padding);
  search_info->kmers_v.resize(db_getsequencecount());
  search_info->kmers = search_info->kmers_v.data();
  search_info->hit_count = 0;
  search_info->uh = unique_init();
  search_info->s = search16_init(opt_match,
                                 opt_mismatch,
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
  search_info->m = minheap_init(tophits);
}


auto query_exit(struct searchinfo_s * search_info) -> void
{
  search16_exit(search_info->s);
  unique_exit(search_info->uh);
  minheap_exit(search_info->m);

  search_info->qsequence = nullptr;
  search_info->hits = nullptr;
  search_info->kmers = nullptr;
}


auto partition_query(struct chimera_info_s * chimera_info) -> void
{
  auto rest = chimera_info->query_len;
  auto * cursor = chimera_info->query_seq.data();
  for (auto i = 0; i < parts; ++i)
    {
      auto const length = (rest + (parts - i - 1)) / (parts - i);

      auto & search_info = chimera_info->si[i];

      search_info.query_no = chimera_info->query_no;
      search_info.strand = 0;
      search_info.qsize = chimera_info->query_size;
      search_info.query_head_len = chimera_info->query_head_len;
      search_info.query_head = chimera_info->query_head.data();
      search_info.qseqlen = length;
      assert(static_cast<std::size_t>(length) <= search_info.qsequence_v.size());
      std::strncpy(search_info.qsequence_v.data(), cursor, length);
      search_info.qsequence_v[length] = '\0';

      rest -= length;
      cursor = std::next(cursor, length);
    }
}


auto chimera_thread_init(struct chimera_info_s * ci) -> void
{

  for (int i = 0; i < maxparts; ++i)
    {
      query_init(&ci->si[i]);
    }

  ci->s = search16_init(opt_match,
                        opt_mismatch,
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
}


auto chimera_thread_exit(struct chimera_info_s * ci) -> void
{
  search16_exit(ci->s);

  for (auto & a_search_info : ci->si) {
    query_exit(&a_search_info);
  }
}


auto chimera_thread_core(struct chimera_info_s * ci) -> uint64_t
{
  chimera_thread_init(ci);

  std::vector<struct hit> allhits_list(maxcandidates);

  LinearMemoryAligner lma;

  int64_t * scorematrix = lma.scorematrix_create(opt_match, opt_mismatch);

  lma.set_parameters(scorematrix,
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

  while (true)
    {
      /* get next sequence */

      xpthread_mutex_lock(&mutex_input);

      if (opt_uchime_ref != nullptr)
        {
          if (fasta_next(query_fasta_h, (opt_notrunclabels == 0),
                         chrmap_no_change_vector.data()))
            {
              ci->query_head_len = fasta_get_header_length(query_fasta_h);
              ci->query_len = fasta_get_sequence_length(query_fasta_h);
              ci->query_no = fasta_get_seqno(query_fasta_h);
              ci->query_size = fasta_get_abundance(query_fasta_h);

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              /* copy the data locally (query seq, head) */
              std::strcpy(ci->query_head.data(), fasta_get_header(query_fasta_h));
              std::strcpy(ci->query_seq.data(), fasta_get_sequence(query_fasta_h));
            }
          else
            {
              xpthread_mutex_unlock(&mutex_input);
              break; /* end while loop */
            }
        }
      else
        {
          if (seqno < db_getsequencecount())
            {
              ci->query_no = seqno;
              ci->query_head_len = db_getheaderlen(seqno);
              ci->query_len = db_getsequencelen(seqno);
              ci->query_size = db_getabundance(seqno);

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              std::strcpy(ci->query_head.data(), db_getheader(seqno));
              std::strcpy(ci->query_seq.data(), db_getsequence(seqno));
            }
          else
            {
              xpthread_mutex_unlock(&mutex_input);
              break; /* end while loop */
            }
        }

      xpthread_mutex_unlock(&mutex_input);

      auto status = Status::no_parents;

      /* partition query */
      partition_query(ci);

      /* perform searches and collect candidate parents */
      ci->cand_count = 0;
      auto allhits_count = 0;

      if (ci->query_len >= parts)
        {
          for (auto i = 0; i < parts; ++i)
            {
              struct hit * hits = nullptr;
              auto hit_count = 0;
              search_onequery(&ci->si[i], opt_qmask);
              search_joinhits(&ci->si[i], nullptr, & hits, & hit_count);
              for (auto j = 0; j < hit_count; ++j)
                {
                  if (hits[j].accepted)
                    {
                      allhits_list[allhits_count] = hits[j];
                      ++allhits_count;
                    }
                }
              xfree(hits);
            }
        }

      for (int i = 0; i < allhits_count; ++i)
        {
          unsigned int const target = allhits_list[i].target;

          /* skip duplicates */
          auto k = 0;
          for (k = 0; k < ci->cand_count; ++k)
            {
              if (ci->cand_list[k] == target)
                {
                  break;
                }
            }

          if (k == ci->cand_count)
            {
              ci->cand_list[ci->cand_count] = target;
              ++ci->cand_count;
            }

          /* deallocate cigar */
          if (allhits_list[i].nwalignment != nullptr)
            {
              xfree(allhits_list[i].nwalignment);
              allhits_list[i].nwalignment = nullptr;
            }
        }


      /* align full query to each candidate */

      search16_qprep(ci->s, ci->query_seq.data(), ci->query_len);

      search16(ci->s,
               ci->cand_count,
               ci->cand_list.data(),
               ci->snwscore.data(),
               ci->snwalignmentlength.data(),
               ci->snwmatches.data(),
               ci->snwmismatches.data(),
               ci->snwgaps.data(),
               ci->nwcigar.data());

      for (auto i = 0; i < ci->cand_count; ++i)
        {
          int64_t const target = ci->cand_list[i];
          int64_t nwscore = ci->snwscore[i];
          char * nwcigar = nullptr;
          int64_t nwalignmentlength = 0;
          int64_t nwmatches = 0;
          int64_t nwmismatches = 0;
          int64_t nwgaps = 0;

          if (nwscore == std::numeric_limits<short>::max())
            {
              /* In case the SIMD aligner cannot align,
                 perform a new alignment with the
                 linear memory aligner */

              auto * tseq = db_getsequence(target);
              int64_t const tseqlen = db_getsequencelen(target);

              if (ci->nwcigar[i] != nullptr)
                {
                  xfree(ci->nwcigar[i]);
                }

              nwcigar = xstrdup(lma.align(ci->query_seq.data(),
                                          tseq,
                                          ci->query_len,
                                          tseqlen));
              lma.alignstats(nwcigar,
                             ci->query_seq.data(),
                             tseq,
                             & nwscore,
                             & nwalignmentlength,
                             & nwmatches,
                             & nwmismatches,
                             & nwgaps);

              ci->nwcigar[i] = nwcigar;
              ci->nwscore[i] = nwscore;
              ci->nwalignmentlength[i] = nwalignmentlength;
              ci->nwmatches[i] = nwmatches;
              ci->nwmismatches[i] = nwmismatches;
              ci->nwgaps[i] = nwgaps;
            }
          else
            {
              ci->nwscore[i] = ci->snwscore[i];
              ci->nwalignmentlength[i] = ci->snwalignmentlength[i];
              ci->nwmatches[i] = ci->snwmatches[i];
              ci->nwmismatches[i] = ci->snwmismatches[i];
              ci->nwgaps[i] = ci->snwgaps[i];
            }
        }


      /* find the best pair of parents, then compute score for them */

      if (opt_chimeras_denovo != nullptr)
        {
          /* long high-quality reads */
          if (find_best_parents_long(ci) != 0)
            {
              status = eval_parents_long(ci);
            }
          else
            {
              status = Status::no_parents;
            }
        }
      else
        {
          if (find_best_parents(ci) != 0)
            {
              status = eval_parents(ci);
            }
          else
            {
              status = Status::no_parents;
            }
        }

      /* output results */

      xpthread_mutex_lock(&mutex_output);

      ++total_count;
      total_abundance += ci->query_size;

      if (status == Status::chimeric)
        {
          ++chimera_count;
          chimera_abundance += ci->query_size;

          if (opt_chimeras != nullptr)
            {
              fasta_print_general(fp_chimeras,
                                  nullptr,
                                  ci->query_seq.data(),
                                  ci->query_len,
                                  ci->query_head.data(),
                                  ci->query_head_len,
                                  ci->query_size,
                                  chimera_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( (opt_uchime_ref != nullptr) ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h);

            }
        }

      if (status == Status::suspicious)
        {
          ++borderline_count;
          borderline_abundance += ci->query_size;

          if (opt_borderline != nullptr)
            {
              fasta_print_general(fp_borderline,
                                  nullptr,
                                  ci->query_seq.data(),
                                  ci->query_len,
                                  ci->query_head.data(),
                                  ci->query_head_len,
                                  ci->query_size,
                                  borderline_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( (opt_uchime_ref != nullptr) ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h);

            }
        }

      if (status < Status::suspicious)
        {
          ++nonchimera_count;
          nonchimera_abundance += ci->query_size;

          /* output no parents, no chimeras */
          if ((status < Status::low_score) and (opt_uchimeout != nullptr))
            {
              fprintf(fp_uchimeout, "0.0000\t");

              header_fprint_strip(fp_uchimeout,
                                  ci->query_head.data(),
                                  ci->query_head_len,
                                  opt_xsize,
                                  opt_xee,
                                  opt_xlength);

              if (opt_uchimeout5 != 0)
                {
                  fprintf(fp_uchimeout,
                          "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
                }
              else
                {
                  fprintf(fp_uchimeout,
                          "\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
                }
            }

          if (opt_nonchimeras != nullptr)
            {
              fasta_print_general(fp_nonchimeras,
                                  nullptr,
                                  ci->query_seq.data(),
                                  ci->query_len,
                                  ci->query_head.data(),
                                  ci->query_head_len,
                                  ci->query_size,
                                  nonchimera_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( (opt_uchime_ref != nullptr) ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h);
            }
        }

      if (status < Status::suspicious)
        {
          /* uchime_denovo: add non-chimeras to db */
          if ((opt_uchime_denovo != nullptr) or (opt_uchime2_denovo != nullptr) or (opt_uchime3_denovo != nullptr) or (opt_chimeras_denovo != nullptr))
            {
              dbindex_addsequence(seqno, opt_qmask);
            }
        }

      for (auto i = 0; i < ci->cand_count; ++i)
        {
          if (ci->nwcigar[i] != nullptr)
            {
              xfree(ci->nwcigar[i]);
            }
        }

      if (opt_uchime_ref != nullptr)
        {
          progress = fasta_get_position(query_fasta_h);
        }
      else
        {
          progress += db_getsequencelen(seqno);
        }

      progress_update(progress);

      ++seqno;

      xpthread_mutex_unlock(&mutex_output);
    }

  chimera_thread_exit(ci);

  xfree(scorematrix);

  return 0;
}


auto chimera_thread_worker(void * vp) -> void *
{
  return (void *) chimera_thread_core(cia + (int64_t) vp);
}


auto chimera_threads_run() -> void
{
  xpthread_attr_init(&attr);
  xpthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* create worker threads */
  for (int64_t t = 0; t < opt_threads; ++t)
    {
      xpthread_create(pthread + t, & attr,
                      chimera_thread_worker, (void *) t);
    }

  /* finish worker threads */
  for (int64_t t = 0; t < opt_threads; ++t)
    {
      xpthread_join(pthread[t], nullptr);
    }

  xpthread_attr_destroy(&attr);
}


auto open_chimera_file(std::FILE ** output_stream, char * name) -> void
{
  if (name != nullptr)
    {
      *output_stream = fopen_output(name);
      if (*output_stream == nullptr)
        {
          fatal("Unable to open file %s for writing", name);
        }
    }
  else
    {
      *output_stream = nullptr;
    }
}


auto close_chimera_file(std::FILE * output_stream) -> void
{
  if (output_stream != nullptr)
    {
      std::fclose(output_stream);
    }
}


auto chimera(struct Parameters const & parameters) -> void
{
  open_chimera_file(&fp_chimeras, opt_chimeras);
  open_chimera_file(&fp_nonchimeras, opt_nonchimeras);
  open_chimera_file(&fp_borderline, opt_borderline);

  if (parameters.opt_chimeras_denovo != nullptr)
    {
      open_chimera_file(&fp_uchimealns, opt_alnout);
      open_chimera_file(&fp_uchimeout, opt_tabbedout);
    }
  else
    {
      open_chimera_file(&fp_uchimealns, opt_uchimealns);
      open_chimera_file(&fp_uchimeout, opt_uchimeout);
    }


  /* override any options the user might have set */
  opt_maxaccepts = few;
  opt_maxrejects = rejects;
  opt_id = chimera_id;

  if (parameters.opt_strand)
    {
      fatal("Only --strand plus is allowed with uchime_ref.");
    }

  if (parameters.opt_uchime_ref == nullptr)
    {
      opt_self = 1;
      opt_selfid = 1;
      opt_threads = 1;
      opt_maxsizeratio = 1.0 / opt_abskew;
    }

  tophits = opt_maxaccepts + opt_maxrejects;

  uint64_t progress_total = 0;
  chimera_count = 0;
  nonchimera_count = 0;
  progress = 0;
  seqno = 0;

  /* prepare threads */
  std::vector<pthread_t> pthread_v(opt_threads);
  pthread = pthread_v.data();
  std::vector<struct chimera_info_s> cia_v(opt_threads);
  cia = cia_v.data();

  /* init mutexes for input and output */
  xpthread_mutex_init(&mutex_input, nullptr);
  xpthread_mutex_init(&mutex_output, nullptr);

  char * denovo_dbname = nullptr;

  /* prepare queries / database */
  if (parameters.opt_uchime_ref != nullptr)
    {
      /* check if the reference database may be an UDB file */

      auto const is_udb = udb_detect_isudb(opt_db);

      if (is_udb)
        {
          udb_read(opt_db, true, true);
        }
      else
        {
          db_read(opt_db, 0);
          if (opt_dbmask == MASK_DUST)
            {
              dust_all();
            }
          else if ((opt_dbmask == MASK_SOFT) and (opt_hardmask != 0))
            {
              hardmask_all();
            }
          dbindex_prepare(1, opt_dbmask);
          dbindex_addallsequences(opt_dbmask);
        }

      query_fasta_h = fasta_open(parameters.opt_uchime_ref);
      progress_total = fasta_get_size(query_fasta_h);
    }
  else
    {

      if (parameters.opt_uchime_denovo != nullptr)
        {
          denovo_dbname = parameters.opt_uchime_denovo;
        }
      else if (parameters.opt_uchime2_denovo != nullptr)
        {
          denovo_dbname = parameters.opt_uchime2_denovo;
        }
      else if (parameters.opt_uchime3_denovo != nullptr)
        {
          denovo_dbname = parameters.opt_uchime3_denovo;
        }
      else if (parameters.opt_chimeras_denovo != nullptr)
        {
          denovo_dbname = parameters.opt_chimeras_denovo;
        }
      else {
        fatal("Internal error");
      }

      db_read(denovo_dbname, 0);

      if (parameters.opt_qmask == MASK_DUST)
        {
          dust_all();
        }
      else if ((parameters.opt_qmask == MASK_SOFT) and (opt_hardmask != 0))
        {
          hardmask_all();
        }

      db_sortbyabundance();
      dbindex_prepare(1, parameters.opt_qmask);
      progress_total = db_getnucleotidecount();
    }

  if (parameters.opt_log != nullptr)
    {
      if ((parameters.opt_uchime_ref != nullptr) or (parameters.opt_uchime_denovo != nullptr))
        {
          fprintf(fp_log, "%8.2f  minh\n", opt_minh);
        }
      auto const is_a_uchime_command = (parameters.opt_uchime_ref != nullptr) or
        (parameters.opt_uchime_denovo != nullptr) or
        (parameters.opt_uchime2_denovo != nullptr) or
        (parameters.opt_uchime3_denovo != nullptr);
      if (is_a_uchime_command)
        {
          fprintf(fp_log, "%8.2f  xn\n", opt_xn);
          fprintf(fp_log, "%8.2f  dn\n", opt_dn);
          fprintf(fp_log, "%8.2f  xa\n", 1.0);
        }

      if ((parameters.opt_uchime_ref != nullptr) or (parameters.opt_uchime_denovo != nullptr))
        {
          fprintf(fp_log, "%8.2f  mindiv\n", opt_mindiv);
        }

      fprintf(fp_log, "%8.2f  id\n", opt_id);

      if (is_a_uchime_command)
        {
          fprintf(fp_log, "%8d  maxp\n", 2);
        }

      fprintf(fp_log, "\n");
    }


  progress_init("Detecting chimeras", progress_total);

  chimera_threads_run();

  progress_done();

  if (not parameters.opt_quiet)
    {
      if (total_count > 0)
        {
          if (parameters.opt_chimeras_denovo != nullptr)
            {
              fprintf(stderr,
                      "Found %d (%.1f%%) chimeras and "
                      "%d (%.1f%%) non-chimeras "
                      "in %u unique sequences.\n",
                      chimera_count,
                      100.0 * chimera_count / total_count,
                      nonchimera_count,
                      100.0 * nonchimera_count / total_count,
                      total_count);
            }
          else
            {
              fprintf(stderr,
                      "Found %d (%.1f%%) chimeras, "
                      "%d (%.1f%%) non-chimeras,\n"
                      "and %d (%.1f%%) borderline sequences "
                      "in %u unique sequences.\n",
                      chimera_count,
                      100.0 * chimera_count / total_count,
                      nonchimera_count,
                      100.0 * nonchimera_count / total_count,
                      borderline_count,
                      100.0 * borderline_count / total_count,
                      total_count);
            }
        }
      else
        {
          if (parameters.opt_chimeras_denovo != nullptr)
            {
              fprintf(stderr,
                      "Found %d chimeras and "
                      "%d non-chimeras "
                      "in %u unique sequences.\n",
                      chimera_count,
                      nonchimera_count,
                      total_count);
            }
          else
            {
              fprintf(stderr,
                      "Found %d chimeras, "
                      "%d non-chimeras,\n"
                      "and %d borderline sequences "
                      "in %u unique sequences.\n",
                      chimera_count,
                      nonchimera_count,
                      borderline_count,
                      total_count);
            }
        }

      if (total_abundance > 0)
        {
          if (parameters.opt_chimeras_denovo != nullptr)
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " (%.1f%%) chimeras and "
                      "%" PRId64 " (%.1f%%) non-chimeras "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      100.0 * chimera_abundance / total_abundance,
                      nonchimera_abundance,
                      100.0 * nonchimera_abundance / total_abundance,
                      total_abundance);
            }
          else
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " (%.1f%%) chimeras, "
                      "%" PRId64 " (%.1f%%) non-chimeras,\n"
                      "and %" PRId64 " (%.1f%%) borderline sequences "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      100.0 * chimera_abundance / total_abundance,
                      nonchimera_abundance,
                      100.0 * nonchimera_abundance / total_abundance,
                      borderline_abundance,
                      100.0 * borderline_abundance / total_abundance,
                      total_abundance);
            }
        }
      else
        {
          if (parameters.opt_chimeras_denovo != nullptr)
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " chimeras, "
                      "%" PRId64 " non-chimeras "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      nonchimera_abundance,
                      total_abundance);
            }
          else
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " chimeras, "
                      "%" PRId64 " non-chimeras,\n"
                      "and %" PRId64 " borderline sequences "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      nonchimera_abundance,
                      borderline_abundance,
                      total_abundance);
            }
        }
    }

  if (parameters.opt_log != nullptr)
    {
      if (parameters.opt_uchime_ref != nullptr)
        {
          fprintf(fp_log, "%s", parameters.opt_uchime_ref);
        }
      else
        {
          fprintf(fp_log, "%s", denovo_dbname);
        }

      if (seqno > 0)
        {
          fprintf(fp_log, ": %d/%u chimeras (%.1f%%)\n",
                  chimera_count,
                  seqno,
                  100.0 * chimera_count / seqno);
        }
      else
        {
          fprintf(fp_log, ": %d/%u chimeras\n",
                  chimera_count,
                  seqno);
        }
    }


  if (parameters.opt_uchime_ref != nullptr)
    {
      fasta_close(query_fasta_h);
    }

  dbindex_free();
  db_free();

  xpthread_mutex_destroy(&mutex_output);
  xpthread_mutex_destroy(&mutex_input);

  close_chimera_file(fp_borderline);
  close_chimera_file(fp_uchimeout);
  close_chimera_file(fp_uchimealns);
  close_chimera_file(fp_nonchimeras);
  close_chimera_file(fp_chimeras);

  show_rusage();
}
