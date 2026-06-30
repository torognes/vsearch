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
#include "attributes.h"
#include "chimera.h"
#include "dbindex.h"
#include "linmemalign.h"
#include "mask.h"
#include "minheap.h"
#include "udb.h"
#include "unique.h"
#include "utils/cigar.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/span.hpp"
#include "utils/threads.hpp"
#include <algorithm>  // std::copy, std::fill, std::fill_n, std::max, std::max_element, std::min, std::transform
#include <array>
#include <cassert>
#include <cctype>  // std::tolower
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdlib>  // std::qsort
#include <cstdio>  // std::FILE, std::fprintf, std::sscanf
#include <cstring>  // std::strlen, std::strcpy
#include <iterator>  // std::next
#include <limits>
#include <memory>
#include <mutex>  // std::mutex, std::lock_guard, std::unique_lock
#include <numeric>  // std::accumulate
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
constexpr auto maxparts = 100;
constexpr auto window = 32;
constexpr auto few = 4;
constexpr auto maxcandidates = few * maxparts;
constexpr auto rejects = 16;
constexpr auto chimera_id = 0.55;
static int tophits;
static fastx_handle query_fasta_h;

/* global data protected by mutex_output (file output and the global
   stats counters). mutex_input is not here: it only serializes input
   reading on the CLI path and is owned as a local by chimera_threads_run. */
static std::mutex mutex_output;
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
  int part_alloc = 0; /* the longest query part allocated memory for */

  int query_no = 0;
  std::vector<char> query_head;
  int query_head_len = 0;
  int64_t query_size = 0;
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
  std::vector<char *> nwcigar = std::vector<char *>(maxcandidates);  // this is a test

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
  std::vector<bool> ignore;

  double best_h = 0;

  int parts = 0;  /* number of query parts for chimera detection */

  /* API result fields — populated by eval_parents when result_out is set */
  struct chimera_result_s * result_out = nullptr;

  /* API per-thread working state — initialized by chimera_detect_init,
     reused across chimera_detect_single calls */
  std::vector<struct hit> api_allhits_list;
  std::unique_ptr<LinearMemoryAligner> api_lma_ptr;
};


static struct chimera_info_s * cia;


enum struct Status : unsigned char {
  no_parents,   // (0) non-chimeric
  no_alignment, // (1) score < 0, non-chimeric
  low_score,    // (2) score < minh, non-chimeric
  suspicious,   // (3) score >= minh, not available with uchime2_denovo and uchime3_denovo
  chimeric      // (4) score >= minh && divdiff >= opt_mindiv && ...
};


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

}  // end of anonymous namespace


auto realloc_arrays(struct chimera_info_s * chimera_info) -> void
{
  if (opt_chimeras_denovo != nullptr)
    {
      if (opt_chimeras_parts == 0) {
        chimera_info->parts = (chimera_info->query_len + 99) / 100;
      }
      else {
        chimera_info->parts = opt_chimeras_parts;
      }
      if (chimera_info->parts < 2) {
        chimera_info->parts = 2;
      }
      else if (chimera_info->parts > maxparts) {
        chimera_info->parts = maxparts;
      }
    }
  else
    {
      /* default for uchime, uchime2, and uchime3 */
      chimera_info->parts = 4;
    }

  const int maxhlen = std::max(chimera_info->query_head_len, 1);
  if (maxhlen > chimera_info->head_alloc)
    {
      chimera_info->query_head.resize(static_cast<size_t>(maxhlen) + 1);
    }
  chimera_info->head_alloc = std::max(chimera_info->head_alloc, maxhlen);

  /* realloc arrays based on query length */

  const int maxqlen = std::max(chimera_info->query_len, 1);
  auto const max_2x2_size = static_cast<size_t>(maxcandidates) * static_cast<size_t>(maxqlen);

  if (maxqlen > chimera_info->query_alloc)
    {
      chimera_info->query_alloc = maxqlen;

      chimera_info->query_seq.resize(static_cast<size_t>(maxqlen) + 1);

      chimera_info->maxi.resize(static_cast<size_t>(maxqlen) + 1);
      chimera_info->maxsmooth.resize(static_cast<size_t>(maxqlen));
      chimera_info->match.resize(max_2x2_size);
      chimera_info->insert.resize(max_2x2_size);
      chimera_info->smooth.resize(max_2x2_size);

      chimera_info->scan_p.resize(static_cast<size_t>(maxqlen) + 1);
      chimera_info->scan_q.resize(static_cast<size_t>(maxqlen) + 1);

      const int64_t maxalnlen = static_cast<int64_t>(maxqlen) + (2 * static_cast<int64_t>(db_getlongestsequence()));
      chimera_info->paln.resize(maxparents);
      for (auto & a_parent_alignment : chimera_info->paln) {
        a_parent_alignment.resize(static_cast<size_t>(maxalnlen) + 1);
      }
      chimera_info->qaln.resize(static_cast<size_t>(maxalnlen) + 1);
      chimera_info->diffs.resize(static_cast<size_t>(maxalnlen) + 1);
      chimera_info->votes.resize(static_cast<size_t>(maxalnlen) + 1);
      chimera_info->model.resize(static_cast<size_t>(maxalnlen) + 1);
      chimera_info->ignore.resize(static_cast<size_t>(maxalnlen) + 1);
    }

  // resize query parts if longer than earlier, minimum 100
  const int maxpartlen =
    std::max((maxqlen + chimera_info->parts - 1) / chimera_info->parts, 100);
  if (maxpartlen > chimera_info->part_alloc)
    {
      for (auto & query_info: chimera_info->si)
        {
          query_info.qsequence_v.resize(static_cast<size_t>(maxpartlen) + 1);
          query_info.qsequence = query_info.qsequence_v.data();
        }
      chimera_info->part_alloc = maxpartlen;
    }
}


auto reset_matches(struct chimera_info_s * a_chimera_info) -> void {
  // refactoring: initialization to zero? (useless), or reset to zero??
  std::fill(a_chimera_info->match.begin(), a_chimera_info->match.end(), 0);
  std::fill(a_chimera_info->insert.begin(), a_chimera_info->insert.end(), 0);
}


auto find_matches(struct chimera_info_s * chimera_info) -> void
{
  /* find the positions with matches for each potential parent */
  /* also note the positions with inserts in front */

  auto const & qseq = chimera_info->query_seq;

  for (auto i = 0; i < chimera_info->cand_count; ++i)
    {
      auto const * tseq = db_getsequence(chimera_info->cand_list[static_cast<size_t>(i)]);

      auto qpos = 0;
      auto tpos = 0;

      auto * cigar_start = chimera_info->nwcigar[static_cast<size_t>(i)];
      auto const cigar_length = std::strlen(cigar_start);
      auto const cigar_pairs = parse_cigar_string(Span<char>{cigar_start, cigar_length});

      for (auto const & a_pair: cigar_pairs) {
        auto const operation = a_pair.first;
        auto const runlength = a_pair.second;
        switch (operation) {
        case Operation::match:
          for (auto j = 0; j < runlength; ++j)
            {
              if ((map_4bit(qseq[static_cast<size_t>(qpos)]) &
                   map_4bit(tseq[tpos])) != 0U)
                {
                  chimera_info->match[static_cast<size_t>((i * chimera_info->query_len) + qpos)] = 1;
                }
              ++qpos;
              ++tpos;
            }
          break;

        case Operation::insertion:
          chimera_info->insert[static_cast<size_t>((i * chimera_info->query_len) + qpos)] = static_cast<int>(runlength);
          tpos += static_cast<int>(runlength);
          break;

        case Operation::deletion:
          qpos += static_cast<int>(runlength);
          break;
        }
      }
    }
}


struct parents_info_s
{
  int cand = -1;
  int start = -1;
  int len = 0;
};


auto compare_positions(const void * a, const void * b) -> int
{
  const int lhs = static_cast<parents_info_s const *>(a)->start;
  const int rhs = static_cast<parents_info_s const *>(b)->start;

  if (lhs < rhs) {
    return -1;
  }
  if (lhs > rhs) {
    return +1;
  }
  return 0;
}


auto scan_matches(struct chimera_info_s * ci,
                  int const * matches,
                  int const len,
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
    p[static_cast<size_t>(i) + 1] = p[static_cast<size_t>(i)] + ((matches[i] != 0) ? score_match : score_mismatch);
  }

  q[static_cast<size_t>(len)] = p[static_cast<size_t>(len)];
  for (auto i = len - 1; i >= 0; --i) {
    q[static_cast<size_t>(i)] = std::max(q[static_cast<size_t>(i) + 1], p[static_cast<size_t>(i)]);
  }

  auto best_i = 0;
  auto best_d = -1;
  auto best_c = -1.0;
  auto i = 1;
  auto j = 1;
  while (j <= len)
    {
      auto const c = q[static_cast<size_t>(j)] - p[static_cast<size_t>(i - 1)];
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

  std::vector<struct parents_info_s> best_parents(maxparents);
  std::vector<bool> position_used(static_cast<size_t>(ci->query_len), false);

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
          int j = 0;
          while (j < ci->query_len)
            {
              int start = j;
              int len = 0;
              while ((j < ci->query_len) &&
                     (not position_used[static_cast<size_t>(j)]) &&
                     ((len == 0) or (ci->insert[static_cast<size_t>((i * ci->query_len) + j)] == 0)))
                {
                  ++len;
                  ++j;
                }
              if (len > best_len)
                {
                  int scan_best_start = 0;
                  int scan_best_len = 0;
                  if (scan_matches(ci,
                                   &ci->match[static_cast<size_t>((i * ci->query_len) + start)],
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
          best_parents[static_cast<size_t>(f)].cand = best_cand;
          best_parents[static_cast<size_t>(f)].start = best_start;
          best_parents[static_cast<size_t>(f)].len = best_len;
          ++parents_found;

#if 0
          if (f == 0)
            std::printf("\n");
          std::printf("Best parents long: %d %d %d %d %s %s\n",
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
              position_used[static_cast<size_t>(j)] = true;
            }
          pos_remaining -= best_len;
        }
      else {
        break;
      }
    }

  /* sort parents by position (skip when empty: qsort requires a
     non-null pointer even for zero elements) */
  if (parents_found > 0)
    {
      std::qsort(best_parents.data(),
                 static_cast<size_t>(parents_found),
                 sizeof(struct parents_info_s),
                 compare_positions);
    }

  ci->parents_found = parents_found;

  for (int f = 0; f < parents_found; ++f)
    {
      ci->best_parents[static_cast<size_t>(f)] = best_parents[static_cast<size_t>(f)].cand;
      ci->best_start[static_cast<size_t>(f)] = best_parents[static_cast<size_t>(f)].start;
      ci->best_len[static_cast<size_t>(f)] = best_parents[static_cast<size_t>(f)].len;
    }

#if 0
  if (pos_remaining == 0)
    std::printf("Fully covered!\n");
  else
    std::printf("Not covered completely (%d).\n", pos_remaining);
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
      best_parent_cand[static_cast<size_t>(f)] = -1;
      ci->best_parents[static_cast<size_t>(f)] = -1;
    }

  std::vector<bool> cand_selected(static_cast<size_t>(ci->cand_count), false);

  for (int f = 0; f < 2; ++f)
    {
      if (f > 0)
        {
          /* for all parents except the first */

          /* wipe out matches for all candidates in positions
             covered by the previous parent */

          for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
            {
              int const z = (best_parent_cand[static_cast<size_t>(f - 1)] * ci->query_len) + qpos;
              if (ci->smooth[static_cast<size_t>(z)] == ci->maxsmooth[static_cast<size_t>(qpos)])
                {
                  for (int i = qpos + 1 - window; i <= qpos; ++i)
                    {
                      for (int j = 0; j < ci->cand_count; ++j)
                        {
                          ci->match[static_cast<size_t>((j * ci->query_len) + i)] = 0;
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
          if (not cand_selected[static_cast<size_t>(i)])
            {
              int sum = 0;
              for (int qpos = 0; qpos < ci->query_len; ++qpos)
                {
                  size_t const z = (static_cast<size_t>(i) * static_cast<size_t>(ci->query_len)) + static_cast<size_t>(qpos);
                  sum += ci->match[z];
                  if (qpos >= window)
                    {
                      sum -= ci->match[z - static_cast<size_t>(window)];
                    }
                  if (qpos >= window - 1)
                    {
                      ci->smooth[z] = sum;
                      ci->maxsmooth[static_cast<size_t>(qpos)] = std::max(ci->smooth[z], ci->maxsmooth[static_cast<size_t>(qpos)]);
                    }
                }
            }
        }


      /* find parent with the most wins */

      std::vector<int> wins(static_cast<size_t>(ci->cand_count), 0);

      for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
        {
          if (ci->maxsmooth[static_cast<size_t>(qpos)] != 0)
            {
              for (int i = 0; i < ci->cand_count; ++i)
                {
                  if (not cand_selected[static_cast<size_t>(i)])
                    {
                      size_t const z = (static_cast<size_t>(i) * static_cast<size_t>(ci->query_len)) + static_cast<size_t>(qpos);
                      if (ci->smooth[z] == ci->maxsmooth[static_cast<size_t>(qpos)])
                        {
                          ++wins[static_cast<size_t>(i)];
                        }
                    }
                }
            }
        }

      /* select best parent based on most wins */

      int maxwins = 0;
      for (int i = 0; i < ci->cand_count; ++i)
        {
          int const w = wins[static_cast<size_t>(i)];
          if (w > maxwins)
            {
              maxwins = w;
              best_parent_cand[static_cast<size_t>(f)] = i;
            }
        }

      /* terminate loop if no parent found */

      if (best_parent_cand[static_cast<size_t>(f)] < 0) {
        break;
      }

#if 0
      std::printf("Query %d: Best parent (%d) candidate: %d. Wins: %d\n",
             ci->query_no, f, best_parent_cand[f], maxwins);
#endif

      ci->best_parents[static_cast<size_t>(f)] = best_parent_cand[static_cast<size_t>(f)];
      cand_selected[static_cast<size_t>(best_parent_cand[static_cast<size_t>(f)])] = true;
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


auto fill_max_alignment_length(struct chimera_info_s * chimera_info) -> void
{
  /* find max insertions in front of each position in the query sequence */

  std::fill(chimera_info->maxi.begin(), chimera_info->maxi.end(), 0);

  auto const count = static_cast<size_t>(chimera_info->parents_found);
  assert(count <= chimera_info->best_parents.size());
  auto const best_parents_view = Span<int>{chimera_info->best_parents.data(), count};
  for (auto const best_parent : best_parents_view) {
    auto pos = 0LL;
    auto * cigar_start = chimera_info->nwcigar[static_cast<size_t>(best_parent)];
    auto const cigar_length = std::strlen(cigar_start);
    auto const cigar_pairs = parse_cigar_string(Span<char>{cigar_start, cigar_length});

    for (auto const & a_pair: cigar_pairs) {
      auto const operation = a_pair.first;
      auto const runlength = a_pair.second;
      switch (operation) {
      case Operation::match:
      case Operation::deletion:
        pos += runlength;
        break;

      case Operation::insertion:
        assert(runlength <= std::numeric_limits<int>::max());
        chimera_info->maxi[static_cast<size_t>(pos)] = std::max(static_cast<int>(runlength), chimera_info->maxi[static_cast<size_t>(pos)]);
        break;
      }
    }
  }
}


auto fill_alignment_parents(struct chimera_info_s * ci) -> void
{
  /* fill in alignment strings for the parents */

  for (int i = 0; i < ci->parents_found; ++i)
    {
      auto & alignment = ci->paln[static_cast<size_t>(i)];
      int const cand = ci->best_parents[static_cast<size_t>(i)];
      int const target_seqno = static_cast<int>(ci->cand_list[static_cast<size_t>(cand)]);
      char const * target_seq = db_getsequence(static_cast<uint64_t>(target_seqno));

      auto is_inserted = false;
      int qpos = 0;
      int tpos = 0;
      int alnpos = 0;

      auto * cigar_start = ci->nwcigar[static_cast<size_t>(cand)];
      auto const cigar_length = std::strlen(cigar_start);
      auto const cigar_pairs = parse_cigar_string(Span<char>{cigar_start, cigar_length});
      for (auto const & a_pair: cigar_pairs) {
        auto const operation = a_pair.first;
        auto const runlength = a_pair.second;
        switch (operation) {
        case Operation::insertion:
          for (int j = 0; j < ci->maxi[static_cast<size_t>(qpos)]; ++j)
            {
              if (j < runlength)
                {
                  alignment[static_cast<size_t>(alnpos)] = map_uppercase(target_seq[tpos]);
                  ++tpos;
                  ++alnpos;
                }
              else
                {
                  alignment[static_cast<size_t>(alnpos)] = '-';
                  ++alnpos;
                }
            }
          is_inserted = true;
          break;

        case Operation::match:
        case Operation::deletion:
          for (int j = 0; j < runlength; ++j)
            {
              if (not is_inserted)
                {
                  std::fill_n(&alignment[static_cast<size_t>(alnpos)], ci->maxi[static_cast<size_t>(qpos)], '-');
                  alnpos += ci->maxi[static_cast<size_t>(qpos)];
                }

              if (operation == Operation::match)
                {
                  alignment[static_cast<size_t>(alnpos)] = map_uppercase(target_seq[tpos]);
                  ++tpos;
                  ++alnpos;
                }
              else
                {
                  alignment[static_cast<size_t>(alnpos)] = '-';
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
          std::fill_n(&alignment[static_cast<size_t>(alnpos)], ci->maxi[static_cast<size_t>(qpos)], '-');
          alnpos += ci->maxi[static_cast<size_t>(qpos)];
        }

      /* end of sequence string */
      alignment[static_cast<size_t>(alnpos)] = '\0';
    }
}


auto fill_in_alignment_string_for_query(struct chimera_info_s * chimera_info) -> void {
  auto alnpos = 0;
  auto qpos = 0;
  for (auto const nucleotide: chimera_info->query_seq) {
    // add insertion (if any):
    auto const insert_length = chimera_info->maxi[static_cast<size_t>(qpos)];
    std::fill_n(&chimera_info->qaln[static_cast<size_t>(alnpos)], insert_length, '-');
    alnpos += insert_length;

    // add (mis-)matching position:
    chimera_info->qaln[static_cast<size_t>(alnpos)] = map_uppercase(nucleotide);
    ++alnpos;
    ++qpos;
  }
  // add terminal gap (if any):
  auto const insert_length = chimera_info->maxi[static_cast<size_t>(chimera_info->query_len)];
  std::fill_n(&chimera_info->qaln[static_cast<size_t>(alnpos)], insert_length, '-');
  alnpos += insert_length;
  chimera_info->qaln[static_cast<size_t>(alnpos)] = '\0';
}


auto fill_in_model_string_for_query(struct chimera_info_s * chimera_info) -> void {
  int nth_parent = 0;
  auto alnpos = 0;
  for (int qpos = 0; qpos < chimera_info->query_len; ++qpos)
    {
      if (qpos >= (chimera_info->best_start[static_cast<size_t>(nth_parent)] + chimera_info->best_len[static_cast<size_t>(nth_parent)])) {
        ++nth_parent;
      }
      // add insertion (if any):
      auto const insert_length = chimera_info->maxi[static_cast<size_t>(qpos)];
      std::fill_n(&chimera_info->model[static_cast<size_t>(alnpos)], insert_length, static_cast<char>('A' + nth_parent));
      alnpos += insert_length;

      // add (mis-)matching position:
      chimera_info->model[static_cast<size_t>(alnpos)] = static_cast<char>('A' + nth_parent);
      ++alnpos;
    }
  // add terminal gap (if any):
  auto const insert_length = chimera_info->maxi[static_cast<size_t>(chimera_info->query_len)];
  std::fill_n(&chimera_info->model[static_cast<size_t>(alnpos)], insert_length, static_cast<char>('A' + nth_parent));
  alnpos += insert_length;
  chimera_info->model[static_cast<size_t>(alnpos)] = '\0';
}


auto count_matches_with_parents(struct chimera_info_s const * chimera_info,
                                int const alignment_length) -> std::array<int, maxparents> {
  std::array<int, maxparents> matches {{}};

  for (auto i = 0; i < alignment_length; ++i)
    {
      auto const qsym = map_4bit(chimera_info->qaln[static_cast<size_t>(i)]);

      for (auto f = 0; f < chimera_info->parents_found; ++f)
        {
          auto const psym = map_4bit(chimera_info->paln[static_cast<size_t>(f)][static_cast<size_t>(i)]);
          if (qsym == psym) {
            ++matches[static_cast<size_t>(f)];
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
                   std::vector<unsigned char> const & psym,
                   unsigned char const qsym) -> char {
  auto const all_defined = (qsym != 0U) and
    std::all_of(psym.begin(),
                psym.end(),
                [](unsigned char const symbol) -> bool{ return symbol != 0U; });

  char diff = ' ';

  if (not all_defined) { return diff; }

  auto z = 0;
  for (auto f = 0; f < ci->parents_found; ++f) {
    if (psym[static_cast<size_t>(f)] == qsym) {
      diff = static_cast<char>('A' + f);
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

  std::vector<unsigned char> psym;
  psym.reserve(maxparents);

  for (int i = 0; i < alnlen; ++i)
    {
      auto const qsym = map_4bit(ci->qaln[static_cast<size_t>(i)]);
      for (int f = 0; f < ci->parents_found; ++f) {
        psym.emplace_back(map_4bit(ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i)]));
      }

      /* lower case parent symbols that differ from query */

      for (int f = 0; f < ci->parents_found; ++f) {
        if ((psym[static_cast<size_t>(f)] != 0U) and (psym[static_cast<size_t>(f)] != qsym)) {
          ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i)] = static_cast<char>(std::tolower(ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i)]));
        }
      }

      /* compute diffs */
      ci->diffs[static_cast<size_t>(i)] = compute_diffs(ci, psym, qsym);
      psym.clear();
    }

  ci->diffs[static_cast<size_t>(alnlen)] = '\0';


  auto const match_QP = count_matches_with_parents(ci, alnlen);

  int const seqno_a = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[0])]);
  int const seqno_b = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[1])]);
  int const seqno_c = ci->parents_found > 2 ? static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[2])]) : -1;

  auto const QP = compute_global_similarities_with_parents(match_QP, alnlen);
  auto const QT = *std::max_element(QP.begin(), QP.end());

  double const QA = QP[0];
  double const QB = QP[1];
  double const QC = ci->parents_found > 2 ? QP[2] : 0.00;
  double const QM = 100.00;
  double const divfrac = 100.00 * (QM - QT) / QT;  // divergence of the model with the closest parent

  /* Populate API result struct if requested */
  if (ci->result_out != nullptr)
    {
      auto * r = ci->result_out;
      r->score = 99.9999;  /* chimeras_denovo always reports chimeric */
      std::snprintf(r->query_label, sizeof(r->query_label), "%.*s",
                    ci->query_head_len, ci->query_head.data());
      std::snprintf(r->parent_a_label, sizeof(r->parent_a_label), "%.*s",
                    static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))), db_getheader(static_cast<uint64_t>(seqno_a)));
      std::snprintf(r->parent_b_label, sizeof(r->parent_b_label), "%.*s",
                    static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))), db_getheader(static_cast<uint64_t>(seqno_b)));
      /* closest parent = max of QA, QB */
      if (QA >= QB)
        {
          std::snprintf(r->closest_parent_label, sizeof(r->closest_parent_label),
                        "%.*s", static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))), db_getheader(static_cast<uint64_t>(seqno_a)));
        }
      else
        {
          std::snprintf(r->closest_parent_label, sizeof(r->closest_parent_label),
                        "%.*s", static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))), db_getheader(static_cast<uint64_t>(seqno_b)));
        }
      r->id_query_model = QM;
      r->id_query_a = QA;
      r->id_query_b = QB;
      r->id_a_b = QC;  /* AB not computed in long path; use QC as proxy */
      r->id_query_top = QT;
      r->left_yes = 0;
      r->left_no = 0;
      r->left_abstain = 0;
      r->right_yes = 0;
      r->right_no = 0;
      r->right_abstain = 0;
      r->divergence = divfrac;
      r->flag = 'Y';  /* eval_parents_long is always chimeric */
    }

  std::lock_guard<std::mutex> const output_lock(mutex_output);

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
          int const parent_seqno = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[static_cast<size_t>(f)])]);
          std::fprintf(fp_uchimealns, "\nParent%c (%5" PRIu64 " nt) ",
                       'A' + f,
                       db_getsequencelen(static_cast<uint64_t>(parent_seqno)));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(static_cast<uint64_t>(parent_seqno)),
                              static_cast<int>(db_getheaderlen(static_cast<uint64_t>(parent_seqno))),
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
              if (ci->qaln[static_cast<size_t>(i + j)] != '-')
                {
                  ++qnt;
                }

              for (int f = 0; f < ci->parents_found; ++f) {
                if (ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i + j)] != '-')
                  {
                    ++pnt[static_cast<size_t>(f)];
                  }
              }
            }

          std::fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                  qpos + 1, w, &ci->qaln[static_cast<size_t>(i)], qpos + qnt);

          for (int f = 0; f < ci->parents_found; ++f)
            {
              std::fprintf(fp_uchimealns, "%c %5d %.*s %d\n",
                      'A' + f,
                      ppos[static_cast<size_t>(f)] + 1, w, &ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i)], ppos[static_cast<size_t>(f)] + pnt[static_cast<size_t>(f)]);
            }

          std::fprintf(fp_uchimealns, "Diffs   %.*s\n", w, &ci->diffs[static_cast<size_t>(i)]);
          std::fprintf(fp_uchimealns, "Model   %.*s\n", w, &ci->model[static_cast<size_t>(i)]);
          std::fprintf(fp_uchimealns, "\n");

          rest -= width;
          qpos += qnt;
          for (int f = 0; f < ci->parents_found; ++f) {
            ppos[static_cast<size_t>(f)] += pnt[static_cast<size_t>(f)];
          }
        }

      std::fprintf(fp_uchimealns, "Ids.  QA %.2f%%, QB %.2f%%, QC %.2f%%, "
              "QT %.2f%%, QModel %.2f%%, Div. %+.2f%%\n",
              QA, QB, QC, QT, QM, divfrac);
    }

  if (opt_tabbedout != nullptr)
    {
      std::fprintf(fp_uchimeout, "%.4f\t", 99.9999);

      header_fprint_strip(fp_uchimeout,
                          ci->query_head.data(),
                          ci->query_head_len,
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      std::fprintf(fp_uchimeout, "\t");
      header_fprint_strip(fp_uchimeout,
                          db_getheader(static_cast<uint64_t>(seqno_a)),
                          static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      std::fprintf(fp_uchimeout, "\t");
      header_fprint_strip(fp_uchimeout,
                          db_getheader(static_cast<uint64_t>(seqno_b)),
                          static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      std::fprintf(fp_uchimeout, "\t");
      if (seqno_c >= 0)
        {
          header_fprint_strip(fp_uchimeout,
                              db_getheader(static_cast<uint64_t>(seqno_c)),
                              static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_c))),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
        }
      else
        {
          std::fprintf(fp_uchimeout, "*");
        }
      std::fprintf(fp_uchimeout, "\t");

      std::fprintf(fp_uchimeout,
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
      for (int j = 0; j < ci->maxi[static_cast<size_t>(i)]; ++j)
        {
          *q = '-';
          ++q;
        }
      *q = map_uppercase(ci->query_seq[static_cast<size_t>(qpos)]);
      ++qpos;
      ++q;
    }
  for (int j = 0; j < ci->maxi[static_cast<size_t>(ci->query_len)]; ++j)
    {
      *q = '-';
      ++q;
    }
  *q = 0;

  /* mark positions to ignore in voting */
  std::fill(ci->ignore.begin(), ci->ignore.end(), false);

  for (int i = 0; i < alnlen; ++i)
    {
      auto const qsym  = map_4bit(ci->qaln[static_cast<size_t>(i)]);
      auto const p1sym = map_4bit(ci->paln[0][static_cast<size_t>(i)]);
      auto const p2sym = map_4bit(ci->paln[1][static_cast<size_t>(i)]);

      /* ignore gap positions and those next to the gap */
      if ((qsym == 0U) or (p1sym == 0U) or (p2sym == 0U))
        {
          ci->ignore[static_cast<size_t>(i)] = true;
          if (i > 0)
            {
              ci->ignore[static_cast<size_t>(i - 1)] = true;
            }
          if (i < alnlen - 1)
            {
              ci->ignore[static_cast<size_t>(i + 1)] = true;
            }
        }

      /* ignore ambiguous symbols */
      if (is_ambiguous_4bit(qsym) or
          is_ambiguous_4bit(p1sym) or
          is_ambiguous_4bit(p2sym))
        {
          ci->ignore[static_cast<size_t>(i)] = true;
        }

      /* lower case parent symbols that differ from query */

      if ((p1sym != 0U) and (p1sym != qsym))
        {
          ci->paln[0][static_cast<size_t>(i)] = static_cast<char>(std::tolower(ci->paln[0][static_cast<size_t>(i)]));
        }

      if ((p2sym != 0U) and (p2sym != qsym))
        {
          ci->paln[1][static_cast<size_t>(i)] = static_cast<char>(std::tolower(ci->paln[1][static_cast<size_t>(i)]));
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

      ci->diffs[static_cast<size_t>(i)] = diff;
    }

  ci->diffs[static_cast<size_t>(alnlen)] = '\0';

  /* compute score */

  int sumA = 0;
  int sumB = 0;
  int sumN = 0;

  // refactoring: extract to function, use a struct to pass results
  // std::transform(ci->diffs.begin(),
  //                std::next(ci->diffs.begin(), alnlen),
  //                ci->ignore.begin(),
  //                [&sumA, &sumB, &sumN](char const diff, bool const is_ignored) -> void {
  //                         if (is_ignored) { return; }
  //                         if (diff == 'A') {
  //                             ++sumA;
  //                            }
  //                         else if (diff == 'B') {
  //                             ++sumB;
  //                            }
  //                         else if (diff != ' ') {
  //                             ++sumN;
  //                         }
  //                         return;
  //                 }
  //                );

  for (auto i = 0; i < alnlen; ++i)
    {
      if (ci->ignore[static_cast<size_t>(i)]) { continue; }
      auto const diff = ci->diffs[static_cast<size_t>(i)];

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
      if (not ci->ignore[static_cast<size_t>(i)])
        {
          char const diff = ci->diffs[static_cast<size_t>(i)];
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
                  left_h = left_y / ((opt_xn * (left_n + opt_dn)) + left_a);
                  right_h = right_y / ((opt_xn * (right_n + opt_dn)) + right_a);
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

                  left_h = left_n / ((opt_xn * (left_y + opt_dn)) + left_a);
                  right_h = right_n / ((opt_xn * (right_y + opt_dn)) + right_a);
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
              char const diff = ci->diffs[static_cast<size_t>(i)];
              if (diff == 'A')
                {
                  ci->diffs[static_cast<size_t>(i)] = 'B';
                }
              else if (diff == 'B')
                {
                  ci->diffs[static_cast<size_t>(i)] = 'A';
                }
            }
        }

      /* fill in votes and model */

      for (int i = 0; i < alnlen; ++i)
        {
          char const m = i <= best_i ? 'A' : 'B';
          ci->model[static_cast<size_t>(i)] = m;

          char v = ' ';
          if (not ci->ignore[static_cast<size_t>(i)])
            {
              char const d = ci->diffs[static_cast<size_t>(i)];

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
          ci->votes[static_cast<size_t>(i)] = v;

          /* lower case diffs for no votes */
          if (v == '!')
            {
              ci->diffs[static_cast<size_t>(i)] = static_cast<char>(std::tolower(ci->diffs[static_cast<size_t>(i)]));
            }
        }

      /* fill in crossover region */

      for (int i = best_i + 1; i < alnlen; ++i)
        {
          if ((ci->diffs[static_cast<size_t>(i)] == ' ') or (ci->diffs[static_cast<size_t>(i)] == 'A'))
            {
              ci->model[static_cast<size_t>(i)] = 'x';
            }
          else
            {
              break;
            }
        }

      ci->votes[static_cast<size_t>(alnlen)] = 0;
      ci->model[static_cast<size_t>(alnlen)] = 0;

      /* count matches */

      auto const index_a = best_is_reverse ? 1U : 0U;
      auto const index_b = best_is_reverse ? 0U : 1U;

      int match_QA = 0;
      int match_QB = 0;
      int match_AB = 0;
      int match_QM = 0;
      int cols = 0;

      for (auto i = 0; i < alnlen; i++)
        {
          if (not ci->ignore[static_cast<size_t>(i)])
            {
              ++cols;

              auto const qsym = map_4bit(ci->qaln[static_cast<size_t>(i)]);
              auto const asym = map_4bit(ci->paln[index_a][static_cast<size_t>(i)]);
              auto const bsym = map_4bit(ci->paln[index_b][static_cast<size_t>(i)]);
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

      int const seqno_a = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[index_a])]);
      int const seqno_b = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[index_b])]);

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

      /* Populate API result struct if requested */
      if (ci->result_out != nullptr)
        {
          auto * r = ci->result_out;
          r->score = best_h;
          std::snprintf(r->query_label, sizeof(r->query_label), "%.*s",
                        ci->query_head_len, ci->query_head.data());
          std::snprintf(r->parent_a_label, sizeof(r->parent_a_label), "%.*s",
                        static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))), db_getheader(static_cast<uint64_t>(seqno_a)));
          std::snprintf(r->parent_b_label, sizeof(r->parent_b_label), "%.*s",
                        static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))), db_getheader(static_cast<uint64_t>(seqno_b)));
          if (QA >= QB)
            {
              std::snprintf(r->closest_parent_label, sizeof(r->closest_parent_label),
                            "%.*s", static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))), db_getheader(static_cast<uint64_t>(seqno_a)));
            }
          else
            {
              std::snprintf(r->closest_parent_label, sizeof(r->closest_parent_label),
                            "%.*s", static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))), db_getheader(static_cast<uint64_t>(seqno_b)));
            }
          r->id_query_model = QM;
          r->id_query_a = QA;
          r->id_query_b = QB;
          r->id_a_b = AB;
          r->id_query_top = QT;
          r->left_yes = best_left_y;
          r->left_no = best_left_n;
          r->left_abstain = best_left_a;
          r->right_yes = best_right_y;
          r->right_no = best_right_n;
          r->right_abstain = best_right_a;
          r->divergence = divdiff;
          r->flag = (status == Status::chimeric) ? 'Y' :
                    (status == Status::low_score ? 'N' : '?');
        }

      /* print alignment */

      std::unique_lock<std::mutex> output_lock(mutex_output);

      if ((opt_uchimealns != nullptr) and (status == Status::chimeric))
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

          std::fprintf(fp_uchimealns, "\nParentA (%5" PRIu64 " nt) ",
                  db_getsequencelen(static_cast<uint64_t>(seqno_a)));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(static_cast<uint64_t>(seqno_a)),
                              static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);

          std::fprintf(fp_uchimealns, "\nParentB (%5" PRIu64 " nt) ",
                  db_getsequencelen(static_cast<uint64_t>(seqno_b)));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(static_cast<uint64_t>(seqno_b)),
                              static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          std::fprintf(fp_uchimealns, "\n\n");

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
                  if (ci->qaln[static_cast<size_t>(i + j)] != '-')
                    {
                      ++qnt;
                    }
                  if (ci->paln[0][static_cast<size_t>(i + j)] != '-')
                    {
                      ++p1nt;
                    }
                  if (ci->paln[1][static_cast<size_t>(i + j)] != '-')
                    {
                      ++p2nt;
                    }
                }

              if (not best_is_reverse)
                {
                  std::fprintf(fp_uchimealns, "A %5d %.*s %d\n",
                          p1pos + 1, w, &ci->paln[0][static_cast<size_t>(i)], p1pos + p1nt);
                  std::fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos + 1, w, &ci->qaln[static_cast<size_t>(i)], qpos + qnt);
                  std::fprintf(fp_uchimealns, "B %5d %.*s %d\n",
                          p2pos + 1, w, &ci->paln[1][static_cast<size_t>(i)], p2pos + p2nt);
                }
              else
                {
                  std::fprintf(fp_uchimealns, "A %5d %.*s %d\n",
                          p2pos + 1, w, &ci->paln[1][static_cast<size_t>(i)], p2pos + p2nt);
                  std::fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos + 1, w, &ci->qaln[static_cast<size_t>(i)], qpos + qnt);
                  std::fprintf(fp_uchimealns, "B %5d %.*s %d\n",
                          p1pos + 1, w, &ci->paln[0][static_cast<size_t>(i)], p1pos + p1nt);
                }

              std::fprintf(fp_uchimealns, "Diffs   %.*s\n", w, &ci->diffs[static_cast<size_t>(i)]);
              std::fprintf(fp_uchimealns, "Votes   %.*s\n", w, &ci->votes[static_cast<size_t>(i)]);
              std::fprintf(fp_uchimealns, "Model   %.*s\n", w, &ci->model[static_cast<size_t>(i)]);
              std::fprintf(fp_uchimealns, "\n");

              qpos += qnt;
              p1pos += p1nt;
              p2pos += p2nt;
              rest -= width;
            }

          std::fprintf(fp_uchimealns, "Ids.  QA %.1f%%, QB %.1f%%, AB %.1f%%, "
                  "QModel %.1f%%, Div. %+.1f%%\n",
                  QA, QB, AB, QM, divfrac);

          std::fprintf(fp_uchimealns, "Diffs Left %d: N %d, A %d, Y %d (%.1f%%); "
                  "Right %d: N %d, A %d, Y %d (%.1f%%), Score %.4f\n",
                  sumL, best_left_n, best_left_a, best_left_y,
                  100.0 * best_left_y / sumL,
                  sumR, best_right_n, best_right_a, best_right_y,
                  100.0 * best_right_y / sumR,
                  best_h);
        }

      if (opt_uchimeout != nullptr)
        {
          std::fprintf(fp_uchimeout, "%.4f\t", best_h);

          header_fprint_strip(fp_uchimeout,
                              ci->query_head.data(),
                              ci->query_head_len,
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          std::fprintf(fp_uchimeout, "\t");
          header_fprint_strip(fp_uchimeout,
                              db_getheader(static_cast<uint64_t>(seqno_a)),
                              static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          std::fprintf(fp_uchimeout, "\t");
          header_fprint_strip(fp_uchimeout,
                              db_getheader(static_cast<uint64_t>(seqno_b)),
                              static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          std::fprintf(fp_uchimeout, "\t");

          if (opt_uchimeout5 == 0)
            {
              if (QA >= QB)
                {
                  header_fprint_strip(fp_uchimeout,
                                      db_getheader(static_cast<uint64_t>(seqno_a)),
                                      static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_a))),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                }
              else
                {
                  header_fprint_strip(fp_uchimeout,
                                      db_getheader(static_cast<uint64_t>(seqno_b)),
                                      static_cast<int>(db_getheaderlen(static_cast<uint64_t>(seqno_b))),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                }
              std::fprintf(fp_uchimeout, "\t");
            }

          std::fprintf(fp_uchimeout,
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
      output_lock.unlock();
    }

  return status;
}


auto query_init(struct searchinfo_s * search_info) -> void
{
  static constexpr auto overflow_padding = 16U;  // 16 * sizeof(short) = 32 bytes
  search_info->hits_v.resize(static_cast<size_t>(tophits));
  search_info->hits = search_info->hits_v.data();
  search_info->kmers_v.reserve(db_getsequencecount() + overflow_padding);
  search_info->kmers_v.resize(db_getsequencecount());
  search_info->kmers = search_info->kmers_v.data();
  search_info->hit_count = 0;
  search_info->uh = unique_init();
  search_info->s = search16_init(static_cast<CELL>(opt_match),
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
  for (auto i = 0; i < chimera_info->parts; ++i)
    {
      auto const length =
        (rest + (chimera_info->parts - i - 1)) / (chimera_info->parts - i);

      auto & search_info = chimera_info->si[static_cast<size_t>(i)];

      search_info.query_no = chimera_info->query_no;
      search_info.strand = 0;
      search_info.qsize = chimera_info->query_size;
      search_info.query_head_len = chimera_info->query_head_len;
      search_info.query_head = chimera_info->query_head.data();
      search_info.qseqlen = length;
      assert(static_cast<std::size_t>(length) <= search_info.qsequence_v.size());
      std::copy(cursor, std::next(cursor, length), search_info.qsequence_v.begin());
      search_info.qsequence_v[static_cast<size_t>(length)] = '\0';

      rest -= length;
      cursor = std::next(cursor, length);
    }
}


auto chimera_thread_init(struct chimera_info_s * ci) -> void
{

  for (int i = 0; i < maxparts; ++i)
    {
      query_init(&ci->si[static_cast<size_t>(i)]);
    }

  ci->s = search16_init(static_cast<CELL>(opt_match),
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


auto chimera_thread_exit(struct chimera_info_s * ci) -> void
{
  search16_exit(ci->s);

  for (auto & a_search_info : ci->si) {
    query_exit(&a_search_info);
  }
}


/* Process a single query that has already been loaded into ci.
   Shared by chimera_thread_core (CLI) and chimera_detect_single (API).
   ci->query_seq, query_head, query_len, query_head_len, query_size must
   be populated. allhits_list must be pre-allocated to maxcandidates.
   lma is the per-thread linear memory aligner (fallback for SIMD overflow). */
static auto chimera_process_query(struct chimera_info_s * ci,
                                  std::vector<struct hit> & allhits_list,
                                  LinearMemoryAligner & lma) -> Status
{
  /* partition query */
  partition_query(ci);

  /* perform searches and collect candidate parents */
  ci->cand_count = 0;
  ci->best_h = 0.0;
  auto allhits_count = 0;

  if (ci->query_len >= ci->parts)
    {
      std::vector<struct hit> hits;
      for (auto i = 0; i < ci->parts; ++i)
        {
          search_onequery(&ci->si[static_cast<size_t>(i)], static_cast<int>(opt_qmask));
          search_joinhits(&ci->si[static_cast<size_t>(i)], nullptr, hits);
          for (auto & hit : hits) {
            if (hit.accepted and allhits_count < maxcandidates)
              {
                allhits_list[static_cast<size_t>(allhits_count)] = hit;
                ++allhits_count;
              }
            else
              {
                // Unallocate alignments for weak hits
                if (hit.nwalignment)
                  {
                    xfree(hit.nwalignment);
                    hit.nwalignment = nullptr;
                  }
              }
          }
          hits.clear();
        }
    }

  for (auto i = 0; i < allhits_count; ++i)
    {
      unsigned int const target = static_cast<unsigned int>(allhits_list[static_cast<size_t>(i)].target);

      /* skip duplicates */
      auto k = 0;
      for (k = 0; k < ci->cand_count; ++k)
        {
          if (ci->cand_list[static_cast<size_t>(k)] == target)
            {
              break;
            }
        }

      if (k == ci->cand_count)
        {
          ci->cand_list[static_cast<size_t>(ci->cand_count)] = target;
          ++ci->cand_count;
        }

      /* deallocate cigar */
      if (allhits_list[static_cast<size_t>(i)].nwalignment != nullptr)
        {
          xfree(allhits_list[static_cast<size_t>(i)].nwalignment);
          allhits_list[static_cast<size_t>(i)].nwalignment = nullptr;
        }
    }


  /* align full query to each candidate */

  search16_qprep(ci->s, ci->query_seq.data(), ci->query_len);

  search16(ci->s,
           static_cast<unsigned int>(ci->cand_count),
           ci->cand_list.data(),
           ci->snwscore.data(),
           ci->snwalignmentlength.data(),
           ci->snwmatches.data(),
           ci->snwmismatches.data(),
           ci->snwgaps.data(),
           ci->nwcigar.data());

  for (auto i = 0; i < ci->cand_count; ++i)
    {
      int64_t const target = ci->cand_list[static_cast<size_t>(i)];
      int64_t nwscore = ci->snwscore[static_cast<size_t>(i)];
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

          auto * tseq = db_getsequence(static_cast<uint64_t>(target));
          int64_t const tseqlen = static_cast<int64_t>(db_getsequencelen(static_cast<uint64_t>(target)));

          if (ci->nwcigar[static_cast<size_t>(i)] != nullptr)
            {
              xfree(ci->nwcigar[static_cast<size_t>(i)]);
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

          ci->nwcigar[static_cast<size_t>(i)] = nwcigar;
          ci->nwscore[static_cast<size_t>(i)] = nwscore;
          ci->nwalignmentlength[static_cast<size_t>(i)] = nwalignmentlength;
          ci->nwmatches[static_cast<size_t>(i)] = nwmatches;
          ci->nwmismatches[static_cast<size_t>(i)] = nwmismatches;
          ci->nwgaps[static_cast<size_t>(i)] = nwgaps;
        }
      else
        {
          ci->nwscore[static_cast<size_t>(i)] = ci->snwscore[static_cast<size_t>(i)];
          ci->nwalignmentlength[static_cast<size_t>(i)] = ci->snwalignmentlength[static_cast<size_t>(i)];
          ci->nwmatches[static_cast<size_t>(i)] = ci->snwmatches[static_cast<size_t>(i)];
          ci->nwmismatches[static_cast<size_t>(i)] = ci->snwmismatches[static_cast<size_t>(i)];
          ci->nwgaps[static_cast<size_t>(i)] = ci->snwgaps[static_cast<size_t>(i)];
        }
    }


  /* find the best pair of parents, then compute score for them */

  if (opt_chimeras_denovo != nullptr)
    {
      /* long high-quality reads */
      if (find_best_parents_long(ci) != 0)
        {
          return eval_parents_long(ci);
        }
      else
        {
          return Status::no_parents;
        }
    }
  else
    {
      if (find_best_parents(ci) != 0)
        {
          return eval_parents(ci);
        }
      else
        {
          return Status::no_parents;
        }
    }
}


auto chimera_thread_core(struct chimera_info_s * ci,
                         std::mutex & mutex_input) -> uint64_t
{
  chimera_thread_init(ci);

  std::vector<struct hit> allhits_list(maxcandidates);

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

  while (true)
    {
      /* get next sequence */

      std::unique_lock<std::mutex> input_lock(mutex_input);

      /* Query-file progress position. Read here, under input_lock, into a
         worker-local: fasta_get_position() reads the shared query handle,
         which another worker advances in fasta_next() under the same lock.
         Reading it later under mutex_output (as before) raced that writer. */
      uint64_t query_position = 0;

      if (opt_uchime_ref != nullptr)
        {
          if (fasta_next(query_fasta_h, (opt_notrunclabels == 0),
                         chrmap_no_change_vector.data()))
            {
              ci->query_head_len = static_cast<int>(fasta_get_header_length(query_fasta_h));
              ci->query_len = static_cast<int>(fasta_get_sequence_length(query_fasta_h));
              ci->query_no = static_cast<int>(fasta_get_seqno(query_fasta_h));
              ci->query_size = fasta_get_abundance(query_fasta_h);

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              /* copy the data locally (query seq, head) */
              std::strcpy(ci->query_head.data(), fasta_get_header(query_fasta_h));
              std::strcpy(ci->query_seq.data(), fasta_get_sequence(query_fasta_h));
              query_position = fasta_get_position(query_fasta_h);
            }
          else
            {
              break; /* end while loop; input_lock released by RAII */
            }
        }
      else
        {
          if (seqno < db_getsequencecount())
            {
              ci->query_no = static_cast<int>(seqno);
              ci->query_head_len = static_cast<int>(db_getheaderlen(seqno));
              ci->query_len = static_cast<int>(db_getsequencelen(seqno));
              ci->query_size = static_cast<int64_t>(db_getabundance(seqno));

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              std::strcpy(ci->query_head.data(), db_getheader(seqno));
              std::strcpy(ci->query_seq.data(), db_getsequence(seqno));
            }
          else
            {
              break; /* end while loop; input_lock released by RAII */
            }
        }

      input_lock.unlock();

      auto const status = chimera_process_query(ci, allhits_list, lma);

      /* output results */

      std::lock_guard<std::mutex> const output_lock(mutex_output);

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
                                  static_cast<uint64_t>(ci->query_size),
                                  chimera_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( (opt_uchime_ref != nullptr) ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h,
                                  0);

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
                                  static_cast<uint64_t>(ci->query_size),
                                  borderline_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( (opt_uchime_ref != nullptr) ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h,
                                  0);

            }
        }

      if (status < Status::suspicious)
        {
          ++nonchimera_count;
          nonchimera_abundance += ci->query_size;

          /* output no parents, no chimeras */
          if ((status < Status::low_score) and (opt_uchimeout != nullptr))
            {
              std::fprintf(fp_uchimeout, "%.4f\t", ci->best_h);

              header_fprint_strip(fp_uchimeout,
                                  ci->query_head.data(),
                                  ci->query_head_len,
                                  opt_xsize,
                                  opt_xee,
                                  opt_xlength);

              if (opt_uchimeout5 != 0)
                {
                  std::fprintf(fp_uchimeout,
                          "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
                }
              else
                {
                  std::fprintf(fp_uchimeout,
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
                                  static_cast<uint64_t>(ci->query_size),
                                  nonchimera_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( (opt_uchime_ref != nullptr) ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h,
                                  0);
            }
        }

      if (status < Status::suspicious)
        {
          /* uchime_denovo: add non-chimeras to db */
          if ((opt_uchime_denovo != nullptr) or (opt_uchime2_denovo != nullptr) or (opt_uchime3_denovo != nullptr) or (opt_chimeras_denovo != nullptr))
            {
              dbindex_addsequence(seqno, static_cast<int>(opt_qmask));
            }
        }

      for (auto i = 0; i < ci->cand_count; ++i)
        {
          if (ci->nwcigar[static_cast<size_t>(i)] != nullptr)
            {
              xfree(ci->nwcigar[static_cast<size_t>(i)]);
            }
        }

      if (opt_uchime_ref != nullptr)
        {
          progress = query_position;
        }
      else
        {
          progress += db_getsequencelen(seqno);
        }

      progress_update(progress);

      ++seqno;
    }

  chimera_thread_exit(ci);


  return 0;
}


auto chimera_threads_run() -> void
{
  /* mutex_input serializes input reading among the CLI workers; it is
     owned here rather than at file scope (the API path does not use it). */
  std::mutex mutex_input;

  /* run the worker pool; each worker processes queries until the input
     is exhausted. chimera_thread_core returns a value that the previous
     pthread_join already discarded, so it is ignored here too. */
  ThreadRunner threadrunner(static_cast<std::size_t>(opt_threads),
                            [&mutex_input](uint64_t nth_thread) {
                              chimera_thread_core(cia + nth_thread, mutex_input);
                            });
  threadrunner.run();
}


auto open_chimera_file(std::FILE ** output_stream, char const * name) -> void
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
      fclose_output(output_stream);
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

  tophits = static_cast<int>(opt_maxaccepts + opt_maxrejects);

  uint64_t progress_total = 0;
  chimera_count = 0;
  nonchimera_count = 0;
  progress = 0;
  seqno = 0;

  /* prepare per-thread chimera detection state */
  std::vector<struct chimera_info_s> cia_v(static_cast<size_t>(opt_threads));
  cia = cia_v.data();

  char const * denovo_dbname = nullptr;

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
          dbindex_prepare(1, static_cast<int>(opt_dbmask));
          dbindex_addallsequences(static_cast<int>(opt_dbmask));
        }

      query_fasta_h = fasta_open(parameters.opt_uchime_ref);
      progress_total = fasta_get_size(query_fasta_h);

      /* The query file is parsed inside the worker threads
         (chimera_thread_core). Defer parse errors so a malformed query
         stops the pool cooperatively instead of calling fatal()/std::exit()
         from a worker while siblings are writing output (CC3); reported
         after the pool joins, below, from the main thread. */
      query_fasta_h->defer_errors = true;
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
      dbindex_prepare(1, static_cast<int>(parameters.opt_qmask));
      progress_total = db_getnucleotidecount();
    }

  if (parameters.opt_log != nullptr)
    {
      if ((parameters.opt_uchime_ref != nullptr) or (parameters.opt_uchime_denovo != nullptr))
        {
          std::fprintf(fp_log, "%8.2f  minh\n", opt_minh);
        }
      auto const is_a_uchime_command = (parameters.opt_uchime_ref != nullptr) or
        (parameters.opt_uchime_denovo != nullptr) or
        (parameters.opt_uchime2_denovo != nullptr) or
        (parameters.opt_uchime3_denovo != nullptr);
      if (is_a_uchime_command)
        {
          std::fprintf(fp_log, "%8.2f  xn\n", opt_xn);
          std::fprintf(fp_log, "%8.2f  dn\n", opt_dn);
          std::fprintf(fp_log, "%8.2f  xa\n", 1.0);
        }

      if ((parameters.opt_uchime_ref != nullptr) or (parameters.opt_uchime_denovo != nullptr))
        {
          std::fprintf(fp_log, "%8.2f  mindiv\n", opt_mindiv);
        }

      std::fprintf(fp_log, "%8.2f  id\n", opt_id);

      if (is_a_uchime_command)
        {
          std::fprintf(fp_log, "%8d  maxp\n", 2);
        }

      std::fprintf(fp_log, "\n");
    }


  progress_init("Detecting chimeras", progress_total);

  chimera_threads_run();

  progress_done();

  /* all workers joined; report a deferred query parse error (CC3, uchime_ref
     only) from the main thread so it does not race a worker's output */
  if ((parameters.opt_uchime_ref != nullptr) and fastx_get_error(query_fasta_h))
    {
      fatal("%s", fastx_get_errmsg(query_fasta_h));
    }

  if (not parameters.opt_quiet)
    {
      if (total_count > 0)
        {
          if (parameters.opt_chimeras_denovo != nullptr)
            {
              std::fprintf(stderr,
                      "Found %d (%.1f%%) chimeras and "
                      "%d (%.1f%%) non-chimeras "
                      "in %d unique sequences.\n",
                      chimera_count,
                      100.0 * chimera_count / total_count,
                      nonchimera_count,
                      100.0 * nonchimera_count / total_count,
                      total_count);
            }
          else
            {
              std::fprintf(stderr,
                      "Found %d (%.1f%%) chimeras, "
                      "%d (%.1f%%) non-chimeras,\n"
                      "and %d (%.1f%%) borderline sequences "
                      "in %d unique sequences.\n",
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
              std::fprintf(stderr,
                      "Found %d chimeras and "
                      "%d non-chimeras "
                      "in %d unique sequences.\n",
                      chimera_count,
                      nonchimera_count,
                      total_count);
            }
          else
            {
              std::fprintf(stderr,
                      "Found %d chimeras, "
                      "%d non-chimeras,\n"
                      "and %d borderline sequences "
                      "in %d unique sequences.\n",
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
              std::fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " (%.1f%%) chimeras and "
                      "%" PRId64 " (%.1f%%) non-chimeras "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      100.0 * static_cast<double>(chimera_abundance) / static_cast<double>(total_abundance),
                      nonchimera_abundance,
                      100.0 * static_cast<double>(nonchimera_abundance) / static_cast<double>(total_abundance),
                      total_abundance);
            }
          else
            {
              std::fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " (%.1f%%) chimeras, "
                      "%" PRId64 " (%.1f%%) non-chimeras,\n"
                      "and %" PRId64 " (%.1f%%) borderline sequences "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      100.0 * static_cast<double>(chimera_abundance) / static_cast<double>(total_abundance),
                      nonchimera_abundance,
                      100.0 * static_cast<double>(nonchimera_abundance) / static_cast<double>(total_abundance),
                      borderline_abundance,
                      100.0 * static_cast<double>(borderline_abundance) / static_cast<double>(total_abundance),
                      total_abundance);
            }
        }
      else
        {
          if (parameters.opt_chimeras_denovo != nullptr)
            {
              std::fprintf(stderr,
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
              std::fprintf(stderr,
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
          std::fprintf(fp_log, "%s", parameters.opt_uchime_ref);
        }
      else
        {
          std::fprintf(fp_log, "%s", denovo_dbname);
        }

      if (seqno > 0)
        {
          std::fprintf(fp_log, ": %d/%u chimeras (%.1f%%)\n",
                  chimera_count,
                  seqno,
                  100.0 * chimera_count / seqno);
        }
      else
        {
          std::fprintf(fp_log, ": %d/%u chimeras\n",
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

  close_chimera_file(fp_borderline);
  close_chimera_file(fp_uchimeout);
  close_chimera_file(fp_uchimealns);
  close_chimera_file(fp_nonchimeras);
  close_chimera_file(fp_chimeras);

  show_rusage();
}


/* === Library API implementation === */

auto chimera_info_alloc() -> struct chimera_info_s *
{
  return new chimera_info_s {};
}

auto chimera_info_free(struct chimera_info_s * ci) -> void
{
  if (ci != nullptr)
    {
      delete ci;
    }
}

auto chimera_session_init() -> void
{
  /* Session-level initialization for chimera detection.
     Sets search-shaping globals that chimera() normally sets but
     the library API path bypasses, and initializes static mutexes
     used by eval_parents/eval_parents_long.

     Call once after DB is loaded and indexed, before creating
     per-thread chimera_info_s handles.

     Assumes global opt_* variables are already set (at minimum:
     opt_match, opt_mismatch, opt_gap_*, opt_minh, opt_xn, opt_dn,
     opt_mindiv, opt_mindiffs, opt_abskew, opt_dbmask, opt_qmask,
     opt_wordlength). */

  /* Override search parameters to chimera detection defaults.
     These must match what chimera() sets in the CLI path (lines 2311-2329).
     opt_weak_id defaults to 10.0 (sentinel) — clamp to opt_id so the
     acceptance check in search_acceptable_aligned doesn't reject everything.
     The CLI does this in command dispatch after option parsing. */
  opt_maxaccepts = few;
  opt_maxrejects = rejects;
  opt_id = chimera_id;
  if (opt_weak_id > opt_id)
    {
      opt_weak_id = opt_id;
    }
  tophits = static_cast<int>(opt_maxaccepts + opt_maxrejects);

  /* For denovo mode, set opt_self/opt_selfid so sequences don't match
     themselves as candidate parents, and set opt_maxsizeratio for
     abundance skew filtering. Matches chimera() CLI setup. */
  if (opt_uchime_ref == nullptr)
    {
      opt_self = 1;
      opt_selfid = 1;
      opt_maxsizeratio = 1.0 / opt_abskew;
    }

  /* eval_parents and eval_parents_long acquire mutex_output for file
     output. The API path does not write files, but the lock is still
     taken; mutex_output is a std::mutex and needs no initialization. */
}


auto chimera_session_cleanup() -> void
{
  /* nothing to release: mutex_output is a std::mutex with automatic
     lifetime. Kept as a stable API symbol. */
}


auto chimera_detect_thread_init(struct chimera_info_s * ci) -> void
{
  /* Per-thread initialization: SIMD aligners, k-mer finders, working
     buffers. Safe to call concurrently for different ci instances.
     Requires chimera_session_init() to have been called first. */
  assert(tophits > 0
         && "chimera_session_init() must be called before chimera_detect_thread_init()");

  chimera_thread_init(ci);

  /* Allocate per-thread working state for chimera_process_query.
     These mirror the locals in chimera_thread_core but persist
     across calls to chimera_detect_single. */
  ci->api_allhits_list.resize(maxcandidates);

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
  ci->api_lma_ptr.reset(new LinearMemoryAligner(scoring));
}


auto chimera_detect_init(struct chimera_info_s * ci) -> void
{
  /* Convenience wrapper: session init + per-thread init in one call.
     Use for single-threaded detection (one chimera_info_s per session).
     For multi-threaded detection, call chimera_session_init() once then
     chimera_detect_thread_init() per thread. */
  chimera_session_init();
  chimera_detect_thread_init(ci);
}

auto chimera_detect_single(struct chimera_info_s * ci,
                           const char * query_seq,
                           const char * query_head,
                           int query_len,
                           int64_t query_size,
                           struct chimera_result_s * result) -> int
{
  /* Populate query in the chimera_info_s.
     ci is per-thread state — must NOT be shared across threads. */
  ci->query_no = 0;
  ci->query_head_len = static_cast<int>(std::strlen(query_head));
  ci->query_len = query_len;
  ci->query_size = query_size;

  realloc_arrays(ci);

  std::strcpy(ci->query_head.data(), query_head);
  std::strcpy(ci->query_seq.data(), query_seq);

  /* Clear result. Non-chimeric results will have only query_label and
     flag='N' populated; all other fields remain zero. */
  *result = {};
  ci->result_out = result;

  /* Use the SAME processing code as the CLI path */
  auto const status = chimera_process_query(ci, ci->api_allhits_list,
                                            *ci->api_lma_ptr);

  if (status == Status::no_parents)
    {
      /* Populate result for no-parents case */
      std::snprintf(result->query_label, sizeof(result->query_label), "%.*s",
                    ci->query_head_len, ci->query_head.data());
      result->flag = 'N';
    }

  /* Free CIGAR strings from this detection */
  for (auto i = 0; i < ci->cand_count; ++i)
    {
      if (ci->nwcigar[static_cast<size_t>(i)] != nullptr)
        {
          xfree(ci->nwcigar[static_cast<size_t>(i)]);
          ci->nwcigar[static_cast<size_t>(i)] = nullptr;
        }
    }

  ci->result_out = nullptr;
  return 0;
}

auto chimera_detect_thread_cleanup(struct chimera_info_s * ci) -> void
{
  /* Per-thread cleanup: frees all resources allocated by
     chimera_detect_thread_init (SIMD aligners, unique k-mer finders,
     minheaps, CIGAR strings, linear memory aligner). */
  for (auto & p : ci->nwcigar)
    {
      if (p != nullptr)
        {
          xfree(p);
          p = nullptr;
        }
    }
  chimera_thread_exit(ci);

  /* Release API working state */
  ci->api_lma_ptr.reset();
  ci->api_allhits_list.clear();
  ci->api_allhits_list.shrink_to_fit();
}


auto chimera_detect_cleanup(struct chimera_info_s * ci) -> void
{
  /* Convenience wrapper: per-thread cleanup + session cleanup in one call.
     Use for single-threaded detection (one chimera_info_s per session).
     For multi-threaded detection, call chimera_detect_thread_cleanup()
     per thread, then chimera_session_cleanup() once. */
  chimera_detect_thread_cleanup(ci);
  chimera_session_cleanup();
}


/* === Batch chimera detection API === */


struct chimera_batch_context_s {
  const char ** query_seqs;
  const char ** query_heads;
  const int * query_lens;
  const int64_t * query_sizes;
  int query_count;
  struct chimera_result_s * results;

  /* per-thread chimera state arrays (sized to opt_threads) */
  struct chimera_info_s ** ci_array;

  /* work-stealing counter */
  std::mutex mutex;
  int next_query;
};


static auto chimera_batch_worker_fn(struct chimera_batch_context_s & ctx,
                                    uint64_t tid) -> void
{
  struct chimera_info_s * ci = ctx.ci_array[tid];

  while (true)
    {
      int qi {0};
      {
        std::lock_guard<std::mutex> const lock(ctx.mutex);
        qi = ctx.next_query++;
      }

      if (qi >= ctx.query_count)
        {
          break;
        }

      chimera_detect_single(ci,
                            ctx.query_seqs[qi],
                            ctx.query_heads[qi],
                            ctx.query_lens[qi],
                            ctx.query_sizes[qi],
                            &ctx.results[qi]);
    }
}


auto chimera_detect_batch(const char ** query_seqs,
                          const char ** query_heads,
                          const int * query_lens,
                          const int64_t * query_sizes,
                          int query_count,
                          struct chimera_result_s * results) -> void
{
  if (query_count <= 0)
    {
      return;
    }

  int const nthreads = std::max(1, static_cast<int>(opt_threads));

  /* Save globals that chimera_session_init will clobber */
  int64_t const saved_maxaccepts = opt_maxaccepts;
  int64_t const saved_maxrejects = opt_maxrejects;
  double const saved_id = opt_id;
  double const saved_weak_id = opt_weak_id;
  bool const saved_self = opt_self;
  bool const saved_selfid = opt_selfid;
  double const saved_maxsizeratio = opt_maxsizeratio;

  /* Session-level init (sets search-shaping globals) */
  chimera_session_init();

  /* Allocate per-thread chimera state */
  struct chimera_batch_context_s ctx;
  ctx.query_seqs = query_seqs;
  ctx.query_heads = query_heads;
  ctx.query_lens = query_lens;
  ctx.query_sizes = query_sizes;
  ctx.query_count = query_count;
  ctx.results = results;
  ctx.next_query = 0;

  ctx.ci_array = static_cast<struct chimera_info_s **>(
    xmalloc(static_cast<size_t>(nthreads) * sizeof(struct chimera_info_s *)));

  for (int t = 0; t < nthreads; t++)
    {
      ctx.ci_array[t] = chimera_info_alloc();
      chimera_detect_thread_init(ctx.ci_array[t]);
    }

  /* run all queries through the worker pool (work-stealing on next_query) */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(nthreads),
                              [&ctx](uint64_t tid) {
                                chimera_batch_worker_fn(ctx, tid);
                              });
    threadrunner.run();
  }

  /* Cleanup per-thread state */
  for (int t = 0; t < nthreads; t++)
    {
      chimera_detect_thread_cleanup(ctx.ci_array[t]);
      chimera_info_free(ctx.ci_array[t]);
    }
  xfree(ctx.ci_array);

  /* Session-level cleanup */
  chimera_session_cleanup();

  /* Restore globals that chimera_session_init clobbered */
  opt_maxaccepts = saved_maxaccepts;
  opt_maxrejects = saved_maxrejects;
  opt_id = saved_id;
  opt_weak_id = saved_weak_id;
  opt_self = saved_self;
  opt_selfid = saved_selfid;
  opt_maxsizeratio = saved_maxsizeratio;
}
