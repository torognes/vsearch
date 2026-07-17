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
#include "arch/increment_counters.hpp"  // increment_counters_from_bitmap*
#include "core/align_simd.hpp"
#include "core/db.hpp"  // Database
#include "core/dbindex.hpp"
#include "core/linmemalign.hpp"
#include "core/minheap.hpp"
#include "core/otutable.hpp"
#include "core/searchcore.hpp"  // struct hit, struct searchinfo_s
#include "core/unique.hpp"
#include "os/system.hpp"  // xfree
#include "utils/make_unique.hpp"  // make_unique
#include "utils/seqcmp.hpp"
#include "utils/span.hpp"
#include "utils/string_alloc.hpp"
#include <algorithm>  // std::count_if, std::min, std::max
#include <array>
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::sscanf, std::size_t
#include <cstdlib>  // std::qsort
#include <cstring>  // std::strlen, std::memset, std::strcmp
#include <limits>
#include <vector>


// Deleters for the opaque per-query handles owned by searchinfo_s via
// unique_ptr (declared in searchcore.hpp). unique_ptr only invokes these on a
// non-null pointer, and the *_exit functions free already-allocated buffers, so
// they never fatal() — safe to run during unwinding (noexcept).
auto uhandle_deleter::operator()(uhandle_s * handle) const noexcept -> void { unique_exit(handle); }
auto s16info_deleter::operator()(s16info_s * handle) const noexcept -> void { search16_exit(handle); }
auto minheap_deleter::operator()(minheap_s * handle) const noexcept -> void { minheap_exit(handle); }


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto make_hits_span(struct searchinfo_s const * search_info) -> Span<struct hit> {
    assert(search_info != nullptr);
    assert(search_info->hit_count >= 0);
    auto const length = static_cast<std::size_t>(search_info->hit_count);
    return Span<struct hit>{search_info->hits, length};
  }


  auto count_number_of_hits_to_keep(struct searchinfo_s const * search_info) -> std::size_t {
    if (search_info == nullptr) { return std::size_t{0}; }
    auto const hits = make_hits_span(search_info);
    return static_cast<std::size_t>(std::count_if(hits.cbegin(), hits.cend(),
                                                  [](struct hit const & hit) -> bool {
                                                    return hit.accepted or hit.weak;
                                                  }));
  }


  auto copy_over_hits_to_be_kept(std::vector<struct hit> & hits,
                                 struct searchinfo_s const * search_info) -> void {
    if (search_info == nullptr) { return; }
    for (auto const & hit : make_hits_span(search_info)) {
      if (hit.accepted or hit.weak) {
        hits.emplace_back(hit);
      }
    }
  }


  auto free_rejected_alignments(struct searchinfo_s const * search_info) -> void {
    if (search_info == nullptr) { return; }
    for (auto & hit : make_hits_span(search_info)) {
      if (not (hit.accepted or hit.weak) and hit.aligned) {
        hit.nwalignment.clear();  // std::string; drop the rejected alignment
      }
    }
  }

}  // end of anonymous namespace



/* per thread data */

inline auto hit_compare_byid_typed(struct hit const * lhs, struct hit const * rhs) -> int
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


inline auto hit_compare_bysize_typed(struct hit const * lhs, struct hit const * rhs, struct Database const & db) -> int
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

  auto const lhs_abundance = db.getabundance(static_cast<uint64_t>(lhs->target));
  auto const rhs_abundance = db.getabundance(static_cast<uint64_t>(rhs->target));
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
  return hit_compare_byid_typed(static_cast<struct hit const *>(lhs), static_cast<struct hit const *>(rhs));
}


auto search_enough_kmers(struct searchinfo_s const & searchinfo,
                         unsigned int const count) -> bool
{
  struct Parameters const & parameters = *searchinfo.parameters;
  return (count >= parameters.opt_minwordmatches) or (count >= searchinfo.kmersamplecount);
}


auto search_topscores(struct searchinfo_s * searchinfo) -> void
{
  struct Parameters const & parameters = *searchinfo->parameters;
  /*
    Count the kmer hits in each database sequence and
    make a sorted list of a given number (th)
    of the database sequences with the highest number of matching kmers.
    These are stored in the min heap array.
  */

  /* count kmer hits in the database sequences */
  unsigned int const indexed_count = searchinfo->dbindex->getcount();

  /* zero counts */
  std::memset(searchinfo->kmers, 0, indexed_count * sizeof(count_t));

  minheap_clear(searchinfo->m.get());

  for (auto i = 0U; i < searchinfo->kmersamplecount; i++)
    {
      auto const kmer = searchinfo->kmersample[i];
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
              /* Saturate at INT16_MAX (32767) rather than letting the
                 unsigned-short counter wrap at 65536. The SIMD bitmap path
                 (increment_counters_from_bitmap*) increments these counters
                 with signed saturation and so caps at 32767; matching that
                 here keeps every counter in [0, 32767], where the two paths
                 agree and neither can wrap a high-overlap target's count back
                 to ~0 and silently drop it from the candidate set (the cap is
                 far above any realistic minwordmatches). */
              count_t & counter = searchinfo->kmers[list[j]];
              if (counter < INT16_MAX) { ++counter; }
            }
        }
    }

  auto const minmatches = std::min(static_cast<unsigned int>(parameters.opt_minwordmatches), searchinfo->kmersamplecount);

  for (auto i = 0U; i < indexed_count; i++)
    {
      auto const count = searchinfo->kmers[i];
      if (count >= minmatches)
        {
          auto const seqno = searchinfo->dbindex->getmapping(i);
          unsigned int const length = static_cast<unsigned int>(searchinfo->db->getsequencelen(seqno));

          elem_t novel;
          novel.count = count;
          novel.seqno = seqno;
          novel.length = length;

          minheap_add(searchinfo->m.get(), & novel);
        }
    }

  minheap_sort(searchinfo->m.get());
}


auto align_trim(struct hit * hit, struct Parameters const & parameters) -> void
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

  auto const * const cigar = hit->nwalignment.c_str();
  auto const * p = cigar;
  auto op = '\0';
  int64_t run = 0;
  if (*p != 0)
    {
      run = 1;
      auto scanlength = 0;
      std::sscanf(p, "%" PRId64 "%n", &run, &scanlength);
      op = *(p + scanlength);
      if (op != 'M')
        {
          hit->trim_aln_left = 1 + scanlength;
          if (op == 'D')
            {
              hit->trim_q_left = static_cast<int>(run);
            }
          else
            {
              hit->trim_t_left = static_cast<int>(run);
            }
        }
    }

  /* right trim alignment */

  auto const * e = cigar + hit->nwalignment.size();
  if (e > cigar)
    {
      p = e - 1;
      op = *p;
      if (op != 'M')
        {
          while ((p > cigar) and (*(p - 1) <= '9'))
            {
              --p;
            }
          run = 1;
          std::sscanf(p, "%" PRId64, &run);
          hit->trim_aln_right = static_cast<int>(e - p);
          if (op == 'D')
            {
              hit->trim_q_right = static_cast<int>(run);
            }
          else
            {
              hit->trim_t_right = static_cast<int>(run);
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

  switch (parameters.opt_iddef)
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


namespace {
  // 128-bit unsigned helper type; __extension__ silences -Wpedantic in C++11.
  __extension__ typedef unsigned __int128 uint128_t;

  /* Compare an abundance against `ratio * reference` without the precision
     loss that converting 64-bit abundances to double incurs above 2^53.
     `ratio` is a non-negative size-ratio threshold (e.g. 1/abskew, or a user
     --minsizeratio / --maxsizeratio); `value` and `reference` are non-negative
     abundances.  Returns the sign of (value - ratio*reference): negative if
     value < ratio*reference, zero if equal, positive if value > ratio*reference.
     The double `ratio` is used at its exact stored (dyadic) value, so results
     match the old double comparison wherever that was exact and stay exact for
     abundances beyond 2^53. */
  auto abundance_ratio_cmp(int64_t const value,
                           double const ratio,
                           int64_t const reference) -> int
  {
    if ((reference <= 0) or (ratio <= 0.0))
      {
        return (value > 0) ? 1 : 0;   // ratio * reference == 0
      }
    if (not std::isfinite(ratio))
      {
        return -1;                    // ratio * reference == +inf > value
      }

    /* Below 2^53 both abundances convert to double exactly, so the historical
       double comparison is itself exact; keep it on that path so long-standing
       boundary behaviour for non-dyadic ratios (e.g. 1/abskew when abskew is
       9) is preserved.  Only switch to exact 128-bit arithmetic when an
       abundance exceeds 2^53, where double would silently drop integer bits. */
    constexpr int64_t exact_double_limit = (int64_t{1} << 53);
    if ((value < exact_double_limit) and (reference < exact_double_limit))
      {
        double const product = ratio * static_cast<double>(reference);
        double const value_d = static_cast<double>(value);
        if (value_d < product) { return -1; }
        if (value_d > product) { return 1; }
        return 0;
      }

    /* Exact decomposition of the stored double: ratio = mantissa * 2^exponent,
       with mantissa an integer in [2^52, 2^53). */
    int exponent = 0;
    auto const mantissa =
      static_cast<int64_t>(std::ldexp(std::frexp(ratio, &exponent), 53));
    exponent -= 53;                   // ratio = mantissa * 2^exponent

    /* Compare value against mantissa * 2^exponent * reference in 128 bits,
       applying the power-of-two shift one bit at a time with an overflow guard
       so neither side can leave the 128-bit range (the guard also terminates
       early for extreme ratios such as the maximum-double default). */
    auto lhs = static_cast<uint128_t>(static_cast<uint64_t>(value));
    auto rhs = static_cast<uint128_t>(static_cast<uint64_t>(mantissa))
             * static_cast<uint128_t>(static_cast<uint64_t>(reference));

    for (int shift = exponent; shift > 0; --shift)
      {
        if ((rhs >> 126) != 0) { return -1; }  // ratio*reference >> value
        rhs <<= 1;
      }
    for (int shift = exponent; shift < 0; ++shift)
      {
        if ((lhs >> 126) != 0) { return 1; }   // value >> ratio*reference
        lhs <<= 1;
      }

    if (lhs < rhs) { return -1; }
    if (lhs > rhs) { return 1; }
    return 0;
  }
}  // anonymous namespace


auto search_acceptable_unaligned(struct searchinfo_s const & searchinfo,
                                 int const target) -> bool
{
  struct Parameters const & parameters = *searchinfo.parameters;
  /* opt_maxsizeratio, opt_self and opt_selfid are read through `parameters`: the
     chimera path threads a detection copy that overrides them via si->parameters,
     so the engine sees the detection values without a mutated global (E1). */
  /* consider whether a hit satisfies accepted criteria before alignment */

  // true: needs further consideration
  // false: reject

  auto const target_seqno = static_cast<uint64_t>(target);
  auto const * qseq = searchinfo.qsequence;
  auto const * dlabel = searchinfo.db->getheader(target_seqno);
  auto const * dseq = searchinfo.db->getsequence(target_seqno);
  int64_t const dseqlen = static_cast<int64_t>(searchinfo.db->getsequencelen(target_seqno));
  int64_t const tsize = static_cast<int64_t>(searchinfo.db->getabundance(target_seqno));

  return (
          /* maxqsize */
          (searchinfo.qsize <= parameters.opt_maxqsize)
          and
          /* mintsize */
          (tsize >= parameters.opt_mintsize)
          and
          /* minsizeratio */
          (abundance_ratio_cmp(searchinfo.qsize, parameters.opt_minsizeratio, tsize) >= 0)
          and
          /* maxsizeratio */
          (abundance_ratio_cmp(searchinfo.qsize, parameters.opt_maxsizeratio, tsize) <= 0)
          and
          /* minqt */
          (searchinfo.qseqlen >= parameters.opt_minqt * static_cast<double>(dseqlen))
          and
          /* maxqt */
          (searchinfo.qseqlen <= parameters.opt_maxqt * static_cast<double>(dseqlen))
          and
          /* minsl */
          (searchinfo.qseqlen < dseqlen ?
           searchinfo.qseqlen >= parameters.opt_minsl * static_cast<double>(dseqlen) :
           static_cast<double>(dseqlen) >= parameters.opt_minsl * searchinfo.qseqlen)
          and
          /* maxsl */
          (searchinfo.qseqlen < dseqlen ?
           searchinfo.qseqlen <= parameters.opt_maxsl * static_cast<double>(dseqlen) :
           static_cast<double>(dseqlen) <= parameters.opt_maxsl * searchinfo.qseqlen)
          and
          /* idprefix */
          ((searchinfo.qseqlen >= parameters.opt_idprefix) and
           (dseqlen >= parameters.opt_idprefix) and
           (seqcmp(qseq, dseq, parameters.opt_idprefix) == 0))
          and
          /* idsuffix */
          ((searchinfo.qseqlen >= parameters.opt_idsuffix) and
           (dseqlen >= parameters.opt_idsuffix) and
           (seqcmp(qseq + searchinfo.qseqlen - parameters.opt_idsuffix,
                    dseq + dseqlen - parameters.opt_idsuffix,
                    parameters.opt_idsuffix) == 0))
          and
          /* self */
          ((parameters.opt_self == 0) or (std::strcmp(searchinfo.query_head, dlabel) != 0))
          and
          /* selfid */
          ((parameters.opt_selfid == 0) or
           (searchinfo.qseqlen != dseqlen) or
           (seqcmp(qseq, dseq, searchinfo.qseqlen) != 0))
          );
}


namespace {
  /* Does the alignment use a gap of a class whose '*' (infinite) sentinel is
     set? An infinite gap-open penalty forbids that class outright; an infinite
     gap-extension penalty forbids gaps of that class longer than one. The gap
     classification mirrors LinearMemoryAligner::alignstats: an 'I' cigar op is
     a query gap, a 'D' op a target gap; the first cigar op is a left-terminal
     gap, the last is right-terminal, any other is interior. Callers guard this
     with parameters.opt_gap_penalty_has_infinite to skip the scan on the common
     finite-penalty path. */
  auto alignment_uses_forbidden_gap(char const * cigar,
                                    struct Parameters const & parameters) -> bool
  {
    if (cigar == nullptr) { return false; }
    auto const * cursor = cigar;
    bool first_op = true;
    while (*cursor != '\0')
      {
        int64_t run = 1;
        int scanned = 0;
        std::sscanf(cursor, "%" PRId64 "%n", &run, &scanned);
        cursor += scanned;
        char const operation = *cursor;
        ++cursor;
        if ((operation == 'I') or (operation == 'D'))
          {
            bool const is_query = (operation == 'I');
            bool const is_left = first_op;
            bool const is_right = (*cursor == '\0');
            bool const open_infinite = is_query
              ? (is_left  ? parameters.opt_gap_open_query_left_infinite
                 : is_right ? parameters.opt_gap_open_query_right_infinite
                 : parameters.opt_gap_open_query_interior_infinite)
              : (is_left  ? parameters.opt_gap_open_target_left_infinite
                 : is_right ? parameters.opt_gap_open_target_right_infinite
                 : parameters.opt_gap_open_target_interior_infinite);
            bool const extension_infinite = is_query
              ? (is_left  ? parameters.opt_gap_extension_query_left_infinite
                 : is_right ? parameters.opt_gap_extension_query_right_infinite
                 : parameters.opt_gap_extension_query_interior_infinite)
              : (is_left  ? parameters.opt_gap_extension_target_left_infinite
                 : is_right ? parameters.opt_gap_extension_target_right_infinite
                 : parameters.opt_gap_extension_target_interior_infinite);
            if (open_infinite) { return true; }
            if (extension_infinite and (run > 1)) { return true; }
          }
        first_op = false;
      }
    return false;
  }
}  // anonymous namespace


auto search_acceptable_aligned(struct searchinfo_s const & searchinfo,
                               struct hit * hit) -> bool
{
  struct Parameters const & parameters = *searchinfo.parameters;
  /* opt_weak_id and opt_id are read through `parameters`: the chimera path threads
     a detection copy that overrides them via si->parameters, so the engine sees
     the detection values without a mutated global (E1). */
  if (/* weak_id */
      (hit->id >= 100.0 * parameters.opt_weak_id) and
      /* maxsubs */
      (hit->mismatches <= parameters.opt_maxsubs) and
      /* maxgaps */
      (hit->internal_gaps <= parameters.opt_maxgaps) and
      /* '*' infinite gap penalties forbid a whole gap class (open) or gaps
         longer than one (extension); reject an alignment that used one */
      ((not parameters.opt_gap_penalty_has_infinite) or
       (not alignment_uses_forbidden_gap(hit->nwalignment.c_str(), parameters))) and
      /* mincols */
      (hit->internal_alignmentlength >= parameters.opt_mincols) and
      /* leftjust */
      ((parameters.opt_leftjust == 0) or (hit->trim_q_left +
                           hit->trim_t_left == 0)) and
      /* rightjust */
      ((parameters.opt_rightjust == 0) or (hit->trim_q_right +
                            hit->trim_t_right == 0)) and
      /* query_cov */
      (hit->matches + hit->mismatches >= parameters.opt_query_cov * searchinfo.qseqlen) and
      /* target_cov */
      (hit->matches + hit->mismatches >=
       parameters.opt_target_cov * static_cast<double>(searchinfo.db->getsequencelen(static_cast<uint64_t>(hit->target)))) and
      /* maxid */
      (hit->id <= 100.0 * parameters.opt_maxid) and
      /* mid */
      (100.0 * hit->matches / (hit->matches + hit->mismatches) >= parameters.opt_mid) and
      /* maxdiffs */
      (hit->mismatches + hit->internal_indels <= parameters.opt_maxdiffs))
    {
      if (parameters.opt_cluster_unoise != nullptr)
        {
          const auto mismatches = hit->mismatches;
          auto const skew = 1.0 * static_cast<double>(searchinfo.qsize) / static_cast<double>(searchinfo.db->getabundance(static_cast<uint64_t>(hit->target)));
          auto const beta = 1.0 / std::pow(2, (1.0 * parameters.opt_unoise_alpha * mismatches) + 1);

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

      if (hit->id >= 100.0 * parameters.opt_id)
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
  /* opt_maxaccepts/opt_maxrejects are read through searchinfo->parameters: each
     caller threads a copy carrying its adjustment (search/cluster the seqcount
     clamp, chimera the detection defaults), so no global is mutated (E1). */
  /* compute global alignment */

  std::array<unsigned int, MAXDELAYED> target_list {{}};
  std::array<CELL, MAXDELAYED> nwscore_list {{}};
  std::array<unsigned short, MAXDELAYED> nwalignmentlength_list {{}};
  std::array<unsigned short, MAXDELAYED> nwmatches_list {{}};
  std::array<unsigned short, MAXDELAYED> nwmismatches_list {{}};
  std::array<unsigned short, MAXDELAYED> nwgaps_list {{}};
  std::array<char *, MAXDELAYED> nwcigar_list {{}};

  unsigned int target_count = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++)
    {
      struct hit const * hit = searchinfo->hits + x;
      if (not hit->rejected)
        {
          target_list[target_count++] = static_cast<unsigned int>(hit->target);
        }
    }

  if (target_count != 0)
    {
      search16(searchinfo->s.get(),
               target_count,
               target_list.data(),
               nwscore_list.data(),
               nwalignmentlength_list.data(),
               nwmatches_list.data(),
               nwmismatches_list.data(),
               nwgaps_list.data(),
               nwcigar_list.data(),
               *searchinfo->db);
    }

  unsigned int i = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++)
    {
      /* maxrejects or maxaccepts reached - ignore remaining hits */
      if ((searchinfo->rejects < searchinfo->parameters->opt_maxrejects) and (searchinfo->accepts < searchinfo->parameters->opt_maxaccepts))
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

              int64_t const dseqlen = static_cast<int64_t>(searchinfo->db->getsequencelen(static_cast<uint64_t>(target)));

              if (nwscore == std::numeric_limits<short>::max())
                {
                  /* In case the SIMD aligner cannot align,
                     perform a new alignment with the
                     linear memory aligner */

                  char const * dseq = searchinfo->db->getsequence(static_cast<uint64_t>(target));

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
              if (nwcigar != nullptr)  // search16 may leave the cigar null
                {
                  hit->nwalignment = nwcigar;  // std::string copies the cigar
                  xfree(nwcigar);              // free the owned char* (xstrdup'd or from nwcigar_list[i])
                }
              else
                {
                  hit->nwalignment.clear();
                }
              hit->nwscore = static_cast<int>(nwscore);
              hit->nwdiff = static_cast<int>(nwalignmentlength - nwmatches);
              hit->nwgaps = static_cast<int>(nwgaps);
              hit->nwindels = static_cast<int>(nwalignmentlength - nwmatches - nwmismatches);
              hit->nwalignmentlength = static_cast<int>(nwalignmentlength);
              hit->nwid = 100.0 * static_cast<double>(nwalignmentlength - hit->nwdiff) /
                static_cast<double>(nwalignmentlength);
              hit->matches = static_cast<int>(nwalignmentlength - hit->nwdiff);
              hit->mismatches = hit->nwdiff - hit->nwindels;

              /* trim alignment and compute numbers excluding terminal gaps */
              align_trim(hit, *searchinfo->parameters);

              /* test accept/reject criteria after alignment */
              if (search_acceptable_aligned(*searchinfo, hit))
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


auto search_onequery(struct searchinfo_s * searchinfo, Masking const seqmask) -> void
{
  /* opt_maxaccepts/opt_maxrejects are read through searchinfo->parameters: each
     caller threads a copy carrying its adjustment (search/cluster the seqcount
     clamp, chimera the detection defaults), so no global is mutated (E1). Query
     kmers are extracted at searchinfo->dbindex->wordlength, the effective index width. */
  searchinfo->hit_count = 0;

  search16_qprep(searchinfo->s.get(), searchinfo->qsequence, searchinfo->qseqlen);

  struct Scoring scoring = scoring_from_options(*searchinfo->parameters);


  searchinfo->lma = make_unique<LinearMemoryAligner>(scoring);


  /* extract unique kmer samples from query*/
  unique_count(searchinfo->uh.get(), static_cast<int>(searchinfo->dbindex->wordlength),
               searchinfo->qseqlen, searchinfo->qsequence,
               &searchinfo->kmersamplecount, &searchinfo->kmersample, seqmask);

  /* find database sequences with the most kmer hits */
  search_topscores(searchinfo);

  /* analyse targets with the highest number of kmer hits */
  searchinfo->accepts = 0;
  searchinfo->rejects = 0;
  searchinfo->finalized = 0;

  int delayed = 0;

  while ((searchinfo->finalized + delayed < searchinfo->parameters->opt_maxaccepts + searchinfo->parameters->opt_maxrejects - 1) and
         (searchinfo->rejects < searchinfo->parameters->opt_maxrejects) and
         (searchinfo->accepts < searchinfo->parameters->opt_maxaccepts) and
         (not minheap_isempty(searchinfo->m.get())))
    {
      elem_t const e = minheap_poplast(searchinfo->m.get());

      struct hit * hit = searchinfo->hits + searchinfo->hit_count;

      hit->target = static_cast<int>(e.seqno);
      hit->count = e.count;
      hit->strand = searchinfo->strand;
      hit->rejected = false;
      hit->accepted = false;
      hit->aligned = false;
      hit->weak = false;
      hit->nwalignment.clear();

      /* Test some accept/reject criteria before alignment */
      if (search_acceptable_unaligned(*searchinfo, static_cast<int>(e.seqno)))
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

  searchinfo->lma.reset();  // frees the aligner (also freed by ~searchinfo_s on unwind)
}


auto search_findbest2_byid(struct searchinfo_s const * si_p,
                           struct searchinfo_s const * si_m) -> struct hit *
{
  struct Parameters const & parameters = *si_p->parameters;
  struct hit * best = nullptr;

  for (int i = 0; i < si_p->hit_count; i++)
    {
      if ((best == nullptr) or (hit_compare_byid_typed(si_p->hits + i, best) < 0))
        {
          best = si_p->hits + i;
        }
    }

  if (parameters.opt_strand)
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


auto search_findbest2_bysize(struct searchinfo_s const * si_p,
                             struct searchinfo_s const * si_m) -> struct hit *
{
  struct Parameters const & parameters = *si_p->parameters;
  struct hit * best = nullptr;

  for (int i = 0; i < si_p->hit_count; i++)
    {
      if ((best == nullptr) or (hit_compare_bysize_typed(si_p->hits + i, best, *si_p->db) < 0))
        {
          best = si_p->hits + i;
        }
    }

  if (parameters.opt_strand)
    {
      for (int i = 0; i < si_m->hit_count; i++)
        {
          if ((best == nullptr) or (hit_compare_bysize_typed(si_m->hits + i, best, *si_p->db) < 0))
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


auto search_joinhits(struct searchinfo_s const * si_plus,
                     struct searchinfo_s const * si_minus,
                     std::vector<struct hit> & hits) -> void
{
  /* join and sort accepted and weak hits from both strands */
  /* free the remaining alignments */

  auto const counter = count_number_of_hits_to_keep(si_plus) + count_number_of_hits_to_keep(si_minus);

  /* allocate new array of hits */
  hits.reserve(counter);

  copy_over_hits_to_be_kept(hits, si_plus);
  copy_over_hits_to_be_kept(hits, si_minus);

  free_rejected_alignments(si_plus);
  free_rejected_alignments(si_minus);

  /* last, sort the hits. std::sort (not std::qsort) because struct hit now owns
     a std::string cigar and is no longer trivially copyable — qsort's bitwise
     element moves would corrupt it. */
  std::sort(hits.begin(), hits.end(),
            [](struct hit const & lhs, struct hit const & rhs) -> bool {
              return hit_compare_byid_typed(&lhs, &rhs) < 0;
            });
}
