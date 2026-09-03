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

#pragma once

#include "core/linmemalign.hpp"
#include "core/mask.hpp"  // Masking
#include "core/minheap.hpp"  // Minheap
#include "core/unique.hpp"  // Uniquer
#include "utils/span.hpp"  // Span
#include "utils/view.hpp"  // View
#include <algorithm>  // std::find_if
#include <array>
#include <cassert>
#include <cstddef>  // std::size_t
#include <iterator>  // std::distance
#include <memory>  // std::unique_ptr
#include <string>  // std::string
#include <vector>


struct s16info_s;

/* Deleter so searchinfo_s can own the opaque per-query SIMD aligner handle via
   unique_ptr and free it on unwind (a fatal() thrown in a library session). Its
   operator() is defined in searchcore.cpp, where the matching search16_exit
   function is visible. Stateless, so it adds no size to the unique_ptr. */
struct s16info_deleter { auto operator()(s16info_s * handle) const noexcept -> void; };

/* the number of alignments that can be delayed */
constexpr auto MAXDELAYED = 8U;

/* Default minimum number of word matches for word lengths 3-15 */
constexpr std::array<int, 16> minwordmatches_defaults =
  {{ -1, -1, -1, 18, 17, 16, 15, 14, 12, 11, 10,  9,  8,  7,  5,  3 }};

struct hit
{
  int target;
  int strand;

  /* candidate info */
  unsigned int count;     /* number of unique kmers shared with query */

  bool accepted;          /* is it accepted? */
  bool rejected;          /* is it rejected? */
  bool aligned;           /* has this hit been aligned */
  bool weak;              /* weak hits are aligned with id > weak_id */

  /* info about global alignment, including terminal gaps */

  int nwscore;           /* alignment score */
  int nwdiff;            /* indels and mismatches in global alignment */
  int nwgaps;            /* gaps in global alignment */
  int nwindels;          /* indels in global alignment */
  int nwalignmentlength; /* length of global alignment */
  double nwid;           /* percent identity of global alignment */
  std::string nwalignment; /* alignment string (cigar) of global alignment (owned) */
  int matches;
  int mismatches;

  /* info about alignment excluding terminal gaps */

  int internal_alignmentlength;
  int internal_gaps;
  int internal_indels;
  int trim_q_left;
  int trim_q_right;
  int trim_t_left;
  int trim_t_right;
  int trim_aln_left;
  int trim_aln_right;

  /* more info */

  double id;             /* identity used for ranking */
  double id0;
  double id1;
  double id2;
  double id3;
  double id4;

  int shortest;          /* length of shortest of query and target */
  int longest;           /* length of longest of query and target */
};

/* The number of differences in an alignment: mismatching columns plus
   gapped columns, terminal gaps excluded. The quantity --maxdiffs is
   compared against, and the diffs userfield. Equivalently alnlen - ids. */
inline auto difference_count(struct hit const & hit) noexcept -> int {
  return hit.mismatches + hit.internal_indels;
}

/* Percentage of identity over the columns holding two nucleotides, gaps of
   every kind excluded. The quantity --mid is compared against, and the mid
   userfield. Zero letter pairs yields 0.0 rather than a NaN: no letter pair
   means no identity to measure. */
inline auto letter_pair_identity(struct hit const & hit) noexcept -> double {
  auto const letter_pairs = hit.matches + hit.mismatches;
  return (letter_pairs > 0) ? 100.0 * hit.matches / letter_pairs : 0.0;
}

/* type of kmer hit counter element remember possibility of overflow */
using count_t = unsigned short;

struct searchinfo_s
{
  int query_no = 0;                 /* query number, zero-based */
  int strand = 0;                   /* strand of query being analysed */
  int64_t qsize = 0;                    /* query abundance */
  std::vector<char> query_head_v;  /* owned header storage (the copying paths) */
  View<char> query_head {nullptr, 0};  /* query header: a view into query_head_v, the
                                          database, or a caller-owned buffer */
  std::vector<char> qsequence_v;  /* vector of query sequence chars, grown to
                                       the longest query seen plus
                                       buffer_headroom (the copying paths) */
  Span<char> qsequence;          /* query sequence (length == query length):
                                       a span over qsequence_v, the database, or
                                       a caller-owned buffer */
  Span<char> full_qsequence;     /* the whole query when qsequence holds only a
                                       part of it (the chimera path partitions
                                       each query); empty everywhere else. Read
                                       by the --selfid gate, which compares the
                                       full-length query against the candidate */
  /* the kmers sampled from the query. A view into whatever produced them: the
     Uniquer's own list_ buffer for a plain search, or -- in the sintax
     bootstrap -- a caller-owned subset array that must outlive the search it
     is set for */
  View<unsigned int> kmersample;
  /* one kmer-match counter per indexed database sequence. Sized by the
     per-thread init, which also reserves headroom past the logical end for the
     SIMD counter stores; readers take a Span over it once, rather than caching
     a second pointer that a reallocation could leave stale */
  std::vector<count_t> kmers_v;
  /* the hit buffer, sized once by the per-thread init to the worst case and
     reused across queries; hit_count is the live fill level, so the hits of the
     query at hand are the first hit_count elements and no more */
  std::vector<struct hit> hits_v;
  int hit_count = 0;
  Uniquer uh {};  /* unique kmer finder instance (owned) */
  std::unique_ptr<s16info_s, s16info_deleter> s;   /* SIMD aligner instance (owned) */
  struct nwinfo_s * nw = nullptr;         /* NW aligner instance */
  std::unique_ptr<LinearMemoryAligner> lma;        /* Linear memory aligner instance (owned) */
  int accepts = 0;                  /* number of accepts */
  int rejects = 0;                  /* number of rejects */
  Minheap m;                     /* min heap with the top kmer db seqs (owned) */
  int finalized = 0;
  /* run configuration, set by the per-thread init at each call site (E1
     shared-infra phase). A pointer (default null) so searchinfo_s stays
     default-constructible for the arrays/vectors that hold it; the searchcore
     functions that take a searchinfo_s read parameters->opt_* through it. The
     pointee is the caller's canonical Parameters and must outlive the si. */
  struct Parameters const * parameters = nullptr;
  /* the k-mer index this query is searched against, set by the per-thread init
     at each call site beside parameters (dbindex migration). A pointer (default
     null) so searchinfo_s stays default-constructible; the searchcore functions
     that take a searchinfo_s read the index through it. The pointee is the
     caller's Dbindex and must outlive the si. */
  struct Dbindex const * dbindex = nullptr;
  /* the in-memory sequence database this query is searched against, set by the
     per-thread init at each call site beside dbindex (Database migration). A
     pointer (default null) so searchinfo_s stays default-constructible; the
     searchcore functions that take a searchinfo_s read the sequences through it.
     The pointee is the caller's Database and must outlive the si. */
  struct Database const * db = nullptr;
  /* whether search_acceptable_aligned() applies the UNOISE skew/beta rule on
     top of the ordinary accept criteria. Only --cluster_unoise does; it used
     to be read as "opt_cluster_unoise is non-null", which a library caller
     never sets, so an engine driven through the library read as not-UNOISE no
     matter what it was asked to do. Set by the per-thread init at each call
     site beside parameters; false everywhere else, which is what every other
     search path wants. */
  bool unoise_acceptance = false;
};

/* The hits of the query at hand: the live prefix of the reused hit buffer.
   Everything that reads or edits a query's hits goes through one of these two,
   so that no caller has to pair hits_v with hit_count by hand -- the pairing
   that the raw `struct hit * hits` member used to invite. Use the View form
   unless the hits are to be modified. */
inline auto make_hits_span(struct searchinfo_s * const search_info) -> Span<struct hit> {
  assert(search_info != nullptr);
  assert(search_info->hit_count >= 0);
  auto const length = static_cast<std::size_t>(search_info->hit_count);
  return make_span(search_info->hits_v).first(length);
}

inline auto make_hits_view(struct searchinfo_s const * const search_info) -> View<struct hit> {
  assert(search_info != nullptr);
  assert(search_info->hit_count >= 0);
  auto const length = static_cast<std::size_t>(search_info->hit_count);
  return make_view(search_info->hits_v).first(length);
}

/* The hits to report for this query: all of them, or -- under --top_hits_only
   -- only the leading run tied with the best identity. --top_hits_only does not
   filter the hits, it truncates them at the first one that is not tied with the
   best, so the answer is a prefix and therefore a View.

   Seven writers used to spell this as a break in the middle of their loop.
   find_if stops at the first hit that fails the test, exactly as the break did,
   so the result is unchanged even for an unsorted list. */
inline auto top_hits(View<struct hit> const hits, bool const top_hits_only)
  -> View<struct hit> {
  if ((not top_hits_only) or hits.empty()) { return hits; }
  auto const best_id = hits.front().id;
  auto const * const first_worse =
    std::find_if(hits.begin(), hits.end(),
                 [best_id](struct hit const & candidate) -> bool {
                   return candidate.id < best_id;
                 });
  return hits.first(static_cast<std::size_t>(std::distance(hits.begin(), first_worse)));
}


auto search_topscores(struct searchinfo_s * searchinfo) -> void;

auto search_onequery(struct searchinfo_s * searchinfo, Masking seqmask) -> void;

/* both return a mutable pointer into si_p's or si_m's hit buffer -- the caller
   moves the winning alignment string out of it -- so neither can take a
   pointer-to-const searchinfo_s */
/* si_m is the reverse strand and is a Span, not a pointer, because it is
   absent unless --strand both is in effect: an empty span says so, where the
   null pointer these took before said nothing the callee could check. The
   caller's own opt_strand and the span's emptiness are the same fact, and the
   assert inside says so. */
auto search_findbest2_byid(struct searchinfo_s & si_p,
                           Span<struct searchinfo_s> si_m) -> struct hit *;

auto search_findbest2_bysize(struct searchinfo_s & si_p,
                             Span<struct searchinfo_s> si_m) -> struct hit *;

auto search_acceptable_unaligned(struct searchinfo_s const & searchinfo,
                                 int target) -> bool;

auto search_acceptable_aligned(struct searchinfo_s const & searchinfo,
                               struct hit & hit) -> bool;

auto align_trim(struct hit & hit, struct Parameters const & parameters) -> void;

/* copies the accepted and weak hits of both strands into `hits`, then drops the
   alignment strings of the ones it did not copy -- hence the mutable si */
auto search_joinhits(struct searchinfo_s * si_p,
                     struct searchinfo_s * si_m,
                     std::vector<struct hit> & hits) -> void;

/* The k-mer evidence for one candidate target: how many distinct k-mers it
   shares with the query, and how many it holds at all -- its own ceiling on
   `shared`, and so on the threshold it can be asked to meet. Zero means the
   candidate contributes no k-mer, hence shares none. A struct, because two
   adjacent unsigned ints would be swappable at the call site. */
struct KmerEvidence
{
  unsigned int shared;
  unsigned int target_kmers;
};

auto search_enough_kmers(struct searchinfo_s const & searchinfo,
                         struct KmerEvidence const & evidence) -> bool;
