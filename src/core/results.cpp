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
#include "core/attributes.hpp"
#include "core/db.hpp"  // Database
#include "core/fasta.hpp"  // fasta_print_general
#include "core/searchcore.hpp"  // struct hit, top_hits
#include "core/showalign.hpp"
#include "core/tax.hpp"  // TaxLevel, tax_levels, tax_split
#include "utils/cigar.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/view.hpp"
#include "utils/taxonomic_fields.h"
#include "utils/sequence_digest.hpp"
#include "utils/print_record.hpp"  // OutputRecord, fprint
#include "utils/print_view.hpp"  // fprint
#include "utils/prog_id.hpp"  // PROG_NAME, PROG_VERSION
#include <algorithm>  // std::equal, std::max
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <iterator>  // std::next
#include <string>  // std::string, std::to_string


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto check_if_perfect_match(char const * opt_cluster_fast,
                              struct hit const & hit) -> bool {
    if (opt_cluster_fast != nullptr) {
      /* cluster_fast */
      /* use '=' for identical sequences, ignoring terminal gaps */
      return (hit.matches == hit.internal_alignmentlength);
    }
    /* cluster_size, cluster_smallmem, cluster_unoise */
    /* usearch_global, search_exact, allpairs_global */
    /* use '=' for strictly identical sequences */
    return (hit.matches == hit.nwalignmentlength);
  }

  /* The name of one taxonomic rank, as a window into the header it was parsed
     out of (see TaxLevel in core/tax.hpp). A zero-length rank is one the header
     does not carry, and yields an empty view. */
  auto rank_name(View<char> const header, TaxLevel const & rank) -> View<char> {
    return header.subspan(static_cast<std::size_t>(rank.start),
                          static_cast<std::size_t>(rank.length));
  }

}  // end of anonymous namespace


auto results_show_fastapairs_one(std::FILE * output_handle,
                                 struct hit const & hit,
                                 View<char> const query_head,
                                 View<char> const qsequence,
                                 View<char> const qsequence_rc,
                                 struct Database const & db,
                                 struct Parameters const & parameters) -> void
{
  /* http://www.drive5.com/usearch/manual/fastapairs.html */

  auto const query = (hit.strand != 0) ? qsequence_rc : qsequence;
  auto const qrow = get_alignment_qrow(query,
                                 make_view(hit.nwalignment),
                                 hit.nwalignmentlength);
  fasta_print_general(output_handle,
                      nullptr,
                      View<char>{&qrow[static_cast<std::size_t>(hit.trim_q_left + hit.trim_t_left)],
                                 static_cast<std::size_t>(hit.internal_alignmentlength)},
                      query_head,
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0,
                      0,
                      parameters);

  auto const target = static_cast<uint64_t>(hit.target);
  auto const trow = get_alignment_trow(db.sequence_view(target),
                                 make_view(hit.nwalignment),
                                 hit.nwalignmentlength);
  fasta_print_general(output_handle,
                      nullptr,
                      View<char>{&trow[static_cast<std::size_t>(hit.trim_q_left + hit.trim_t_left)],
                                 static_cast<std::size_t>(hit.internal_alignmentlength)},
                      db.header_view(target),
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0,
                      0,
                      parameters);

  fprint(output_handle, '\n');
}


auto results_show_qsegout_one(std::FILE * output_handle,
                              struct hit const & hit,
                              View<char> const query_head,
                              View<char> const qsequence,
                              View<char> const qsequence_rc,
                              struct Parameters const & parameters) -> void
{
  auto const query = (hit.strand != 0) ? qsequence_rc : qsequence;
  auto const qseglen = static_cast<int64_t>(query.size())
    - hit.trim_q_left - hit.trim_q_right;
  auto const qseg = query.subspan(static_cast<std::size_t>(hit.trim_q_left),
                                  static_cast<std::size_t>(qseglen));

  fasta_print_general(output_handle,
                      nullptr,
                      qseg,
                      query_head,
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0,
                      0,
                      parameters);
}


auto results_show_tsegout_one(std::FILE * output_handle,
                              struct hit const & hit,
                              struct Database const & db,
                              struct Parameters const & parameters) -> void
{
  auto const target = static_cast<uint64_t>(hit.target);
  auto const target_sequence = db.sequence_view(target);
  auto const tseglen = static_cast<int64_t>(target_sequence.size())
    - hit.trim_t_left - hit.trim_t_right;
  auto const tseg = target_sequence.subspan(static_cast<std::size_t>(hit.trim_t_left),
                                            static_cast<std::size_t>(tseglen));

  fasta_print_general(output_handle,
                      nullptr,
                      tseg,
                      db.header_view(target),
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0,
                      0,
                      parameters);
}


auto results_show_blast6out_one(std::FILE * output_handle,
                                struct hit const * hit,
                                View<char> const query_head,
                                int64_t const qseqlen,
                                struct Database const & db) -> void
{

  /*
    http://www.drive5.com/usearch/manual/blast6out.html

    query label
    target label
    percent identity
    alignment length
    number of mismatches
    number of gap opens
    1-based position of start in query
    1-based position of end in query
    1-based position of start in target
    1-based position of end in target
    E-value
    bit score

    Note that USEARCH shows 13 fields when there is no hit,
    but only 12 when there is a hit. Fixed in VSEARCH.
  */

  if (hit == nullptr) {
    fprint(output_handle, query_head);
    fprint(output_handle, "\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n");
    return;
  }
  // if 'hit->strand' then 'minus strand' else 'plus strand'
  auto const target = static_cast<uint64_t>(hit->target);
  int const qstart = (hit->strand != 0) ? static_cast<int>(qseqlen) : 1;
  int const qend = (hit->strand != 0) ? 1 : static_cast<int>(qseqlen);

  OutputRecord record {output_handle};
  fprint(record, query_head);
  fprint(record, '\t');
  fprint(record, db.header_view(target));
  fprint(record, '\t');
  std::fprintf(record.stream(), "%.1f", hit->id);
  fprint(record, '\t');
  fprint_integer(record, hit->internal_alignmentlength);
  fprint(record, '\t');
  fprint_integer(record, hit->mismatches);
  fprint(record, '\t');
  fprint_integer(record, hit->internal_gaps);
  fprint(record, '\t');
  fprint_integer(record, qstart);
  fprint(record, '\t');
  fprint_integer(record, qend);
  fprint(record, '\t');
  fprint_integer(record, 1);
  fprint(record, '\t');
  fprint_integer(record, db.getsequencelen(target));
  fprint(record, '\t');
  fprint_integer(record, -1);
  fprint(record, '\t');
  fprint_integer(record, 0);
  fprint(record, '\n');
}


auto results_show_uc_one(std::FILE * output_handle,
                         struct hit const * hit,
                         View<char> const query_head,
                         int64_t const qseqlen,
                         int const clusterno,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void
{
  /*
    http://www.drive5.com/usearch/manual/ucout.html

    Columns:
    H/N
    cluster no (0-based) (target sequence no)
    sequence length (query)
    percent identity
    strand: + or -
    0
    0
    compressed alignment, e.g. 9I92M14D, or "=" if perfect alignment
    query label
    target label
  */

  if (hit == nullptr) {
    fprint(output_handle, "N\t*\t*\t*\t.\t*\t*\t*\t");
    fprint(output_handle, query_head);
    fprint(output_handle, "\t*\n");
    return;
  }

  auto const is_perfect_match = check_if_perfect_match(parameters.opt_cluster_fast, *hit);

  OutputRecord record {output_handle};
  fprint(record, "H\t");
  fprint_integer(record, clusterno);
  fprint(record, '\t');
  fprint_integer(record, qseqlen);
  fprint(record, '\t');
  std::fprintf(record.stream(), "%.1f", hit->id);
  fprint(record, '\t');
  fprint(record, (hit->strand != 0) ? '-' : '+');
  fprint(record, "\t0\t0\t");
  if (is_perfect_match) {
    fprint(record, "=");
  } else {
    fprint(record, make_view(hit->nwalignment));
  }
  fprint(record, '\t');
  auto const target = static_cast<uint64_t>(hit->target);
  /* the two header fields come from a helper that writes to the stream, so the
     record flushes around them: three fwrites for this line instead of one,
     still against the 24 stdio calls the unbuffered form needs */
  header_fprint_strip(record.stream(),
                      query_head,
                      parameters.opt_xsize,
                      parameters.opt_xee,
                      parameters.opt_xlength);
  fprint(record, '\t');
  header_fprint_strip(record.stream(),
                      db.header_view(target),
                      parameters.opt_xsize,
                      parameters.opt_xee,
                      parameters.opt_xlength);
  fprint(record, '\n');
}


namespace {
/* One --userout field, written as itself. Lifted verbatim out of the loop in
   results_show_userout_one, whose only remaining job is the tab that separates
   the fields -- which is what its index was for.

   The three View<char> parameters are the same three, in the same order, that
   the enclosing writer takes. `parameters` is deliberately not among them: none
   of the cases reads it. */
auto print_userfield(std::FILE * output_handle,
                     int const field,
                     struct hit const * const hit,
                     View<char> const query_head,
                     View<char> const qsequence,
                     View<char> const qsequence_rc,
                     struct Database const & db) -> void
{
  auto const qseqlen = static_cast<int64_t>(qsequence.size());

  View<char> tsequence;
  int64_t tseqlen = 0;
  View<char> t_head;

  if (hit != nullptr)
    {
      auto const target = static_cast<uint64_t>(hit->target);
      tsequence = db.sequence_view(target);
      tseqlen = static_cast<int64_t>(tsequence.size());
      t_head = db.header_view(target);
    }

  switch (field)
    {
    case 0: /* query */
      fprint(output_handle, query_head);
      break;
    case 1: /* target */
      if (hit != nullptr) { fprint(output_handle, t_head); }
      else { fprint(output_handle, '*'); }
      break;
    case 2: /* evalue */
      fprint(output_handle, "-1");
      break;
    case 3: /* id */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? hit->id : 0.0);
      break;
    case 4: /* pctpv */
      std::fprintf(output_handle, "%.1f", ((hit != nullptr) and (hit->internal_alignmentlength > 0)) ? 100.0 * hit->matches / hit->internal_alignmentlength : 0.0);
      break;
    case 5: /* pctgaps */
      std::fprintf(output_handle, "%.1f", ((hit != nullptr) and (hit->internal_alignmentlength > 0)) ? 100.0 * hit->internal_indels / hit->internal_alignmentlength : 0.0);
      break;
    case 6: /* pairs */
      fprint_integer(output_handle, (hit != nullptr) ? hit->matches + hit->mismatches : 0);
      break;
    case 7: /* gaps */
      fprint_integer(output_handle, (hit != nullptr) ? hit->internal_indels : 0);
      break;
    case 8: /* qlo */
      fprint_integer(output_handle, (hit != nullptr) ? ((hit->strand != 0) ? qseqlen : 1) : 0);
      break;
    case 9: /* qhi */
      fprint_integer(output_handle, (hit != nullptr) ? ((hit->strand != 0) ? 1 : qseqlen) : 0);
      break;
    case 10: /* tlo */
      fprint_integer(output_handle, (hit != nullptr) ? 1 : 0);
      break;
    case 11: /* thi */
      fprint_integer(output_handle, tseqlen);
      break;
    case 12: /* pv */
      fprint_integer(output_handle, (hit != nullptr) ? hit->matches : 0);
      break;
    case 13: /* ql */
      fprint_integer(output_handle, qseqlen);
      break;
    case 14: /* tl */
      fprint_integer(output_handle, (hit != nullptr) ? tseqlen : 0);
      break;
    case 15: /* qs */
      fprint_integer(output_handle, qseqlen);
      break;
    case 16: /* ts */
      fprint_integer(output_handle, (hit != nullptr) ? tseqlen : 0);
      break;
    case 17: /* alnlen */
      fprint_integer(output_handle, (hit != nullptr) ? hit->internal_alignmentlength : 0);
      break;
    case 18: /* opens */
      fprint_integer(output_handle, (hit != nullptr) ? hit->internal_gaps : 0);
      break;
    case 19: /* exts */
      fprint_integer(output_handle, (hit != nullptr) ? hit->internal_indels - hit->internal_gaps : 0);
      break;
    case 20: /* raw */
      fprint_integer(output_handle, (hit != nullptr) ? hit->nwscore : 0);
      break;
    case 21: /* bits */
      fprint_integer(output_handle, 0);
      break;
    case 22: /* aln */
      if (hit != nullptr)
        {
          print_uncompressed_cigar(output_handle, make_view(hit->nwalignment));
        }
      break;
    case 23: /* caln */
      if (hit != nullptr)
        {
          fprint(output_handle, make_view(hit->nwalignment));
        }
      break;
    case 24: /* qstrand */
      if (hit != nullptr)
        {
          fprint(output_handle, (hit->strand != 0) ? '-' : '+');
        }
      break;
    case 25: /* tstrand */
      if (hit != nullptr)
        {
          fprint(output_handle, '+');
        }
      break;
    case 26: /* qrow */
      if (hit != nullptr)
        {
          auto const query = (hit->strand != 0) ? qsequence_rc : qsequence;
          auto const qrow = get_alignment_qrow(query,
                                         make_view(hit->nwalignment),
                                         hit->nwalignmentlength);
          fprint(output_handle,
                 View<char>{&qrow[static_cast<std::size_t>(hit->trim_q_left + hit->trim_t_left)],
                            static_cast<std::size_t>(hit->internal_alignmentlength)});
        }
      break;
    case 27: /* trow */
      if (hit != nullptr)
        {
          auto const trow = get_alignment_trow(tsequence,
                                         make_view(hit->nwalignment),
                                         hit->nwalignmentlength);
          fprint(output_handle,
                 View<char>{&trow[static_cast<std::size_t>(hit->trim_q_left + hit->trim_t_left)],
                            static_cast<std::size_t>(hit->internal_alignmentlength)});
        }
      break;
    case 28: /* qframe */
    case 29: /* tframe */
      fprint(output_handle, "+0");
      break;
    case 30: /* mism */
      fprint_integer(output_handle, (hit != nullptr) ? hit->mismatches : 0);
      break;
    case 31: /* ids */
      fprint_integer(output_handle, (hit != nullptr) ? hit->matches : 0);
      break;
    case 32: /* qcov */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? 100.0 * (hit->matches + hit->mismatches) / static_cast<double>(qseqlen) : 0.0);
      break;
    case 33: /* tcov */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? 100.0 * (hit->matches + hit->mismatches) / static_cast<double>(tseqlen) : 0.0);
      break;
    case 34: /* id0 */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? hit->id0 : 0.0);
      break;
    case 35: /* id1 */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? hit->id1 : 0.0);
      break;
    case 36: /* id2 */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? hit->id2 : 0.0);
      break;
    case 37: /* id3 */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? hit->id3 : 0.0);
      break;
    case 38: /* id4 */
      std::fprintf(output_handle, "%.1f", (hit != nullptr) ? hit->id4 : 0.0);
      break;

      /* new internal alignment coordinates */

    case 39: /* qilo */
      fprint_integer(output_handle, (hit != nullptr) ? hit->trim_q_left + 1 : 0);
      break;
    case 40: /* qihi */
      fprint_integer(output_handle, (hit != nullptr) ? qseqlen - hit->trim_q_right : 0);
      break;
    case 41: /* tilo */
      fprint_integer(output_handle, (hit != nullptr) ? hit->trim_t_left + 1 : 0);
      break;
    case 42: /* tihi */
      fprint_integer(output_handle, (hit != nullptr) ? tseqlen - hit->trim_t_right : 0);
      break;
    default:
      /* userfields_requested only ever holds validated indices (0..42),
         so this is unreachable today. It guards against a userfields_names
         entry being added or reordered in utils/userfields.cpp without a matching
         case here — the positional coupling would otherwise print nothing
         silently (E2). */
      fatal("Internal error: unknown userfield index in results_show_userout_one");
    }
}
}  // anonymous namespace


auto results_show_userout_one(std::FILE * output_handle, struct hit const * hit,
                              View<char> const query_head,
                              View<char> const qsequence,
                              View<char> const qsequence_rc,
                              struct Database const & db,
                              struct Parameters const & parameters) -> void
{

  /*
    http://drive5.com/usearch/manual/userout.html
    qlo, qhi, tlo, thi and raw are given more meaningful values here
  */

  auto const & userfields_requested = parameters.opt_userfields;

  /* The requested fields, tab-separated. The first field carries no leading
     tab and the rest do, which is the whole of what the index this loop used to
     carry was for; drop(1) is empty rather than out of range when a single
     field was requested, but front() above it needs the guard. */
  if (not userfields_requested.empty())
    {
      print_userfield(output_handle, userfields_requested.front(),
                      hit, query_head, qsequence, qsequence_rc, db);
    }

  for (auto const field : make_view(userfields_requested).drop(1))
    {
      fprint(output_handle, '\t');
      print_userfield(output_handle, field,
                      hit, query_head, qsequence, qsequence_rc, db);
    }

  fprint(output_handle, '\n');
}


auto results_show_lcaout(std::FILE * output_handle,
                         View<struct hit> const hits,
                         View<char> const query_head,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void
{
  /* Output last common ancestor (LCA) of the hits,
     in a similar way to the Sintax command */

  /* Use a modified Boyer-Moore majority voting algorithm at each taxonomic
     level to find the most common name at each level */

  fprint(output_handle, query_head);
  fprint(output_handle, '\t');

  if (hits.empty()) {
    fprint(output_handle, '\n');
    return;
  }

  constexpr auto levels = static_cast<std::size_t>(tax_levels);
  std::array<int, tax_levels> votes {{}};
  std::array<int, tax_levels> cand;
  cand.fill(-1);
  /* the ranks of the candidate at each level, as offsets into that candidate's
     own header (see TaxLevel in core/tax.hpp); one array of records where two
     parallel start/length arrays used to sit side by side, indexed in lockstep */
  std::array<std::array<TaxLevel, tax_levels>, tax_levels> cand_levels {{}};
  std::array<int, tax_levels> level_match {{}};

  auto const top = top_hits(hits, parameters.opt_top_hits_only != 0);

  for (auto const & hit : top)
    {
      int const seqno = hit.target;
      std::array<TaxLevel, tax_levels> new_levels {{}};
      tax_split(seqno, new_levels, db);
      auto const hit_header = db.header_view(static_cast<uint64_t>(seqno));

      for (std::size_t k = 0; k < levels; ++k)
        {
          if (votes[k] == 0)
            {
              cand[k] = seqno;
              votes[k] = 1;
              cand_levels[k] = new_levels;
            }
          else
            {
              /* Does this hit agree with the candidate at level k on every rank
                 down to k? The two sets of offsets index two different headers,
                 so both headers are part of the predicate -- and the two
                 arguments are interchangeable, equality being symmetric. */
              auto const cand_header = db.header_view(static_cast<uint64_t>(cand[k]));
              auto const match =
                std::equal(cand_levels[k].begin(),
                           std::next(cand_levels[k].begin(), static_cast<std::ptrdiff_t>(k) + 1),
                           new_levels.begin(),
                           [&](TaxLevel const & cand_rank, TaxLevel const & hit_rank) -> bool {
                             return rank_name(cand_header, cand_rank)
                               == rank_name(hit_header, hit_rank);
                           });
              if (match)
                {
                  ++votes[k];
                }
              else
                {
                  --votes[k];
                }
            }
        }
    }

  /* count actual matches to the candidate at each level */

  for (auto const & hit : top)
    {
      auto const seqno = hit.target;
      std::array<TaxLevel, tax_levels> new_levels {{}};
      tax_split(seqno, new_levels, db);
      auto const hit_header = db.header_view(static_cast<uint64_t>(seqno));

      for (std::size_t k = 0; k < levels; ++k)
        {
          auto const cand_header = db.header_view(static_cast<uint64_t>(cand[k]));
          auto const match =
            std::equal(cand_levels[k].begin(),
                       std::next(cand_levels[k].begin(), static_cast<std::ptrdiff_t>(k) + 1),
                       new_levels.begin(),
                       [&](TaxLevel const & cand_rank, TaxLevel const & hit_rank) -> bool {
                         return rank_name(cand_header, cand_rank)
                           == rank_name(hit_header, hit_rank);
                       });
          if (match)
            {
              ++level_match[k];
            }
        }
    }

  /* output results */

  if (top.empty()) {
    fprint(output_handle, '\n');
    return;
  }
  auto comma = false;
  for (std::size_t j = 0; j < levels; ++j)
    {
      if (1.0 * level_match[j] / static_cast<double>(top.size()) < parameters.opt_lca_cutoff)
        {
          break;
        }

      if (cand_levels[j][j].length > 0)
        {
          std::fputs((comma ? "," : ""), output_handle);
          fprint(output_handle, static_cast<char>(taxonomic_fields[j]));
          fprint(output_handle, ':');
          fprint(output_handle,
                 rank_name(db.header_view(static_cast<uint64_t>(cand[j])),
                           cand_levels[j][j]));
          comma = true;
        }
    }

  fprint(output_handle, '\n');
}


auto results_show_alnout(std::FILE * output_handle,
                         View<struct hit> const hits,
                         View<char> const query_head,
                         View<char> const qsequence,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void
{
  /* http://drive5.com/usearch/manual/alnout.html */


  if (hits.empty()) {
    if (parameters.opt_output_no_hits != 0) {
      fprint(output_handle, '\n');
      fprint(output_handle, "Query >");
      fprint(output_handle, query_head);
      fprint(output_handle, "\nNo hits\n");
    }
    return;
  }

  auto const qseqlen = static_cast<int64_t>(qsequence.size());

  fprint(output_handle, '\n');

  fprint(output_handle, "Query >");
  fprint(output_handle, query_head);
  fprint(output_handle, "\n %Id   TLen  Target\n");

  auto const top = top_hits(hits, parameters.opt_top_hits_only != 0);

  for (auto const & hit : top)
    {
      auto const target = static_cast<uint64_t>(hit.target);
      OutputRecord record {output_handle};
      std::fprintf(record.stream(), "%3.0f", hit.id);
      fprint(record, "% ");
      fprint_integer(record, db.getsequencelen(target), 6);
      fprint(record, "  ");
      fprint(record, db.header_view(target));
      fprint(record, '\n');
    }

  for (auto const & hit : top)
    {
      fprint(output_handle, '\n');


      auto const target = static_cast<uint64_t>(hit.target);
      auto const dseqlen = static_cast<int64_t>(db.getsequencelen(target));

      /* the width of the wider of the two lengths, in decimal digits. The
         digits view carries it, so no std::string is built to be measured and
         thrown away -- and on libstdc++ <= 10 std::to_string is itself a
         vsnprintf call with a format string. One buffer serves both: the first
         view is dead before the second is taken. */
      decimal::Buffer length_digits {};
      auto const qlenlen = decimal::to_decimal(length_digits, qseqlen).size();
      auto const tlenlen = decimal::to_decimal(length_digits, dseqlen).size();
      auto const numwidth = static_cast<int>(std::max(qlenlen, tlenlen));

      {
        OutputRecord record {output_handle};
        fprint(record, " Query ");
        fprint_integer(record, qseqlen, static_cast<std::size_t>(numwidth));
        fprint(record, "nt >");
        fprint(record, query_head);
        fprint(record, "\nTarget ");
        fprint_integer(record, dseqlen, static_cast<std::size_t>(numwidth));
        fprint(record, "nt >");
        fprint(record, db.header_view(target));
        fprint(record, '\n');
      }

      int64_t const rowlen = (parameters.opt_rowlen == 0) ? (qseqlen + dseqlen) : parameters.opt_rowlen;

      align_show(output_handle,
                 qsequence,
                 hit.trim_q_left,
                 "Qry",
                 db.sequence_view(target),
                 hit.trim_t_left,
                 "Tgt",
                 View<char>{std::next(hit.nwalignment.c_str(), hit.trim_aln_left),
                            hit.nwalignment.size()
                            - static_cast<std::size_t>(hit.trim_aln_left)
                            - static_cast<std::size_t>(hit.trim_aln_right)},
                 numwidth,
                 3,
                 rowlen,
                 hit.strand,
                 parameters);

      OutputRecord record {output_handle};
      fprint(record, '\n');
      fprint_integer(record, hit.internal_alignmentlength);
      fprint(record, " cols, ");
      fprint_integer(record, hit.matches);
      fprint(record, " ids (");
      std::fprintf(record.stream(), "%3.1f", hit.id);
      fprint(record, "%), ");
      fprint_integer(record, hit.internal_indels);
      fprint(record, " gaps (");
      std::fprintf(record.stream(), "%3.1f", hit.internal_alignmentlength > 0 ?
              100.0 * hit.internal_indels / hit.internal_alignmentlength :
              0.0);
      fprint(record, "%)\n");
    }
}


namespace {
auto build_sam_strings(View<char> const alignment,
                       View<char> const queryseq,
                       View<char> const targetseq,
                       std::string & cigar,
                       std::string & md) -> void
{
  /*
    convert cigar to sam format:
    add "1" to operations without run length
    flip direction of indels in cigar string

    build MD-string with substitutions
  */

  cigar.clear();
  md.clear();

  auto const queryseqlen = static_cast<int64_t>(queryseq.size());
  auto const targetseqlen = static_cast<int64_t>(targetseq.size());

  /* the alignment view is over a std::string, so the NUL that read_runlength()
     needs below sits just past its end */
  auto const * p = alignment.begin();
  auto const * const e = alignment.end();

  /* indices into the two sequence views, hence std::size_t; the bound checks
     below promote them to int64_t so that qpos + run cannot overflow */
  std::size_t qpos = 0;
  std::size_t tpos = 0;

  auto matched = 0;
  auto flag = false; /* 1: MD string ends with a number */

  while (p < e)
    {
      /* read_runlength() rather than find_runlength_of_leftmost_operation():
         the latter clamps to 1, which would hide the malformed run lengths the
         guard below reports */
      char const * next_operation = nullptr;
      auto const scanned_run = read_runlength(p, next_operation);
      /* no run-length number: the cigar convention is an implicit run of 1 */
      auto const run = (next_operation == p) ? 1LL : scanned_run;
      p = next_operation;
      auto const op = *p;
      p = std::next(p);

      /*
        Guard against a CIGAR whose run-lengths do not sum to the query
        and target lengths: walking qpos/tpos past the sequence ends would
        read out of bounds. Well-formed input from the in-tree aligner
        never trips this; a corrupted or externally supplied CIGAR could.
        'M' advances both positions, 'D' only the query, 'I' only the
        target, so each op is bounded against exactly what it reads.
      */
      if (run < 0)
        {
          fatal("Invalid CIGAR string: negative run length");
        }

      switch (op)
        {
        case 'M':
          if ((static_cast<int64_t>(qpos) + run > queryseqlen) or
              (static_cast<int64_t>(tpos) + run > targetseqlen))
            {
              fatal("Invalid CIGAR string: run length exceeds sequence bounds");
            }
          cigar += std::to_string(run);
          cigar += 'M';

          for (auto i = 0LL; i < run; ++i)
            {
              if (is_same_4bit(queryseq[qpos], targetseq[tpos]))
                {
                  ++matched;
                }
              else
                {
                  if (not flag)
                    {
                      md += std::to_string(matched);
                      matched = 0;
                      flag = true;
                    }

                  md += targetseq[tpos];
                  flag = false;
                }
              ++qpos;
              ++tpos;
            }

          break;

        case 'D':
          if (static_cast<int64_t>(qpos) + run > queryseqlen)
            {
              fatal("Invalid CIGAR string: run length exceeds sequence bounds");
            }
          cigar += std::to_string(run);
          cigar += 'I';
          qpos += static_cast<std::size_t>(run);
          break;

        case 'I':
          if (static_cast<int64_t>(tpos) + run > targetseqlen)
            {
              fatal("Invalid CIGAR string: run length exceeds sequence bounds");
            }
          cigar += std::to_string(run);
          cigar += 'D';

          if (not flag)
            {
              md += std::to_string(matched);
              matched = 0;
              flag = true;
            }

          md += '^';
          for (auto i = 0LL; i < run; ++i)
            {
              md += targetseq[tpos];
              ++tpos;
            }
          flag = false;
          break;
        default:
          fatal("Invalid CIGAR string: unknown operation");
        }
    }

  if (not flag)
    {
      md += std::to_string(matched);
    }
}
}  // anonymous namespace

auto results_show_samheader(std::FILE * output_handle,
                            char const * dbname,
                            struct Database const & db,
                            struct Parameters const & parameters) -> void
{
  if ((parameters.opt_samout != nullptr) and parameters.opt_samheader)
    {
      fprint(output_handle, "@HD\tVN:1.0\tSO:unsorted\tGO:query\n");

      std::array<char, len_hex_dig_md5> md5hex;
      for (uint64_t i = 0; i < db.getsequencecount(); ++i)
        {
          get_hex_seq_digest_md5(md5hex, db.sequence_view(i));
          fprint(output_handle, "@SQ\tSN:");
          fprint(output_handle, db.header_view(i));
          fprint(output_handle, "\tLN:");
          fprint_integer(output_handle, db.getsequencelen(i));
          fprint(output_handle, "\tM5:");
          std::fputs(md5hex.data(), output_handle);
          fprint(output_handle, "\tUR:file:");
          std::fputs(dbname, output_handle);
          fprint(output_handle, '\n');
        }

      fprint(output_handle, "@PG\tID:");
      std::fputs(PROG_NAME, output_handle);
      fprint(output_handle, "\tVN:");
      std::fputs(PROG_VERSION, output_handle);
      fprint(output_handle, "\tCL:");
      fprint(output_handle, make_view(parameters.command_line));
      fprint(output_handle, '\n');
    }
}


auto results_show_samout(std::FILE * output_handle,
                         View<struct hit> const hits,
                         View<char> const query_head,
                         View<char> const qsequence,
                         View<char> const qsequence_rc,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void
{
  /*
    SAM format output

    http://samtools.github.io/hts-specs/SAMv1.pdf
    http://www.drive5.com/usearch/manual/sam_files.html
    http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output
    http://davetang.org/muse/2011/01/28/perl-and-sam/

    1: qname, query template name
    2: flag, bitwise flag (12 bits)
    (0x004=unmapped, 0x010=rev strand, 0x100 sec. alignment)
    3: rname, reference sequence name
    4: pos, 1-based leftmost mapping position (1)
    5: mapq, mapping quality (255)
    6: cigar, cigar string (MID)
    7: rnext, ref name of next/paired read (*)
    8: pnest, position of next/paired read (0)
    9: tlen, obs template length (target length)
    10: seq, segment of sequence
    11: qual, ascii of phred based quality+33 (*)
    12: optional tags (tag:type:value)

    Optional tags AS, XN, XM, XO, XG, NM, MD and YT used in usearch8.

    Usearch8:

    AS:i:? alignment score (i.e percent identity)
    XN:i:? next best alignment score (always 0?)
    XM:i:? number of mismatches
    XO:i:? number of gap opens (excluding terminal gaps)
    XG:i:? number of gap extensions (excluding terminal gaps)
    NM:i:? edit distance (sum of XM and XG)
    MD:Z:? variant string
    YT:Z:UU string representing alignment type

  */


  if (hits.empty()) {
    if (parameters.opt_output_no_hits != 0) {
      fprint(output_handle, query_head);
      fprint(output_handle, '\t');
      fprint_integer(output_handle, 0x04);
      fprint(output_handle, '\t');
      fprint(output_handle, "*");
      fprint(output_handle, '\t');
      fprint_integer(output_handle, static_cast<uint64_t>(0));
      fprint(output_handle, '\t');
      fprint_integer(output_handle, 255);
      fprint(output_handle, '\t');
      fprint(output_handle, "*");
      fprint(output_handle, '\t');
      fprint(output_handle, "*");
      fprint(output_handle, '\t');
      fprint_integer(output_handle, static_cast<uint64_t>(0));
      fprint(output_handle, '\t');
      fprint_integer(output_handle, static_cast<uint64_t>(0));
      fprint(output_handle, '\t');
      fprint(output_handle, qsequence);
      fprint(output_handle, "\t*\n");
    }
    return;
  }

  auto const top = top_hits(hits, parameters.opt_top_hits_only != 0);

  /* indexed rather than a range-for: the SAM flag marks every hit after the
     first as a secondary alignment (0x100), so the position is part of the
     output */
  for (std::size_t t = 0; t < top.size(); ++t)
    {
      auto const & hit = top[t];

      std::string cigar;
      std::string md;

      auto const target = static_cast<uint64_t>(hit.target);
      auto const query = (hit.strand != 0) ? qsequence_rc : qsequence;
      build_sam_strings(make_view(hit.nwalignment),
                        query,
                        db.sequence_view(target),
                        cigar,
                        md);

      OutputRecord record {output_handle};
      fprint(record, query_head);
      fprint(record, '\t');
      fprint_integer(record, (0x10 * hit.strand) | (t > 0 ? 0x100 : 0));
      fprint(record, '\t');
      fprint(record, db.header_view(target));
      fprint(record, '\t');
      fprint_integer(record, static_cast<uint64_t>(1));
      fprint(record, '\t');
      fprint_integer(record, 255);
      fprint(record, '\t');
      fprint(record, make_view(cigar));
      fprint(record, '\t');
      fprint(record, "*");
      fprint(record, '\t');
      fprint_integer(record, static_cast<uint64_t>(0));
      fprint(record, '\t');
      fprint_integer(record, static_cast<uint64_t>(0));
      fprint(record, '\t');
      fprint(record, query);
      fprint(record, '\t');
      fprint(record, "*");
      fprint(record, "\tAS:i:");
      std::fprintf(record.stream(), "%.0f", hit.id);
      fprint(record, "\tXN:i:");
      fprint_integer(record, 0);
      fprint(record, "\tXM:i:");
      fprint_integer(record, hit.mismatches);
      fprint(record, "\tXO:i:");
      fprint_integer(record, hit.internal_gaps);
      fprint(record, "\tXG:i:");
      fprint_integer(record, hit.internal_indels);
      fprint(record, "\tNM:i:");
      fprint_integer(record, hit.mismatches + hit.internal_indels);
      fprint(record, "\tMD:Z:");
      fprint(record, make_view(md));
      fprint(record, "\tYT:Z:");
      fprint(record, "UU");
      fprint(record, '\n');
    }
}
