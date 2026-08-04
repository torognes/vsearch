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
#include "core/searchcore.hpp"  // struct hit
#include "core/showalign.hpp"
#include "core/tax.hpp"
#include "utils/cigar.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/view.hpp"
#include "utils/taxonomic_fields.h"
#include "utils/sequence_digest.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/prog_id.hpp"  // PROG_NAME, PROG_VERSION
#include <algorithm>  // std::max
#include <array>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <iterator>  // std::next
#include <string>  // std::string, std::to_string


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto check_if_perfect_match(char const * opt_cluster_fast,
                              struct hit const * hits) -> bool {
    if (opt_cluster_fast != nullptr) {
      /* cluster_fast */
      /* use '=' for identical sequences, ignoring terminal gaps */
      return (hits->matches == hits->internal_alignmentlength);
    }
    /* cluster_size, cluster_smallmem, cluster_unoise */
    /* usearch_global, search_exact, allpairs_global */
    /* use '=' for strictly identical sequences */
    return (hits->matches == hits->nwalignmentlength);
  }

}  // end of anonymous namespace


auto results_show_fastapairs_one(std::FILE * output_handle,
                                 struct hit const * hits,
                                 View<char> const query_head,
                                 View<char> const qsequence,
                                 View<char> const qsequence_rc,
                                 struct Database const & db,
                                 struct Parameters const & parameters) -> void
{
  /* http://www.drive5.com/usearch/manual/fastapairs.html */

  if (hits == nullptr) {
    return;
  }

  auto const query = (hits->strand != 0) ? qsequence_rc : qsequence;
  auto const qrow = get_alignment_qrow(query,
                                 View<char>{hits->nwalignment.c_str(), hits->nwalignment.size()},
                                 hits->nwalignmentlength);
  fasta_print_general(output_handle,
                      nullptr,
                      View<char>{&qrow[static_cast<std::size_t>(hits->trim_q_left + hits->trim_t_left)],
                                 static_cast<std::size_t>(hits->internal_alignmentlength)},
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

  auto const target = static_cast<uint64_t>(hits->target);
  auto const trow = get_alignment_trow(db.sequence_view(target),
                                 View<char>{hits->nwalignment.c_str(), hits->nwalignment.size()},
                                 hits->nwalignmentlength);
  fasta_print_general(output_handle,
                      nullptr,
                      View<char>{&trow[static_cast<std::size_t>(hits->trim_q_left + hits->trim_t_left)],
                                 static_cast<std::size_t>(hits->internal_alignmentlength)},
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
                              struct hit const * hits,
                              View<char> const query_head,
                              View<char> const qsequence,
                              View<char> const qsequence_rc,
                              struct Parameters const & parameters) -> void
{
  if (hits == nullptr) {
    return;
  }

  auto const query = (hits->strand != 0) ? qsequence_rc : qsequence;
  auto const qseglen = static_cast<int64_t>(query.size())
    - hits->trim_q_left - hits->trim_q_right;
  auto const qseg = query.subspan(static_cast<std::size_t>(hits->trim_q_left),
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
                              struct hit const * hits,
                              struct Database const & db,
                              struct Parameters const & parameters) -> void
{
  if (hits == nullptr) {
    return;
  }
  auto const target = static_cast<uint64_t>(hits->target);
  auto const target_sequence = db.sequence_view(target);
  auto const tseglen = static_cast<int64_t>(target_sequence.size())
    - hits->trim_t_left - hits->trim_t_right;
  auto const tseg = target_sequence.subspan(static_cast<std::size_t>(hits->trim_t_left),
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
                                struct hit const * hits,
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

  if (hits == nullptr) {
    fprint(output_handle, query_head);
    fprint(output_handle, "\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n");
    return;
  }
  // if 'hp->strand' then 'minus strand' else 'plus strand'
  auto const target = static_cast<uint64_t>(hits->target);
  int const qstart = (hits->strand != 0) ? static_cast<int>(qseqlen) : 1;
  int const qend = (hits->strand != 0) ? 1 : static_cast<int>(qseqlen);

  fprint(output_handle, query_head);
  fprint(output_handle, '\t');
  fprint(output_handle, db.header_view(target));
  fprint(output_handle, '\t');
  std::fprintf(output_handle, "%.1f", hits->id);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, hits->internal_alignmentlength);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, hits->mismatches);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, hits->internal_gaps);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, qstart);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, qend);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, 1);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, db.getsequencelen(target));
  fprint(output_handle, '\t');
  fprint_integer(output_handle, -1);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, 0);
  fprint(output_handle, '\n');
}


auto results_show_uc_one(std::FILE * output_handle,
                         struct hit const * hits,
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

  if (hits == nullptr) {
    fprint(output_handle, "N\t*\t*\t*\t.\t*\t*\t*\t");
    fprint(output_handle, query_head);
    fprint(output_handle, "\t*\n");
    return;
  }

  auto const is_perfect_match = check_if_perfect_match(parameters.opt_cluster_fast, hits);

  fprint(output_handle, "H\t");
  fprint_integer(output_handle, clusterno);
  fprint(output_handle, '\t');
  fprint_integer(output_handle, qseqlen);
  fprint(output_handle, '\t');
  std::fprintf(output_handle, "%.1f", hits->id);
  fprint(output_handle, '\t');
  fprint(output_handle, (hits->strand != 0) ? '-' : '+');
  fprint(output_handle, "\t0\t0\t");
  std::fputs(is_perfect_match ? "=" : hits->nwalignment.c_str(), output_handle);
  fprint(output_handle, '\t');
  auto const target = static_cast<uint64_t>(hits->target);
  header_fprint_strip(output_handle,
                      query_head,
                      parameters.opt_xsize,
                      parameters.opt_xee,
                      parameters.opt_xlength);
  fprint(output_handle, '\t');
  header_fprint_strip(output_handle,
                      db.header_view(target),
                      parameters.opt_xsize,
                      parameters.opt_xee,
                      parameters.opt_xlength);
  fprint(output_handle, '\n');
}


auto results_show_userout_one(std::FILE * output_handle, struct hit const * hits,
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

  auto const qseqlen = static_cast<int64_t>(qsequence.size());
  auto const & userfields_requested = parameters.opt_userfields;

  for (std::size_t c = 0; c < userfields_requested.size(); ++c)
    {
      if (c != 0)
        {
          fprint(output_handle, '\t');
        }

      auto const field = userfields_requested[c];

      View<char> tsequence;
      int64_t tseqlen = 0;
      char const * t_head = nullptr;

      if (hits != nullptr)
        {
          auto const target = static_cast<uint64_t>(hits->target);
          tsequence = db.sequence_view(target);
          tseqlen = static_cast<int64_t>(db.getsequencelen(target));
          t_head = db.getheader(target);
        }


      switch (field)
        {
        case 0: /* query */
          fprint(output_handle, query_head);
          break;
        case 1: /* target */
          std::fputs((hits != nullptr) ? t_head : "*", output_handle);
          break;
        case 2: /* evalue */
          fprint(output_handle, "-1");
          break;
        case 3: /* id */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? hits->id : 0.0);
          break;
        case 4: /* pctpv */
          std::fprintf(output_handle, "%.1f", ((hits != nullptr) and (hits->internal_alignmentlength > 0)) ? 100.0 * hits->matches / hits->internal_alignmentlength : 0.0);
          break;
        case 5: /* pctgaps */
          std::fprintf(output_handle, "%.1f", ((hits != nullptr) and (hits->internal_alignmentlength > 0)) ? 100.0 * hits->internal_indels / hits->internal_alignmentlength : 0.0);
          break;
        case 6: /* pairs */
          fprint_integer(output_handle, (hits != nullptr) ? hits->matches + hits->mismatches : 0);
          break;
        case 7: /* gaps */
          fprint_integer(output_handle, (hits != nullptr) ? hits->internal_indels : 0);
          break;
        case 8: /* qlo */
          fprint_integer(output_handle, (hits != nullptr) ? ((hits->strand != 0) ? qseqlen : 1) : 0);
          break;
        case 9: /* qhi */
          fprint_integer(output_handle, (hits != nullptr) ? ((hits->strand != 0) ? 1 : qseqlen) : 0);
          break;
        case 10: /* tlo */
          fprint_integer(output_handle, (hits != nullptr) ? 1 : 0);
          break;
        case 11: /* thi */
          fprint_integer(output_handle, tseqlen);
          break;
        case 12: /* pv */
          fprint_integer(output_handle, (hits != nullptr) ? hits->matches : 0);
          break;
        case 13: /* ql */
          fprint_integer(output_handle, qseqlen);
          break;
        case 14: /* tl */
          fprint_integer(output_handle, (hits != nullptr) ? tseqlen : 0);
          break;
        case 15: /* qs */
          fprint_integer(output_handle, qseqlen);
          break;
        case 16: /* ts */
          fprint_integer(output_handle, (hits != nullptr) ? tseqlen : 0);
          break;
        case 17: /* alnlen */
          fprint_integer(output_handle, (hits != nullptr) ? hits->internal_alignmentlength : 0);
          break;
        case 18: /* opens */
          fprint_integer(output_handle, (hits != nullptr) ? hits->internal_gaps : 0);
          break;
        case 19: /* exts */
          fprint_integer(output_handle, (hits != nullptr) ? hits->internal_indels - hits->internal_gaps : 0);
          break;
        case 20: /* raw */
          fprint_integer(output_handle, (hits != nullptr) ? hits->nwscore : 0);
          break;
        case 21: /* bits */
          fprint_integer(output_handle, 0);
          break;
        case 22: /* aln */
          if (hits != nullptr)
            {
              print_uncompressed_cigar(output_handle, View<char>{hits->nwalignment.c_str(), hits->nwalignment.size()});
            }
          break;
        case 23: /* caln */
          if (hits != nullptr)
            {
              fprint(output_handle, View<char>{hits->nwalignment.data(), hits->nwalignment.size()});
            }
          break;
        case 24: /* qstrand */
          if (hits != nullptr)
            {
              fprint(output_handle, (hits->strand != 0) ? '-' : '+');
            }
          break;
        case 25: /* tstrand */
          if (hits != nullptr)
            {
              fprint(output_handle, '+');
            }
          break;
        case 26: /* qrow */
          if (hits != nullptr)
            {
              auto const query = (hits->strand != 0) ? qsequence_rc : qsequence;
              auto const qrow = get_alignment_qrow(query,
                                             View<char>{hits->nwalignment.c_str(), hits->nwalignment.size()},
                                             hits->nwalignmentlength);
              fprint(output_handle,
                     View<char>{&qrow[static_cast<std::size_t>(hits->trim_q_left + hits->trim_t_left)],
                                static_cast<std::size_t>(hits->internal_alignmentlength)});
            }
          break;
        case 27: /* trow */
          if (hits != nullptr)
            {
              auto const trow = get_alignment_trow(tsequence,
                                             View<char>{hits->nwalignment.c_str(), hits->nwalignment.size()},
                                             hits->nwalignmentlength);
              fprint(output_handle,
                     View<char>{&trow[static_cast<std::size_t>(hits->trim_q_left + hits->trim_t_left)],
                                static_cast<std::size_t>(hits->internal_alignmentlength)});
            }
          break;
        case 28: /* qframe */
        case 29: /* tframe */
          fprint(output_handle, "+0");
          break;
        case 30: /* mism */
          fprint_integer(output_handle, (hits != nullptr) ? hits->mismatches : 0);
          break;
        case 31: /* ids */
          fprint_integer(output_handle, (hits != nullptr) ? hits->matches : 0);
          break;
        case 32: /* qcov */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? 100.0 * (hits->matches + hits->mismatches) / static_cast<double>(qseqlen) : 0.0);
          break;
        case 33: /* tcov */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? 100.0 * (hits->matches + hits->mismatches) / static_cast<double>(tseqlen) : 0.0);
          break;
        case 34: /* id0 */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? hits->id0 : 0.0);
          break;
        case 35: /* id1 */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? hits->id1 : 0.0);
          break;
        case 36: /* id2 */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? hits->id2 : 0.0);
          break;
        case 37: /* id3 */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? hits->id3 : 0.0);
          break;
        case 38: /* id4 */
          std::fprintf(output_handle, "%.1f", (hits != nullptr) ? hits->id4 : 0.0);
          break;

          /* new internal alignment coordinates */

        case 39: /* qilo */
          fprint_integer(output_handle, (hits != nullptr) ? hits->trim_q_left + 1 : 0);
          break;
        case 40: /* qihi */
          fprint_integer(output_handle, (hits != nullptr) ? qseqlen - hits->trim_q_right : 0);
          break;
        case 41: /* tilo */
          fprint_integer(output_handle, (hits != nullptr) ? hits->trim_t_left + 1 : 0);
          break;
        case 42: /* tihi */
          fprint_integer(output_handle, (hits != nullptr) ? tseqlen - hits->trim_t_right : 0);
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
  fprint(output_handle, '\n');
}


auto results_show_lcaout(std::FILE * output_handle,
                         struct hit const * hits,
                         int const hitcount,
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

  if (hitcount == 0) {
    fprint(output_handle, '\n');
    return;
  }

  constexpr auto levels = static_cast<std::size_t>(tax_levels);
  std::array<int, tax_levels> votes {{}};
  std::array<int, tax_levels> cand;
  cand.fill(-1);
  std::array<std::array<int, tax_levels>, tax_levels> cand_level_start {{}};
  std::array<std::array<int, tax_levels>, tax_levels> cand_level_len {{}};
  std::array<int, tax_levels> level_match {{}};

  auto const top_hit_id = hits[0].id;
  auto tophitcount = 0;

  for (auto t = 0; t < hitcount; ++t)
    {
      struct hit const * hp = hits + t;

      if ((parameters.opt_top_hits_only != 0) and (hp->id < top_hit_id))
        {
          break;
        }

      ++tophitcount;

      int const seqno = hp->target;
      std::array<int, tax_levels> new_level_start {{}};  // refactoring: std::array<struct a_level{.start, .length}, tax_levels>
      std::array<int, tax_levels> new_level_len {{}};
      tax_split(seqno, new_level_start.data(), new_level_len.data(), db);

      for (std::size_t k = 0; k < levels; ++k)
        {
          if (votes[k] == 0)
            {
              cand[k] = seqno;
              votes[k] = 1;
              for (std::size_t j = 0; j < levels; ++j)
                {
                  cand_level_start[k][j] = new_level_start[j];
                  cand_level_len[k][j] = new_level_len[j];
                }
            }
          else
            {
              auto match = true;
              for (std::size_t j = 0; j <= k; ++j)
                {
                  auto const cand_level = db.header_view(static_cast<uint64_t>(cand[k]))
                    .subspan(static_cast<std::size_t>(cand_level_start[k][j]),
                             static_cast<std::size_t>(cand_level_len[k][j]));
                  auto const query_level = db.header_view(static_cast<uint64_t>(seqno))
                    .subspan(static_cast<std::size_t>(new_level_start[j]),
                             static_cast<std::size_t>(new_level_len[j]));
                  if (cand_level != query_level)
                    {
                      match = false;
                      break;
                    }
                }
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

  for (auto t = 0; t < tophitcount; ++t)
    {
      auto const seqno = hits[t].target;
      std::array<int, tax_levels> new_level_start {{}};
      std::array<int, tax_levels> new_level_len {{}};
      tax_split(seqno, new_level_start.data(), new_level_len.data(), db);

      for (std::size_t k = 0; k < levels; ++k)
        {
          auto match = true;
          for (std::size_t j = 0; j <= k; ++j)
            {
              auto const cand_level = db.header_view(static_cast<uint64_t>(cand[k]))
                .subspan(static_cast<std::size_t>(cand_level_start[k][j]),
                         static_cast<std::size_t>(cand_level_len[k][j]));
              auto const query_level = db.header_view(static_cast<uint64_t>(seqno))
                .subspan(static_cast<std::size_t>(new_level_start[j]),
                         static_cast<std::size_t>(new_level_len[j]));
              if (cand_level != query_level)
                {
                  match = false;
                  break;
                }
            }
          if (match)
            {
              ++level_match[k];
            }
        }
    }

  /* output results */

  if (tophitcount == 0) {
    fprint(output_handle, '\n');
    return;
  }
  auto comma = false;
  for (std::size_t j = 0; j < levels; ++j)
    {
      if (1.0 * level_match[j] / tophitcount < parameters.opt_lca_cutoff)
        {
          break;
        }

      if (cand_level_len[j][j] > 0)
        {
          std::fputs((comma ? "," : ""), output_handle);
          fprint(output_handle, static_cast<char>(taxonomic_fields[j]));
          fprint(output_handle, ':');
          fprint(output_handle,
                 db.header_view(static_cast<uint64_t>(cand[j]))
                   .subspan(static_cast<std::size_t>(cand_level_start[j][j]),
                            static_cast<std::size_t>(cand_level_len[j][j])));
          comma = true;
        }
    }

  fprint(output_handle, '\n');
}


auto results_show_alnout(std::FILE * output_handle,
                         struct hit const * hits,
                         int const hitcount,
                         View<char> const query_head,
                         View<char> const qsequence,
                         struct Database const & db,
                         struct Parameters const & parameters) -> void
{
  /* http://drive5.com/usearch/manual/alnout.html */


  if (hitcount == 0) {
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

  auto const top_hit_id = hits[0].id;

  for (auto t = 0; t < hitcount; ++t)
    {
      auto const * hp = hits + t;

      if ((parameters.opt_top_hits_only != 0) and (hp->id < top_hit_id))
        {
          break;
        }

      auto const target = static_cast<uint64_t>(hp->target);
      std::fprintf(output_handle, "%3.0f", hp->id);
      fprint(output_handle, "% ");
      fprint_integer(output_handle, db.getsequencelen(target), 6);
      fprint(output_handle, "  ");
      fprint(output_handle, db.header_view(target));
      fprint(output_handle, '\n');
    }

  for (auto t = 0; t < hitcount; ++t)
    {
      auto const * hp = hits + t;

      if ((parameters.opt_top_hits_only != 0) and (hp->id < top_hit_id))
        {
          break;
        }

      fprint(output_handle, '\n');


      auto const target = static_cast<uint64_t>(hp->target);
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

      fprint(output_handle, " Query ");
      fprint_integer(output_handle, qseqlen, static_cast<std::size_t>(numwidth));
      fprint(output_handle, "nt >");
      fprint(output_handle, query_head);
      fprint(output_handle, "\nTarget ");
      fprint_integer(output_handle, dseqlen, static_cast<std::size_t>(numwidth));
      fprint(output_handle, "nt >");
      fprint(output_handle, db.header_view(target));
      fprint(output_handle, '\n');

      int64_t const rowlen = (parameters.opt_rowlen == 0) ? (qseqlen + dseqlen) : parameters.opt_rowlen;

      align_show(output_handle,
                 qsequence,
                 hp->trim_q_left,
                 "Qry",
                 db.sequence_view(target),
                 hp->trim_t_left,
                 "Tgt",
                 View<char>{std::next(hp->nwalignment.c_str(), hp->trim_aln_left),
                            hp->nwalignment.size()
                            - static_cast<std::size_t>(hp->trim_aln_left)
                            - static_cast<std::size_t>(hp->trim_aln_right)},
                 numwidth,
                 3,
                 rowlen,
                 hp->strand,
                 parameters);

      fprint(output_handle, '\n');
      fprint_integer(output_handle, hp->internal_alignmentlength);
      fprint(output_handle, " cols, ");
      fprint_integer(output_handle, hp->matches);
      fprint(output_handle, " ids (");
      std::fprintf(output_handle, "%3.1f", hp->id);
      fprint(output_handle, "%), ");
      fprint_integer(output_handle, hp->internal_indels);
      fprint(output_handle, " gaps (");
      std::fprintf(output_handle, "%3.1f", hp->internal_alignmentlength > 0 ?
              100.0 * hp->internal_indels / hp->internal_alignmentlength :
              0.0);
      fprint(output_handle, "%)\n");
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
      auto const scanned_run = read_runlength(p, &next_operation);
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
          get_hex_seq_digest_md5(Span<char>{md5hex.data(), md5hex.size()},
                                 db.sequence_view(i));
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
      fprint(output_handle, View<char>{parameters.command_line.data(), parameters.command_line.size()});
      fprint(output_handle, '\n');
    }
}


auto results_show_samout(std::FILE * output_handle,
                         struct hit const * hits,
                         int const hitcount,
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


  if (hitcount == 0) {
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

  auto const top_hit_id = hits[0].id;

  for (auto t = 0; t < hitcount; ++t)
    {
      auto const * hp = hits + t;

      if ((parameters.opt_top_hits_only != 0) and (hp->id < top_hit_id))
        {
          break;
        }


      std::string cigar;
      std::string md;

      auto const target = static_cast<uint64_t>(hp->target);
      auto const query = (hp->strand != 0) ? qsequence_rc : qsequence;
      build_sam_strings(View<char>{hp->nwalignment.c_str(), hp->nwalignment.size()},
                        query,
                        db.sequence_view(target),
                        cigar,
                        md);

      fprint(output_handle, query_head);
      fprint(output_handle, '\t');
      fprint_integer(output_handle, (0x10 * hp->strand) | (t > 0 ? 0x100 : 0));
      fprint(output_handle, '\t');
      fprint(output_handle, db.header_view(target));
      fprint(output_handle, '\t');
      fprint_integer(output_handle, static_cast<uint64_t>(1));
      fprint(output_handle, '\t');
      fprint_integer(output_handle, 255);
      fprint(output_handle, '\t');
      fprint(output_handle, View<char>{cigar.data(), cigar.size()});
      fprint(output_handle, '\t');
      fprint(output_handle, "*");
      fprint(output_handle, '\t');
      fprint_integer(output_handle, static_cast<uint64_t>(0));
      fprint(output_handle, '\t');
      fprint_integer(output_handle, static_cast<uint64_t>(0));
      fprint(output_handle, '\t');
      fprint(output_handle, query);
      fprint(output_handle, '\t');
      fprint(output_handle, "*");
      fprint(output_handle, "\tAS:i:");
      std::fprintf(output_handle, "%.0f", hp->id);
      fprint(output_handle, "\tXN:i:");
      fprint_integer(output_handle, 0);
      fprint(output_handle, "\tXM:i:");
      fprint_integer(output_handle, hp->mismatches);
      fprint(output_handle, "\tXO:i:");
      fprint_integer(output_handle, hp->internal_gaps);
      fprint(output_handle, "\tXG:i:");
      fprint_integer(output_handle, hp->internal_indels);
      fprint(output_handle, "\tNM:i:");
      fprint_integer(output_handle, hp->mismatches + hp->internal_indels);
      fprint(output_handle, "\tMD:Z:");
      fprint(output_handle, View<char>{md.data(), md.size()});
      fprint(output_handle, "\tYT:Z:");
      fprint(output_handle, "UU");
      fprint(output_handle, '\n');
    }
}
