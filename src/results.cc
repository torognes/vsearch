/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "attributes.h"
#include "showalign.h"
#include "tax.h"
#include "userfields.h"
#include <algorithm>  // std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::snprintf, std::sscanf
#include <cstring>  // std::strlen, std::strncmp


auto results_show_fastapairs_one(std::FILE * output_handle,
                                 struct hit * hits,
                                 char * query_head,
                                 char * qsequence,
                                 char * qsequence_rc) -> void
{
  /* http://www.drive5.com/usearch/manual/fastapairs.html */

  if (hits == nullptr) {
    return;
  }

  auto * qrow = align_getrow(hits->strand ? qsequence_rc : qsequence,
                             hits->nwalignment,
                             hits->nwalignmentlength,
                             0);
  fasta_print_general(output_handle,
                      nullptr,
                      qrow + hits->trim_q_left + hits->trim_t_left,
                      hits->internal_alignmentlength,
                      query_head,
                      strlen(query_head),
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0);
  xfree(qrow);

  auto * trow = align_getrow(db_getsequence(hits->target),
                             hits->nwalignment,
                             hits->nwalignmentlength,
                             1);
  fasta_print_general(output_handle,
                      nullptr,
                      trow + hits->trim_q_left + hits->trim_t_left,
                      hits->internal_alignmentlength,
                      db_getheader(hits->target),
                      db_getheaderlen(hits->target),
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0);
  xfree(trow);

  fprintf(output_handle, "\n");
}


auto results_show_qsegout_one(std::FILE * output_handle,
                              struct hit * hits,
                              char * query_head,
                              char * qsequence,
                              int64_t qseqlen,
                              char * qsequence_rc) -> void
{
  if (hits == nullptr) {
    return;
  }

  char * qseg = (hits->strand ? qsequence_rc : qsequence) + hits->trim_q_left;
  int const qseglen = qseqlen - hits->trim_q_left - hits->trim_q_right;

  fasta_print_general(output_handle,
                      nullptr,
                      qseg,
                      qseglen,
                      query_head,
                      strlen(query_head),
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0);
}


auto results_show_tsegout_one(std::FILE * output_handle,
                              struct hit * hits) -> void
{
  if (hits == nullptr) {
    return;
  }
  auto * tseg = db_getsequence(hits->target) + hits->trim_t_left;
  int const tseglen = db_getsequencelen(hits->target) - hits->trim_t_left - hits->trim_t_right;

  fasta_print_general(output_handle,
                      nullptr,
                      tseg,
                      tseglen,
                      db_getheader(hits->target),
                      db_getheaderlen(hits->target),
                      0,
                      0,
                      -1.0,
                      -1,
                      -1,
                      nullptr,
                      0.0);
}


auto results_show_blast6out_one(std::FILE * output_handle,
                                struct hit * hits,
                                char * query_head,
                                int64_t qseqlen) -> void
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
    fprintf(output_handle, "%s\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n", query_head);
    return;
  }
  // if 'hp->strand' then 'minus strand' else 'plus strand'
  const int qstart = hits->strand ? qseqlen : 1;
  const int qend = hits->strand ? 1 : qseqlen;

  fprintf(output_handle,
          "%s\t%s\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%" PRIu64 "\t%d\t%d\n",
          query_head,
          db_getheader(hits->target),
          hits->id,
          hits->internal_alignmentlength,
          hits->mismatches,
          hits->internal_gaps,
          qstart,
          qend,
          1,
          db_getsequencelen(hits->target),
          -1,
          0);
}


auto results_show_uc_one(std::FILE * output_handle,
                         struct hit * hits,
                         char * query_head,
                         int64_t qseqlen,
                         int clusterno) -> void
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

  if (hits != nullptr)
    {
      auto perfect = false;

      if (opt_cluster_fast)
        {
          /* cluster_fast */
          /* use = for identical sequences ignoring terminal gaps */
          perfect = (hits->matches == hits->internal_alignmentlength);
        }
      else
        {
          /* cluster_size, cluster_smallmem, cluster_unoise */
          /* usearch_global, search_exact, allpairs_global */
          /* use = for strictly identical sequences */
          perfect = (hits->matches == hits->nwalignmentlength);
        }

      fprintf(output_handle,
              "H\t%d\t%" PRId64 "\t%.1f\t%c\t0\t0\t%s\t",
              clusterno,
              qseqlen,
              hits->id,
              hits->strand ? '-' : '+',
              perfect ? "=" : hits->nwalignment);
      header_fprint_strip(output_handle,
                          query_head,
                          strlen(query_head),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(output_handle, "\t");
      header_fprint_strip(output_handle,
                          db_getheader(hits->target),
                          db_getheaderlen(hits->target),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(output_handle, "\n");
    }
  else
    {
      fprintf(output_handle, "N\t*\t*\t*\t.\t*\t*\t*\t%s\t*\n", query_head);
    }
}


auto results_show_userout_one(std::FILE * output_handle, struct hit * hits,
                              char * query_head,
                              char * qsequence, int64_t qseqlen,
                              char * qsequence_rc) -> void
{

  /*
    http://drive5.com/usearch/manual/userout.html
    qlo, qhi, tlo, thi and raw are given more meaningful values here
  */

  for (auto c = 0; c < userfields_requested_count; c++)
    {
      if (c != 0)
        {
          fprintf(output_handle, "\t");
        }

      auto const field = userfields_requested[c];

      char * tsequence = nullptr;
      int64_t tseqlen = 0;
      char * t_head = nullptr;

      if (hits != nullptr)
        {
          tsequence = db_getsequence(hits->target);
          tseqlen = db_getsequencelen(hits->target);
          t_head = db_getheader(hits->target);
        }

      char * qrow = nullptr;
      char * trow = nullptr;

      switch (field)
        {
        case 0: /* query */
          fprintf(output_handle, "%s", query_head);
          break;
        case 1: /* target */
          fprintf(output_handle, "%s", hits ? t_head : "*");
          break;
        case 2: /* evalue */
          fprintf(output_handle, "-1");
          break;
        case 3: /* id */
          fprintf(output_handle, "%.1f", hits ? hits->id : 0.0);
          break;
        case 4: /* pctpv */
          fprintf(output_handle, "%.1f", (hits && (hits->internal_alignmentlength > 0)) ? 100.0 * hits->matches / hits->internal_alignmentlength : 0.0);
          break;
        case 5: /* pctgaps */
          fprintf(output_handle, "%.1f", (hits && (hits->internal_alignmentlength > 0)) ? 100.0 * hits->internal_indels / hits->internal_alignmentlength : 0.0);
          break;
        case 6: /* pairs */
          fprintf(output_handle, "%d", hits ? hits->matches + hits->mismatches : 0);
          break;
        case 7: /* gaps */
          fprintf(output_handle, "%d", hits ? hits->internal_indels : 0);
          break;
        case 8: /* qlo */
          fprintf(output_handle, "%" PRId64, hits ? (hits->strand ? qseqlen : 1) : 0);
          break;
        case 9: /* qhi */
          fprintf(output_handle, "%" PRId64, hits ? (hits->strand ? 1 : qseqlen) : 0);
          break;
        case 10: /* tlo */
          fprintf(output_handle, "%d", hits ? 1 : 0);
          break;
        case 11: /* thi */
          fprintf(output_handle, "%" PRId64, tseqlen);
          break;
        case 12: /* pv */
          fprintf(output_handle, "%d", hits ? hits->matches : 0);
          break;
        case 13: /* ql */
          fprintf(output_handle, "%" PRId64, qseqlen);
          break;
        case 14: /* tl */
          fprintf(output_handle, "%" PRId64, hits ? tseqlen : 0);
          break;
        case 15: /* qs */
          fprintf(output_handle, "%" PRId64, qseqlen);
          break;
        case 16: /* ts */
          fprintf(output_handle, "%" PRId64, hits ? tseqlen : 0);
          break;
        case 17: /* alnlen */
          fprintf(output_handle, "%d", hits ? hits->internal_alignmentlength : 0);
          break;
        case 18: /* opens */
          fprintf(output_handle, "%d", hits ? hits->internal_gaps : 0);
          break;
        case 19: /* exts */
          fprintf(output_handle, "%d", hits ? hits->internal_indels - hits->internal_gaps : 0);
          break;
        case 20: /* raw */
          fprintf(output_handle, "%d", hits ? hits->nwscore : 0);
          break;
        case 21: /* bits */
          fprintf(output_handle, "%d", 0);
          break;
        case 22: /* aln */
          if (hits)
            {
              align_fprint_uncompressed_alignment(output_handle, hits->nwalignment);
            }
          break;
        case 23: /* caln */
          if (hits)
            {
              fprintf(output_handle, "%s", hits->nwalignment);
            }
          break;
        case 24: /* qstrand */
          if (hits)
            {
              fprintf(output_handle, "%c", hits->strand ? '-' : '+');
            }
          break;
        case 25: /* tstrand */
          if (hits)
            {
              fprintf(output_handle, "%c", '+');
            }
          break;
        case 26: /* qrow */
          if (hits)
            {
              qrow = align_getrow(hits->strand ? qsequence_rc : qsequence,
                                  hits->nwalignment,
                                  hits->nwalignmentlength,
                                  0);
              fprintf(output_handle, "%.*s",
                      hits->internal_alignmentlength,
                      qrow + hits->trim_q_left + hits->trim_t_left);
              xfree(qrow);
            }
          break;
        case 27: /* trow */
          if (hits)
            {
              trow = align_getrow(tsequence,
                                  hits->nwalignment,
                                  hits->nwalignmentlength,
                                  1);
              fprintf(output_handle, "%.*s",
                      hits->internal_alignmentlength,
                      trow + hits->trim_q_left + hits->trim_t_left);
              xfree(trow);
            }
          break;
        case 28: /* qframe */
          fprintf(output_handle, "+0");
          break;
        case 29: /* tframe */
          fprintf(output_handle, "+0");
          break;
        case 30: /* mism */
          fprintf(output_handle, "%d", hits ? hits->mismatches : 0);
          break;
        case 31: /* ids */
          fprintf(output_handle, "%d", hits ? hits->matches : 0);
          break;
        case 32: /* qcov */
          fprintf(output_handle, "%.1f",
                  hits ? 100.0 * (hits->matches + hits->mismatches) / qseqlen : 0.0);
          break;
        case 33: /* tcov */
          fprintf(output_handle, "%.1f",
                  hits ? 100.0 * (hits->matches + hits->mismatches) / tseqlen : 0.0);
          break;
        case 34: /* id0 */
          fprintf(output_handle, "%.1f", hits ? hits->id0 : 0.0);
          break;
        case 35: /* id1 */
          fprintf(output_handle, "%.1f", hits ? hits->id1 : 0.0);
          break;
        case 36: /* id2 */
          fprintf(output_handle, "%.1f", hits ? hits->id2 : 0.0);
          break;
        case 37: /* id3 */
          fprintf(output_handle, "%.1f", hits ? hits->id3 : 0.0);
          break;
        case 38: /* id4 */
          fprintf(output_handle, "%.1f", hits ? hits->id4 : 0.0);
          break;

          /* new internal alignment coordinates */

        case 39: /* qilo */
          fprintf(output_handle, "%d", hits ? hits->trim_q_left + 1 : 0);
          break;
        case 40: /* qihi */
          fprintf(output_handle, "%" PRId64, hits ? qseqlen - hits->trim_q_right : 0);
          break;
        case 41: /* tilo */
          fprintf(output_handle, "%d", hits ? hits->trim_t_left + 1 : 0);
          break;
        case 42: /* tihi */
          fprintf(output_handle, "%" PRId64, hits ? tseqlen - hits->trim_t_right : 0);
          break;
        }
    }
  fprintf(output_handle, "\n");
}


auto results_show_lcaout(std::FILE * output_handle,
                         struct hit * hits,
                         int hitcount,
                         char * query_head) -> void
{
  /* Output last common ancestor (LCA) of the hits,
     in a similar way to the Sintax command */

  /* Use a modified Boyer-Moore majority voting algorithm at each taxonomic
     level to find the most common name at each level */

  fprintf(output_handle, "%s\t", query_head);

  int votes[tax_levels];
  int cand[tax_levels];
  int cand_level_start[tax_levels][tax_levels];
  int cand_level_len[tax_levels][tax_levels];
  int level_match[tax_levels];

  for (auto k = 0; k < tax_levels; k++)
    {
      votes[k] = 0;
      cand[k] = -1;
      level_match[k] = 0;
    }

  auto const top_hit_id = hits[0].id;
  auto tophitcount = 0;

  for (auto t = 0; t < hitcount; t++)
    {
      struct hit * hp = hits + t;

      if (opt_top_hits_only && (hp->id < top_hit_id))
        {
          break;
        }

      tophitcount++;

      int const seqno = hp->target;
      int new_level_start[tax_levels];
      int new_level_len[tax_levels];
      tax_split(seqno, new_level_start, new_level_len);

      for (auto k = 0; k < tax_levels; k++)
        {
          if (votes[k] == 0)
            {
              cand[k] = seqno;
              votes[k] = 1;
              for (auto j = 0; j < tax_levels; j++)
                {
                  cand_level_start[k][j] = new_level_start[j];
                  cand_level_len[k][j] = new_level_len[j];
                }
            }
          else
            {
              auto match = true;
              for (auto j = 0; j <= k; j++)
                {
                  if ((new_level_len[j] != cand_level_len[k][j]) ||
                      (strncmp(db_getheader(cand[k]) + cand_level_start[k][j],
                               db_getheader(seqno) + new_level_start[j],
                               new_level_len[j]) != 0))
                    {
                      match = false;
                      break;
                    }
                }
              if (match)
                {
                  votes[k]++;
                }
              else
                {
                  votes[k]--;
                }
            }
        }
    }

  /* count actual matches to the candidate at each level */

  for (auto t = 0; t < tophitcount; t++)
    {
      auto const seqno = hits[t].target;
      int new_level_start[tax_levels];
      int new_level_len[tax_levels];
      tax_split(seqno, new_level_start, new_level_len);

      for (auto k = 0; k < tax_levels; k++)
        {
          auto match = true;
          for (auto j = 0; j <= k; j++)
            {
              if ((new_level_len[j] != cand_level_len[k][j]) ||
                  (strncmp(db_getheader(cand[k]) + cand_level_start[k][j],
                           db_getheader(seqno) + new_level_start[j],
                           new_level_len[j]) != 0))
                {
                  match = false;
                  break;
                }
            }
          if (match)
            {
              level_match[k]++;
            }
        }
    }

  /* output results */

  if (tophitcount > 0)
    {
      auto comma = false;
      for (auto j = 0; j < tax_levels; j++)
        {
          if (1.0 * level_match[j] / tophitcount < opt_lca_cutoff)
            {
              break;
            }

          if (cand_level_len[j][j] > 0)
            {
              fprintf(output_handle,
                      "%s%c:%.*s",
                      (comma ? "," : ""),
                      tax_letters[j],
                      cand_level_len[j][j],
                      db_getheader(cand[j]) + cand_level_start[j][j]);
              comma = true;
            }
        }
    }

  fprintf(output_handle, "\n");
}


auto results_show_alnout(std::FILE * output_handle,
                         struct hit * hits,
                         int hitcount,
                         char * query_head,
                         char * qsequence,
                         int64_t qseqlen) -> void
{
  /* http://drive5.com/usearch/manual/alnout.html */

  if (hitcount > 0)
    {
      fprintf(output_handle, "\n");

      fprintf(output_handle,"Query >%s\n", query_head);
      fprintf(output_handle," %%Id   TLen  Target\n");

      auto const top_hit_id = hits[0].id;

      for (auto t = 0; t < hitcount; t++)
        {
          auto * hp = hits + t;

          if (opt_top_hits_only && (hp->id < top_hit_id))
            {
              break;
            }

          fprintf(output_handle,"%3.0f%% %6" PRIu64 "  %s\n",
                  hp->id,
                  db_getsequencelen(hp->target),
                  db_getheader(hp->target));
        }

      for (auto t = 0; t < hitcount; t++)
        {
          auto * hp = hits + t;

          if (opt_top_hits_only && (hp->id < top_hit_id))
            {
              break;
            }

          fprintf(output_handle,"\n");


          auto * dseq = db_getsequence(hp->target);
          int64_t const dseqlen = db_getsequencelen(hp->target);

          auto const qlenlen = snprintf(nullptr, 0, "%" PRId64, qseqlen);
          auto const tlenlen = snprintf(nullptr, 0, "%" PRId64, dseqlen);
          auto const numwidth = std::max(qlenlen, tlenlen);

          fprintf(output_handle," Query %*" PRId64 "nt >%s\n", numwidth,
                  qseqlen, query_head);
          fprintf(output_handle,"Target %*" PRId64 "nt >%s\n", numwidth,
                  dseqlen, db_getheader(hp->target));

          int const rowlen = opt_rowlen == 0 ? qseqlen + dseqlen : opt_rowlen;

          align_show(output_handle,
                     qsequence,
                     qseqlen,
                     hp->trim_q_left,
                     "Qry",
                     dseq,
                     dseqlen,
                     hp->trim_t_left,
                     "Tgt",
                     hp->nwalignment + hp->trim_aln_left,
                     strlen(hp->nwalignment)
                     - hp->trim_aln_left - hp->trim_aln_right,
                     numwidth,
                     3,
                     rowlen,
                     hp->strand);

          fprintf(output_handle, "\n%d cols, %d ids (%3.1f%%), %d gaps (%3.1f%%)\n",
                  hp->internal_alignmentlength,
                  hp->matches,
                  hp->id,
                  hp->internal_indels,
                  hp->internal_alignmentlength > 0 ?
                  100.0 * hp->internal_indels / hp->internal_alignmentlength :
                  0.0);

#if 0
          fprintf(output_handle, "%d kmers, %d score, %d gap opens. %s %s %d %d %d %d %d\n",
                  hp->count, hp->nwscore, hp->nwgaps,
                  hp->accepted ? "accepted" : "not accepted",
                  hp->nwalignment, hp->nwalignmentlength,
                  hp->trim_q_left, hp->trim_q_right,
                  hp->trim_t_left, hp->trim_t_right
                  );
#endif
        }
    }
  else if (opt_output_no_hits)
    {
      fprintf(output_handle, "\n");
      fprintf(output_handle,"Query >%s\n", query_head);
      fprintf(output_handle,"No hits\n");
    }
}


auto inline nucleotide_equal(char lhs, char rhs) -> bool
{
  return chrmap_4bit[(int) lhs] == chrmap_4bit[(int) rhs];
}


auto build_sam_strings(char * alignment,
                       char * queryseq,
                       char * targetseq,
                       xstring * cigar,
                       xstring * md) -> void
{
  /*
    convert cigar to sam format:
    add "1" to operations without run length
    flip direction of indels in cigar string

    build MD-string with substitutions
  */

  cigar->empty();
  md->empty();

  char * p = alignment;
  char * e = p + strlen(p);

  auto qpos = 0;
  auto tpos = 0;

  auto matched = 0;
  auto flag = false; /* 1: MD string ends with a number */

  while(p < e)
    {
      auto run = 1;
      auto scanned = 0;
      sscanf(p, "%d%n", &run, &scanned);
      p += scanned;
      char const op = *p++;

      switch (op)
        {
        case 'M':
          cigar->add_d(run);
          cigar->add_c('M');

          for (int i = 0; i < run; i++)
            {
              if (nucleotide_equal(queryseq[qpos], targetseq[tpos]))
                {
                  matched++;
                }
              else
                {
                  if (! flag)
                    {
                      md->add_d(matched);
                      matched = 0;
                      flag = true;
                    }

                  md->add_c(targetseq[tpos]);
                  flag = false;
                }
              qpos++;
              tpos++;
            }

          break;

        case 'D':
          cigar->add_d(run);
          cigar->add_c('I');
          qpos += run;
          break;

        case 'I':
          cigar->add_d(run);
          cigar->add_c('D');

          if (! flag)
            {
              md->add_d(matched);
              matched = 0;
              flag = true;
            }

          md->add_c('^');
          for (auto i = 0; i < run; i++)
            {
              md->add_c(targetseq[tpos++]);
            }
          flag = false;
          break;
        }
    }

  if (! flag)
    {
      md->add_d(matched);
      matched = 0;
      flag = true;
    }
}

auto results_show_samheader(std::FILE * output_handle,
                            char * cmdline,
                            char * dbname) -> void
{
  if (opt_samout && opt_samheader)
    {
      fprintf(output_handle, "@HD\tVN:1.0\tSO:unsorted\tGO:query\n");

      for (uint64_t i = 0; i < db_getsequencecount(); i++)
        {
          char md5hex[len_hex_dig_md5];
          get_hex_seq_digest_md5(md5hex,
                                 db_getsequence(i),
                                 db_getsequencelen(i));
          fprintf(output_handle,
                  "@SQ\tSN:%s\tLN:%" PRIu64 "\tM5:%s\tUR:file:%s\n",
                  db_getheader(i),
                  db_getsequencelen(i),
                  md5hex,
                  dbname);
        }

      fprintf(output_handle,
              "@PG\tID:%s\tVN:%s\tCL:%s\n",
              PROG_NAME,
              PROG_VERSION,
              cmdline);
    }
}

auto results_show_samout(std::FILE * output_handle,
                         struct hit * hits,
                         int hitcount,
                         char * query_head,
                         char * qsequence,
                         char * qsequence_rc) -> void
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

  if (hitcount > 0)
    {
      auto const top_hit_id = hits[0].id;

      for (auto t = 0; t < hitcount; t++)
        {
          auto * hp = hits + t;

          if (opt_top_hits_only && (hp->id < top_hit_id))
            {
              break;
            }

          /*

           */

          xstring cigar;
          xstring md;

          build_sam_strings(hp->nwalignment,
                            hp->strand ? qsequence_rc : qsequence,
                            db_getsequence(hp->target),
                            &cigar,
                            &md);

          fprintf(output_handle,
                  "%s\t%u\t%s\t%" PRIu64
                  "\t%u\t%s\t%s\t%" PRIu64
                  "\t%" PRIu64
                  "\t%s\t%s\t"
                  "AS:i:%.0f\tXN:i:%d\tXM:i:%d\tXO:i:%d\t"
                  "XG:i:%d\tNM:i:%d\tMD:Z:%s\tYT:Z:%s\n",
                  query_head,
                  (0x10 * hp->strand) | (t > 0 ? 0x100 : 0),
                  db_getheader(hp->target),
                  (uint64_t) 1,
                  255,
                  cigar.get_string(),
                  "*",
                  (uint64_t) 0,
                  (uint64_t) 0,
                  hp->strand ? qsequence_rc : qsequence,
                  "*",
                  hp->id,
                  0,
                  hp->mismatches,
                  hp->internal_gaps,
                  hp->internal_indels,
                  hp->mismatches + hp->internal_indels,
                  md.get_string(),
                  "UU");
        }
    }
  else if (opt_output_no_hits)
    {
      fprintf(output_handle,
              "%s\t%u\t%s\t%" PRIu64 "\t%u\t%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t%s\t%s\n",
              query_head,
              0x04,
              "*",
              (uint64_t) 0,
              255,
              "*",
              "*",
              (uint64_t) 0,
              (uint64_t) 0,
              qsequence,
              "*");
    }
}
