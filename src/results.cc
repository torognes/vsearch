/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

void results_show_fastapairs_one(FILE * fp,
                                 struct hit * hp,
                                 char * query_head,
                                 char * qsequence,
                                 int64_t qseqlen,
                                 char * rc)
{
  /* http://www.drive5.com/usearch/manual/fastapairs.html */
  
  if (hp)
    {
      char * qrow = align_getrow(hp->strand ? rc : qsequence,
                                 hp->nwalignment,
                                 hp->nwalignmentlength,
                                 0);
      fasta_print_general(fp,
                          0,
                          qrow + hp->trim_q_left + hp->trim_t_left,
                          hp->internal_alignmentlength,
                          query_head,
                          strlen(query_head),
                          0,
                          0,
                          -1,
                          -1,
                          0,
                          0.0);
      xfree(qrow);

      char * trow = align_getrow(db_getsequence(hp->target),
                                 hp->nwalignment,
                                 hp->nwalignmentlength,
                                 1);
      fasta_print_general(fp,
                          0,
                          trow + hp->trim_q_left + hp->trim_t_left,
                          hp->internal_alignmentlength,
                          db_getheader(hp->target),
                          db_getheaderlen(hp->target),
                          0,
                          0,
                          -1,
                          -1,
                          0,
                          0.0);
      xfree(trow);

      fprintf(fp, "\n");
    }
}


void results_show_blast6out_one(FILE * fp,
                                struct hit * hp,
                                char * query_head,
                                char * qsequence, 
                                int64_t qseqlen,
                                char * rc)
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
  
  if (hp)
    {
      int qstart, qend;
      
      if (hp->strand)
        {
          /* minus strand */
          qstart = qseqlen;
          qend = 1;
        }
      else
        {
          /* plus strand */
          qstart = 1;
          qend = qseqlen;
        }
      
      fprintf(fp,
              "%s\t%s\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%" PRIu64 "\t%d\t%d\n",
              query_head,
              db_getheader(hp->target),
              hp->id,
              hp->internal_alignmentlength,
              hp->mismatches,
              hp->internal_gaps,
              qstart,
              qend,
              1,
              db_getsequencelen(hp->target),
              -1,
              0);
    }
  else
    {
        fprintf(fp, "%s\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n", query_head);
    }
}

void results_show_uc_one(FILE * fp,
                         struct hit * hp,
                         char * query_head,
                         char * qsequence,
                         int64_t qseqlen,
                         char * rc,
                         int clusterno)
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
    compressed alignment, e.g. 9I92M14D, or "=" if prefect alignment
    query label
    target label
  */

  if (hp)
    {
      bool perfect = (hp->matches == qseqlen) &&
        ((uint64_t)(qseqlen) == db_getsequencelen(hp->target));

      fprintf(fp,
              "H\t%d\t%" PRId64 "\t%.1f\t%c\t0\t0\t%s\t%s\t%s\n",
              clusterno,
              qseqlen,
              hp->id,
              hp->strand ? '-' : '+',
              perfect ? "=" : hp->nwalignment,
              query_head,
              db_getheader(hp->target));
    }
  else
    fprintf(fp, "N\t*\t*\t*\t.\t*\t*\t*\t%s\t*\n", query_head);
}

void results_show_userout_one(FILE * fp, struct hit * hp,
                              char * query_head,
                              char * qsequence, int64_t qseqlen,
                              char * rc)
{

  /*
    http://drive5.com/usearch/manual/userout.html
    qlo, qhi, tlo, thi and raw are given more meaningful values here
  */

  for (int c = 0; c < userfields_requested_count; c++)
    {
      if (c)
        fprintf(fp, "\t");

      int field = userfields_requested[c];
          
      char * tsequence = 0;
      int64_t tseqlen = 0;
      char * t_head = 0;

      if (hp)
        {
          tsequence = db_getsequence(hp->target);
          tseqlen = db_getsequencelen(hp->target);
          t_head = db_getheader(hp->target);
        }

      char * qrow;
      char * trow;

      switch (field)
        {
        case 0: /* query */
          fprintf(fp, "%s", query_head);
          break;
        case 1: /* target */
          fprintf(fp, "%s", hp ? t_head : "*");
          break;
        case 2: /* evalue */
          fprintf(fp, "-1");
          break;
        case 3: /* id */
          fprintf(fp, "%.1f", hp ? hp->id : 0.0);
          break;
        case 4: /* pctpv */
          fprintf(fp, "%.1f", (hp && (hp->internal_alignmentlength > 0)) ? 100.0 * hp->matches / hp->internal_alignmentlength : 0.0);
          break;
        case 5: /* pctgaps */
          fprintf(fp, "%.1f", (hp && (hp->internal_alignmentlength > 0)) ? 100.0 * hp->internal_indels / hp->internal_alignmentlength : 0.0);
          break;
        case 6: /* pairs */
          fprintf(fp, "%d", hp ? hp->matches + hp->mismatches : 0);
          break;
        case 7: /* gaps */
          fprintf(fp, "%d", hp ? hp->internal_indels : 0);
          break;
        case 8: /* qlo */
          fprintf(fp, "%" PRId64, hp ? (hp->strand ? qseqlen : 1) : 0);
          break;
        case 9: /* qhi */
          fprintf(fp, "%" PRId64, hp ? (hp->strand ? 1 : qseqlen) : 0);
          break;
        case 10: /* tlo */
          fprintf(fp, "%d", hp ? 1 : 0);
          break;
        case 11: /* thi */
          fprintf(fp, "%" PRId64, tseqlen);
          break;
        case 12: /* pv */
          fprintf(fp, "%d", hp ? hp->matches : 0);
          break;
        case 13: /* ql */
          fprintf(fp, "%" PRId64, qseqlen);
          break;
        case 14: /* tl */
          fprintf(fp, "%" PRId64, hp ? tseqlen : 0);
          break;
        case 15: /* qs */
          fprintf(fp, "%" PRId64, qseqlen);
          break;
        case 16: /* ts */
          fprintf(fp, "%" PRId64, hp ? tseqlen : 0);
          break;
        case 17: /* alnlen */
          fprintf(fp, "%d", hp ? hp->internal_alignmentlength : 0);
          break;
        case 18: /* opens */
          fprintf(fp, "%d", hp ? hp->internal_gaps : 0);
          break;
        case 19: /* exts */
          fprintf(fp, "%d", hp ? hp->internal_indels - hp->internal_gaps : 0);
          break;
        case 20: /* raw */
          fprintf(fp, "%d", hp ? hp->nwscore : 0);
          break;
        case 21: /* bits */
          fprintf(fp, "%d", 0);
          break;
        case 22: /* aln */ 
          if (hp)
            align_fprint_uncompressed_alignment(fp, hp->nwalignment);
          break;
        case 23: /* caln */
          if (hp)
            fprintf(fp, "%s", hp->nwalignment);
          break;
        case 24: /* qstrand */
          if (hp)
            fprintf(fp, "%c", hp->strand ? '-' : '+');
          break;
        case 25: /* tstrand */
          if (hp)
            fprintf(fp, "%c", '+');
          break;
        case 26: /* qrow */
          if (hp)
            {
              qrow = align_getrow(hp->strand ? rc : qsequence,
                                  hp->nwalignment,
                                  hp->nwalignmentlength,
                                  0);
              fprintf(fp, "%.*s",
                      (int)(hp->internal_alignmentlength),
                      qrow + hp->trim_q_left + hp->trim_t_left);
              xfree(qrow);
            }
          break;
        case 27: /* trow */
          if (hp)
            {
              trow = align_getrow(tsequence,
                                  hp->nwalignment,
                                  hp->nwalignmentlength,
                                  1);
              fprintf(fp, "%.*s",
                      (int)(hp->internal_alignmentlength),
                      trow + hp->trim_q_left + hp->trim_t_left);
              xfree(trow);
            }
          break;
        case 28: /* qframe */
          fprintf(fp, "+0");
          break;
        case 29: /* tframe */
          fprintf(fp, "+0");
          break;
        case 30: /* mism */
          fprintf(fp, "%d", hp ? hp->mismatches : 0);
          break;
        case 31: /* ids */
          fprintf(fp, "%d", hp ? hp->matches : 0);
          break;
        case 32: /* qcov */
          fprintf(fp, "%.1f",
                  hp ? 100.0 * (hp->matches + hp->mismatches) / qseqlen : 0.0);
          break;
        case 33: /* tcov */
          fprintf(fp, "%.1f",
                  hp ? 100.0 * (hp->matches + hp->mismatches) / tseqlen : 0.0);
          break;
        case 34: /* id0 */
          fprintf(fp, "%.1f", hp ? hp->id0 : 0.0);
          break;
        case 35: /* id1 */
          fprintf(fp, "%.1f", hp ? hp->id1 : 0.0);
          break;
        case 36: /* id2 */
          fprintf(fp, "%.1f", hp ? hp->id2 : 0.0);
          break;
        case 37: /* id3 */
          fprintf(fp, "%.1f", hp ? hp->id3 : 0.0);
          break;
        case 38: /* id4 */
          fprintf(fp, "%.1f", hp ? hp->id4 : 0.0);
          break;

          /* new internal alignment coordinates */

        case 39: /* qilo */
          fprintf(fp, "%d", hp ? hp->trim_q_left + 1 : 0);
          break;
        case 40: /* qihi */
          fprintf(fp, "%" PRId64, hp ? qseqlen - hp->trim_q_right : 0);
          break;
        case 41: /* tilo */
          fprintf(fp, "%d", hp ? hp->trim_t_left + 1 : 0);
          break;
        case 42: /* tihi */
          fprintf(fp, "%" PRId64, hp ? tseqlen - hp->trim_t_right : 0);
          break;
        }
    }
  fprintf(fp, "\n");
}

void results_show_alnout(FILE * fp,
                         struct hit * hits,
                         int hitcount,
                         char * query_head,
                         char * qsequence,
                         int64_t qseqlen,
                         char * rc)
{
  /* http://drive5.com/usearch/manual/alnout.html */

  if (hitcount)
    {
      fprintf(fp, "\n");

      fprintf(fp,"Query >%s\n", query_head);
      fprintf(fp," %%Id   TLen  Target\n");
      
      double top_hit_id = hits[0].id;

      for(int t = 0; t < hitcount; t++)
        {
          struct hit * hp = hits + t;
          
          if (opt_top_hits_only && (hp->id < top_hit_id))
            break;

          fprintf(fp,"%3.0f%% %6" PRIu64 "  %s\n",
                  hp->id,
                  db_getsequencelen(hp->target),
                  db_getheader(hp->target));
        }

      for(int t = 0; t < hitcount; t++)
        {
          struct hit * hp = hits + t;
          
          if (opt_top_hits_only && (hp->id < top_hit_id))
            break;

          fprintf(fp,"\n");
          

          char * dseq = db_getsequence(hp->target);
          int64_t dseqlen = db_getsequencelen(hp->target);
          
          int qlenlen = snprintf(0, 0, "%" PRId64, qseqlen);
          int tlenlen = snprintf(0, 0, "%" PRId64, dseqlen);
          int numwidth = MAX(qlenlen, tlenlen);
          
          fprintf(fp," Query %*" PRId64 "nt >%s\n", numwidth,
                  qseqlen, query_head);
          fprintf(fp,"Target %*" PRId64 "nt >%s\n", numwidth,
                  dseqlen, db_getheader(hp->target));
          
          int rowlen = opt_rowlen == 0 ? qseqlen+dseqlen : opt_rowlen;

          align_show(fp,
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
          
          fprintf(fp, "\n%d cols, %d ids (%3.1f%%), %d gaps (%3.1f%%)\n",
                  hp->internal_alignmentlength,
                  hp->matches,
                  hp->id,
                  hp->internal_indels,
                  hp->internal_alignmentlength > 0 ?
                  100.0 * hp->internal_indels / hp->internal_alignmentlength :
                  0.0);

#if 0
          fprintf(fp, "%d kmers, %d score, %d gap opens. %s %s %d %d %d %d %d\n",
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
      fprintf(fp, "\n");
      fprintf(fp,"Query >%s\n", query_head);
      fprintf(fp,"No hits\n");
    }
}

bool inline nucleotide_equal(char a, char b)
{
  return chrmap_4bit[(int)a] == chrmap_4bit[(int)b];
}

void build_sam_strings(char * alignment,
                       char * queryseq,
                       char * targetseq,
                       xstring * cigar,
                       xstring * md)
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

  int qpos = 0;
  int tpos = 0;

  int matched = 0;
  bool flag = 0; /* 1: MD string ends with a number */

  while(p < e)
    {
      int run = 1;
      int scanned = 0;
      sscanf(p, "%d%n", & run, & scanned);
      p += scanned;
      char op = *p++;

      switch (op)
        {
        case 'M':
          cigar->add_d(run);
          cigar->add_c('M');

          for(int i=0; i<run; i++)
            {
              if (nucleotide_equal(queryseq[qpos], targetseq[tpos]))
                matched++;
              else
                {
                  if (!flag)
                    {
                      md->add_d(matched);
                      matched = 0;
                      flag = 1;
                    }

                  md->add_c(targetseq[tpos]);
                  flag = 0;
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

          if (!flag)
            {
              md->add_d(matched);
              matched = 0;
              flag = 1;
            }
          
          md->add_c('^');
          for(int i=0; i<run; i++)
            md->add_c(targetseq[tpos++]);
          flag = 0;
          break;
        }
    }

  if (!flag)
    {
      md->add_d(matched);
      matched = 0;
      flag = 1;
    }
}

void results_show_samheader(FILE * fp,
                            char * cmdline,
                            char * dbname)
{
  if (opt_samout && opt_samheader)
    {
      fprintf(fp, "@HD\tVN:1.0\tSO:unsorted\tGO:query\n");
      
      for(uint64_t i=0; i<db_getsequencecount(); i++)
        {
          char md5hex[LEN_HEX_DIG_MD5];
          get_hex_seq_digest_md5(md5hex,
                                 db_getsequence(i),
                                 db_getsequencelen(i));
          fprintf(fp,
                  "@SQ\tSN:%s\tLN:%" PRIu64 "\tM5:%s\tUR:file:%s\n",
                  db_getheader(i),
                  db_getsequencelen(i),
                  md5hex,
                  dbname);
        }

      fprintf(fp,
              "@PG\tID:%s\tVN:%s\tCL:%s\n",
              PROG_NAME,
              PROG_VERSION,
              cmdline);
    }
}

void results_show_samout(FILE * fp,
                         struct hit * hits,
                         int hitcount,
                         char * query_head,
                         char * qsequence,
                         int64_t qseqlen,
                         char * rc)
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
      double top_hit_id = hits[0].id;
      
      for(int t = 0; t < hitcount; t++)
        {
          struct hit * hp = hits + t;
          
          if (opt_top_hits_only && (hp->id < top_hit_id))
            break;
          
          /*

          */

          xstring cigar;
          xstring md;

          build_sam_strings(hp->nwalignment,
                            hp->strand ? rc : qsequence,
                            db_getsequence(hp->target),
                            & cigar,
                            & md);

          fprintf(fp,
                  "%s\t%u\t%s\t%" PRIu64
                  "\t%u\t%s\t%s\t%" PRIu64
                  "\t%" PRIu64
                  "\t%s\t%s\t"
                  "AS:i:%.0f\tXN:i:%d\tXM:i:%d\tXO:i:%d\t"
                  "XG:i:%d\tNM:i:%d\tMD:Z:%s\tYT:Z:%s\n",
                  query_head,
                  0x10 * hp->strand | (t>0 ? 0x100 : 0),
                  db_getheader(hp->target),
                  (uint64_t) 1,
                  255,
                  cigar.get_string(),
                  "*",
                  (uint64_t) 0,
                  (uint64_t) 0,
                  hp->strand ? rc : qsequence,
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
      fprintf(fp,
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
