/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "pvalue.h"

/* Use default PEAR scoring method and statistics */
const int score_method = 2;
const double pvalue_level = 0.01;
const double alpha = 1.0;
const double beta = -1.0;

/* 
   TODO:
   - Parallelize with pthreads
   - Check matching labels of forward and reverse reads (/1 and /2)
*/

/* static variables, needs lock */

static FILE * fp_fastqout = 0;
static FILE * fp_fastaout = 0;
static FILE * fp_fastqout_notmerged_fwd = 0;
static FILE * fp_fastqout_notmerged_rev = 0;
static FILE * fp_fastaout_notmerged_fwd = 0;
static FILE * fp_fastaout_notmerged_rev = 0;
static FILE * fp_eetabbedout = 0;

static long merged = 0;
static long notmerged = 0;
static long total = 0;

FILE * fileopenw(char * filename)
{
  FILE * fp = 0;
  fp = fopen(filename, "w");
  if (!fp)
    fatal("Unable to open file for writing (%s)", filename);
  return fp;
}

int get_qual(char q)
{
  int qual = q - opt_fastq_ascii;
  char msg[200];

  if (qual < opt_fastq_qmin)
    {
      snprintf(msg, 200, "FASTQ quality value (%d) below qmin (%ld)",
               qual, opt_fastq_qmin);
      fatal(msg);
    }
  else if (qual > opt_fastq_qmax)
    {
      snprintf(msg, 200, "FASTQ quality value (%d) above qmax (%ld)",
               qual, opt_fastq_qmax);
      fatal(msg);
    }
  return qual;
}

static char merge_qual_same[128][128];
static char merge_qual_diff[128][128];
static double match_oes[128][128];
static double match_score[128][128];
static double mism_oes[128][128];
static double mism_score[128][128];
static double q2p[128];
static double ambig_score;
static double ambig_oes;

double q_to_p(double q)
{
  int x = q - opt_fastq_ascii;
  if (x < 2)
    return 0.75;
  else
    return exp10(-x/10.0);
}

void precompute_qual()
{
  /* Precompute tables of scores etc */

  ambig_oes = alpha * 0.25 + beta * 0.75;

  switch (score_method)
    {
    case 1:
      ambig_score = alpha * 0.25 + beta * 0.75;
      break;
    case 2:
      ambig_score = beta * 0.75;
      break;
    case 3:
      ambig_score = beta;
      break;
    }

  for (int x = 33; x < 126; x++)
    {
      double px = q_to_p(x);
      q2p[x] = px;

      for (int y = 33; y < 126; y++)
        {
          double py = q_to_p(y);

          double p, q, oes, score;
          
          /* Quality score equations from Edgar & Flyvbjerg (2015) */

          /* Match */
          p = px * py / 3.0 / (1.0 - px - py + 4.0 * px * py / 3.0);
          q = opt_fastq_ascii + MIN(round(-10.0*log10(p)), opt_fastq_qmaxout);
          merge_qual_same[x][y] = q;
          
          /* Mismatch, x is highest quality */
          p = px * (1.0 - py / 3.0) / (px + py - 4.0 * px * py / 3.0);
          q = opt_fastq_ascii + MIN(round(-10.0*log10(p)), opt_fastq_qmaxout);
          merge_qual_diff[x][y] = q;
          
          /* Score weights from PEAR */

          /* Match */

          /* probability that they really are similar,
             given that they look similar and have
             error probabilites of px and py, resp. */

          p = 1.0 - px - py + px * py * 4.0 / 3.0;

          oes = alpha * p + beta * (1.0 - p);

          switch (score_method)
            {
            case 1:
              score = alpha * p + beta * (1.0 - p);
              break;
            case 2:
              score = alpha * p;
              break;
            case 3:
              score = alpha;
              break;
            }
          
          match_oes[x][y] = oes;
          match_score[x][y] = score;

          /* Mismatch */

          /* Probability that they really are different,
             given that they look different and have
             error probabilities of px and py, resp. */

          p = 1.0 - px / 3.0 - py / 3.0 + px * py * 4.0 / 9.0;
          
          oes = (beta - alpha) * p + alpha;
          
          switch (score_method)
            {
            case 1:
              score = alpha * (1.0 - p) + beta * p;
              break;
            case 2:
              score = beta * p;
              break;
            case 3:
              score = beta;
              break;
            }

          mism_oes[x][y] = oes;
          mism_score[x][y] = score;
        }
    }
}

void merge_sym(char * sym,       char * qual,
               char fwd_sym,     char rev_sym,
               char fwd_qual,    char rev_qual)
{
  if (rev_sym == 'N')
    {
      * sym = fwd_sym;
      * qual = fwd_qual;
    }
  else if (fwd_sym == 'N')
    {
      * sym = rev_sym;
      * qual = rev_qual;
    }
  else if (fwd_sym == rev_sym)
    {
      /* agreement */
      * sym = fwd_sym;
      * qual = merge_qual_same[(unsigned)fwd_qual][(unsigned)rev_qual];
    }
  else
    {
      /* disagreement */
      if (fwd_qual > rev_qual)
        {
          * sym = fwd_sym;
          * qual = merge_qual_diff[(unsigned)fwd_qual][(unsigned)rev_qual];
        }
      else
        {
          * sym = rev_sym;
          * qual = merge_qual_diff[(unsigned)rev_qual][(unsigned)fwd_qual];
        }
    }
}

static char * merged_sequence = 0;
static char * merged_quality = 0;
static char * merged_header = 0;

static size_t merged_sequence_alloc = 0;
static size_t merged_header_alloc = 0;

void keep(char * fwd_header,   char * rev_header,
          char * fwd_sequence, char * rev_sequence,
          char * fwd_quality,  char * rev_quality,
          long fwd_trunc,      long rev_trunc,
          long offset)
{
  /* The offset is the distance between the (truncated) 3' ends of the two
     sequences */

  long rev_3prime_overhang = offset > fwd_trunc ? offset - fwd_trunc : 0;
  long fwd_5prime_overhang = fwd_trunc > offset ? fwd_trunc - offset : 0;
  long mergelen = fwd_trunc + rev_trunc - offset;

  if (opt_label_suffix)
    {
      size_t need = strlen(fwd_header) + strlen(opt_label_suffix) + 1;
      if (need > merged_header_alloc)
        {
          merged_header_alloc = need;
          merged_header = (char*)xrealloc(merged_header, merged_header_alloc);
        }
    }
  
  if ((size_t)(mergelen + 1) > merged_sequence_alloc)
    {
      merged_sequence_alloc = mergelen + 1;
      merged_sequence = (char*)xrealloc(merged_sequence,merged_sequence_alloc);
      merged_quality  = (char*)xrealloc(merged_quality,merged_sequence_alloc);
    }

  char * merged_sequence = (char *) xmalloc(mergelen + 1);
  char * merged_quality = (char *) xmalloc(mergelen + 1);

  double ee = 0.0;

  long fwd_pos = 0;
  long rev_pos = rev_trunc - 1 + fwd_5prime_overhang - rev_3prime_overhang;

  long fwd_errors = 0;
  long rev_errors = 0;

  for(long i = 0; i < mergelen; i++)
    {
      bool has_fwd = 0;
      if ((fwd_pos >= 0) && (fwd_pos < fwd_trunc))
        has_fwd = 1;

      bool has_rev = 0;
      if ((rev_pos >= 0) && (rev_pos < rev_trunc))
        has_rev = 1;

      char sym;
      char qual;
      
      if (has_fwd && has_rev)
        {
          char fwd_sym = fwd_sequence[fwd_pos];
          char rev_sym = chrmap_complement[(int)(rev_sequence[rev_pos])];
          char fwd_qual = fwd_quality[fwd_pos];
          char rev_qual = rev_quality[rev_pos];

          merge_sym(& sym,
                    & qual,
                    fwd_qual < 2 ? 'N' : fwd_sym,
                    rev_qual < 2 ? 'N' : rev_sym,
                    fwd_qual,
                    rev_qual);

          if (sym != fwd_sym)
            fwd_errors++;
          if (sym != rev_sym)
            rev_errors++;
        }
      else if (has_fwd)
        {
          sym = fwd_sequence[fwd_pos];
          qual = fwd_quality[fwd_pos];
        }
      else
        {
          sym = chrmap_complement[(int)(rev_sequence[rev_pos])];
          qual = rev_quality[rev_pos];
        }

      merged_sequence[i] = sym;
      merged_quality[i] = qual;

      ee += q2p[(unsigned)qual];

      fwd_pos++;
      rev_pos--;
    }

  merged_sequence[mergelen] = 0;
  merged_quality[mergelen] = 0;

  if (ee <= opt_fastq_maxee)
    {
      merged++;

      if (opt_label_suffix)
        (void) sprintf(merged_header, "%s%s", fwd_header, opt_label_suffix);
      
      if (opt_fastqout)
        {
          if (opt_fastq_eeout)
            fastq_print_with_ee(fp_fastqout,
                                opt_label_suffix ? merged_header : fwd_header,
                                merged_sequence,
                                merged_quality,
                                ee);
          else
            fastq_print(fp_fastqout,
                        opt_label_suffix ? merged_header : fwd_header,
                        merged_sequence,
                        merged_quality);
        }

      if (opt_fastaout)
        fasta_print(fp_fastaout,
                    opt_label_suffix ? merged_header : fwd_header,
                    merged_sequence,
                    mergelen);

      if (opt_eetabbedout)
        {
          double ee_fwd = 0.0;
          for(int i=0; i<fwd_trunc; i++)
            ee_fwd += q2p[(unsigned)fwd_quality[i]];
          
          double ee_rev = 0.0;
          for(int i=0; i<rev_trunc; i++)
            ee_rev += q2p[(unsigned)rev_quality[rev_trunc-i-1]];
          
          fprintf(fp_eetabbedout, "%.2lf\t%.2lf\t%ld\t%ld\n",
                  ee_fwd, ee_rev, fwd_errors, rev_errors);
        }
    }
  else
    {
      notmerged++;
    }

  free(merged_sequence);
  merged_sequence = 0;
  free(merged_quality);
  merged_quality = 0;
}

void discard(char * fwd_header,   char * rev_header,
             char * fwd_sequence, char * rev_sequence,
             char * fwd_quality,  char * rev_quality,
             long fwd_length,     long rev_length)
{
  notmerged++;

  if (opt_fastqout_notmerged_fwd)
    fastq_print(fp_fastqout_notmerged_fwd,
                fwd_header,
                fwd_sequence,
                fwd_quality);

  if (opt_fastqout_notmerged_rev)
    fastq_print(fp_fastqout_notmerged_rev,
                rev_header,
                rev_sequence,
                rev_quality);

  if (opt_fastaout_notmerged_fwd)
    fasta_print(fp_fastaout_notmerged_fwd,
                fwd_header,
                fwd_sequence,
                fwd_length);

  if (opt_fastaout_notmerged_rev)
    fasta_print(fp_fastaout_notmerged_rev,
                rev_header,
                rev_sequence,
                rev_length);
}


double overlap_score(char * fwd_sequence, char * rev_sequence,
                     char * fwd_quality,  char * rev_quality,
                     long fwd_pos_start, long rev_pos_start,
                     long overlap,
                     long * pdiffs)
{
  long fwd_pos = fwd_pos_start;
  long rev_pos = rev_pos_start;
  
  double score = 0.0;
  long diffs = 0;

  for (long j=0; j < overlap; j++)
    {
      char fwd_sym = fwd_sequence[fwd_pos];
      char rev_sym = chrmap_complement[(int)(rev_sequence[rev_pos])];
      unsigned int fwd_qual = fwd_quality[fwd_pos];
      unsigned int rev_qual = rev_quality[rev_pos];
      fwd_pos--;
      rev_pos++;
      
      if (fwd_sym == rev_sym)
        {
          score += match_score[fwd_qual][rev_qual];
        }
      else
        {
          diffs++;
          score += mism_score[fwd_qual][rev_qual];
        }
    }

  *pdiffs = diffs;
  return score;
}

double overlap_oes(char * fwd_sequence, char * rev_sequence,
                   char * fwd_quality,  char * rev_quality,
                   long fwd_pos_start, long rev_pos_start,
                   long overlap)
{
  long fwd_pos = fwd_pos_start;
  long rev_pos = rev_pos_start;
  double oes = 0.0;

  for (long j=0; j < overlap; j++)
    {
      char fwd_sym = fwd_sequence[fwd_pos];
      char rev_sym = chrmap_complement[(int)(rev_sequence[rev_pos])];
      unsigned int fwd_qual = fwd_quality[fwd_pos];
      unsigned int rev_qual = rev_quality[rev_pos];
      fwd_pos--;
      rev_pos++;

      if (fwd_sym == rev_sym)
        oes += match_oes[fwd_qual][rev_qual];
      else
        oes += mism_oes[fwd_qual][rev_qual];
    }

  return oes;
}

long merge(char * fwd_sequence, char * rev_sequence,
           char * fwd_quality,  char * rev_quality,
           long fwd_trunc, long rev_trunc)
{
  long i1 = opt_fastq_minovlen;

  i1 = MAX(i1, fwd_trunc + rev_trunc - opt_fastq_maxmergelen);

  long i2 = opt_fastq_allowmergestagger ?
    (fwd_trunc + rev_trunc - opt_fastq_minovlen) :
    fwd_trunc;

  i2 = MIN(i2, fwd_trunc + rev_trunc - opt_fastq_minmergelen);

  long best_i = 0;
  double best_score = 0.0;
  double best_oes = 0.0;

  for(long i = i1; i <= i2; i++)
    {
      long fwd_3prime_overhang = i > rev_trunc ? i - rev_trunc : 0;
      long rev_3prime_overhang = i > fwd_trunc ? i - fwd_trunc : 0;
      long overlap = i - fwd_3prime_overhang - rev_3prime_overhang;
      long diffs = 0;
      long fwd_pos_start = fwd_trunc - fwd_3prime_overhang - 1;
      long rev_pos_start = rev_trunc - rev_3prime_overhang - overlap;
      
      double score = overlap_score(fwd_sequence, rev_sequence,
                                   fwd_quality, rev_quality,
                                   fwd_pos_start, rev_pos_start,
                                   overlap,
                                   & diffs);
      
      if ((diffs <= opt_fastq_maxdiffs) && (score > best_score))
        {
          best_score = score;
          best_i = i;
        }
    }

  /* compute oes for best alignment */

  if (best_score > 0.0)
    {
      long i = best_i;
      long fwd_3prime_overhang = i > rev_trunc ? i - rev_trunc : 0;
      long rev_3prime_overhang = i > fwd_trunc ? i - fwd_trunc : 0;
      long overlap = i - fwd_3prime_overhang - rev_3prime_overhang;
      long fwd_pos_start = fwd_trunc - fwd_3prime_overhang - 1;
      long rev_pos_start = rev_trunc - rev_3prime_overhang - overlap;
      
      best_oes = overlap_oes(fwd_sequence, rev_sequence,
                             fwd_quality, rev_quality,
                             fwd_pos_start, rev_pos_start,
                             overlap);
    }

  /* Statistical test as in PEAR */

  const int basefreqpct = (int) (0.25 * 100.0);

  int overlapregion = opt_fastq_minovlen;
  if (overlapregion > 99)
    overlapregion = 99;
  if (overlapregion < 1)
    overlapregion = 1;

  double cutoff;

  if (pvalue_level == 1.0)
    {
      cutoff = DBL_MIN;
    }
  else if (pvalue_level == 0.01)
    {
      cutoff = precomp2_01[overlapregion][basefreqpct];
    }
  else if (pvalue_level == 0.05)
    {
      cutoff = precomp2_05[overlapregion][basefreqpct];
    }
  else if (pvalue_level == 0.001)
    {
      cutoff = precomp2_001[overlapregion][basefreqpct];
    }
  else
    {
      cutoff = precomp2_0001[overlapregion][basefreqpct];
    }

  if (best_oes <= cutoff)
    best_i = 0;

  return best_i;
}

void fastq_mergepairs()
{
  /* open input files */
  
  fastq_handle fastq_fwd = fastq_open(opt_fastq_mergepairs);
  fastq_handle fastq_rev = fastq_open(opt_reverse);

  /* open output files */

  if (opt_fastqout)
    fp_fastqout = fileopenw(opt_fastqout);
  if (opt_fastaout)
    fp_fastaout = fileopenw(opt_fastaout);
  if (opt_fastqout_notmerged_fwd)
    fp_fastqout_notmerged_fwd = fileopenw(opt_fastqout_notmerged_fwd);
  if (opt_fastqout_notmerged_rev)
    fp_fastqout_notmerged_rev = fileopenw(opt_fastqout_notmerged_rev);
  if (opt_fastaout_notmerged_fwd)
    fp_fastaout_notmerged_fwd = fileopenw(opt_fastaout_notmerged_fwd);
  if (opt_fastaout_notmerged_rev)
    fp_fastaout_notmerged_rev = fileopenw(opt_fastaout_notmerged_rev);
  if (opt_eetabbedout)
    fp_eetabbedout = fileopenw(opt_eetabbedout);

  /* precompute merged quality values */

  precompute_qual();

  /* init progress */

  unsigned long filesize = fastq_get_size(fastq_fwd);
  progress_init("Merging reads", filesize);

  /* start loop */
  
  while(fastq_next(fastq_fwd, 1, chrmap_upcase))
    {
      if (! fastq_next(fastq_rev, 1, chrmap_upcase))
        fatal("More forward reads than reverse reads");

      /* TODO: Check that labels match: label/1 and label/2 */

      total++;

      long fwd_length = fastq_get_sequence_length(fastq_fwd);
      long rev_length = fastq_get_sequence_length(fastq_rev);
      
      char * fwd_header = fastq_get_header(fastq_fwd);
      char * rev_header = fastq_get_header(fastq_rev);

      char * fwd_sequence = fastq_get_sequence(fastq_fwd);
      char * rev_sequence = fastq_get_sequence(fastq_rev);

      char * fwd_quality = fastq_get_quality(fastq_fwd);
      char * rev_quality = fastq_get_quality(fastq_rev);

      long fwd_trunc = fwd_length;
      long rev_trunc = rev_length;
      
      bool skip = 0;

      /* check length */

      if ((fwd_length < opt_fastq_minlen) ||
          (rev_length < opt_fastq_minlen))
        skip = 1;

      /* truncate sequences by quality */

      if (!skip)
        {
          for (long i = 0; i < fwd_length; i++)
            if (get_qual(fwd_quality[i]) <= opt_fastq_truncqual)
              {
                fwd_trunc = i;
                break;
              }
          if (fwd_trunc < opt_fastq_minlen)
            skip = 1;
        }

      if (!skip)
        {          
          for (long i = 0; i < rev_length; i++)
            if (get_qual(rev_quality[i]) <= opt_fastq_truncqual)
              {
                rev_trunc = i;
                break;
              }
          if (rev_trunc < opt_fastq_minlen)
            skip = 1;
        }
      
      /* count n's */

      /* replace quality of N's by zero */

      if (!skip)
        {
          long fwd_ncount = 0;
          for (long i = 0; i < fwd_trunc; i++)
            if (fwd_sequence[i] == 'N')
              {
                fwd_quality[i] = opt_fastq_ascii;
                fwd_ncount++;
              }
          if (fwd_ncount > opt_fastq_maxns)
            skip = 1;
        }

      if (!skip)
        {
          long rev_ncount = 0;
          for (long i = 0; i < rev_trunc; i++)
            if (rev_sequence[i] == 'N')
              {
                rev_quality[i] = opt_fastq_ascii;
                rev_ncount++;
              }
          if (rev_ncount > opt_fastq_maxns)
            skip = 1;
        }
      
      long offset = 0;
      
      if (!skip)
        {
          offset = merge(fwd_sequence, rev_sequence,
                         fwd_quality,  rev_quality,
                         fwd_trunc,    rev_trunc);
        }

      if (offset)
        {
          keep(fwd_header,   rev_header,
               fwd_sequence, rev_sequence,
               fwd_quality,  rev_quality,
               fwd_trunc,    rev_trunc,
               offset);
        }
      else
        discard(fwd_header,   rev_header,
                fwd_sequence, rev_sequence,
                fwd_quality,  rev_quality,
                fwd_length,   rev_length);

      progress_update(fastq_get_position(fastq_fwd));
    }
  
  progress_done();
  
  if (fastq_next(fastq_rev, 1, chrmap_upcase))
    fatal("More reverse reads than forward reads");

  fprintf(stderr,
          "%lu read pairs total\n",
          total);

  fprintf(stderr,
          "%lu pairs merged (%.1lf%%) and %lu pairs not merged (%.1lf%%).\n",
          merged,
          100.0 * merged / total,
          notmerged,
          100.0 * notmerged / total);

  if (opt_eetabbedout)
    fclose(fp_eetabbedout);
  if (opt_fastaout_notmerged_rev)
    fclose(fp_fastaout_notmerged_rev);
  if (opt_fastaout_notmerged_fwd)
    fclose(fp_fastaout_notmerged_fwd);
  if (opt_fastqout_notmerged_rev)
    fclose(fp_fastqout_notmerged_rev);
  if (opt_fastqout_notmerged_fwd)
    fclose(fp_fastqout_notmerged_fwd);
  if (opt_fastaout)
    fclose(fp_fastaout);
  if (opt_fastqout)
    fclose(fp_fastqout);

  fastq_close(fastq_rev);
  fastq_rev = 0;
  fastq_close(fastq_fwd);
  fastq_fwd = 0;

  if (merged_sequence)
    free(merged_sequence);
  if (merged_quality)
    free(merged_quality);
  if (merged_header)
    free(merged_header);
}
