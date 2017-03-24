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

#define INPUTCHUNKSIZE 10000

// scores in bits

const double merge_minscore =    16.0;
const double gap_penalty    =     5.7;

/* fragment length distribution */

const double mean = 450.0;
const double sigma = 30.0;

/* static variables */

static FILE * fp_fastqout = 0;
static FILE * fp_fastaout = 0;
static FILE * fp_fastqout_notmerged_fwd = 0;
static FILE * fp_fastqout_notmerged_rev = 0;
static FILE * fp_fastaout_notmerged_fwd = 0;
static FILE * fp_fastaout_notmerged_rev = 0;
static FILE * fp_eetabbedout = 0;
static fastx_handle fastq_fwd;
static fastx_handle fastq_rev;
static int64_t merged = 0;
static int64_t notmerged = 0;
static int64_t total = 0;
static double sum_squared_fragment_length = 0.0;
static double sum_fragment_length = 0.0;
static pthread_t * pthread;
static pthread_attr_t attr;
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_progress;
static pthread_mutex_t mutex_counters;
static pthread_mutex_t mutex_output;
static pthread_cond_t cond_output;
static int64_t output_next;
static char merge_qual_same[128][128];
static char merge_qual_diff[128][128];
static double match_score[128][128];
static double mism_score[128][128];
static double q2p[128];

static double sum_ee_fwd = 0.0;
static double sum_ee_rev = 0.0;
static double sum_ee_merged = 0.0;
static uint64_t sum_errors_fwd = 0.0;
static uint64_t sum_errors_rev = 0.0;

/* reasons for not merging:
   - input seq too short (after truncation)
   - input seq too long
   - too many Ns in input
   - no significant overlap found
   - too many differences (maxdiff)
   - staggered
   - indels in overlap region
   - potential repeats in overlap region / multiple overlaps
   - merged sequence too long
   - merged sequence too short
*/


typedef struct merge_data_s
{
  char * fwd_header;
  char * rev_header;
  char * fwd_sequence;
  char * rev_sequence;
  char * fwd_quality;
  char * rev_quality;
  int64_t header_alloc;
  int64_t seq_alloc;
  int64_t fwd_length;
  int64_t rev_length;
  int64_t fwd_trunc;
  int64_t rev_trunc;
  int64_t pair_no;
  char * merged_header;
  char * merged_sequence;
  char * merged_quality;
  int64_t merged_length;
  int64_t merged_header_alloc;
  int64_t merged_seq_alloc;
  double ee_merged;
  double ee_fwd;
  double ee_rev;
  int64_t fwd_errors;
  int64_t rev_errors;
  int64_t offset;
  bool merged;
  char * cigar;
  int64_t cigar_length;
  int64_t fwd_merge_start;
  int64_t fwd_merge_end;
  int64_t rev_merge_start;
  int64_t rev_merge_end;
} merge_data_t;

FILE * fileopenw(char * filename)
{
  FILE * fp = 0;
  fp = fopen(filename, "w");
  if (!fp)
    fatal("Unable to open file for writing (%s)", filename);
  return fp;
}

inline int get_qual(char q)
{
  int qual = q - opt_fastq_ascii;
  char msg[200];

  if (qual < opt_fastq_qmin)
    {
      snprintf(msg, 200, "FASTQ quality value (%d) below qmin (%" PRId64 ")",
               qual, opt_fastq_qmin);
      fatal(msg);
    }
  else if (qual > opt_fastq_qmax)
    {
      snprintf(msg, 200, "FASTQ quality value (%d) above qmax (%" PRId64 ")",
               qual, opt_fastq_qmax);
      fatal(msg);
    }
  return qual;
}

inline double q_to_p(int q)
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

  for (int x = 33; x <= 126; x++)
    {
      double px = q_to_p(x);
      q2p[x] = px;

      for (int y = 33; y <= 126; y++)
        {
          double py = q_to_p(y);

          double p, q;

          /* Quality score equations from Edgar & Flyvbjerg (2015) */

          /* Match */
          p = px * py / 3.0 / (1.0 - px - py + 4.0 * px * py / 3.0);
          q = opt_fastq_ascii + MIN(round(-10.0*log10(p)), opt_fastq_qmaxout);
          merge_qual_same[x][y] = q;

          /* Mismatch, x is highest quality */
          p = px * (1.0 - py / 3.0) / (px + py - 4.0 * px * py / 3.0);
          q = opt_fastq_ascii + MIN(round(-10.0*log10(p)), opt_fastq_qmaxout);
          merge_qual_diff[x][y] = q;

          /* Match */

          /* observed match,
             probability that they truly are identical,
             if error probabilites of px and py, resp. */

          p = 1.0 - px - py + px * py * 4.0 / 3.0;
          //          match_score[x][y] = log2(p/0.25);
          match_score[x][y] = log2(p/0.25);

          /* Mismatch */

          /* observed mismatch,
             probability that they truly are identical,
             if error probabilites of px and py, resp. */
          
          /* Error: p = (px + py) / 3.0 - px * py * 4.0 / 9.0; */
          mism_score[x][y] = log2((1-p)/0.75);
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

void keep(merge_data_t * ip)
{
  merged++;
  sum_fragment_length += ip->merged_length;
  sum_squared_fragment_length += ip->merged_length * ip->merged_length;

  sum_ee_merged += ip->ee_merged;
  sum_ee_fwd += ip->ee_fwd;
  sum_ee_rev += ip->ee_rev;
  sum_errors_fwd += ip->fwd_errors;
  sum_errors_rev += ip->rev_errors;

  if (opt_fastqout)
    {
      if (opt_fastq_eeout || opt_eeout)
        fastq_print_with_ee(fp_fastqout,
                            ip->merged_header,
                            ip->merged_sequence,
                            ip->merged_quality,
                            ip->ee_merged);
      else
        fastq_print(fp_fastqout,
                    ip->merged_header,
                    ip->merged_sequence,
                    ip->merged_quality);
    }

  if (opt_fastaout)
    {
      if (opt_fastq_eeout || opt_eeout)
        fasta_print_ee(fp_fastaout,
                       ip->merged_header,
                       ip->merged_sequence,
                       strlen(ip->merged_sequence),
                       ip->ee_merged);
      else
        fasta_print(fp_fastaout,
                    ip->merged_header,
                    ip->merged_sequence,
                    strlen(ip->merged_sequence));
    }

  if (opt_eetabbedout)
    fprintf(fp_eetabbedout, "%.2lf\t%.2lf\t%" PRId64 "\t%" PRId64 "\n",
            ip->ee_fwd, ip->ee_rev, ip->fwd_errors, ip->rev_errors);
}

void discard(merge_data_t * ip)
{
  notmerged++;

  if (opt_fastqout_notmerged_fwd)
    fastq_print(fp_fastqout_notmerged_fwd,
                ip->fwd_header,
                ip->fwd_sequence,
                ip->fwd_quality);

  if (opt_fastqout_notmerged_rev)
    fastq_print(fp_fastqout_notmerged_rev,
                ip->rev_header,
                ip->rev_sequence,
                ip->rev_quality);

  if (opt_fastaout_notmerged_fwd)
    fasta_print(fp_fastaout_notmerged_fwd,
                ip->fwd_header,
                ip->fwd_sequence,
                ip->fwd_length);

  if (opt_fastaout_notmerged_rev)
    fasta_print(fp_fastaout_notmerged_rev,
                ip->rev_header,
                ip->rev_sequence,
                ip->rev_length);
}

void merge(merge_data_t * ip)
{
  /* find offset at left side of merged region */

  int64_t insertions = 0;
  int64_t deletions = 0;
  char * p = ip->cigar;
  while (char c = *p++)
    {
      if (c == 'I')
        insertions++;
      if (c == 'D')
        deletions++;
    }

  /* The offset is the distance (overlap) between
     the (truncated) 3' ends of the two sequences */

  /* adjust offset for position in forward sequence */

  int64_t adjusted_offset = ip->offset + deletions - insertions;
  
  /* length of 5' overhang of the forward sequence not merged
     with the reverse sequence */

  int64_t fwd_5prime_overhang = ip->fwd_trunc > adjusted_offset ?
    ip->fwd_trunc - adjusted_offset : 0;

  ip->ee_merged = 0.0;
  ip->ee_fwd = 0.0;
  ip->ee_rev = 0.0;
  ip->fwd_errors = 0;
  ip->rev_errors = 0;

  char sym, qual;
  char fwd_sym, fwd_qual, rev_sym, rev_qual;
  int64_t fwd_pos, rev_pos, merged_pos;
  double ee;

  merged_pos = 0;

  // 5' overhang in forward sequence

  fwd_pos = 0;
  
  while(fwd_pos < fwd_5prime_overhang)
    {
      sym = ip->fwd_sequence[fwd_pos];
      qual = ip->fwd_quality[fwd_pos];
      
      ip->merged_sequence[merged_pos] = sym;
      ip->merged_quality[merged_pos] = qual;

      ee = q2p[(unsigned)qual];
      ip->ee_merged += ee;
      ip->ee_fwd += ee;

      fwd_pos++;
      merged_pos++;
    }
  
  // merged region
  
  // find start

  int64_t rev_3prime_overhang = adjusted_offset > ip->fwd_trunc ?
    adjusted_offset - ip->fwd_trunc : 0;

  rev_pos = ip->rev_trunc - 1 - rev_3prime_overhang;

  p = ip->cigar + ip->cigar_length;

  //  while ((fwd_pos < ip->fwd_trunc) && (rev_pos >= 0))

  while (p > ip->cigar)
    {
      char c = *--p;

      switch (c)
        {
        case 'M':
          
          fwd_sym = ip->fwd_sequence[fwd_pos];
          rev_sym = chrmap_complement[(int)(ip->rev_sequence[rev_pos])];
          fwd_qual = ip->fwd_quality[fwd_pos];
          rev_qual = ip->rev_quality[rev_pos];
          
          merge_sym(& sym,
                    & qual,
                    fwd_qual < 2 ? 'N' : fwd_sym,
                    rev_qual < 2 ? 'N' : rev_sym,
                    fwd_qual,
                    rev_qual);
          
          if (sym != fwd_sym)
            ip->fwd_errors++;
          if (sym != rev_sym)
            ip->rev_errors++;
          
          ip->merged_sequence[merged_pos] = sym;
          ip->merged_quality[merged_pos] = qual;
          ip->ee_merged += q2p[(unsigned)qual];
          ip->ee_fwd += q2p[(unsigned)fwd_qual];
          ip->ee_rev += q2p[(unsigned)rev_qual];
          
          fwd_pos++;
          rev_pos--;
          merged_pos++;
          break;

        case 'D':

          /* base in forward read only */

          qual = ip->fwd_quality[fwd_pos];
          ee = q2p[(unsigned)qual];

          // printf("Merge with D, qual %d ee %f\n", qual-33, ee);

          // best: ee < 0.003 & always -> 1333
          // always, never -> 1341 
          // never, always -> 1354
          // always, always -> 1338
          // never, never -> 1348
          
          if (1)
            {
              sym = ip->fwd_sequence[fwd_pos];
              
              ip->merged_sequence[merged_pos] = sym;
              ip->merged_quality[merged_pos] = opt_fastq_ascii + 2; // qual;
              merged_pos++;
              
              ip->ee_merged += ee;
              ip->ee_fwd += ee;
            }

          fwd_pos++;

          break;
          
        case 'I':
          
          /* base in reverse read only */

          qual = ip->rev_quality[rev_pos];
          ee = q2p[(unsigned)qual];

          //          printf("Merge with I, qual %d ee %f\n", qual-33, ee);
          
          if (0)
            {
              sym = chrmap_complement[(int)(ip->rev_sequence[rev_pos])];
              
              ip->merged_sequence[merged_pos] = sym;
              ip->merged_quality[merged_pos] = qual;
              merged_pos++;
              
              ip->ee_merged += ee;
              ip->ee_rev += ee;
              
            }

          rev_pos--;

          break;
        }
    }

  // 5' overhang in reverse sequence

  while (rev_pos >= 0)
    {
      sym = chrmap_complement[(int)(ip->rev_sequence[rev_pos])];
      qual = ip->rev_quality[rev_pos];
      
      ip->merged_sequence[merged_pos] = sym;
      ip->merged_quality[merged_pos] = qual;
      merged_pos++;

      ee = q2p[(unsigned)qual];
      ip->ee_merged += ee;
      ip->ee_rev += ee;

      rev_pos--;
    }

  int64_t mergelen = merged_pos;
  ip->merged_length = mergelen;
  
  ip->merged_sequence[mergelen] = 0;
  ip->merged_quality[mergelen] = 0;

  if (ip->ee_merged <= opt_fastq_maxee)
    {
      if (opt_label_suffix)
        (void) sprintf(ip->merged_header, "%s%s",
                       ip->fwd_header, opt_label_suffix);
      else
        strcpy(ip->merged_header, ip->fwd_header);

      ip->merged = 1;
    }
}

double normal_dist_prob_density(double x, double mean, double stdev)
{
  return exp(-0.5 * (x - mean) * (x - mean) / (stdev * stdev)) /
    sqrt(2.0 * M_PI * stdev * stdev);
}


static int cigarcount = 0;

int64_t optimize(merge_data_t * ip)
{
  // i = position in forward sequence
  // j = position in reverse sequence
  // m = length of forward sequence (after truncation)
  // n = length of reverse sequence (after truncation)

  int64_t m = ip->fwd_trunc;
  int64_t n = ip->rev_trunc;

  double h[m];
  char trace[m][n];
  int gapped[m];

  for (int i = 0; i < m; i++)
    {
      h[i] = 0.0;
      gapped[i] = 0;
    }

  double best_score = 0.0;
  int64_t best_i = 0;
  
  int64_t j_min = n > m ? n - m : 0;

  int hits = 0;

  for (int64_t j = n-1; j >= j_min; j--)
    {
      char rev_sym = chrmap_complement[(int)(ip->rev_sequence[j])];
      unsigned int rev_qual = ip->rev_quality[j];

      int64_t i_min = m > n ? m - j - 1 : n - j - 1;
      
      double v, f;

      if (i_min == 0)
        {
          v = - (n-1-j) * gap_penalty;
        }
      else
        {
          v = h[i_min - 1];
        }

      f = - (n-j+1) * gap_penalty;

      f = - 1e20; /* illegal gap unless staggered */

      for (int64_t i = i_min; i < m; i++)
        {
          double temp = h[i];

          double e = temp - gap_penalty;

          char fwd_sym = ip->fwd_sequence[i];
          unsigned int fwd_qual = ip->fwd_quality[i];
          
          double match =
            fwd_sym == rev_sym ?
            match_score[fwd_qual][rev_qual] :
            mism_score[fwd_qual][rev_qual];
          
          double g = v + match;
          
          char dir;
          if (g >= e)
            if (g >= f)
              {
                v = g; dir = 'M';
              }
            else
              {
                v = f; dir = 'D';
              }
          else
            if (e >= f)
              {
                v = e; dir = 'I';
              }
            else
              {
                v = f; dir = 'D';
              }
          
          h[i] = v;
          trace[i][j] = dir;

          if (dir != 'M')
            gapped[i-(n-1-j)] = 1;

          f = v - gap_penalty;

          v = temp;
        }
      
      double final = h[m-1];

#if 0
      bool has_gap = gapped[m-1];

      if (has_gap)
        final = -1000.0;
#endif
      
      double len = m + j; // m + n - (n-j) = m+j

      double overlap = m + n - len;

      /* correct for deviation from normal distribution */

      /* 2x250bp MiSeq len 300 stdev 30, no gaps: 0 or 1 err on 100X cov. */
      /* 2x250bp MiSeq len 450 stdev 30, no gaps: 21845 errors on 100X cov. */
      /* lacking merge of 21589 reads, 205 wrong insert, 51 fp (of 928300) */
      /* generally, too much overlap */

      double prob = normal_dist_prob_density(len, mean, sigma);
      double prob_max = normal_dist_prob_density(mean, mean, sigma);
      
      double adjust = log2(prob/prob_max);

#if 0
      final += adjust;
#endif
      
#if 0
      final = overlap / ((2*overlap - final) * (2*overlap - final) + 1); 
#endif

#if 0
      final = final / overlap;
#endif

      if (final >= merge_minscore)
        {
#if 0
          fprintf(stderr, 
                  "len=%.0f, prob=%f, adjust=%f, final=%f\n",
                  len, prob, adjust, final);
#endif
          hits++;
        }


      if (final > best_score)
        {
          best_score = final;
          best_i = n - j;
        }
    }
  
#if 0
  if (hits > 1)
    fprintf(stderr, "Seq: %s Hits: %d\n", ip->fwd_header, hits);
#endif

  if (best_score >= merge_minscore)
    {

      ip->fwd_merge_end = m - 1;
      ip->rev_merge_end = n - best_i;

#if 1
      char * cigar = ip->cigar;
      char * p = cigar;

      int64_t i = m - 1;
      int64_t j = n - best_i;
      
      char c = 0;

      int qual = -1;

      while ((i>=0) && (j<n))
        {
          c = trace[i][j];
          *p++ = c;
          
          switch(c)
            {
            case 'M':
              i--;
              j++;
              break;
            case 'I':
              qual = ip->rev_quality[j];
              j++;
              break;
            case 'D':
              qual = ip->fwd_quality[i];
              i--;
              break;
            }
        }
      *p = 0;
      ip->cigar_length = p - cigar;

      if (c == 'M')
        {
          ip->fwd_merge_start = i+1;
          ip->rev_merge_start = j-1;
        }
      else if (c == 'I')
        {
          ip->fwd_merge_start = i;
          ip->rev_merge_start = j-1;
        }
      else if (c == 'D')
        {
          ip->fwd_merge_start = i+1;
          ip->rev_merge_start = j;
        }

      if (qual >= 0)
        return 0; // Do not align sequences with gaps

#if 0
      if (qual >= 0)
        {
          fprintf(stderr, "R1: %s CIGAR (%d): ", ip->fwd_header, ++cigarcount);

          char last = 0;
          int64_t run = 0;
          while (p > cigar)
            {
              char op = *--p;
              if (op == last)
                run++;
              else
                {
                  if (run > 1)
                    fprintf(stderr, "%c%" PRIi64, last, run);
                  else
                    fprintf(stderr, "%c", last);
                  run = 1;
                  last = op;
                }
            }
          if (run > 1)
            fprintf(stderr, "%c%" PRIi64, last, run);
          else if (run)
            fprintf(stderr, "%c", last);
          fprintf(stderr, " Qual: %d", qual - (int)opt_fastq_ascii);
          fprintf(stderr, "\n");
        }
#endif

#endif      
      return best_i;
    }
  else
    return 0;

#if 0


  int64_t i1 = 1;

  i1 = MAX(i1, ip->fwd_trunc + ip->rev_trunc - opt_fastq_maxmergelen);

  int64_t i2 = opt_fastq_allowmergestagger ?
    (ip->fwd_trunc + ip->rev_trunc - opt_fastq_minovlen) :
    ip->fwd_trunc;

  i2 = MIN(i2, ip->fwd_trunc + ip->rev_trunc - opt_fastq_minmergelen);

  double best_score = 0.0;
  int64_t best_i = 0;
  int64_t best_diffs = 0;

  for(int64_t i = i1; i <= i2; i++)
    {
      int64_t fwd_3prime_overhang = i > ip->rev_trunc ? i - ip->rev_trunc : 0;
      int64_t rev_3prime_overhang = i > ip->fwd_trunc ? i - ip->fwd_trunc : 0;
      int64_t overlap = i - fwd_3prime_overhang - rev_3prime_overhang;
      int64_t fwd_pos_start = ip->fwd_trunc - fwd_3prime_overhang - 1;
      int64_t rev_pos_start = ip->rev_trunc - rev_3prime_overhang - overlap;

      int64_t fwd_pos = fwd_pos_start;
      int64_t rev_pos = rev_pos_start;
      double score = 0.0;
      
      int64_t diffs = 0;
      
      for (int64_t j=0; j < overlap; j++)
        {
          char fwd_sym = ip->fwd_sequence[fwd_pos];
          char rev_sym = chrmap_complement[(int)(ip->rev_sequence[rev_pos])];
          
          unsigned int fwd_qual = ip->fwd_quality[fwd_pos];
          unsigned int rev_qual = ip->rev_quality[rev_pos];
          fwd_pos--;
          rev_pos++;
          
          if (fwd_sym == rev_sym)
            score += match_score[fwd_qual][rev_qual];
          else
            {
              score += mism_score[fwd_qual][rev_qual];
              diffs++;
              
              if (score < merge_giveup)
                break;
            }
        }
      
      if (score > best_score)
        {
          best_score = score;
          best_i = i;
          best_diffs = diffs;
        }
    }
  
  if ((best_score >= merge_minscore) && (best_diffs <= opt_fastq_maxdiffs))
    return best_i;
  else
    return 0;

#endif

}

void process(merge_data_t * ip)
{
  ip->merged = 0;

  bool skip = 0;

  /* check length */

  if ((ip->fwd_length < opt_fastq_minlen) ||
      (ip->rev_length < opt_fastq_minlen) ||
      (ip->fwd_length > opt_fastq_maxlen) ||
      (ip->rev_length > opt_fastq_maxlen))
    skip = 1;

  /* truncate sequences by quality */

  int64_t fwd_trunc = ip->fwd_length;

  if (!skip)
    {
      for (int64_t i = 0; i < ip->fwd_length; i++)
        if (get_qual(ip->fwd_quality[i]) <= opt_fastq_truncqual)
          {
            fwd_trunc = i;
            break;
          }
      if (fwd_trunc < opt_fastq_minlen)
        skip = 1;
    }

  ip->fwd_trunc = fwd_trunc;

  int64_t rev_trunc = ip->rev_length;

  if (!skip)
    {
      for (int64_t i = 0; i < ip->rev_length; i++)
        if (get_qual(ip->rev_quality[i]) <= opt_fastq_truncqual)
          {
            rev_trunc = i;
            break;
          }
      if (rev_trunc < opt_fastq_minlen)
        skip = 1;
    }

  ip->rev_trunc = rev_trunc;

  /* count n's */

  /* replace quality of N's by zero */

  if (!skip)
    {
      int64_t fwd_ncount = 0;
      for (int64_t i = 0; i < fwd_trunc; i++)
        if (ip->fwd_sequence[i] == 'N')
          {
            ip->fwd_quality[i] = opt_fastq_ascii;
            fwd_ncount++;
          }
      if (fwd_ncount > opt_fastq_maxns)
        skip = 1;
    }

  if (!skip)
    {
      int64_t rev_ncount = 0;
      for (int64_t i = 0; i < rev_trunc; i++)
        if (ip->rev_sequence[i] == 'N')
          {
            ip->rev_quality[i] = opt_fastq_ascii;
            rev_ncount++;
          }
      if (rev_ncount > opt_fastq_maxns)
        skip = 1;
    }

  ip->offset = 0;

  if (!skip)
    ip->offset = optimize(ip);

  if (ip->offset >= opt_fastq_minovlen)
    merge(ip);
}

bool read_pair(merge_data_t * ip)
{
  int64_t suffix_len = opt_label_suffix ? strlen(opt_label_suffix) : 0;

  if (fastq_next(fastq_fwd, 0, chrmap_upcase))
    {
      if (! fastq_next(fastq_rev, 0, chrmap_upcase))
        fatal("More forward reads than reverse reads");

      /* allocate more memory if necessary */

      int64_t fwd_header_len = fastq_get_header_length(fastq_fwd);
      int64_t rev_header_len = fastq_get_header_length(fastq_rev);
      int64_t header_needed = MAX(fwd_header_len, rev_header_len) + 1;

      if (header_needed > ip->header_alloc)
        {
          ip->header_alloc = header_needed;
          ip->fwd_header = (char*) xrealloc(ip->fwd_header, header_needed);
          ip->rev_header = (char*) xrealloc(ip->rev_header, header_needed);
        }

      ip->fwd_length = fastq_get_sequence_length(fastq_fwd);
      ip->rev_length = fastq_get_sequence_length(fastq_rev);
      int64_t seq_needed = MAX(ip->fwd_length, ip->rev_length) + 1;

      if (seq_needed > ip->seq_alloc)
        {
          ip->seq_alloc = seq_needed;
          ip->fwd_sequence = (char*) xrealloc(ip->fwd_sequence, seq_needed);
          ip->rev_sequence = (char*) xrealloc(ip->rev_sequence, seq_needed);
          ip->fwd_quality  = (char*) xrealloc(ip->fwd_quality,  seq_needed);
          ip->rev_quality  = (char*) xrealloc(ip->rev_quality,  seq_needed);
        }

      int64_t merged_seq_needed = ip->fwd_length + ip->rev_length + 1;

      if (merged_seq_needed > ip->merged_seq_alloc)
        {
          ip->merged_seq_alloc = merged_seq_needed;
          ip->merged_sequence = (char*) xrealloc(ip->merged_sequence,
                                                 merged_seq_needed);
          ip->merged_quality = (char*) xrealloc(ip->merged_quality,
                                                merged_seq_needed);
          ip->cigar = (char*) xrealloc(ip->cigar, merged_seq_needed);
        }

      int64_t merged_header_needed = fwd_header_len + suffix_len + 1;

      if (merged_header_needed > ip->merged_header_alloc)
        {
          ip->merged_header_alloc = merged_header_needed;
          ip->merged_header = (char*) xrealloc(ip->merged_header,
                                               merged_header_needed);
        }

      /* make local copies of the seq, header and qual */

      strcpy(ip->fwd_header,   fastq_get_header(fastq_fwd));
      strcpy(ip->rev_header,   fastq_get_header(fastq_rev));
      strcpy(ip->fwd_sequence, fastq_get_sequence(fastq_fwd));
      strcpy(ip->rev_sequence, fastq_get_sequence(fastq_rev));
      strcpy(ip->fwd_quality,  fastq_get_quality(fastq_fwd));
      strcpy(ip->rev_quality,  fastq_get_quality(fastq_rev));

      ip->merged_header[0] = 0;
      ip->merged_sequence[0] = 0;
      ip->merged_quality[0] = 0;
      ip->merged = 0;
      ip->pair_no = total++;
      ip->cigar[0] = 0;

      return 1;
    }
  else
    return 0;
}

void pair()
{
  merge_data_t * merge_buffer = (merge_data_t*)
    xmalloc(INPUTCHUNKSIZE * sizeof(merge_data_s));

  for(int64_t i=0; i<INPUTCHUNKSIZE; i++)
    {
      merge_data_s * ip = merge_buffer + i;
      ip->fwd_header = 0;
      ip->rev_header = 0;
      ip->fwd_sequence = 0;
      ip->rev_sequence = 0;
      ip->fwd_quality = 0;
      ip->rev_quality = 0;
      ip->header_alloc = 0;
      ip->seq_alloc = 0;
      ip->fwd_length = 0;
      ip->rev_length = 0;
      ip->pair_no = 0;
      ip->cigar = 0;
    }

  bool more = 1;
  while (more)
    {
      pthread_mutex_lock(&mutex_input);
      progress_update(fastq_get_position(fastq_fwd));
      int64_t pairs = 0;
      while (more && (pairs < INPUTCHUNKSIZE))
        {
          more = read_pair(merge_buffer + pairs);
          if (more)
            pairs++;
        }
      pthread_mutex_unlock(&mutex_input);

      for(int64_t i=0; i<pairs; i++)
        process(merge_buffer + i);

      pthread_mutex_lock(&mutex_output);
      while(output_next < merge_buffer[0].pair_no)
        pthread_cond_wait(&cond_output, &mutex_output);

      for(int64_t i=0; i<pairs; i++)
        {
          merge_data_s * ip = merge_buffer + i;
          if (ip->merged)
            keep(ip);
          else
            discard(ip);
          output_next++;
        }

      pthread_cond_broadcast(&cond_output);
      pthread_mutex_unlock(&mutex_output);
    }

  for(int64_t i=0; i<INPUTCHUNKSIZE; i++)
    {
      merge_data_s * ip = merge_buffer + i;
      if (ip->fwd_header)
        xfree(ip->fwd_header);
      if (ip->rev_header)
        xfree(ip->rev_header);
      if (ip->fwd_sequence)
        xfree(ip->fwd_sequence);
      if (ip->rev_sequence)
        xfree(ip->rev_sequence);
      if (ip->fwd_quality)
        xfree(ip->fwd_quality);
      if (ip->rev_quality)
        xfree(ip->rev_quality);
      if (ip->cigar)
        xfree(ip->cigar);
    }

  xfree(merge_buffer);
}

void * pair_worker(void * vp)
{
  int64_t t = (int64_t) vp;
  (void) t;
  pair();
  return 0;
}

void pair_all()
{
  pthread_mutex_init(&mutex_input, NULL);
  pthread_mutex_init(&mutex_counters, NULL);
  pthread_mutex_init(&mutex_progress, NULL);

  pthread_mutex_init(&mutex_output, NULL);
  pthread_cond_init(&cond_output, 0);
  output_next = 0;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  pthread = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));

  for(int t=0; t<opt_threads; t++)
    if (pthread_create(pthread+t, &attr, pair_worker, (void*)(int64_t)t))
      fatal("Cannot create thread");

  for(int t=0; t<opt_threads; t++)
    if (pthread_join(pthread[t], NULL))
      fatal("Cannot join thread");

  xfree(pthread);

  pthread_attr_destroy(&attr);

  pthread_cond_destroy(&cond_output);
  pthread_mutex_destroy(&mutex_output);

  pthread_mutex_destroy(&mutex_input);
  pthread_mutex_destroy(&mutex_counters);
  pthread_mutex_destroy(&mutex_progress);
}

void fastq_mergepairs()
{
  /* open input files */

  fastq_fwd = fastq_open(opt_fastq_mergepairs);
  fastq_rev = fastq_open(opt_reverse);

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

  /* main */

  uint64_t filesize = fastq_get_size(fastq_fwd);
  progress_init("Merging reads", filesize);

  pair_all();

  progress_done();

  if (fastq_next(fastq_rev, 1, chrmap_upcase))
    fatal("More reverse reads than forward reads");

  fprintf(stderr,
          "%10" PRIu64 "  Pairs\n",
          total);

  fprintf(stderr,
          "%10" PRIu64 "  Merged (%.1lf%%)\n",
          merged,
          100.0 * merged / total);

  fprintf(stderr,
          "%10" PRIu64 "  Not merged (%.1lf%%)\n",
          notmerged,
          100.0 * notmerged / total);

  double mean = sum_fragment_length / merged;
  double stdev = sqrt((sum_squared_fragment_length
                       - 2.0 * mean * sum_fragment_length
                       + mean * mean * merged)
                      / (merged + 0.0));
  
  fprintf(stderr,
          "%10.2f  Mean fragment length\n",
          mean);
  
  fprintf(stderr,
          "%10.2f  Standard deviation of fragment length\n",
          stdev);

  fprintf(stderr,
          "%10.2f  Mean expected error in forward sequences\n",
          sum_ee_fwd / merged);

  fprintf(stderr,
          "%10.2f  Mean expected error in reverse sequences\n",
          sum_ee_rev / merged);
  
  fprintf(stderr,
          "%10.2f  Mean expected error in merged sequences\n",
          sum_ee_merged / merged);

  fprintf(stderr,
          "%10.2f  Mean observed errors in merged region of forward sequences\n",
          1.0 * sum_errors_fwd / merged);

  fprintf(stderr,
          "%10.2f  Mean observed errors in merged region of reverse sequences\n",
          1.0 * sum_errors_rev / merged);

  fprintf(stderr,
          "%10.2f  Mean observed errors in merged region\n",
          1.0 * (sum_errors_fwd + sum_errors_rev) / merged);

  
  /* clean up */

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
}
