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

/* chunk constants */

static const int chunk_size = 500; /* read pairs per chunk */
static const int chunk_factor = 2; /* chunks per thread */

/* scores in bits */

static const int k                        = 5;
static const double merge_minscore        = 16.0;
static const double merge_dropmax         = 16.0;
static const int merge_mindiagcount       = 4;
static const int merge_minrepeatdiagcount = 12;
static const double merge_mismatchmax     = -4.0;

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

static double sum_read_length = 0.0;
static double sum_squared_fragment_length = 0.0;
static double sum_fragment_length = 0.0;

static pthread_t * pthread;
static pthread_attr_t attr;

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

static uint64_t failed_undefined = 0;
static uint64_t failed_minlen = 0;
static uint64_t failed_maxlen = 0;
static uint64_t failed_maxns = 0;
static uint64_t failed_minovlen = 0;
static uint64_t failed_maxdiffs = 0;
static uint64_t failed_staggered = 0;
static uint64_t failed_indel = 0;
static uint64_t failed_repeat = 0;
static uint64_t failed_minmergelen = 0;
static uint64_t failed_maxmergelen = 0;
static uint64_t failed_maxee = 0;
static uint64_t failed_minscore = 0;
static uint64_t failed_nokmers = 0;

/* reasons for not merging:
   - undefined
   - ok
   - input seq too short (after truncation)
   - input seq too long
   - too many Ns in input
   - overlap too short
   - too many differences (maxdiff)
   - staggered
   - indels in overlap region
   - potential repeats in overlap region / multiple overlaps
   - merged sequence too short
   - merged sequence too long
   - expected error too high
   - alignment score too low, insignificant, potential indel
   - too few kmers on same diag found
*/

enum reason_enum
  {
    undefined,
    ok,
    minlen,
    maxlen,
    maxns,
    minovlen,
    maxdiffs,
    staggered,
    indel,
    repeat,
    minmergelen,
    maxmergelen,
    maxee,
    minscore,
    nokmers
  };

enum state_enum
  {
    empty,
    filled,
    inprogress,
    processed
  };

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
  reason_enum reason;
  state_enum state;
} merge_data_t;


typedef struct chunk_s
{
  int size; /* size of merge_data = number of pairs of reads */
  state_enum state; /* state of chunk: empty, read, processed */
  merge_data_t * merge_data; /* data for merging */
} chunk_t;

static chunk_t * chunks; /* pointer to array of chunks */

static int chunk_count;
static int chunk_read_next;
static int chunk_process_next;
static int chunk_write_next;
static bool finished_reading = false;
static bool finished_all = false;
static int pairs_read = 0;
static int pairs_written = 0;

static pthread_mutex_t mutex_chunks;
static pthread_cond_t cond_chunks;


FILE * fileopenw(char * filename)
{
  FILE * fp = 0;
  fp = fopen_output(filename);
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

          /*
            observed match,
            p = probability that they truly are identical,
            given error probabilites of px and py, resp.
          */

          // Given two initially identical aligned bases, and
          // the error probabilities px and py,
          // what is the probability of observing a match (or a mismatch)?

          p = 1.0 - px - py + px * py * 4.0 / 3.0;
          match_score[x][y] = log2(p/0.25);

          // Use a minimum mismatch penalty

          mism_score[x][y] = MIN(log2((1.0-p)/0.75), merge_mismatchmax);
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
      fastq_print_general(fp_fastqout,
                          ip->merged_sequence,
                          ip->merged_length,
                          ip->merged_header,
                          strlen(ip->merged_header),
                          ip->merged_quality,
                          0,
                          merged,
                          (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                          ip->ee_merged);
    }

  if (opt_fastaout)
    {
      fasta_print_general(fp_fastaout,
                          0,
                          ip->merged_sequence,
                          ip->merged_length,
                          ip->merged_header,
                          strlen(ip->merged_header),
                          0,
                          merged,
                          -1,
                          -1,
                          (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                          ip->ee_merged);
    }

  if (opt_eetabbedout)
    fprintf(fp_eetabbedout, "%.2lf\t%.2lf\t%" PRId64 "\t%" PRId64 "\n",
            ip->ee_fwd, ip->ee_rev, ip->fwd_errors, ip->rev_errors);
}

void discard(merge_data_t * ip)
{
  switch(ip->reason)
    {
    case undefined:
      failed_undefined++;
      break;

    case ok:
      break;

    case minlen:
      failed_minlen++;
      break;

    case maxlen:
      failed_maxlen++;
      break;

    case maxns:
      failed_maxns++;
      break;

    case minovlen:
      failed_minovlen++;
      break;

    case maxdiffs:
      failed_maxdiffs++;
      break;

    case staggered:
      failed_staggered++;
      break;

    case indel:
      failed_indel++;
      break;

    case repeat:
      failed_repeat++;
      break;

    case minmergelen:
      failed_minmergelen++;
      break;

    case maxmergelen:
      failed_maxmergelen++;
      break;

    case maxee:
      failed_maxee++;
      break;

    case minscore:
      failed_minscore++;
      break;

    case nokmers:
      failed_nokmers++;
      break;
    }

  notmerged++;

  if (opt_fastqout_notmerged_fwd)
    fastq_print_general(fp_fastqout_notmerged_fwd,
                        ip->fwd_sequence,
                        ip->fwd_length,
                        ip->fwd_header,
                        strlen(ip->fwd_header),
                        ip->fwd_quality,
                        0,
                        notmerged,
                        0, 0.0);
                        
  if (opt_fastqout_notmerged_rev)
    fastq_print_general(fp_fastqout_notmerged_rev,
                        ip->rev_sequence,
                        ip->rev_length,
                        ip->rev_header,
                        strlen(ip->rev_header),
                        ip->rev_quality,
                        0,
                        notmerged,
                        0, 0.0);

  if (opt_fastaout_notmerged_fwd)
    fasta_print_general(fp_fastaout_notmerged_fwd,
                        0,
                        ip->fwd_sequence,
                        ip->fwd_length,
                        ip->fwd_header,
                        strlen(ip->fwd_header),
                        0,
                        notmerged,
                        -1, -1,
                        0, 0.0);
                        
  if (opt_fastaout_notmerged_rev)
    fasta_print_general(fp_fastaout_notmerged_rev,
                        0,
                        ip->rev_sequence,
                        ip->rev_length,
                        ip->rev_header,
                        strlen(ip->rev_header),
                        0,
                        notmerged,
                        -1, -1,
                        0, 0.0);
}

void merge(merge_data_t * ip)
{
  /* length of 5' overhang of the forward sequence not merged
     with the reverse sequence */

  int64_t fwd_5prime_overhang = ip->fwd_trunc > ip->offset ?
    ip->fwd_trunc - ip->offset : 0;

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

  // Merged region

  int64_t rev_3prime_overhang = ip->offset > ip->fwd_trunc ?
    ip->offset - ip->fwd_trunc : 0;

  rev_pos = ip->rev_trunc - 1 - rev_3prime_overhang;

  while ((fwd_pos < ip->fwd_trunc) && (rev_pos >= 0))
    {
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

      ip->reason = ok;
      ip->merged = 1;
    }
  else
    {
      ip->reason = maxee;
    }
}

int64_t optimize(merge_data_t * ip,
                 kh_handle_s * kmerhash)
{
  /* ungapped alignment in each diagonal */

  int64_t i1 = 1;
  int64_t i2 = ip->fwd_trunc + ip->rev_trunc - 1;

  double best_score = 0.0;
  int64_t best_i = 0;
  int64_t best_diffs = 0;

  int hits = 0;

  int kmers = 0;

  int diags[ip->fwd_trunc + ip->rev_trunc];

  kh_insert_kmers(kmerhash, k, ip->fwd_sequence, ip->fwd_trunc);
  kh_find_diagonals(kmerhash, k, ip->rev_sequence, ip->rev_trunc, diags);

  for(int64_t i = i1; i <= i2; i++)
    {
      int diag = ip->rev_trunc + ip->fwd_trunc - i;
      int diagcount = diags[diag];

      if (diagcount >= merge_mindiagcount)
        {
          kmers = 1;

          if (diagcount >= merge_minrepeatdiagcount)
            {
              hits++;
              if (hits > 1)
                break;
            }

          /* for each interesting diagonal */

          int64_t fwd_3prime_overhang
            = i > ip->rev_trunc ? i - ip->rev_trunc : 0;
          int64_t rev_3prime_overhang
            = i > ip->fwd_trunc ? i - ip->fwd_trunc : 0;
          int64_t overlap
            = i - fwd_3prime_overhang - rev_3prime_overhang;
          int64_t fwd_pos_start
            = ip->fwd_trunc - fwd_3prime_overhang - 1;
          int64_t rev_pos_start
            = ip->rev_trunc - rev_3prime_overhang - overlap;

          int64_t fwd_pos = fwd_pos_start;
          int64_t rev_pos = rev_pos_start;
          double score = 0.0;

          int64_t diffs = 0;
          double score_high = 0.0;
          double dropmax = 0.0;

          for (int64_t j=0; j < overlap; j++)
            {
              /* for each pair of bases in the overlap */

              char fwd_sym
                = ip->fwd_sequence[fwd_pos];
              char rev_sym
                = chrmap_complement[(int)(ip->rev_sequence[rev_pos])];

              unsigned int fwd_qual = ip->fwd_quality[fwd_pos];
              unsigned int rev_qual = ip->rev_quality[rev_pos];

              fwd_pos--;
              rev_pos++;

              if (fwd_sym == rev_sym)
                {
                  score += match_score[fwd_qual][rev_qual];
                  if (score > score_high)
                    score_high = score;
                }
              else
                {
                  score += mism_score[fwd_qual][rev_qual];
                  diffs++;
                  if (score < score_high - dropmax)
                    dropmax = score_high - score;
                }
            }

          if (dropmax >= merge_dropmax)
            score = 0.0;

          if (score > best_score)
            {
              best_score = score;
              best_i = i;
              best_diffs = diffs;
            }
        }

    }

  if (hits > 1)
    {
      ip->reason = repeat;
      return 0;
    }

  if ((! opt_fastq_allowmergestagger) && (best_i > ip->fwd_trunc))
    {
      ip->reason = staggered;
      return 0;
    }

  if (best_diffs > opt_fastq_maxdiffs)
    {
      ip->reason = maxdiffs;
      return 0;
    }

  if (kmers == 0)
    {
      ip->reason = nokmers;
      return 0;
    }

  if (best_score < merge_minscore)
    {
      ip->reason = minscore;
      return 0;
    }

  if (best_i < opt_fastq_minovlen)
    {
      ip->reason = minovlen;
      return 0;
    }

  int mergelen = ip->fwd_trunc + ip->rev_trunc - best_i;

  if (mergelen < opt_fastq_minmergelen)
    {
      ip->reason = minmergelen;
      return 0;
    }

  if (mergelen > opt_fastq_maxmergelen)
    {
      ip->reason = maxmergelen;
      return 0;
    }

  return best_i;
}

void process(merge_data_t * ip,
             struct kh_handle_s * kmerhash)
{
  ip->merged = 0;

  bool skip = 0;

  /* check length */

  if ((ip->fwd_length < opt_fastq_minlen) ||
      (ip->rev_length < opt_fastq_minlen))
    {
      ip->reason = minlen;
      skip = 1;
    }

  if ((ip->fwd_length > opt_fastq_maxlen) ||
      (ip->rev_length > opt_fastq_maxlen))
    {
      ip->reason = maxlen;
      skip = 1;
    }

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
        {
          ip->reason = minlen;
          skip = 1;
        }
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
        {
          ip->reason = minlen;
          skip = 1;
        }
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
        {
          ip->reason = maxns;
          skip = 1;
        }
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
        {
          ip->reason = maxns;
          skip = 1;
        }
    }

  ip->offset = 0;

  if (!skip)
    ip->offset = optimize(ip, kmerhash);

  if (ip->offset > 0)
    merge(ip);

  ip->state = processed;
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

      sum_read_length += ip->fwd_length + ip->rev_length;

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

      return 1;
    }
  else
    return 0;
}

void keep_or_discard(merge_data_t * ip)
{
  if (ip->merged)
    keep(ip);
  else
    discard(ip);
}

void init_merge_data(merge_data_t * ip)
{
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
  ip->fwd_trunc = 0;
  ip->rev_trunc = 0;
  ip->pair_no = 0;
  ip->reason = undefined;
  ip->merged_seq_alloc = 0;
  ip->merged_sequence = 0;
  ip->merged_header = 0;
  ip->merged_header_alloc = 0;
  ip->merged_quality = 0;
  ip->merged_length = 0;
}

void free_merge_data(merge_data_t * ip)
{
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

  if (ip->merged_header)
    xfree(ip->merged_header);
  if (ip->merged_sequence)
    xfree(ip->merged_sequence);
  if (ip->merged_quality)
    xfree(ip->merged_quality);
}

inline void chunk_perform_read()
{
  while((!finished_reading) && (chunks[chunk_read_next].state == empty))
    {
      pthread_mutex_unlock(&mutex_chunks);
      progress_update(fastq_get_position(fastq_fwd));
      int r = 0;
      while ((r < chunk_size) &&
             read_pair(chunks[chunk_read_next].merge_data + r))
        r++;
      chunks[chunk_read_next].size = r;
      pthread_mutex_lock(&mutex_chunks);
      pairs_read += r;
      if (r > 0)
        {
          chunks[chunk_read_next].state = filled;
          chunk_read_next = (chunk_read_next + 1) % chunk_count;
        }
      if (r < chunk_size)
        finished_reading = true;
      pthread_cond_broadcast(&cond_chunks);
    }
}

inline void chunk_perform_write()
{
  while (chunks[chunk_write_next].state == processed)
    {
      pthread_mutex_unlock(&mutex_chunks);
      for(int i = 0; i < chunks[chunk_write_next].size; i++)
        keep_or_discard(chunks[chunk_write_next].merge_data + i);
      pthread_mutex_lock(&mutex_chunks);
      pairs_written += chunks[chunk_write_next].size;
      chunks[chunk_write_next].state = empty;
      if (finished_reading && (pairs_written >= pairs_read))
        finished_all = true;
      chunk_write_next = (chunk_write_next + 1) % chunk_count;
      pthread_cond_broadcast(&cond_chunks);
    }
}

inline void chunk_perform_process(struct kh_handle_s * kmerhash)
{
  int chunk_current = chunk_process_next;
  if (chunks[chunk_current].state == filled)
    {
      chunks[chunk_current].state = inprogress;
      chunk_process_next = (chunk_current + 1) % chunk_count;
      pthread_cond_broadcast(&cond_chunks);
      pthread_mutex_unlock(&mutex_chunks);
      for(int i=0; i<chunks[chunk_current].size; i++)
        process(chunks[chunk_current].merge_data + i, kmerhash);
      pthread_mutex_lock(&mutex_chunks);
      chunks[chunk_current].state = processed;
      pthread_cond_broadcast(&cond_chunks);
    }
}

void * pair_worker(void * vp)
{
  /* new */

  int64_t t = (int64_t) vp;

  struct kh_handle_s * kmerhash = kh_init();

  pthread_mutex_lock(&mutex_chunks);

  while (! finished_all)
    {
      if (opt_threads == 1)
        {
          /* One thread does it all */
          chunk_perform_read();
          chunk_perform_process(kmerhash);
          chunk_perform_write();
        }
      else if (opt_threads == 2)
        {
          if (t == 0)
            {
              /* first thread reads and processes */
              while (!
                     (
                      finished_all
                      ||
                      (chunks[chunk_process_next].state == filled)
                      ||
                      ((!finished_reading) &&
                       chunks[chunk_read_next].state == empty)))
                pthread_cond_wait(&cond_chunks, &mutex_chunks);

              chunk_perform_read();
              chunk_perform_process(kmerhash);
            }
          else /* t == 1 */
            {
              /* second thread writes and processes */
              while (!
                     (
                      finished_all
                      ||
                      (chunks[chunk_process_next].state == filled)
                      ||
                      (chunks[chunk_write_next].state == processed)
                      )
                     )
                pthread_cond_wait(&cond_chunks, &mutex_chunks);

              chunk_perform_write();
              chunk_perform_process(kmerhash);
            }
        }
      else
        {
          if (t == 0)
            {
              /* first thread reads and processes */
              while (!
                     (
                      finished_all
                      ||
                      ((!finished_reading) &&
                       (chunks[chunk_read_next].state == empty))
                      ||
                      (chunks[chunk_process_next].state == filled)
                      )
                     )
                pthread_cond_wait(&cond_chunks, &mutex_chunks);

              chunk_perform_read();
              chunk_perform_process(kmerhash);
            }
          else if (t == opt_threads - 1)
            {
              /* last thread writes and processes */
              while (!
                     (
                      finished_all
                      ||
                      (chunks[chunk_write_next].state == processed)
                      ||
                      (chunks[chunk_process_next].state == filled)
                      )
                     )
                pthread_cond_wait(&cond_chunks, &mutex_chunks);

              chunk_perform_write();
              chunk_perform_process(kmerhash);
            }
          else
            {
              /* the other threads are only processing */
              while (!
                     (
                      finished_all
                      ||
                      (chunks[chunk_process_next].state == filled)
                      )
                     )
                pthread_cond_wait(&cond_chunks, &mutex_chunks);

              chunk_perform_process(kmerhash);
            }
        }
    }

  pthread_mutex_unlock(&mutex_chunks);

  kh_exit(kmerhash);

  return 0;
}


void pair_all()
{
  /* prepare chunks */

  chunk_count = chunk_factor * opt_threads;
  chunk_read_next = 0;
  chunk_process_next = 0;
  chunk_write_next = 0;

  chunks = (chunk_t *) xmalloc(chunk_count * sizeof(chunk_t));

  for (int i = 0; i < chunk_count; i++)
    {
      chunks[i].state = empty;
      chunks[i].size = 0;
      chunks[i].merge_data =
        (merge_data_t *) xmalloc(chunk_size * sizeof(merge_data_t));
      for(int64_t j=0; j<chunk_size; j++)
        init_merge_data(chunks[i].merge_data + j);
    }

  pthread_mutex_init(&mutex_chunks, NULL);
  pthread_cond_init(&cond_chunks, 0);

  /* prepare threads */

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));

  for(int t=0; t<opt_threads; t++)
    if (pthread_create(pthread+t, &attr, pair_worker, (void*)(int64_t)t))
      fatal("Cannot create thread");

  /* wait for threads to terminate */

  for(int t=0; t<opt_threads; t++)
    if (pthread_join(pthread[t], NULL))
      fatal("Cannot join thread");

  /* free threads */

  xfree(pthread);
  pthread_attr_destroy(&attr);

  /* free chunks */

  pthread_cond_destroy(&cond_chunks);
  pthread_mutex_destroy(&mutex_chunks);

  for (int i = 0; i < chunk_count; i++)
    {
      for (int j=0; j < chunk_size; j++)
        free_merge_data(chunks[i].merge_data + j);
      xfree(chunks[i].merge_data);
      chunks[i].merge_data = 0;
    }
  xfree(chunks);
  chunks = 0;
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

  double mean_read_length = sum_read_length / (2.0 * pairs_read);

  if (notmerged > 0)
    fprintf(stderr, "\nPairs that failed merging due to various reasons:\n");

  if (failed_undefined)
    fprintf(stderr,
            "%10" PRIu64 "  undefined reason\n",
            failed_undefined);

  if (failed_minlen)
    fprintf(stderr,
            "%10" PRIu64 "  reads too short (after truncation)\n",
            failed_minlen);

  if (failed_maxlen)
    fprintf(stderr,
            "%10" PRIu64 "  reads too long (after truncation)\n",
            failed_maxlen);

  if (failed_maxns)
    fprintf(stderr,
            "%10" PRIu64 "  too many N's\n",
            failed_maxns);

  if (failed_nokmers)
    fprintf(stderr,
            "%10" PRIu64 "  too few kmers found on same diagonal\n",
            failed_nokmers);

  if (failed_repeat)
    fprintf(stderr,
            "%10" PRIu64 "  potential tandem repeat\n",
            failed_repeat);

  if (failed_maxdiffs)
    fprintf(stderr,
            "%10" PRIu64 "  too many differences\n",
            failed_maxdiffs);

  if (failed_minscore)
    fprintf(stderr,
            "%10" PRIu64 "  alignment score too low, or score drop to high\n",
            failed_minscore);

  if (failed_minovlen)
    fprintf(stderr,
            "%10" PRIu64 "  overlap too short\n",
            failed_minovlen);

  if (failed_maxee)
    fprintf(stderr,
            "%10" PRIu64 "  expected error too high\n",
            failed_maxee);

  if (failed_minmergelen)
    fprintf(stderr,
            "%10" PRIu64 "  merged fragment too short\n",
            failed_minmergelen);

  if (failed_maxmergelen)
    fprintf(stderr,
            "%10" PRIu64 "  merged fragment too long\n",
            failed_maxmergelen);

  if (failed_staggered)
    fprintf(stderr,
            "%10" PRIu64 "  staggered read pairs\n",
            failed_staggered);

  if (failed_indel)
    fprintf(stderr,
            "%10" PRIu64 "  indel errors\n",
            failed_indel);

  fprintf(stderr, "\n");

  fprintf(stderr, "Statistics of all reads:\n");

  fprintf(stderr,
          "%10.2f  Mean read length\n",
          mean_read_length);

  fprintf(stderr, "\n");

  fprintf(stderr, "Statistics of merged reads:\n");

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
