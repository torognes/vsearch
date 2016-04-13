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

#define INPUTCHUNKSIZE 10000
#define SCOREMETHOD 2

/* scores */

const double alpha = 4.0;
const double beta = -5.0;

/* static variables */

static FILE * fp_fastqout = 0;
static FILE * fp_fastaout = 0;
static FILE * fp_fastqout_notmerged_fwd = 0;
static FILE * fp_fastqout_notmerged_rev = 0;
static FILE * fp_fastaout_notmerged_fwd = 0;
static FILE * fp_fastaout_notmerged_rev = 0;
static FILE * fp_eetabbedout = 0;
static fastq_handle fastq_fwd;
static fastq_handle fastq_rev;
static long merged = 0;
static long notmerged = 0;
static long total = 0;
static char * merged_sequence = 0;
static char * merged_quality = 0;
static char * merged_header = 0;
static pthread_t * pthread;
static pthread_attr_t attr;
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_progress;
static pthread_mutex_t mutex_counters;
static pthread_mutex_t mutex_output;
static pthread_cond_t cond_output;
static long output_next;
static char merge_qual_same[128][128];
static char merge_qual_diff[128][128];
static double match_score[128][128];
static double mism_score[128][128];
static double q2p[128];
static double ambig_score;

typedef struct merge_data_s
{
  char * fwd_header;
  char * rev_header;
  char * fwd_sequence;
  char * rev_sequence;
  char * rev_seq_comp;
  char * fwd_quality;
  char * rev_quality;
  long header_alloc;
  long seq_alloc;
  long fwd_length;
  long rev_length;
  long fwd_trunc;
  long rev_trunc;
  long pair_no;
  char * merged_header;
  char * merged_sequence;
  char * merged_quality;
  long merged_header_alloc;
  long merged_seq_alloc;
  double ee_merged;
  double ee_fwd;
  double ee_rev;
  long fwd_errors;
  long rev_errors;
  long offset;
  bool merged;
} merge_data_t;

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

double q_to_p(int q)
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

  ambig_score = alpha * 0.25 + beta * 0.75;

  for (int x = 33; x < 126; x++)
    {
      double px = q_to_p(x);
      q2p[x] = px;

      for (int y = 33; y < 126; y++)
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

          /* Score weights from PEAR */

          /* Match */

          /* probability that they really are similar,
             given that they look similar and have
             error probabilites of px and py, resp. */

          p = 1.0 - px - py + px * py * 4.0 / 3.0;

#if SCOREMETHOD == 2
          match_score[x][y] = alpha * p + beta * (1.0 - p);
#else
          match_score[x][y] = alpha * p;
#endif

          /* Mismatch */

          /* Probability that they really are different,
             given that they look different and have
             error probabilities of px and py, resp. */

          p = 1.0 - (px + py) / 3.0 + px * py * 4.0 / 9.0;

#if SCOREMETHOD == 2
          mism_score[x][y] = alpha * (1.0 - p) + beta * p;
#else
          mism_score[x][y] = beta * p;
#endif

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
    fprintf(fp_eetabbedout, "%.2lf\t%.2lf\t%ld\t%ld\n",
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
  /* The offset is the distance between the (truncated) 3' ends of the two
     sequences */

  long rev_3prime_overhang = ip->offset > ip->fwd_trunc ?
    ip->offset - ip->fwd_trunc : 0;
  long fwd_5prime_overhang = ip->fwd_trunc > ip->offset ?
    ip->fwd_trunc - ip->offset : 0;
  long mergelen = ip->fwd_trunc + ip->rev_trunc - ip->offset;

  ip->ee_merged = 0.0;
  ip->ee_fwd = 0.0;
  ip->ee_rev = 0.0;

  long fwd_pos = 0;
  long rev_pos = ip->rev_trunc - 1 +
    fwd_5prime_overhang - rev_3prime_overhang;

  ip->fwd_errors = 0;
  ip->rev_errors = 0;

  for(long i = 0; i < mergelen; i++)
    {
      bool has_fwd = 0;
      if ((fwd_pos >= 0) && (fwd_pos < ip->fwd_trunc))
        has_fwd = 1;

      bool has_rev = 0;
      if ((rev_pos >= 0) && (rev_pos < ip->rev_trunc))
        has_rev = 1;

      char sym;
      char qual;

      if (has_fwd && has_rev)
        {
          char fwd_sym = ip->fwd_sequence[fwd_pos];
          char rev_sym = ip->rev_seq_comp[rev_pos];
          char fwd_qual = ip->fwd_quality[fwd_pos];
          char rev_qual = ip->rev_quality[rev_pos];

          merge_sym(& sym,
                    & qual,
                    fwd_qual < 2 ? 'N' : fwd_sym,
                    rev_qual < 2 ? 'N' : rev_sym,
                    fwd_qual,
                    rev_qual);

          if (sym != fwd_sym)
            (ip->fwd_errors)++;
          if (sym != rev_sym)
            (ip->rev_errors)++;
        }
      else if (has_fwd)
        {
          sym = ip->fwd_sequence[fwd_pos];
          qual = ip->fwd_quality[fwd_pos];
        }
      else
        {
          sym = ip->rev_seq_comp[rev_pos];
          qual = ip->rev_quality[rev_pos];
        }

      ip->merged_sequence[i] = sym;
      ip->merged_quality[i] = qual;
      ip->ee_merged += q2p[(unsigned)qual];

      if (has_fwd)
        ip->ee_fwd += q2p[(unsigned)(ip->fwd_quality[fwd_pos])];
      if (has_rev)
        ip->ee_rev += q2p[(unsigned)(ip->rev_quality[rev_pos])];

      fwd_pos++;
      rev_pos--;
    }

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

double overlap_score(merge_data_t * ip,
                     long fwd_pos_start,
                     long rev_pos_start,
                     long overlap)
{
  long fwd_pos = fwd_pos_start;
  long rev_pos = rev_pos_start;
  double score = 0.0;
  long diffs = 0;

  for (long j=0; j < overlap; j++)
    {
      char fwd_sym = ip->fwd_sequence[fwd_pos];
      char rev_sym = ip->rev_seq_comp[rev_pos];
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
          if (diffs > opt_fastq_maxdiffs)
            return -1000.0;
        }
    }
  return score;
}

long optimize(merge_data_t * ip)
{
  //  long i1 = opt_fastq_minovlen;
  long i1 = 1;

  i1 = MAX(i1, ip->fwd_trunc + ip->rev_trunc - opt_fastq_maxmergelen);

  long i2 = opt_fastq_allowmergestagger ?
    (ip->fwd_trunc + ip->rev_trunc - opt_fastq_minovlen) :
    ip->fwd_trunc;

  i2 = MIN(i2, ip->fwd_trunc + ip->rev_trunc - opt_fastq_minmergelen);

  long best_i = 0;
  double best_score = 0.0;

  for(long i = i1; i <= i2; i++)
    {
      long fwd_3prime_overhang = i > ip->rev_trunc ? i - ip->rev_trunc : 0;
      long rev_3prime_overhang = i > ip->fwd_trunc ? i - ip->fwd_trunc : 0;
      long overlap = i - fwd_3prime_overhang - rev_3prime_overhang;
      long fwd_pos_start = ip->fwd_trunc - fwd_3prime_overhang - 1;
      long rev_pos_start = ip->rev_trunc - rev_3prime_overhang - overlap;

      double score = overlap_score(ip,
                                   fwd_pos_start,
                                   rev_pos_start,
                                   overlap);

      if (score > best_score)
        {
          best_score = score;
          best_i = i;
        }
    }

  return best_i;
}

void process(merge_data_t * ip)
{
  ip->merged = 0;

  bool skip = 0;

  /* check length */

  if ((ip->fwd_length < opt_fastq_minlen) ||
      (ip->rev_length < opt_fastq_minlen))
    skip = 1;

  /* truncate sequences by quality */

  long fwd_trunc = ip->fwd_length;

  if (!skip)
    {
      for (long i = 0; i < ip->fwd_length; i++)
        if (get_qual(ip->fwd_quality[i]) <= opt_fastq_truncqual)
          {
            fwd_trunc = i;
            break;
          }
      if (fwd_trunc < opt_fastq_minlen)
        skip = 1;
    }

  ip->fwd_trunc = fwd_trunc;

  long rev_trunc = ip->rev_length;

  if (!skip)
    {
      for (long i = 0; i < ip->rev_length; i++)
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
      long fwd_ncount = 0;
      for (long i = 0; i < fwd_trunc; i++)
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
      long rev_ncount = 0;
      for (long i = 0; i < rev_trunc; i++)
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

  //  if (ip->offset)
  if (ip->offset >= opt_fastq_minovlen)
    merge(ip);
}

bool read_pair(merge_data_t * ip)
{
  long suffix_len = opt_label_suffix ? strlen(opt_label_suffix) : 0;

  if (fastq_next(fastq_fwd, 0, chrmap_upcase))
    {
      if (! fastq_next(fastq_rev, 0, chrmap_upcase))
        fatal("More forward reads than reverse reads");

      /* allocate more memory if necessary */

      long fwd_header_len = fastq_get_header_length(fastq_fwd);
      long rev_header_len = fastq_get_header_length(fastq_rev);
      long header_needed = MAX(fwd_header_len, rev_header_len) + 1;

      if (header_needed > ip->header_alloc)
        {
          ip->header_alloc = header_needed;
          ip->fwd_header = (char*) xrealloc(ip->fwd_header, header_needed);
          ip->rev_header = (char*) xrealloc(ip->rev_header, header_needed);
        }

      ip->fwd_length = fastq_get_sequence_length(fastq_fwd);
      ip->rev_length = fastq_get_sequence_length(fastq_rev);
      long seq_needed = MAX(ip->fwd_length, ip->rev_length) + 1;

      if (seq_needed > ip->seq_alloc)
        {
          ip->seq_alloc = seq_needed;
          ip->fwd_sequence = (char*) xrealloc(ip->fwd_sequence, seq_needed);
          ip->rev_sequence = (char*) xrealloc(ip->rev_sequence, seq_needed);
          ip->rev_seq_comp = (char*) xrealloc(ip->rev_seq_comp, seq_needed);
          ip->fwd_quality  = (char*) xrealloc(ip->fwd_quality, seq_needed);
          ip->rev_quality  = (char*) xrealloc(ip->rev_quality, seq_needed);

        }

      long merged_seq_needed = ip->fwd_length + ip->rev_length + 1;

      if (merged_seq_needed > ip->merged_seq_alloc)
        {
          ip->merged_seq_alloc = merged_seq_needed;
          ip->merged_sequence = (char*) xrealloc(ip->merged_sequence,
                                                 merged_seq_needed);
          ip->merged_quality = (char*) xrealloc(ip->merged_quality,
                                                merged_seq_needed);
        }

      long merged_header_needed = fwd_header_len + suffix_len + 1;

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

      /* make a complementary sequence of the reverse sequence */

      for(long i=0; i<ip->rev_length; i++)
        ip->rev_seq_comp[i] = chrmap_complement[(int)(ip->rev_sequence[i])];
      ip->rev_seq_comp[ip->rev_length] = 0;

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

void pair()
{
  merge_data_t * merge_buffer = (merge_data_t*)
    xmalloc(INPUTCHUNKSIZE * sizeof(merge_data_s));

  for(long i=0; i<INPUTCHUNKSIZE; i++)
    {
      merge_data_s * ip = merge_buffer + i;
      ip->fwd_header = 0;
      ip->rev_header = 0;
      ip->fwd_sequence = 0;
      ip->rev_sequence = 0;
      ip->rev_seq_comp = 0;
      ip->fwd_quality = 0;
      ip->rev_quality = 0;
      ip->header_alloc = 0;
      ip->seq_alloc = 0;
      ip->fwd_length = 0;
      ip->rev_length = 0;
      ip->pair_no = 0;
    }

  bool more = 1;
  while (more)
    {
      pthread_mutex_lock(&mutex_input);
      progress_update(fastq_get_position(fastq_fwd));
      long pairs = 0;
      while (more && (pairs < INPUTCHUNKSIZE))
        {
          more = read_pair(merge_buffer + pairs);
          if (more)
            pairs++;
        }
      pthread_mutex_unlock(&mutex_input);

      for(long i=0; i<pairs; i++)
        process(merge_buffer + i);

      pthread_mutex_lock(&mutex_output);
      while(output_next < merge_buffer[0].pair_no)
        pthread_cond_wait(&cond_output, &mutex_output);

      for(long i=0; i<pairs; i++)
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

  for(long i=0; i<INPUTCHUNKSIZE; i++)
    {
      merge_data_s * ip = merge_buffer + i;
      if (ip->fwd_header)
        free(ip->fwd_header);
      if (ip->rev_header)
        free(ip->rev_header);
      if (ip->fwd_sequence)
        free(ip->fwd_sequence);
      if (ip->rev_sequence)
        free(ip->rev_sequence);
      if (ip->rev_seq_comp)
        free(ip->rev_seq_comp);
      if (ip->fwd_quality)
        free(ip->fwd_quality);
      if (ip->rev_quality)
        free(ip->rev_quality);
    }

  free(merge_buffer);
}

void * pair_worker(void * vp)
{
  long t = (long) vp;
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
    if (pthread_create(pthread+t, &attr, pair_worker, (void*)(long)t))
      fatal("Cannot create thread");

  for(int t=0; t<opt_threads; t++)
    if (pthread_join(pthread[t], NULL))
      fatal("Cannot join thread");

  free(pthread);

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

  unsigned long filesize = fastq_get_size(fastq_fwd);
  progress_init("Merging reads", filesize);

  pair_all();

  progress_done();

  if (fastq_next(fastq_rev, 1, chrmap_upcase))
    fatal("More reverse reads than forward reads");

  fprintf(stderr,
          "%10lu  Pairs\n",
          total);

  fprintf(stderr,
          "%10lu  Merged (%.1lf%%)\n",
          merged,
          100.0 * merged / total);

  fprintf(stderr,
          "%10lu  Not merged (%.1lf%%)\n",
          notmerged,
          100.0 * notmerged / total);

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

  if (merged_sequence)
    free(merged_sequence);
  if (merged_quality)
    free(merged_quality);
  if (merged_header)
    free(merged_header);
}
