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

int fastq_get_qual(char q)
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

void filter(bool fastq_only, char * filename)
{
  fastx_handle h = fastx_open(filename);

  if (!h)
    fatal("Unrecognized file type (not proper FASTA or FASTQ format)");

  if (fastq_only && ! h->is_fastq)
    fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");

  if ((opt_fastqout || opt_fastqout_discarded) && ! h->is_fastq)
    fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");

  uint64_t filesize = fastx_get_size(h);

  FILE * fp_fastaout = 0;
  FILE * fp_fastqout = 0;
  FILE * fp_fastaout_discarded = 0;
  FILE * fp_fastqout_discarded = 0;

  if (opt_fastaout)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (!fp_fastaout)
        fatal("Unable to open FASTA output file for writing");
    }

  if (opt_fastqout)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (!fp_fastqout)
        fatal("Unable to open FASTQ output file for writing");
    }

  if (opt_fastaout_discarded)
    {
      fp_fastaout_discarded = fopen_output(opt_fastaout_discarded);
      if (!fp_fastaout_discarded)
        fatal("Unable to open FASTA output file for writing");
    }

  if (opt_fastqout_discarded)
    {
      fp_fastqout_discarded = fopen_output(opt_fastqout_discarded);
      if (!fp_fastqout_discarded)
        fatal("Unable to open FASTQ output file for writing");
    }

  uint64_t header_alloc = 0;
  char * header = 0;
  if (opt_relabel)
    {
      header_alloc = strlen(opt_relabel) + 25;
      header = (char*) xmalloc(header_alloc);
    }

  progress_init("Reading input file", filesize);

  int64_t kept = 0;
  int64_t discarded = 0;
  int64_t truncated = 0;

  while(fastx_next(h, 0, chrmap_no_change))
    {
      int64_t length = fastx_get_sequence_length(h);
      char * d = fastx_get_header(h);
      char * p = fastx_get_sequence(h);
      char * q = fastx_get_quality(h);
      int64_t abundance = fastx_get_abundance(h);
      
      /* strip initial part */
      if (opt_fastq_stripleft > 0)
        {
          if (opt_fastq_stripleft < length)
            {
              p += opt_fastq_stripleft;
              q += opt_fastq_stripleft;
              length -= opt_fastq_stripleft;
            }
          else
            {
              p += length;
              q += length;
              length = 0;
            }
        }

      /* strip right end */
      if (opt_fastq_stripright > 0)
        {
          if (opt_fastq_stripright < length)
            length -= opt_fastq_stripright;
          else
            length = 0;
        }

      /* truncate trailing part */
      if (opt_fastq_trunclen >= 0)
        {
          if (length >= opt_fastq_trunclen)
            length = opt_fastq_trunclen;
          else
            length = 0;
        }
      
      /* truncate trailing part, but keep if short */
      if ((opt_fastq_trunclen_keep >= 0) && (length > opt_fastq_trunclen_keep))
        length = opt_fastq_trunclen_keep;

      /* quality and ee truncation */
      double ee = 0.0;
      if (h->is_fastq)
        {
          for (int64_t i = 0; i < length; i++)
            {
              int qual = fastq_get_qual(q[i]);
              ee += exp10(- qual / 10.0);
              
              if ((qual <= opt_fastq_truncqual) ||
                  (ee > opt_fastq_truncee))
                {
                  ee -= exp10(- qual / 10.0);
                  length = i;
                  break;
                }
            }
        }

      /* count n's */
      int64_t ncount = 0;
      for (int64_t i = 0; i < length; i++)
        {
          int pc = p[i];
          if ((pc == 'N') || (pc == 'n'))
            ncount++;
        }

      if ((length >= opt_fastq_minlen) &&
          (length <= opt_fastq_maxlen) &&
          ((opt_fastq_trunclen < 0) || (length >= opt_fastq_trunclen)) &&
          (ncount <= opt_fastq_maxns) &&
          (ee <= opt_fastq_maxee) &&
          ((length == 0) || (ee / length <= opt_fastq_maxee_rate)) &&
          ((opt_minsize == 0) || (abundance >= opt_minsize)) &&
          ((opt_maxsize == 0) || (abundance <= opt_maxsize)))
        {
          /* keep the sequence */

          kept++;

          if ((uint64_t)(length) < fastx_get_sequence_length(h))
            {
              truncated++;
              p[length] = 0;
              if (h->is_fastq)
                q[length] = 0;
            }

          if (opt_fastaout)
            fasta_print_general(fp_fastaout,
                                0,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                abundance,
                                kept,
                                -1,
                                -1,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);

          if (opt_fastqout)
            fastq_print_general(fp_fastqout,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                q,
                                abundance,
                                kept,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);
        }
      else
        {
          /* discard the sequence */

          discarded++;

          if (opt_fastaout_discarded)
            fasta_print_general(fp_fastaout_discarded,
                                0,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                abundance,
                                discarded,
                                -1,
                                -1,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);

          if (opt_fastqout_discarded)
            fastq_print_general(fp_fastqout_discarded,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                q,
                                abundance,
                                discarded,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);
        }

      progress_update(fastx_get_position(h));
    }
  progress_done();

  if (! opt_quiet)
    fprintf(stderr,
            "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
            kept,
            truncated,
            discarded);

  if (opt_log)
    fprintf(fp_log,
            "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
            kept,
            truncated,
            discarded);

  if (header)
    xfree(header);

  if (opt_fastaout)
    fclose(fp_fastaout);
  
  if (opt_fastqout)
    fclose(fp_fastqout);

  if (opt_fastaout_discarded)
    fclose(fp_fastaout_discarded);
  
  if (opt_fastqout_discarded)
    fclose(fp_fastqout_discarded);

  fastx_close(h);
}

void fastq_filter()
{
  filter(1, opt_fastq_filter);
}

void fastx_filter()
{
  filter(0, opt_fastx_filter);
}

void fastq_chars()
{
  uint64_t sequence_chars[256];
  uint64_t quality_chars[256];
  uint64_t tail_chars[256];
  uint64_t total_chars = 0;
  int maxrun[256];

  for(int c=0; c<256; c++)
    {
      sequence_chars[c] = 0;
      quality_chars[c] = 0;
      tail_chars[c] = 0;
      maxrun[c] = 0;
    }

  fastx_handle h = fastq_open(opt_fastq_chars);

  uint64_t filesize = fastq_get_size(h);

  progress_init("Reading FASTQ file", filesize);

  uint64_t seq_count = 0;
  
  int qmin_n = 255, qmax_n = 0;

  while(fastq_next(h, 0, chrmap_upcase))
    {
      int64_t len = fastq_get_sequence_length(h);
      char * p = fastq_get_sequence(h);
      char * q = fastq_get_quality(h);

      seq_count++;
      total_chars += len;

      int run_char = -1;
      int run = 0;
      
      int64_t i = 0;
      while(i<len)
        {
          int pc = *p++;
          int qc = *q++;
          sequence_chars[pc]++;
          quality_chars[qc]++;
          
          if ((pc == 'N') || (pc == 'n'))
            {
              if (qc < qmin_n)
                qmin_n = qc;
              if (qc > qmax_n)
                qmax_n = qc;
            }

          if (pc == run_char)
            {
              run++;
              if (run > maxrun[run_char])
                maxrun[run_char] = run;
            }
          else
            {
              run_char = pc;
              run = 0;
            }
          
          i++;
        }

      if (len >= opt_fastq_tail)
        {
          q = fastq_get_quality(h) + len - 1;
          int tail_char = *q--;
          int tail_len = 1;
          while(*q-- == tail_char)
            {
              tail_len++;
              if (tail_len >= opt_fastq_tail)
                break;
            }
          if (tail_len >= opt_fastq_tail)
            tail_chars[tail_char]++;
        }

      progress_update(fastq_get_position(h));
    }
  progress_done();

  fastq_close(h);
  
  char qmin = 0;
  char qmax = 0;

  for(int c=0; c<=255; c++)
    {
      if (quality_chars[c])
        {
          qmin = c;
          break;
        }
    }

  for(int c=255; c>=0; c--)
    {
      if (quality_chars[c])
        {
          qmax = c;
          break;
        }
    }

  char fastq_ascii, fastq_qmin, fastq_qmax;

  if (qmin < 59)
    fastq_ascii = 33;
  else
    fastq_ascii = 64;

  fastq_qmax = qmax - fastq_ascii;
  fastq_qmin = qmin - fastq_ascii;

  if (!opt_quiet)
    {
      fprintf(stderr, "Read %" PRIu64 " sequences.\n", seq_count);

      fprintf(stderr, "Qmin %d, QMax %d, Range %d\n",
              qmin, qmax, qmax-qmin+1);

      fprintf(stderr, "Guess: -fastq_qmin %d -fastq_qmax %d -fastq_ascii %d\n",
              fastq_qmin, fastq_qmax, fastq_ascii);

      if (fastq_ascii == 64)
        {
          if (qmin < 64)
            fprintf(stderr, "Guess: Solexa format (phred+64)\n");
          else if (qmin < 66)
            fprintf(stderr, "Guess: Illumina 1.3+ format (phred+64)\n");
          else
            fprintf(stderr, "Guess: Illumina 1.5+ format (phred+64)\n");
        }
      else
        {
          if (qmax > 73)
            fprintf(stderr, "Guess: Illumina 1.8+ format (phred+33)\n");
          else
            fprintf(stderr, "Guess: Original Sanger format (phred+33)\n");
        }

      fprintf(stderr, "\n");
      fprintf(stderr, "Letter          N   Freq MaxRun\n");
      fprintf(stderr, "------ ---------- ------ ------\n");

      for(int c=0; c<256; c++)
        {
          if (sequence_chars[c] > 0)
            {
              fprintf(stderr, "     %c %10" PRIu64 " %5.1f%% %6d",
                      c,
                      sequence_chars[c],
                      100.0 * sequence_chars[c] / total_chars,
                      maxrun[c]);
              if ((c == 'N') || (c == 'n'))
                {
                  if (qmin_n < qmax_n)
                    fprintf(stderr, "  Q=%c..%c", qmin_n, qmax_n);
                  else
                    fprintf(stderr, "  Q=%c", qmin_n);
                }
              fprintf(stderr, "\n");
            }
        }

      fprintf(stderr, "\n");
      fprintf(stderr, "Char  ASCII    Freq       Tails\n");
      fprintf(stderr, "----  -----  ------  ----------\n");

      for(int c=qmin; c<=qmax; c++)
        {
          if (quality_chars[c] > 0)
            {
              fprintf(stderr, " '%c'  %5d  %5.1f%%  %10" PRIu64 "\n",
                      c,
                      c,
                      100.0 * quality_chars[c] / total_chars,
                      tail_chars[c]);
            }
        }
    }

  if (opt_log)
    {
      fprintf(fp_log, "Read %" PRIu64 " sequences.\n", seq_count);

      fprintf(fp_log, "Qmin %d, QMax %d, Range %d\n",
              qmin, qmax, qmax-qmin+1);

      fprintf(fp_log, "Guess: -fastq_qmin %d -fastq_qmax %d -fastq_ascii %d\n",
              fastq_qmin, fastq_qmax, fastq_ascii);

      if (fastq_ascii == 64)
        {
          if (qmin < 64)
            fprintf(fp_log, "Guess: Solexa format (phred+64)\n");
          else if (qmin < 66)
            fprintf(fp_log, "Guess: Illumina 1.3+ format (phred+64)\n");
          else
            fprintf(fp_log, "Guess: Illumina 1.5+ format (phred+64)\n");
        }
      else
        {
          if (qmax > 73)
            fprintf(fp_log, "Guess: Illumina 1.8+ format (phred+33)\n");
          else
            fprintf(fp_log, "Guess: Original Sanger format (phred+33)\n");
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "Letter          N   Freq MaxRun\n");
      fprintf(fp_log, "------ ---------- ------ ------\n");

      for(int c=0; c<256; c++)
        {
          if (sequence_chars[c] > 0)
            {
              fprintf(fp_log, "     %c %10" PRIu64 " %5.1f%% %6d",
                      c,
                      sequence_chars[c],
                      100.0 * sequence_chars[c] / total_chars,
                      maxrun[c]);
              if ((c == 'N') || (c == 'n'))
                {
                  if (qmin_n < qmax_n)
                    fprintf(fp_log, "  Q=%c..%c", qmin_n, qmax_n);
                  else
                    fprintf(fp_log, "  Q=%c", qmin_n);
                }
              fprintf(fp_log, "\n");
            }
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "Char  ASCII    Freq       Tails\n");
      fprintf(fp_log, "----  -----  ------  ----------\n");

      for(int c=qmin; c<=qmax; c++)
        {
          if (quality_chars[c] > 0)
            {
              fprintf(fp_log, " '%c'  %5d  %5.1f%%  %10" PRIu64 "\n",
                      c,
                      c,
                      100.0 * quality_chars[c] / total_chars,
                      tail_chars[c]);
            }
        }
    }
}

double q2p(double q)
{
  return exp10(- q / 10.0);
}

void fastq_stats()
{
  fastx_handle h = fastq_open(opt_fastq_stats);

  uint64_t filesize = fastq_get_size(h);

  progress_init("Reading FASTQ file", filesize);

  uint64_t seq_count = 0;
  uint64_t symbols = 0;
  
  int64_t read_length_alloc = 512;

  uint64_t * read_length_table = (uint64_t*) xmalloc(sizeof(uint64_t) * read_length_alloc);
  memset(read_length_table, 0, sizeof(uint64_t) * read_length_alloc);

  uint64_t * qual_length_table = (uint64_t*) xmalloc(sizeof(uint64_t) * read_length_alloc * 256);
  memset(qual_length_table, 0, sizeof(uint64_t) * read_length_alloc * 256);

  uint64_t * ee_length_table = (uint64_t *) xmalloc(sizeof(uint64_t) * read_length_alloc * 4);
  memset(ee_length_table, 0, sizeof(uint64_t) * read_length_alloc * 4);

  uint64_t * q_length_table = (uint64_t *) xmalloc(sizeof(uint64_t) * read_length_alloc * 4);
  memset(q_length_table, 0, sizeof(uint64_t) * read_length_alloc * 4);

  double * sumee_length_table = (double *) xmalloc(sizeof(double) * read_length_alloc);
  memset(sumee_length_table, 0, sizeof(double) * read_length_alloc);

  int64_t len_min = LONG_MAX;
  int64_t len_max = 0;
  
  int qmin = +1000;
  int qmax = -1000;

  uint64_t quality_chars[256];
  for(int c=0; c<256; c++)
    quality_chars[c] = 0;
  
  while(fastq_next(h, 0, chrmap_upcase))
    {
      seq_count++;

      int64_t len = fastq_get_sequence_length(h);
      char * q = fastq_get_quality(h);

      /* update length statistics */

      if (len+1 > read_length_alloc)
        {
          read_length_table = (uint64_t*) xrealloc(read_length_table,
                                              sizeof(uint64_t) * (len+1));
          memset(read_length_table + read_length_alloc, 0, 
                 sizeof(uint64_t) * (len + 1 - read_length_alloc));

          qual_length_table = (uint64_t*) xrealloc(qual_length_table,
                                              sizeof(uint64_t) * (len+1) * 256);
          memset(qual_length_table + 256 * read_length_alloc, 0, 
                 sizeof(uint64_t) * (len + 1 - read_length_alloc) * 256);

          ee_length_table = (uint64_t*) xrealloc(ee_length_table,
                                            sizeof(uint64_t) * (len+1) * 4);
          memset(ee_length_table + 4 * read_length_alloc, 0, 
                 sizeof(uint64_t) * (len + 1 - read_length_alloc) * 4);

          q_length_table = (uint64_t*) xrealloc(q_length_table,
                                           sizeof(uint64_t) * (len+1) * 4);
          memset(q_length_table + 4 * read_length_alloc, 0, 
                 sizeof(uint64_t) * (len + 1 - read_length_alloc) * 4);

          sumee_length_table = (double *) xrealloc(sumee_length_table,
                                                   sizeof(double) * (len+1));
          memset(sumee_length_table + read_length_alloc, 0,
                 sizeof(double) * (len + 1 - read_length_alloc));

          read_length_alloc = len + 1;
        }

      read_length_table[len]++;

      if (len < len_min)
        len_min = len;
      if (len > len_max)
        len_max = len;
      
      /* update quality statistics */
      
      symbols += len;

      double ee_limit[4] = { 1.0, 0.5, 0.25, 0.1 };
      
      double ee = 0.0;
      int qmin_this = 1000;
      for(int64_t i=0; i < len; i++)
        {
          int qc = q[i];

          int qual = qc - opt_fastq_ascii;
          if ((qual < opt_fastq_qmin) || (qual > opt_fastq_qmax))
            {
              char * msg;
              if (xsprintf(& msg,
                           "FASTQ quality value (%d) out of range (%" PRId64 "-%" PRId64 ").\n"
                           "Please adjust the FASTQ quality base character or range with the\n"
                           "--fastq_ascii, --fastq_qmin or --fastq_qmax options. For a complete\n"
                           "diagnosis with suggested values, please run vsearch --fastq_chars file.",
                           qual, opt_fastq_qmin, opt_fastq_qmax) > 0)
                fatal(msg);
              else
                fatal("Out of memory");
              xfree(msg);
            }
          
          quality_chars[qc]++;
          if (qc < qmin)
            qmin = qc;
          if (qc > qmax)
            qmax = qc;

          qual_length_table[256*i + qc]++;

          ee += q2p(qual);
          
          sumee_length_table[i] += ee;

          for(int z=0; z<4; z++)
            {
              if (ee <= ee_limit[z])
                ee_length_table[4*i+z]++;
              else
                break;
            }

          if (qual < qmin_this)
            qmin_this = qual;

          for(int z=0; z<4; z++)
            {
              if (qmin_this > 5*(z+1)) 
                q_length_table[4*i+z]++;
              else
                break;
            }
        }

      progress_update(fastq_get_position(h));
    }
  progress_done();

  /* compute various distributions */

  uint64_t * length_dist = (uint64_t*) xmalloc(sizeof(uint64_t) * (len_max+1));
  int64_t * symb_dist = (int64_t*) xmalloc(sizeof(int64_t) * (len_max+1));

  double * rate_dist = (double*) xmalloc(sizeof(double) * (len_max+1));
  double * avgq_dist = (double*) xmalloc(sizeof(double) * (len_max+1));
  double * avgee_dist = (double*) xmalloc(sizeof(double) * (len_max+1));
  double * avgp_dist = (double*) xmalloc(sizeof(double) * (len_max+1));

  int64_t length_accum = 0;
  int64_t symb_accum = 0;

  for(int64_t i = 0; i <= len_max; i++)
    {
      length_accum += read_length_table[i];
      length_dist[i] = length_accum;

      symb_accum += seq_count - length_accum;
      symb_dist[i] = symb_accum;

      int64_t q = 0;
      int64_t x = 0;
      double e_sum = 0.0;
      for(int c=qmin; c<=qmax; c++)
        {
          int qual = c - opt_fastq_ascii;
          x += qual_length_table[256*i + c];
          q += qual_length_table[256*i + c] * qual;
          e_sum += qual_length_table[256*i + c] * q2p(qual);
        }
      avgq_dist[i] = 1.0 * q / x;
      avgp_dist[i] = e_sum / x;
      avgee_dist[i] = sumee_length_table[i] / x;
      rate_dist[i] = avgee_dist[i] / (i+1);
    }

  if (fp_log)
    {
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Read length distribution\n");
      fprintf(fp_log, "      L           N      Pct   AccPct\n");
      fprintf(fp_log, "-------  ----------  -------  -------\n");
      
      for(int64_t i = len_max; i >= len_min; i--)
        {
          if (read_length_table[i] > 0)
            fprintf(fp_log, "%2s%5" PRId64 "  %10" PRIu64 "   %5.1lf%%   %5.1lf%%\n",
                    (i == len_max ? ">=" : "  "),
                    i,
                    read_length_table[i],
                    read_length_table[i] * 100.0 / seq_count,
                    100.0 * (seq_count - length_dist[i-1]) / seq_count);
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "Q score distribution\n");
      fprintf(fp_log, "ASCII    Q       Pe           N      Pct   AccPct\n");
      fprintf(fp_log, "-----  ---  -------  ----------  -------  -------\n");

      int64_t qual_accum = 0;
      for(int c = qmax ; c >= qmin ; c--)
        {
          if (quality_chars[c] > 0)
            {
              qual_accum += quality_chars[c];
              fprintf(fp_log,
                      "    %c  %3" PRId64 "  %7.5lf  %10" PRIu64 "  %6.1lf%%  %6.1lf%%\n",
                      c,
                      c - opt_fastq_ascii,
                      q2p(c - opt_fastq_ascii),
                      quality_chars[c],
                      100.0 * quality_chars[c] / symbols,
                      100.0 * qual_accum / symbols);
            }
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct\n");
      fprintf(fp_log, "-----  -------  ----  -------  --------  -----  ---------  --------\n");

      for(int64_t i = 2; i <= len_max; i++)
        {
          double PctRecs = 100.0 * (seq_count - length_dist[i-1]) / seq_count;
          double AvgQ = avgq_dist[i-1];
          double AvgP = avgp_dist[i-1];
          double AvgEE = avgee_dist[i-1];
          double Rate = rate_dist[i-1];

          fprintf(fp_log,
                  "%5" PRId64 "  %6.1lf%%  %4.1lf  %7.5lf  %8.6lf  %5.2lf  %9.6lf  %7.3lf%%\n",
                  i,
                  PctRecs,
                  AvgQ,
                  q2p(AvgQ),
                  AvgP,
                  AvgEE,
                  Rate,
                  100.0 * Rate);
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000\n");
      fprintf(fp_log, "-----  -------  -------  -------  -------  -------  -------  -------  -------\n");
      
      for(int64_t i = len_max; i >= 1; i--)
        {
          int64_t read_count[4];
          double read_percentage[4];
          
          for(int z=0; z<4; z++)
            {
              read_count[z] = ee_length_table[4*(i-1)+z];
              read_percentage[z] = 100.0 * read_count[z] / seq_count;
            }
          
          if (read_count[0] > 0)
            {
              fprintf(fp_log,
                      "%5" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  "
                      "%6.2lf%%  %6.2lf%%  %6.2lf%%  %6.2lf%%\n",
                      i,
                      read_count[0], read_count[1],
                      read_count[2], read_count[3],
                      read_percentage[0], read_percentage[1],
                      read_percentage[2], read_percentage[3]);
            }
        }

      
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Truncate at first Q\n");
      fprintf(fp_log, "  Len     Q=5    Q=10    Q=15    Q=20\n");
      fprintf(fp_log, "-----  ------  ------  ------  ------\n");

      for(int64_t i = len_max; i >= len_max/2; i--)
        {
          double read_percentage[4];
          
          for(int z=0; z<4; z++)
            read_percentage[z] = 100.0 * q_length_table[4*(i-1)+z] / seq_count;

          fprintf(fp_log, "%5" PRId64 "  %5.1lf%%  %5.1lf%%  %5.1lf%%  %5.1lf%%\n",
                  i,
                  read_percentage[0], read_percentage[1],
                  read_percentage[2], read_percentage[3]);
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "%10" PRIu64 "  Recs (%.1lfM), 0 too long\n",
              seq_count, seq_count / 1.0e6);
      fprintf(fp_log, "%10.1lf  Avg length\n", 1.0 * symbols / seq_count);
      fprintf(fp_log, "%9.1lfM  Bases\n", symbols / 1.0e6);
    }
  
  xfree(read_length_table);
  xfree(qual_length_table);
  xfree(ee_length_table);
  xfree(q_length_table);
  xfree(sumee_length_table);

  xfree(length_dist);
  xfree(symb_dist);
  xfree(rate_dist);
  xfree(avgq_dist);
  xfree(avgee_dist);
  xfree(avgp_dist);

  fastq_close(h);
  
  fprintf(stderr, "Read %" PRIu64 " sequences.\n", seq_count);
}

void fastx_revcomp()
{
  uint64_t buffer_alloc = 512;
  char * seq_buffer = (char*) xmalloc(buffer_alloc);
  char * qual_buffer = (char*) xmalloc(buffer_alloc);

  uint64_t header_alloc = 512;
  char * header = (char*) xmalloc(header_alloc);

  uint64_t suffix_length = opt_label_suffix ? strlen(opt_label_suffix) : 0;

  fastx_handle h = fastx_open(opt_fastx_revcomp);

  if (!h)
    fatal("Unrecognized file type (not proper FASTA or FASTQ format)");

  if (opt_fastqout && ! h->is_fastq)
    fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");

  uint64_t filesize = fastx_get_size(h);

  FILE * fp_fastaout = 0;
  FILE * fp_fastqout = 0;

  if (opt_fastaout)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (!fp_fastaout)
        fatal("Unable to open FASTA output file for writing");
    }

  if (opt_fastqout)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (!fp_fastqout)
        fatal("Unable to open FASTQ output file for writing");
    }

  if (h->is_fastq)
    progress_init("Reading FASTQ file", filesize);
  else
    progress_init("Reading FASTA file", filesize);

  int count = 0;

  while(fastx_next(h, 0, chrmap_no_change))
    {
      count++;
      
      /* header */
      
      uint64_t hlen = fastx_get_header_length(h);

      if (hlen + suffix_length + 1 > header_alloc)
        {
          header_alloc = hlen + suffix_length + 1;
          header = (char*) xrealloc(header, header_alloc);
        }

      char * d = fastx_get_header(h);

      if (opt_label_suffix)
        snprintf(header, header_alloc, "%s%s", d, opt_label_suffix);
      else
        snprintf(header, header_alloc, "%s", d);


      /* sequence */

      uint64_t length = fastx_get_sequence_length(h);

      if (length + 1 > buffer_alloc)
        {
          buffer_alloc = length + 1;
          seq_buffer = (char *) xrealloc(seq_buffer, buffer_alloc);
          qual_buffer = (char *) xrealloc(qual_buffer, buffer_alloc);
        }

      char * p = fastx_get_sequence(h);
      reverse_complement(seq_buffer, p, length);


      /* quality values */

      char * q = fastx_get_quality(h);

      if (fastx_is_fastq(h))
        {
          /* reverse quality values */
          for(uint64_t i=0; i<length; i++)
            qual_buffer[i] = q[length-1-i];
          qual_buffer[length] = 0;
        }

      if (opt_fastaout)
        fasta_print_general(fp_fastaout,
                            0,
                            seq_buffer,
                            length,
                            header,
                            hlen,
                            0,
                            count,
                            -1, -1, 0, 0.0);

      if (opt_fastqout)
        fastq_print_general(fp_fastqout,
                            seq_buffer,
                            length,
                            header,
                            hlen,
                            qual_buffer,
                            0,
                            count,
                            0, 0.0);
      
      progress_update(fastx_get_position(h));
    }
  progress_done();
      
  if (opt_fastaout)
    fclose(fp_fastaout);
  
  if (opt_fastqout)
    fclose(fp_fastqout);

  fastx_close(h);

  xfree(header);
  xfree(seq_buffer);
  xfree(qual_buffer);
}

void fastq_convert()
{
  fastx_handle h = fastq_open(opt_fastq_convert);

  if (!h)
    fatal("Unable to open FASTQ file");

  uint64_t filesize = fastq_get_size(h);

  FILE * fp_fastqout = 0;

  fp_fastqout = fopen_output(opt_fastqout);
  if (!fp_fastqout)
    fatal("Unable to open FASTQ output file for writing");

  progress_init("Reading FASTQ file", filesize);

  while(fastq_next(h, 0, chrmap_no_change))
    {
      /* header */
      
      char * header = fastq_get_header(h);

      /* sequence */

      uint64_t length = fastq_get_sequence_length(h);
      char * sequence = fastq_get_sequence(h);

      /* convert quality values */
      
      char * quality = fastq_get_quality(h);
      for(uint64_t i=0; i<length; i++)
        {
          int q = quality[i] - opt_fastq_ascii;
          if (q < opt_fastq_qmin)
            {
              fprintf(stderr,
                      "\nFASTQ quality score (%d) below minimum (%" PRId64 ") in entry no %" PRIu64 " starting on line %" PRIu64 "\n", 
                      q,
                      opt_fastq_qmin,
                      fastq_get_seqno(h) + 1,
                      fastq_get_lineno(h));
              fatal("FASTQ quality score too low");
            }
          if (q > opt_fastq_qmax)
            {
              fprintf(stderr,
                      "\nFASTQ quality score (%d) above maximum (%" PRId64 ") in entry no %" PRIu64 " starting on line %" PRIu64 "\n", 
                      q,
                      opt_fastq_qmax,
                      fastq_get_seqno(h) + 1,
                      fastq_get_lineno(h));
              fatal("FASTQ quality score too high");
            }
          if (q < opt_fastq_qminout)
            q = opt_fastq_qminout;
          if (q > opt_fastq_qmaxout)
            q = opt_fastq_qmaxout;
          q += opt_fastq_asciiout;
          if (q < 33)
            q = 33;
          if (q > 126)
            q = 126;
          quality[i] = q;
        }
      quality[length] = 0;

      fastq_print(fp_fastqout, header, sequence, quality);
      
      progress_update(fastq_get_position(h));
    }

  progress_done();

  fclose(fp_fastqout);
  fastq_close(h);
}
