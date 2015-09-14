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

int fastq_get_qual(char q)
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

void fastq_filter()
{
  fastq_handle h = fastq_open(opt_fastq_filter);

  unsigned long filesize = fastq_get_size(h);

  FILE * fp_fastaout = 0;
  FILE * fp_fastqout = 0;
  FILE * fp_fastaout_discarded = 0;
  FILE * fp_fastqout_discarded = 0;

  if (opt_fastaout)
    {
      fp_fastaout = fopen(opt_fastaout, "w");
      if (!fp_fastaout)
        fatal("Unable to open fasta output file for writing");
    }

  if (opt_fastqout)
    {
      fp_fastqout = fopen(opt_fastqout, "w");
      if (!fp_fastqout)
        fatal("Unable to open fastq output file for writing");
    }

  if (opt_fastaout_discarded)
    {
      fp_fastaout_discarded = fopen(opt_fastaout_discarded, "w");
      if (!fp_fastaout_discarded)
        fatal("Unable to open fasta output file for writing");
    }

  if (opt_fastqout_discarded)
    {
      fp_fastqout_discarded = fopen(opt_fastqout_discarded, "w");
      if (!fp_fastqout_discarded)
        fatal("Unable to open fastq output file for writing");
    }

  unsigned long header_alloc = 0;
  char * header = 0;
  if (opt_relabel)
    {
      header_alloc = strlen(opt_relabel) + 25;
      header = (char*) xmalloc(header_alloc);
    }

  progress_init("Reading fastq file", filesize);

  long kept = 0;
  long discarded = 0;
  long truncated = 0;

  char hex_md5[LEN_HEX_DIG_MD5];
  char hex_sha1[LEN_HEX_DIG_SHA1];

  while(fastq_next(h, 0, chrmap_upcase))
    {
      long length = fastq_get_sequence_length(h);
      char * d = fastq_get_header(h);
      char * p = fastq_get_sequence(h);
      char * q = fastq_get_quality(h);

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
      
      /* truncate trailing part */
      if (opt_fastq_trunclen > 0)
        {
          if (length >= opt_fastq_trunclen)
            length = opt_fastq_trunclen;
          else
            length = 0;
        }
      
      /* quality truncation */
      for (long i = 0; i < length; i++)
        {
          int qual = fastq_get_qual(q[i]);

          if (qual <= opt_fastq_truncqual)
            {
              length = i;
              break;
            }
        }

      /* count n's */
      long ncount = 0;
      for (long i = 0; i < length; i++)
        {
          int pc = p[i];
          if ((pc == 'N') || (pc == 'n'))
            ncount++;
        }

      /* compute ee */
      double ee = 0.0;
      for (long i = 0; i < length; i++)
        {
          int qual = fastq_get_qual(q[i]);
          ee += pow(10.0, - qual / 10.0);
        }

      if ((length >= opt_fastq_minlen) &&
          ((opt_fastq_trunclen == 0) || (length >= opt_fastq_trunclen)) &&
          (ee <= opt_fastq_maxee) &&
          (ee / length <= opt_fastq_maxee_rate) &&
          (ncount <= opt_fastq_maxns))
        {
          /* keep the sequence */
          kept++;
          if ((unsigned long)(length) < fastq_get_sequence_length(h))
            {
              truncated++;
              p[length] = 0;
              q[length] = 0;
            }

          if (opt_relabel)
            {
              (void) snprintf(header, header_alloc,
                              "%s%ld", opt_relabel, kept);
              d = header;
            }
          else if (opt_relabel_md5)
            {
              get_hex_seq_digest_md5(hex_md5, p, length);
              d = hex_md5;
            }
          else if (opt_relabel_sha1)
            {
              get_hex_seq_digest_sha1(hex_sha1, p, length);
              d = hex_sha1;
            }

          if (opt_fastaout)
            {
              fprint_fasta_hdr_only(fp_fastaout, d);
              fprint_fasta_seq_only(fp_fastaout, p, length, opt_fasta_width);
            }
          if (opt_fastqout)
            fprint_fastq(fp_fastqout, d, p, q, opt_eeout, ee);
        }
      else
        {
          discarded++;
          p = fastq_get_sequence(h);
          q = fastq_get_quality(h);

          if (opt_relabel)
            {
              (void) snprintf(header, header_alloc, "%s%ld", opt_relabel, discarded);
              d = header;
            }
          else if (opt_relabel_md5)
            {
              get_hex_seq_digest_md5(hex_md5, p, length);
              d = hex_md5;
            }
          else if (opt_relabel_sha1)
            {
              get_hex_seq_digest_sha1(hex_sha1, p, length);
              d = hex_sha1;
            }

          if (opt_fastaout_discarded)
            {
              fprint_fasta_hdr_only(fp_fastaout_discarded, d);
              fprint_fasta_seq_only(fp_fastaout_discarded,
                                    p,
                                    length,
                                    opt_fasta_width);
            }
          if (opt_fastqout_discarded)
            fprint_fastq(fp_fastqout_discarded, d, p, q, opt_eeout, ee);
        }

      progress_update(fastq_get_position(h));
    }
  progress_done();

  fprintf(stderr,
          "%ld sequences kept (of which %ld truncated), %ld sequences discarded.\n",
          kept,
          truncated,
          discarded);

  if (header)
    free(header);

  if (opt_fastaout)
    fclose(fp_fastaout);
  
  if (opt_fastqout)
    fclose(fp_fastqout);

  if (opt_fastaout_discarded)
    fclose(fp_fastaout_discarded);
  
  if (opt_fastqout_discarded)
    fclose(fp_fastqout_discarded);

  fastq_close(h);
}

void fastq_chars()
{
  unsigned long sequence_chars[256];
  unsigned long quality_chars[256];
  unsigned long tail_chars[256];
  unsigned long total_chars = 0;
  int maxrun[256];

  for(int c=0; c<256; c++)
    {
      sequence_chars[c] = 0;
      quality_chars[c] = 0;
      tail_chars[c] = 0;
      maxrun[c] = 0;
    }

  fastq_handle h = fastq_open(opt_fastq_chars);

  unsigned long filesize = fastq_get_size(h);

  progress_init("Reading fastq file", filesize);

  unsigned long seq_count = 0;
  
  int qmin_n = 255, qmax_n = 0;

  while(fastq_next(h, 0, chrmap_upcase))
    {
      long len = fastq_get_sequence_length(h);
      char * p = fastq_get_sequence(h);
      char * q = fastq_get_quality(h);

      seq_count++;
      total_chars += len;

      int run_char = -1;
      int run = 0;
      
      long i = 0;
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
  
  fprintf(stderr, "Read %lu sequences.\n", seq_count);

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

  fprintf(stderr, "Qmin %d, QMax %d, Range %d\n",
          qmin, qmax, qmax-qmin+1);

  fastq_qmax = qmax - fastq_ascii;
  fastq_qmin = qmin - fastq_ascii;
  
  fprintf(stderr, "Guess: -fastq_qmin %d -fastq_qmax %d -fastq_ascii %d\n",
          fastq_qmin, fastq_qmax, fastq_ascii);
  
  if (fastq_ascii == 64)
    {
      if (qmin < 64)
        fprintf(stderr, "Guess: Solexa format\n");
      else if (qmin < 66)
        fprintf(stderr, "Guess: Illumina 1.3+ format\n");
      else
        fprintf(stderr, "Guess: Illumina 1.5+ format\n");
    }
  else
    {
      fprintf(stderr, "Guess: Sanger / Illumina 1.8+ format\n");
    }

  fprintf(stderr, "\n");
  fprintf(stderr, "Letter          N   Freq MaxRun\n");
  fprintf(stderr, "------ ---------- ------ ------\n");
  
  for(int c=0; c<256; c++)
    {
      if (sequence_chars[c] > 0)
        {
          fprintf(stderr, "     %c %10lu %5.1f%% %6d",
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
          fprintf(stderr, " '%c'  %5d  %5.1f%%  %10lu\n",
                  c,
                  c,
                  100.0 * quality_chars[c] / total_chars,
                  tail_chars[c]);
        }
    }
}

double q2p(double q)
{
  return pow(10.0, - q / 10.0);
}

void fastq_stats()
{
  fastq_handle h = fastq_open(opt_fastq_stats);

  unsigned long filesize = fastq_get_size(h);

  progress_init("Reading fastq file", filesize);

  unsigned long seq_count = 0;
  unsigned long symbols = 0;
  
  long read_length_alloc = 512;
  int * read_length_table = (int*) xmalloc(sizeof(int) * read_length_alloc);
  memset(read_length_table, 0, sizeof(int) * read_length_alloc);
  int * qual_length_table = (int*) xmalloc(sizeof(int) * read_length_alloc * 256);
  memset(qual_length_table, 0, sizeof(int) * read_length_alloc * 256);
  int * ee_length_table = (int *) xmalloc(sizeof(int) * read_length_alloc * 4);
  memset(ee_length_table, 0, sizeof(int) * read_length_alloc * 4);
  int * q_length_table = (int *) xmalloc(sizeof(int) * read_length_alloc * 4);
  memset(q_length_table, 0, sizeof(int) * read_length_alloc * 4);

  long len_min = LONG_MAX;
  long len_max = 0;
  
  int qmin = +1000;
  int qmax = -1000;

  unsigned long quality_chars[256];
  for(int c=0; c<256; c++)
    quality_chars[c] = 0;
  
  while(fastq_next(h, 0, chrmap_upcase))
    {
      seq_count++;

      long len = fastq_get_sequence_length(h);
      char * q = fastq_get_quality(h);

      /* update length statistics */

      if (len+1 > read_length_alloc)
        {
          read_length_table = (int*) xrealloc(read_length_table,
                                              sizeof(int) * (len+1));
          memset(read_length_table + read_length_alloc, 0, 
                 sizeof(int) * (len + 1 - read_length_alloc));
          qual_length_table = (int*) xrealloc(qual_length_table,
                                              sizeof(int) * (len+1) * 256);
          memset(qual_length_table + 256 * read_length_alloc, 0, 
                 sizeof(int) * (len + 1 - read_length_alloc) * 256);
          ee_length_table = (int*) xrealloc(ee_length_table,
                                            sizeof(int) * (len+1) * 4);
          memset(ee_length_table + 4 * read_length_alloc, 0, 
                 sizeof(int) * (len + 1 - read_length_alloc) * 4);
          q_length_table = (int*) xrealloc(q_length_table,
                                           sizeof(int) * (len+1) * 4);
          memset(q_length_table + 4 * read_length_alloc, 0, 
                 sizeof(int) * (len + 1 - read_length_alloc) * 4);
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
      for(long i=0; i < len; i++)
        {
          int qc = q[i];

          int qual = qc - opt_fastq_ascii;
          if ((qual < opt_fastq_qmin) || (qual > opt_fastq_qmax))
            fatal("FASTQ quality value out of range");
          
          quality_chars[qc]++;
          if (qc < qmin)
            qmin = qc;
          if (qc > qmax)
            qmax = qc;

          qual_length_table[256*i + qc]++;

          ee += q2p(qual);
          
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

  int * length_dist = (int*) xmalloc(sizeof(int) * (len_max+1));
  long * symb_dist = (long*) xmalloc(sizeof(long) * (len_max+1));

  double * rate_dist = (double*) xmalloc(sizeof(double) * (len_max+1));
  double * avgq_dist = (double*) xmalloc(sizeof(double) * (len_max+1));
  double * avgee_dist = (double*) xmalloc(sizeof(double) * (len_max+1));
  double * avgp_dist = (double*) xmalloc(sizeof(double) * (len_max+1));

  long length_accum = 0;
  long symb_accum = 0;
  double p_sum = 0.0;

  for(long i = 0; i <= len_max; i++)
    {
      length_accum += read_length_table[i];
      length_dist[i] = length_accum;

      symb_accum += seq_count - length_accum;
      symb_dist[i] = symb_accum;

      long q = 0;
      long x = 0;
      double e_sum = 0.0;
      for(int c=qmin; c<=qmax; c++)
        {
          int qual = c - opt_fastq_ascii;
          x += qual_length_table[256*i + c];

          q += qual_length_table[256*i + c] * qual;

          p_sum += qual_length_table[256*i + c] * q2p(qual);
          e_sum += qual_length_table[256*i + c] * q2p(qual);
        }
      avgq_dist[i] = 1.0 * q / x;
      avgp_dist[i] = e_sum / x;
      avgee_dist[i] = 1.0 * p_sum / x;
      rate_dist[i] = 1.0 * p_sum / symb_accum;
    }

  if (fp_log)
    {
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Read length distribution\n");
      fprintf(fp_log, "      L           N      Pct   AccPct\n");
      fprintf(fp_log, "-------  ----------  -------  -------\n");
      
      for(long i = len_max; i >= len_min; i--)
        {
          fprintf(fp_log, "%2s%5ld  %10d   %5.1lf%%   %5.1lf%%\n",
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

      long qual_accum = 0;
      for(int c = qmax ; c >= qmin ; c--)
        {
          if (quality_chars[c] > 0)
            {
              qual_accum += quality_chars[c];
              fprintf(fp_log,
                      "    %c  %3ld  %7.5lf  %10lu  %6.1lf%%  %6.1lf%%\n",
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

      for(long i = 2; i <= len_max; i++)
        {
          double PctRecs = 100.0 * (seq_count - length_dist[i-1]) / seq_count;
          double AvgQ = avgq_dist[i-1];
          double AvgP = avgp_dist[i-1];
          double AvgEE = avgee_dist[i-1];
          double Rate = rate_dist[i-1];

          fprintf(fp_log,
                  "%5ld  %6.1lf%%  %4.1lf  %7.5lf  %8.6lf  %5.2lf  %9.6lf  %7.3lf%%\n",
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
      
      for(long i = len_max+1; i >= 1; i--)
        {
          long read_count[4];
          double read_percentage[4];
          
          for(int z=0; z<4; z++)
            {
              if (i>=2)
                read_count[z] = ee_length_table[4*(i-2)+z];
              else
                read_count[z] = seq_count;
              read_percentage[z] = 100.0 * read_count[z] / seq_count;
            }
          
          fprintf(fp_log,
                  "%5ld  %7ld  %7ld  %7ld  %7ld  "
                  "%6.2lf%%  %6.2lf%%  %6.2lf%%  %6.2lf%%\n",
                  i,
                  read_count[0], read_count[1],
                  read_count[2], read_count[3],
                  read_percentage[0], read_percentage[1],
                  read_percentage[2], read_percentage[3]);
        }

      
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Truncate at first Q\n");
      fprintf(fp_log, "  Len     Q=5    Q=10    Q=15    Q=20\n");
      fprintf(fp_log, "-----  ------  ------  ------  ------\n");

      for(long i = len_max; i >= len_max/2; i--)
        {
          double read_percentage[4];
          
          for(int z=0; z<4; z++)
            read_percentage[z] = 100.0 * q_length_table[4*(i-1)+z] / seq_count;

          fprintf(fp_log, "%5ld  %5.1lf%%  %5.1lf%%  %5.1lf%%  %5.1lf%%\n",
                  i,
                  read_percentage[0], read_percentage[1],
                  read_percentage[2], read_percentage[3]);
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "%10lu  Recs (%.1lfM), 0 too long\n",
              seq_count, seq_count / 1.0e6);
      fprintf(fp_log, "%10.1lf  Avg length\n", 1.0 * symbols / seq_count);
      fprintf(fp_log, "%9.1lfM  Bases\n", symbols / 1.0e6);
    }
  
  free(read_length_table);
  free(qual_length_table);
  free(ee_length_table);
  free(q_length_table);

  fastq_close(h);
  
  fprintf(stderr, "Read %lu sequences.\n", seq_count);
}

void fastx_revcomp()
{
  int filetype = fastx_detect(opt_fastx_revcomp);

  unsigned long buffer_alloc = 512;
  char * seq_buffer = (char*) xmalloc(buffer_alloc);
  char * qual_buffer = (char*) xmalloc(buffer_alloc);

  unsigned long header_alloc = 512;
  char * header = (char*) xmalloc(header_alloc);

  unsigned long suffix_length = opt_label_suffix ? strlen(opt_label_suffix) : 0;

  if (filetype == 1)
    {
      /* fasta */
      fasta_handle h = fasta_open(opt_fastx_revcomp);
      unsigned long filesize = fasta_get_size(h);
      
      FILE * fp_fastaout = 0;
      FILE * fp_fastqout = 0;

      if (opt_fastaout)
        {
          fp_fastaout = fopen(opt_fastaout, "w");
          if (!fp_fastaout)
            fatal("Unable to open fasta output file for writing");
        }

      if (opt_fastqout)
        {
          fp_fastqout = fopen(opt_fastqout, "w");
          if (!fp_fastqout)
            fatal("Unable to open fastq output file for writing");

          fprintf(stderr, "WARNING: Writing FASTQ output without base quality information; using max value\n");
        }

      progress_init("Reading fasta file", filesize);
      
      while(fasta_next(h, 0, chrmap_no_change))
        {
          unsigned long length = fasta_get_sequence_length(h);
          unsigned long hlen = fasta_get_header_length(h);
          char * d = fasta_get_header(h);
          char * p = fasta_get_sequence(h);

          if (hlen + suffix_length + 1 > header_alloc)
            header_alloc = hlen + suffix_length + 1;
          header = (char*) xrealloc(header, header_alloc);

          if (length + 1 > buffer_alloc)
            buffer_alloc = length + 1;
          seq_buffer = (char *) xrealloc(seq_buffer, buffer_alloc);
          qual_buffer = (char *) xrealloc(qual_buffer, buffer_alloc);
          
          if (opt_label_suffix)
            snprintf(header, header_alloc, "%s%s", d, opt_label_suffix);
          else
            snprintf(header, header_alloc, "%s", d);

          reverse_complement(seq_buffer, p, length);
          
          /* set quality values to max */
          for(unsigned long i=0; i<length; i++)
            qual_buffer[i] = opt_fastq_ascii + opt_fastq_qmaxout;
          qual_buffer[length] = 0;

          if (opt_fastaout)
            {
              fprint_fasta_hdr_only(fp_fastaout, header);
              fprint_fasta_seq_only(fp_fastaout,
                                    seq_buffer,
                                    length,
                                    opt_fasta_width);
            }
          if (opt_fastqout)
            fprint_fastq(fp_fastqout,
                         header,
                         seq_buffer,
                         qual_buffer,
                         0,
                         0);
                    
                    
          progress_update(fasta_get_position(h));
        }
      progress_done();
      
      if (opt_fastaout)
        fclose(fp_fastaout);

      if (opt_fastqout)
        fclose(fp_fastqout);

      fasta_close(h);
    }
  else if (filetype == 2)
    {
      /* fastq */
      fastq_handle h = fastq_open(opt_fastx_revcomp);
      unsigned long filesize = fastq_get_size(h);
      
      FILE * fp_fastaout = 0;
      FILE * fp_fastqout = 0;

      if (opt_fastaout)
        {
          fp_fastaout = fopen(opt_fastaout, "w");
          if (!fp_fastaout)
            fatal("Unable to open fasta output file for writing");
        }

      if (opt_fastqout)
        {
          fp_fastqout = fopen(opt_fastqout, "w");
          if (!fp_fastqout)
            fatal("Unable to open fastq output file for writing");
        }

      progress_init("Reading fastq file", filesize);
      
      while(fastq_next(h, 0, chrmap_no_change))
        {
          unsigned long length = fastq_get_sequence_length(h);
          unsigned long hlen = fastq_get_header_length(h);
          char * d = fastq_get_header(h);
          char * p = fastq_get_sequence(h);
          char * q = fastq_get_quality(h);

          if (hlen + suffix_length + 1 > header_alloc)
            header_alloc = hlen + suffix_length + 1;
          header = (char*) xrealloc(header, header_alloc);

          if (length + 1 > buffer_alloc)
            buffer_alloc = length + 1;
          seq_buffer = (char *) xrealloc(seq_buffer, buffer_alloc);
          qual_buffer = (char *) xrealloc(qual_buffer, buffer_alloc);
          
          if (opt_label_suffix)
            snprintf(header, header_alloc, "%s%s", d, opt_label_suffix);
          else
            snprintf(header, header_alloc, "%s", d);

          reverse_complement(seq_buffer, p, length);
          
          /* reverse quality values */
          for(unsigned long i=0; i<length; i++)
            qual_buffer[i] = q[length-1-i];
          qual_buffer[length] = 0;

          if (opt_fastaout)
            {
              fprint_fasta_hdr_only(fp_fastaout, header);
              fprint_fasta_seq_only(fp_fastaout,
                                    seq_buffer,
                                    length,
                                    opt_fasta_width);
            }
          if (opt_fastqout)
            fprint_fastq(fp_fastqout,
                         header,
                         seq_buffer,
                         qual_buffer,
                         0,
                         0);
                    
          progress_update(fastq_get_position(h));
        }
      progress_done();
      
      if (opt_fastaout)
        fclose(fp_fastaout);
  
      if (opt_fastqout)
        fclose(fp_fastqout);

      fastq_close(h);
    }
  else
    fatal("Unable to determine file type.");

  free(header);
  free(seq_buffer);
  free(qual_buffer);
}
