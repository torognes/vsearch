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

double q2p(int q)
{
  return exp10(- q / 10.0);
}

int64_t ee_start(int pos, int resolution)
{
  return pos * (resolution * (pos + 1) + 2) / 2;
}

void fastq_eestats()
{
  fastx_handle h = fastq_open(opt_fastq_eestats);

  uint64_t filesize = fastq_get_size(h);

  FILE * fp_output = 0;

  if (opt_output)
    {
      fp_output = fopen_output(opt_output);
      if (!fp_output)
        fatal("Unable to open output file for writing");
    }

  progress_init("Reading FASTQ file", filesize);

  uint64_t seq_count = 0;
  uint64_t symbols = 0;

  int64_t len_alloc = 10;

  const int resolution = 1000;
  int max_quality = opt_fastq_qmax - opt_fastq_qmin + 1;

  int64_t ee_size = ee_start(len_alloc, resolution);

  uint64_t * read_length_table = (uint64_t*) xmalloc(sizeof(uint64_t) * len_alloc);
  memset(read_length_table, 0, sizeof(uint64_t) * len_alloc);

  uint64_t * qual_length_table = (uint64_t*) xmalloc(sizeof(uint64_t) * len_alloc *
                                           (max_quality+1));
  memset(qual_length_table, 0, sizeof(uint64_t) * len_alloc * (max_quality+1));

  uint64_t * ee_length_table = (uint64_t*) xmalloc(sizeof(uint64_t) * ee_size);
  memset(ee_length_table, 0, sizeof(uint64_t) * ee_size);

  double * sum_ee_length_table = (double*) xmalloc(sizeof(double) * len_alloc);
  memset(sum_ee_length_table, 0, sizeof(double) * len_alloc);

  double * sum_pe_length_table = (double*) xmalloc(sizeof(double) * len_alloc);
  memset(sum_pe_length_table, 0, sizeof(double) * len_alloc);

  int64_t len_min = LONG_MAX;
  int64_t len_max = 0;

  while(fastq_next(h, 0, chrmap_upcase))
    {
      seq_count++;

      int64_t len = fastq_get_sequence_length(h);
      char * q = fastq_get_quality(h);

      /* update length statistics */

      int64_t new_alloc = len + 1;

      if (new_alloc > len_alloc)
        {
          int64_t new_ee_size = ee_start(new_alloc, resolution);

          read_length_table = (uint64_t*) xrealloc(read_length_table,
                                              sizeof(uint64_t) * new_alloc);
          memset(read_length_table + len_alloc, 0,
                 sizeof(uint64_t) * (new_alloc - len_alloc));

          qual_length_table = (uint64_t*) xrealloc(qual_length_table, sizeof(uint64_t) *
                                              new_alloc * (max_quality+1));
          memset(qual_length_table + (max_quality+1) * len_alloc, 0,
                 sizeof(uint64_t) * (new_alloc - len_alloc) * (max_quality+1));

          ee_length_table = (uint64_t*) xrealloc(ee_length_table, sizeof(uint64_t) *
                                            new_ee_size);
          memset(ee_length_table + ee_size, 0,
                 sizeof(uint64_t) * (new_ee_size - ee_size));

          sum_ee_length_table = (double*) xrealloc(sum_ee_length_table,
                                              sizeof(double) * new_alloc);
          memset(sum_ee_length_table + len_alloc, 0,
                 sizeof(double) * (new_alloc - len_alloc));

          sum_pe_length_table = (double*) xrealloc(sum_pe_length_table,
                                              sizeof(double) * new_alloc);
          memset(sum_pe_length_table + len_alloc, 0,
                 sizeof(double) * (new_alloc - len_alloc));

          len_alloc = new_alloc;
          ee_size = new_ee_size;
        }

      if (len < len_min)
        len_min = len;
      if (len > len_max)
        len_max = len;

      /* update quality statistics */

      symbols += len;

      double ee = 0.0;

      for(int64_t i=0; i < len; i++)
        {
          read_length_table[i]++;

          /* quality score */

          int qual = q[i] - opt_fastq_ascii;

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

          if (qual < 0)
            qual = 0;

          qual_length_table[(max_quality+1)*i + qual]++;


          /* Pe */

          double pe = q2p(qual);
          sum_pe_length_table[i] += pe;


          /* expected number of errors */

          ee += pe;

          int64_t e_int = MIN(resolution*(i+1), (int)(resolution * ee));
          ee_length_table[ee_start(i, resolution) + e_int]++;

          sum_ee_length_table[i] += ee;
        }
      progress_update(fastq_get_position(h));
    }
  progress_done();

  fprintf(fp_output,
          "Pos\tRecs\tPctRecs\t"
          "Min_Q\tLow_Q\tMed_Q\tMean_Q\tHi_Q\tMax_Q\t"
          "Min_Pe\tLow_Pe\tMed_Pe\tMean_Pe\tHi_Pe\tMax_Pe\t"
          "Min_EE\tLow_EE\tMed_EE\tMean_EE\tHi_EE\tMax_EE\n");

  for(int64_t i=0; i<len_max; i++)
    {
      int64_t reads = read_length_table[i];
      double pctrecs = 100.0 * reads / seq_count;


      /* q */

      double min_q = -1.0;
      double low_q = -1.0;
      double med_q = -1.0;
      double hi_q  = -1.0;
      double max_q = -1.0;

      double qsum = 0;
      double n = 0;
      for(int q=0; q<=max_quality; q++)
        {
          double x = qual_length_table[(max_quality+1)*i+q];

          if (x > 0)
            {
              qsum += q * x;
              n += x;

              if (min_q<0)
                min_q = q;

              if ((low_q<0) && (n >= 0.25 * reads))
                low_q = q;

              if ((med_q<0) && (n >= 0.50 * reads))
                med_q = q;

              if ((hi_q<0)  && (n >= 0.75 * reads))
                hi_q = q;

              max_q = q;
            }
        }

      double mean_q = 1.0 * qsum / reads;


      /* pe */

      double min_pe = -1.0;
      double low_pe = -1.0;
      double med_pe = -1.0;
      double hi_pe  = -1.0;
      double max_pe = -1.0;

      double pesum = 0;
      n = 0;
      for(int q=max_quality; q>=0; q--)
        {
          double x = qual_length_table[(max_quality+1)*i+q];

          if (x > 0)
            {
              double pe = q2p(q);
              pesum += pe * x;
              n += x;

              if (min_pe<0)
                min_pe = pe;

              if ((low_pe<0) && (n >= 0.25 * reads))
                low_pe = pe;

              if ((med_pe<0) && (n >= 0.50 * reads))
                med_pe = pe;

              if ((hi_pe<0)  && (n >= 0.75 * reads))
                hi_pe = pe;

              max_pe = pe;
            }
        }

      double mean_pe = 1.0 * pesum / reads;


      /* expected errors */

      double min_ee = -1.0;
      double low_ee = -1.0;
      double med_ee = -1.0;
      double hi_ee  = -1.0;
      double max_ee = -1.0;

      int64_t ee_offset = ee_start(i, resolution);
      int64_t max_errors = resolution * (i+1);

      n = 0;
      for(int64_t e=0; e<=max_errors; e++)
        {
          int64_t x = ee_length_table[ee_offset + e];

          if (x > 0)
            {
              n += x;

              if (min_ee<0)
                min_ee = e;

              if ((low_ee<0) && (n >= 0.25 * reads))
                low_ee = e;

              if ((med_ee<0) && (n >= 0.50 * reads))
                med_ee = e;

              if ((hi_ee<0)  && (n >= 0.75 * reads))
                hi_ee = e;

              max_ee = e;
            }
        }

      double mean_ee = sum_ee_length_table[i] / reads;

      min_ee  = (min_ee  + 0.5) / resolution;
      low_ee  = (low_ee  + 0.5) / resolution;
      med_ee  = (med_ee  + 0.5) / resolution;
      hi_ee   = (hi_ee   + 0.5) / resolution;
      max_ee  = (max_ee  + 0.5) / resolution;

      fprintf(fp_output,
              "%" PRId64 "\t%" PRId64 "\t%.1lf"
              "\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf"
              "\t%.2lg\t%.2lg\t%.2lg\t%.2lg\t%.2lg\t%.2lg"
              "\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
              i+1, reads, pctrecs,
              min_q,  low_q,  med_q,  mean_q,  hi_q,  max_q,
              min_pe, low_pe, med_pe, mean_pe, hi_pe, max_pe,
              min_ee, low_ee, med_ee, mean_ee, hi_ee, max_ee);
    }

  xfree(read_length_table);
  xfree(qual_length_table);
  xfree(ee_length_table);
  xfree(sum_ee_length_table);
  xfree(sum_pe_length_table);

  fclose(fp_output);

  fastq_close(h);
}

void fastq_eestats2()
{
  fastx_handle h = fastq_open(opt_fastq_eestats2);

  uint64_t filesize = fastq_get_size(h);

  FILE * fp_output = 0;

  if (opt_output)
    {
      fp_output = fopen_output(opt_output);
      if (!fp_output)
        fatal("Unable to open output file for writing");
    }

  progress_init("Reading FASTQ file", filesize);

  uint64_t seq_count = 0;
  uint64_t symbols = 0;
  uint64_t longest = 0;

  int len_steps = 0;

  uint64_t * count_table = 0;

  while(fastq_next(h, 0, chrmap_upcase))
    {
      seq_count++;

      uint64_t len = fastq_get_sequence_length(h);
      char * q = fastq_get_quality(h);

      /* update length statistics */

      if (len > longest)
        {
          longest = len;
          int new_len_steps = 1 + MAX(0, (MIN(longest, (uint64_t)opt_length_cutoffs_longest) - opt_length_cutoffs_shortest) / opt_length_cutoffs_increment);

          if (new_len_steps > len_steps)
            {
              count_table = (uint64_t *) xrealloc(count_table, sizeof(uint64_t) * new_len_steps * opt_ee_cutoffs_count);
              memset(count_table + len_steps * opt_ee_cutoffs_count, 0, sizeof(uint64_t) * (new_len_steps - len_steps) * opt_ee_cutoffs_count);
              len_steps = new_len_steps;
            }
        }

      /* update quality statistics */

      symbols += len;

      double ee = 0.0;

      for(uint64_t i=0; i < len; i++)
        {
          /* quality score */

          int qual = q[i] - opt_fastq_ascii;

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

          if (qual < 0)
            qual = 0;

          double pe = q2p(qual);

          ee += pe;

          for (int x = 0; x < len_steps; x++)
            {
              uint64_t len_cutoff = opt_length_cutoffs_shortest + x * opt_length_cutoffs_increment;
              if (i+1 == len_cutoff)
                for (int y = 0; y < opt_ee_cutoffs_count; y++)
                  if (ee <= opt_ee_cutoffs_values[y])
                    count_table[x * opt_ee_cutoffs_count + y]++;
            }
        }

      for (int x = 0; x < len_steps; x++)
        {
          uint64_t len_cutoff = opt_length_cutoffs_shortest + x * opt_length_cutoffs_increment;
          if (len < len_cutoff)
            for (int y = 0; y < opt_ee_cutoffs_count; y++)
              if (ee <= opt_ee_cutoffs_values[y])
                count_table[x * opt_ee_cutoffs_count + y]++;
        }

      progress_update(fastq_get_position(h));
    }
  progress_done();

  fprintf(fp_output,
          "%" PRIu64 " reads, max len %" PRIu64 ", avg %.1f\n\n",
          seq_count, longest, 1.0 * symbols / seq_count);

  fprintf(fp_output, "Length");
  for (int y = 0; y < opt_ee_cutoffs_count; y++)
    fprintf(fp_output, "         MaxEE %.2f", opt_ee_cutoffs_values[y]);
  fprintf(fp_output, "\n");
  fprintf(fp_output, "------");
  for (int y = 0; y < opt_ee_cutoffs_count; y++)
    fprintf(fp_output, "   ----------------");
  fprintf(fp_output, "\n");

  for (int x = 0; x < len_steps; x++)
    {
      int len_cutoff = opt_length_cutoffs_shortest + x * opt_length_cutoffs_increment;

      if (len_cutoff > opt_length_cutoffs_longest)
        break;

      fprintf(fp_output, "%6d", len_cutoff);

      for (int y = 0; y < opt_ee_cutoffs_count; y++)
        fprintf(fp_output,
                "   %8" PRIu64 "(%5.1f%%)",
                count_table[x * opt_ee_cutoffs_count + y],
                100.0 * count_table[x * opt_ee_cutoffs_count + y] / seq_count);
      fprintf(fp_output, "\n");
    }

  if (fp_log)
    {
      fprintf(fp_log,
              "%" PRIu64 " reads, max len %" PRIu64 ", avg %.1f\n\n",
              seq_count, longest, 1.0 * symbols / seq_count);

      fprintf(fp_log, "Length");
      for (int y = 0; y < opt_ee_cutoffs_count; y++)
        fprintf(fp_log, "         MaxEE %.2f", opt_ee_cutoffs_values[y]);
      fprintf(fp_log, "\n");
      fprintf(fp_log, "------");
      for (int y = 0; y < opt_ee_cutoffs_count; y++)
        fprintf(fp_log, "   ----------------");
      fprintf(fp_log, "\n");

      for (int x = 0; x < len_steps; x++)
        {
          int len_cutoff = opt_length_cutoffs_shortest + x * opt_length_cutoffs_increment;

          if (len_cutoff > opt_length_cutoffs_longest)
            break;

          fprintf(fp_log, "%6d", len_cutoff);

          for (int y = 0; y < opt_ee_cutoffs_count; y++)
            fprintf(fp_log,
                    "   %8" PRIu64 "(%5.1f%%)",
                    count_table[x * opt_ee_cutoffs_count + y],
                    100.0 * count_table[x * opt_ee_cutoffs_count + y] / seq_count);
          fprintf(fp_log, "\n");
        }
    }

  if (count_table)
    {
      xfree(count_table);
      count_table = 0;
    }

  fclose(fp_output);

  fastq_close(h);
}
