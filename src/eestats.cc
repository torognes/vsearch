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

#include "vsearch5d.h"

double q2p(int q)
{
  return exp10(- q / 10.0);
}

long ee_start(int pos, int resolution)
{
  return pos * (resolution * (pos + 1) + 2) / 2;
}

void fastq_eestats()
{
  fastq_handle h = fastq_open(opt_fastq_eestats);

  unsigned long filesize = fastq_get_size(h);

  FILE * fp_output = 0;

  if (opt_output)
    {
      fp_output = fopen(opt_output, "w");
      if (!fp_output)
        fatal("Unable to open output file for writing");
    }

  progress_init("Reading fastq file", filesize);

  unsigned long seq_count = 0;
  unsigned long symbols = 0;
  
  long len_alloc = 10;

  const int resolution = 1000;
  int max_quality = opt_fastq_qmax - opt_fastq_qmin + 1;
  
  long ee_size = ee_start(len_alloc, resolution);
  
  int * read_length_table = (int*) xmalloc(sizeof(int) * len_alloc);
  memset(read_length_table, 0, sizeof(int) * len_alloc);
  
  int * qual_length_table = (int*) xmalloc(sizeof(int) * len_alloc *
                                           (max_quality+1));
  memset(qual_length_table, 0, sizeof(int) * len_alloc * (max_quality+1));
  
  int * ee_length_table = (int*) xmalloc(sizeof(int) * ee_size);
  memset(ee_length_table, 0, sizeof(int) * ee_size);

  double * sum_ee_length_table = (double*) xmalloc(sizeof(double) * len_alloc);
  memset(sum_ee_length_table, 0, sizeof(double) * len_alloc);
  
  double * sum_pe_length_table = (double*) xmalloc(sizeof(double) * len_alloc);
  memset(sum_pe_length_table, 0, sizeof(double) * len_alloc);
  
  long len_min = LONG_MAX;
  long len_max = 0;

  while(fastq_next(h, 0, chrmap_upcase))
    {
      seq_count++;

      long len = fastq_get_sequence_length(h);
      char * q = fastq_get_quality(h);

      /* update length statistics */

      long new_alloc = len + 1;

      if (new_alloc > len_alloc)
        {
          long new_ee_size = ee_start(new_alloc, resolution);

          read_length_table = (int*) xrealloc(read_length_table,
                                              sizeof(int) * new_alloc);
          memset(read_length_table + len_alloc, 0, 
                 sizeof(int) * (new_alloc - len_alloc));
          
          qual_length_table = (int*) xrealloc(qual_length_table, sizeof(int) *
                                              new_alloc * (max_quality+1));
          memset(qual_length_table + (max_quality+1) * len_alloc, 0, 
                 sizeof(int) * (new_alloc - len_alloc) * (max_quality+1));
          
          ee_length_table = (int*) xrealloc(ee_length_table, sizeof(int) * 
                                            new_ee_size);
          memset(ee_length_table + ee_size, 0,
                 sizeof(int) * (new_ee_size - ee_size));
          
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

      for(long i=0; i < len; i++)
        {
          read_length_table[i]++;

          /* quality score */

          int qual = q[i] - opt_fastq_ascii;
          
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

          if (qual < 0)
            qual = 0;

          qual_length_table[(max_quality+1)*i + qual]++;


          /* Pe */

          double pe = q2p(qual);
          sum_pe_length_table[i] += pe;


          /* expected number of errors */

          ee += pe;

          int e_int = MIN(resolution*(i+1), (int)(resolution * ee));
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

  for(long i=0; i<len_max; i++)
    {
      long reads = read_length_table[i];
      double pctrecs = 100.0 * reads / seq_count;
      
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
          qsum += q * x;
          n += x;
          
          if ((min_q<0) && (x > 0))
            min_q = q;

          if ((low_q<0) && (n >= 0.25 * reads))
            low_q = q;

          if ((med_q<0) && (n >= 0.50 * reads))
            med_q = q;

          if ((hi_q<0)  && (n >= 0.75 * reads))
            hi_q = q;

          if (x > 0)
            max_q = q;
        }

      double mean_q = 1.0 * qsum / reads;


      /* pe */

      double min_pe = q2p(max_q);
      double low_pe = q2p(hi_q);
      double med_pe = q2p(med_q);
      double hi_pe = q2p(low_q);
      double max_pe = q2p(min_q);
      double mean_pe = sum_pe_length_table[i] / reads;


      /* expected errors */

      double min_ee = -1.0;
      double low_ee = -1.0;
      double med_ee = -1.0;
      double hi_ee  = -1.0;
      double max_ee = -1.0;

      long ee_offset = ee_start(i, resolution);
      long max_errors = resolution * (i+1);

      n = 0;
      for(int e=0; e<=max_errors; e++)
        {
          int x = ee_length_table[ee_offset + e];
          n += x;

          if ((min_ee<0) && (x > 0))
            min_ee = e;

          if ((low_ee<0) && (n >= 0.25 * reads))
            low_ee = e;
          
          if ((med_ee<0) && (n >= 0.50 * reads))
            med_ee = e;
          
          if ((hi_ee<0)  && (n >= 0.75 * reads))
            hi_ee = e;

          if (x > 0)
            max_ee = e;
        }

      double mean_ee = sum_ee_length_table[i] / reads;

      min_ee  = (min_ee  + 0.5) / resolution;
      low_ee  = (low_ee  + 0.5) / resolution;
      med_ee  = (med_ee  + 0.5) / resolution;
      hi_ee   = (hi_ee   + 0.5) / resolution;
      max_ee  = (max_ee  + 0.5) / resolution;

      fprintf(fp_output,
              "%ld\t%ld\t%.1lf"
              "\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf"
              "\t%.2lg\t%.2lg\t%.2lg\t%.2lg\t%.2lg\t%.2lg"
              "\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
              i+1, reads, pctrecs,
              min_q,  low_q,  med_q,  mean_q,  hi_q,  max_q,
              min_pe, low_pe, med_pe, mean_pe, hi_pe, max_pe,
              min_ee, low_ee, med_ee, mean_ee, hi_ee, max_ee);
    }
  
  free(read_length_table);
  free(qual_length_table);
  free(ee_length_table);
  free(sum_ee_length_table);
  free(sum_pe_length_table);

  fclose(fp_output);

  fastq_close(h);
}
