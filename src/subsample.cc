/*

  VSEARCH5D: a modified version of VSEARCH

  Copyright (C) 2016, Akifumi S. Tanabe

  Contact: Akifumi S. Tanabe <akifumi.tanabe@gmail.com>

  Original version of VSEARCH
  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri

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

void subsample()
{
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

  db_read(opt_fastx_subsample, 0);
  show_rusage();

  if (fp_fastqout && ! db_is_fastq())
    fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");

  int dbsequencecount = db_getsequencecount();

  unsigned long mass_total = 0;

  if (!opt_sizein)
    mass_total = dbsequencecount;
  else
    for(int i=0; i<dbsequencecount; i++)
      mass_total += db_getabundance(i);
  
  fprintf(stderr, "Got %lu reads from %d amplicons\n",
          mass_total, dbsequencecount);


  int * abundance = (int*) xmalloc(dbsequencecount * sizeof(int));

  for(int i=0; i<dbsequencecount; i++)
    abundance[i] = 0;

  random_init();

  unsigned long n;                              /* number of reads to sample */
  if (opt_sample_size)
    n = opt_sample_size;
  else
    n = mass_total * opt_sample_pct / 100.0;

  if (n > mass_total)
    fatal("Cannot subsample more reads than in the original sample");

  unsigned long x = n;                          /* number of reads left */
  int a = 0;                                    /* amplicon number */
  unsigned long r = 0;                          /* read being checked */
  unsigned long m = 0;                          /* accumulated mass */

  unsigned long mass =                          /* mass of current amplicon */
    opt_sizein ? db_getabundance(0) : 1;
  
  progress_init("Subsampling", mass_total);
  while (x > 0)
    {
      unsigned long random = random_ulong(mass_total - r);

      if (random < x)
        {
          /* selected read r from amplicon a */
          abundance[a]++;
          x--;
        }

      r++;
      m++;
      if (m >= mass)
        {
          /* next amplicon */
          a++;
          mass = opt_sizein ? db_getabundance(a) : 1;
          m = 0;
        }
      progress_update(r);
    }
  progress_done();

  int samples = 0;
  progress_init("Writing output", dbsequencecount);
  for(int i=0; i<dbsequencecount; i++)
    {
      if (abundance[i]>0)
        {
          samples++;

          if (opt_fastaout)
            fasta_print_relabel(fp_fastaout,
                                db_getsequence(i),
                                db_getsequencelen(i),
                                db_getheader(i),
                                db_getheaderlen(i),
                                abundance[i],
                                samples);

          if (opt_fastqout)
            fastq_print_relabel(fp_fastqout,
                                db_getsequence(i),
                                db_getsequencelen(i),
                                db_getheader(i),
                                db_getheaderlen(i),
                                db_getquality(i),
                                abundance[i],
                                samples);
        }
      progress_update(i);
    }
  progress_done();

  free(abundance);

  fprintf(stderr, "Subsampled %lu reads from %d amplicons\n", n, samples);

  db_free();

  if (opt_fastaout)
    fclose(fp_fastaout);

  if (opt_fastqout)
    fclose(fp_fastqout);
}
