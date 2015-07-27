/*
  Copyright (C) 2014-2015 Torbjorn Rognes

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "vsearch.h"

void subsample()
{
  FILE * fp_output = fopen(opt_fastaout, "w");
  if (!fp_output)
    fatal("Unable to open subsampling output file for writing");

  db_read(opt_fastx_subsample, 0);
  show_rusage();

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

  int sampled = 0;
  progress_init("Writing output", dbsequencecount);
  for(int i=0; i<dbsequencecount; i++)
    {
      if (abundance[i]>0)
        {
          if (opt_sizeout)
            db_fprint_fasta_with_size(fp_output, i, abundance[i]);
          else if (opt_xsize)
            db_fprint_fasta_strip_size(fp_output, i);
          else
            db_fprint_fasta(fp_output, i);

          sampled++;
        }
      progress_update(i);
    }
  progress_done();

  free(abundance);

  fprintf(stderr, "Subsampled %lu reads from %d amplicons\n", n, sampled);

  db_free();
  fclose(fp_output);
}
