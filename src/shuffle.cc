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

void shuffle()
{
  FILE * fp_output = fopen(opt_output, "w");
  if (!fp_output)
    fatal("Unable to open shuffle output file for writing");

  db_read(opt_shuffle, 0);
  show_rusage();

  int dbsequencecount = db_getsequencecount();
  int * deck = (int*) xmalloc(dbsequencecount * sizeof(int));

  for(int i=0; i<dbsequencecount; i++)
    deck[i] = i;

  random_init();

  int passed = 0;
  progress_init("Shuffling", dbsequencecount-1);
  for(int i=dbsequencecount-1; i>0; i--)
    {
      /* generate a random number j in the range 0 to i, inclusive */
      int j = random_int(i+1);

      /* exchange elements i and j */
      int t = deck[i];
      deck[i] = deck[j];
      deck[j] = t;

      passed++;
      progress_update(i);
    }
  progress_done();
  show_rusage();

  passed = MIN(dbsequencecount, opt_topn);

  progress_init("Writing output", passed);
  for(int i=0; i<passed; i++)
    {
      db_fprint_fasta(fp_output, deck[i]);
      progress_update(i);
    }
  progress_done();

  show_rusage();

  free(deck);
  db_free();
  fclose(fp_output);
}
