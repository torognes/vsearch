/*
  Copyright (C) 2014 Torbjorn Rognes

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

int * deck;

long random_int(long n)
{
  /*
    Generate a random integer in the range 0 to n-1, inclusive.
    The random() function returns a random number in the range
    0 to 2147483647 (=2^31-1=RAND_MAX), inclusive.
    We should avoid some of the upper generated numbers to
    avoid modulo bias.
  */

  long random_max = RAND_MAX;
  long limit = random_max - (random_max + 1) % n;
  long r = random();
  while (r > limit)
    r = random();
  return r % n;
}

void shuffle()
{
  FILE * fp_output = fopen(opt_output, "w");
  if (!fp_output)
    fatal("Unable to open shuffle output file for writing");

  db_read(opt_shuffle, 0);
  show_rusage();

  int dbsequencecount = db_getsequencecount();
  deck = (int*) xmalloc(dbsequencecount * sizeof(int));

  /* initialize pseudo-random number generator */
  unsigned int seed = opt_seed;
  if (seed == 0)
    {
      int fd = open("/dev/urandom", O_RDONLY);
      if (fd < 0)
        fatal("Unable to open /dev/urandom");
      if (read(fd, & seed, sizeof(seed)) < 0)
        fatal("Unable to read from /dev/urandom");
      close(fd);
    }
  srandom(seed);

  for(int i=0; i<dbsequencecount; i++)
    deck[i] = i;

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
#if 0
      if (opt_relabel)
        {
          if (opt_sizeout)
            fprintf(fp_output, ">%s%d;size=%lu;\n", opt_relabel, i+1, db_getabundance(i));
          else
            fprintf(fp_output, ">%s%d\n", opt_relabel, i+1);

          char * seq = db_getsequence(deck[i]);
          long len = db_getsequencelen(deck[i]);
          fprint_fasta_seq_only(fp_output, seq, len, opt_fasta_width);
        }
      else
        {
#endif
          db_fprint_fasta(fp_output, deck[i]);
#if 0
        }
#endif

      progress_update(i);
    }
  progress_done();
  show_rusage();

  db_free();
  fclose(fp_output);
}
