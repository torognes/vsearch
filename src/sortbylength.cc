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

static struct sortinfo_s
{
  unsigned int length;
  unsigned int size;
  unsigned int seqno;
} * sortinfo;

int sortbylength_compare(const void * a, const void * b)
{
  struct sortinfo_s * x = (struct sortinfo_s *) a;
  struct sortinfo_s * y = (struct sortinfo_s *) b;

  /* longest first, then most abundant, then by label, otherwise keep order */

  if (x->length < y->length)
    return +1;
  else if (x->length > y->length)
    return -1;
  else
    if (x->size < y->size)
      return +1;
    else if (x->size > y->size)
      return -1;
    else
      {
        int r = strcmp(db_getheader(x->seqno), db_getheader(y->seqno));
        if (r != 0)
          return r;
        else
          {
            if (x->seqno < y->seqno)
              return -1;
            else if (x->seqno > y->seqno)
              return +1;
            else
              return 0;
          }
      }
}

void sortbylength()
{
  FILE * fp_output = fopen_output(opt_output);
  if (!fp_output)
    fatal("Unable to open sortbylength output file for writing");

  db_read(opt_sortbylength, 0);
  show_rusage();

  int dbsequencecount = db_getsequencecount();
  sortinfo = (struct sortinfo_s *) xmalloc(dbsequencecount * sizeof(sortinfo_s));

  int passed = 0;

  progress_init("Getting lengths", dbsequencecount);
  for(int i=0; i<dbsequencecount; i++)
    {
      sortinfo[passed].seqno = i;
      sortinfo[passed].length = db_getsequencelen(i);
      sortinfo[passed].size = db_getabundance(i);
      passed++;
      progress_update(i);
    }
  progress_done();
  show_rusage();

  progress_init("Sorting", 100);
  qsort(sortinfo, passed, sizeof(sortinfo_s), sortbylength_compare);
  progress_done();

  double median = 0.0;
  if (passed > 0)
    {
      if (passed % 2)
        median = sortinfo[(passed-1)/2].length;
      else
        median = (sortinfo[(passed/2)-1].length +
                  sortinfo[passed/2].length) / 2.0;
    }

  if (!opt_quiet)
    fprintf(stderr, "Median length: %.0f\n", median);

  if (opt_log)
    fprintf(fp_log, "Median length: %.0f\n", median);

  show_rusage();
  
  passed = MIN(passed, opt_topn);

  progress_init("Writing output", passed);
  for(int i=0; i<passed; i++)
    {
      fasta_print_db_relabel(fp_output, sortinfo[i].seqno, i+1);
      progress_update(i);
    }
  progress_done();
  show_rusage();

  xfree(sortinfo);
  db_free();
  fclose(fp_output);
}
