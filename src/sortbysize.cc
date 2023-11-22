/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <cstdlib>

static struct sortinfo_size_s
{
  unsigned int size;
  unsigned int seqno;
} * sortinfo;

int sortbysize_compare(const void * a, const void * b)
{
  auto * x = (struct sortinfo_size_s *) a;
  auto * y = (struct sortinfo_size_s *) b;

  /* highest abundance first, then by label, otherwise keep order */

  if (x->size < y->size)
    {
      return +1;
    }
  else if (x->size > y->size)
    {
      return -1;
    }
  else
    {
      int r = strcmp(db_getheader(x->seqno), db_getheader(y->seqno));
      if (r != 0)
        {
          return r;
        }
      else
        {
          if (x->seqno < y->seqno)
            {
              return -1;
            }
          else if (x->seqno > y->seqno)
            {
              return +1;
            }
          else
            {
              return 0;
            }
        }
    }
}


[[nodiscard]]
auto find_median_abundance(const int valid_amplicons, const sortinfo_size_s * sortinfo) -> double
{
  // function returns a round value or a value with a remainder of 0.5

  if (valid_amplicons == 0) {
    return 0.0;
  }

  // refactoring C++11: use const& std::vector.size()
  const auto midarray = std::div(valid_amplicons, 2);

  // odd number of valid amplicons
  if (valid_amplicons % 2 != 0)  {
    return sortinfo[midarray.quot].size * 1.0;  // a round value
  }

  // even number of valid amplicons
  // (average of two ints is either round or has a remainder of .5)
  return (sortinfo[midarray.quot - 1].size +
          sortinfo[midarray.quot].size) / 2.0;
}


void sortbysize()
{
  if (!opt_output)
    fatal("FASTA output file for sortbysize must be specified with --output");

  FILE * fp_output = fopen_output(opt_output);
  if (!fp_output)
    {
      fatal("Unable to open sortbysize output file for writing");
    }

  db_read(opt_sortbysize, 0);

  show_rusage();

  int dbsequencecount = db_getsequencecount();

  progress_init("Getting sizes", dbsequencecount);

  // refactoring C++11: use std::vector
  sortinfo = (struct sortinfo_size_s*)
    xmalloc(dbsequencecount * sizeof(sortinfo_size_s));

  int passed = 0;

  for(int i=0; i<dbsequencecount; i++)
    {
      const int64_t size = db_getabundance(i);

      if((size >= opt_minsize) && (size <= opt_maxsize))
        {
          sortinfo[passed].seqno = i;
          sortinfo[passed].size = (unsigned int) size;
          passed++;
        }
      progress_update(i);
    }

  progress_done();

  show_rusage();

  progress_init("Sorting", 100);
  qsort(sortinfo, passed, sizeof(sortinfo_size_s), sortbysize_compare);
  progress_done();

  const double median = find_median_abundance(passed, sortinfo);

  if (! opt_quiet)
    {
      fprintf(stderr, "Median abundance: %.0f\n", median);  // drop remainder
    }

  if (opt_log)
    {
      fprintf(fp_log, "Median abundance: %.0f\n", median);  // drop remainder
    }

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
