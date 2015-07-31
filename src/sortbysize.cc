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

static struct sortinfo_s
{
  unsigned int size;
  unsigned int seqno;
} * sortinfo;

int sortbysize_compare(const void * a, const void * b)
{
  struct sortinfo_s * x = (struct sortinfo_s *) a;
  struct sortinfo_s * y = (struct sortinfo_s *) b;

  /* highest abundance first, then by label, otherwise keep order */

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

void sortbysize()
{
  FILE * fp_output = fopen(opt_output, "w");
  if (!fp_output)
    fatal("Unable to open sortbysize output file for writing");

  db_read(opt_sortbysize, 0);

  show_rusage();

  int dbsequencecount = db_getsequencecount();
  
  progress_init("Getting sizes", dbsequencecount);

  sortinfo = (sortinfo_s*) xmalloc(dbsequencecount * sizeof(sortinfo_s));

  int passed = 0;

  for(int i=0; i<dbsequencecount; i++)
    {
      long size = db_getabundance(i);
      
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
  qsort(sortinfo, passed, sizeof(sortinfo_s), sortbysize_compare);
  progress_done();
  
  double median = 0.0;
  if (passed > 0)
    {
      if (passed % 2)
        median = sortinfo[(passed-1)/2].size;
      else
        median = (sortinfo[(passed/2)-1].size +
                  sortinfo[passed/2].size) / 2.0;
    }

  if (! opt_quiet)
    fprintf(stderr, "Median abundance: %.0f\n", median);

  if (opt_log)
    fprintf(fp_log, "Median abundance: %.0f\n", median);

  show_rusage();
  
  passed = MIN(passed, opt_topn);

  progress_init("Writing output", passed);
  for(int i=0; i<passed; i++)
    {
      unsigned int seqno = sortinfo[i].seqno;
      unsigned int size = sortinfo[i].size;

      if (opt_relabel_sha1 || opt_relabel_md5)
        {
          char * seq = db_getsequence(seqno);
          unsigned int len = db_getsequencelen(seqno);

          fprintf(fp_output, ">");

          if (opt_relabel_sha1)
            fprint_seq_digest_sha1(fp_output, seq, len);
          else
            fprint_seq_digest_md5(fp_output, seq, len);

          if (opt_sizeout)
            fprintf(fp_output, ";size=%u;\n", size);
          else
            fprintf(fp_output, "\n");

          db_fprint_fasta_seq_only(fp_output, seqno);
        }
      else if (opt_relabel)
        {
          if (opt_sizeout)
            fprintf(fp_output, ">%s%d;size=%u;\n", opt_relabel, i+1, size);
          else
            fprintf(fp_output, ">%s%d\n", opt_relabel, i+1);

          db_fprint_fasta_seq_only(fp_output, seqno);
        }
      else
        {
          db_fprint_fasta(fp_output, seqno);
        }
      
      progress_update(i);
    }
  progress_done();
  show_rusage();
  
  free(sortinfo);
  db_free();
  fclose(fp_output);
}
