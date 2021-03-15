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

/* static variables */

FILE * join_fileopenw(char * filename)
{
  FILE * fp = 0;
  fp = fopen_output(filename);
  if (!fp)
    fatal("Unable to open file for writing (%s)", filename);
  return fp;
}

void fastq_join()
{
  FILE * fp_fastqout = 0;
  FILE * fp_fastaout = 0;

  fastx_handle fastq_fwd = 0;
  fastx_handle fastq_rev = 0;

  uint64_t total = 0;

  /* check input and options */

  if (!opt_reverse)
    fatal("No reverse reads file specified with --reverse");

  if ((!opt_fastqout) && (!opt_fastaout))
    fatal("No output files specified");

  char * padgap = 0;
  char * padgapq = 0;

  if (opt_join_padgap)
    padgap = xstrdup(opt_join_padgap);
  else
    padgap = xstrdup("NNNNNNNN");

  uint64_t padlen = strlen(padgap);

  if (opt_join_padgapq)
    padgapq = xstrdup(opt_join_padgapq);
  else
    {
      padgapq = (char *) xmalloc(padlen + 1);
      for(uint64_t i = 0; i < padlen; i++)
        padgapq[i] = 'I';
      padgapq[padlen] = 0;
    }

  if (padlen != strlen(padgapq))
    fatal("Strings given by --join_padgap and --join_padgapq differ in length");

  /* open input files */

  fastq_fwd = fastq_open(opt_fastq_join);
  fastq_rev = fastq_open(opt_reverse);

  /* open output files */

  if (opt_fastqout)
    fp_fastqout = join_fileopenw(opt_fastqout);
  if (opt_fastaout)
    fp_fastaout = join_fileopenw(opt_fastaout);

  /* main */

  uint64_t filesize = fastq_get_size(fastq_fwd);
  progress_init("Joining reads", filesize);

  /* do it */

  total = 0;

  uint64_t alloc = 0;
  uint64_t len = 0;
  char * seq = 0;
  char * qual = 0;

  while(fastq_next(fastq_fwd, 0, chrmap_no_change))
    {
      if (! fastq_next(fastq_rev, 0, chrmap_no_change))
        fatal("More forward reads than reverse reads");

      uint64_t fwd_seq_length = fastq_get_sequence_length(fastq_fwd);
      uint64_t rev_seq_length = fastq_get_sequence_length(fastq_rev);

      /* allocate enough mem */

      uint64_t needed = fwd_seq_length + rev_seq_length + padlen + 1;
      if (alloc < needed)
        {
          seq = (char *) xrealloc(seq, needed);
          qual = (char *) xrealloc(qual, needed);
          alloc = needed;
        }

      /* join them */

      strcpy(seq, fastq_get_sequence(fastq_fwd));
      strcpy(qual, fastq_get_quality(fastq_fwd));
      len = fwd_seq_length;

      strcpy(seq + len, padgap);
      strcpy(qual + len, padgapq);
      len += padlen;

      /* reverse complement reverse read */

      char * rev_seq = fastq_get_sequence(fastq_rev);
      char * rev_qual = fastq_get_quality(fastq_rev);

      for(uint64_t i = 0; i < rev_seq_length; i++)
        {
          uint64_t rev_pos = rev_seq_length - 1 - i;
          seq[len]  = chrmap_complement[(int)(rev_seq[rev_pos])];
          qual[len] = rev_qual[rev_pos];
          len++;
        }
      seq[len] = 0;
      qual[len] = 0;

      /* write output */

      if (opt_fastqout)
        {
          fastq_print_general(fp_fastqout,
                              seq,
                              len,
                              fastq_get_header(fastq_fwd),
                              fastq_get_header_length(fastq_fwd),
                              qual,
                              0,
                              total + 1,
                              -1.0);
        }

      if (opt_fastaout)
        {
          fasta_print_general(fp_fastaout,
                              0,
                              seq,
                              len,
                              fastq_get_header(fastq_fwd),
                              fastq_get_header_length(fastq_fwd),
                              0,
                              total + 1,
                              -1.0,
                              -1,
                              -1,
                              0,
                              0);
        }

      total++;
      progress_update(fastq_get_position(fastq_fwd));
    }

  progress_done();

  if (fastq_next(fastq_rev, 0, chrmap_no_change))
    fatal("More reverse reads than forward reads");

  fprintf(stderr,
          "%" PRIu64 " pairs joined\n",
          total);

  /* clean up */

  if (opt_fastaout)
    fclose(fp_fastaout);
  if (opt_fastqout)
    fclose(fp_fastqout);

  fastq_close(fastq_rev);
  fastq_rev = 0;
  fastq_close(fastq_fwd);
  fastq_fwd = 0;

  if (seq)
    xfree(seq);
  if (qual)
    xfree(qual);
  xfree(padgap);
  xfree(padgapq);
}
