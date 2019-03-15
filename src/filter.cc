/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2019, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

inline int fastq_get_qual(char q)
{
  int qual = q - opt_fastq_ascii;
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
  return qual;
}

struct analysis_res
{
  bool discarded;
  bool truncated;
  int start;
  int length;
  double ee;
};

struct analysis_res analyse(fastx_handle h)
{
  struct analysis_res res = { false, false, 0, 0, -1.0 };
  res.length = fastx_get_sequence_length(h);
  int64_t old_length = res.length;

  /* strip left (5') end */
  if (opt_fastq_stripleft < res.length)
    {
      res.start += opt_fastq_stripleft;
      res.length -= opt_fastq_stripleft;
    }
  else
    {
      res.start = res.length;
      res.length = 0;
    }

  /* strip right (3') end */
  if (opt_fastq_stripright < res.length)
    res.length -= opt_fastq_stripright;
  else
    res.length = 0;

  /* truncate trailing (3') part */
  if (opt_fastq_trunclen >= 0)
    if (res.length > opt_fastq_trunclen)
      res.length = opt_fastq_trunclen;

  /* truncate trailing (3') part, but keep if short */
  if (opt_fastq_trunclen_keep >= 0)
    if (res.length > opt_fastq_trunclen_keep)
      res.length = opt_fastq_trunclen_keep;

  if (h->is_fastq)
    {
      /* truncate by quality and expected errors (ee) */
      res.ee = 0.0;
      char * q = fastx_get_quality(h) + res.start;
      for (int64_t i = 0; i < res.length; i++)
        {
          int qual = fastq_get_qual(q[i]);
          double e = exp10(-0.1 * qual);
          res.ee += e;

          if ((qual <= opt_fastq_truncqual) ||
              (res.ee > opt_fastq_truncee))
            {
              res.ee -= e;
              res.length = i;
              break;
            }
        }

      /* filter by expected errors (ee) */
      if (res.ee > opt_fastq_maxee)
        res.discarded = true;
      if ((res.length > 0) && (res.ee / res.length > opt_fastq_maxee_rate))
        res.discarded = true;
    }

  /* filter by length */
  if ((opt_fastq_trunclen >= 0) && (res.length < opt_fastq_trunclen))
    res.discarded = true;
  if (res.length < opt_fastq_minlen)
    res.discarded = true;
  if (res.length > opt_fastq_maxlen)
    res.discarded = true;

  /* filter by n's */
  int64_t ncount = 0;
  char * p = fastx_get_sequence(h) + res.start;
  for (int64_t i = 0; i < res.length; i++)
    {
      int pc = p[i];
      if ((pc == 'N') || (pc == 'n'))
        ncount++;
    }
  if (ncount > opt_fastq_maxns)
    res.discarded = true;

  /* filter by abundance */
  int64_t abundance = fastx_get_abundance(h);
  if (abundance < opt_minsize)
    res.discarded = true;
  if (abundance > opt_maxsize)
    res.discarded = true;

  res.truncated = res.length < old_length;

  return res;
}

void filter(bool fastq_only, char * filename)
{
  if ((!opt_fastqout) && (!opt_fastaout) &&
      (!opt_fastqout_discarded) && (!opt_fastaout_discarded) &&
      (!opt_fastqout_rev) && (!opt_fastaout_rev) &&
      (!opt_fastqout_discarded_rev) && (!opt_fastaout_discarded_rev))
    fatal("No output files specified");

  fastx_handle h1 = 0;
  fastx_handle h2 = 0;

  h1 = fastx_open(filename);

  if (!h1)
    fatal("Unrecognized file type (not proper FASTA or FASTQ format)");

  if (fastq_only && ! h1->is_fastq)
    fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");

  if ((opt_fastqout || opt_fastqout_discarded) && ! h1->is_fastq)
    fatal("Cannot write FASTQ output with FASTA input file (no quality scores)");

  uint64_t filesize = fastx_get_size(h1);

  if (opt_reverse)
    {
      h2 = fastx_open(opt_reverse);

      if (!h2)
        fatal("Unrecognized file type (not proper FASTA or FASTQ format) for reverse reads");

      if (fastq_only && ! h2->is_fastq)
        fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");

      if ((opt_fastqout_rev || opt_fastqout_discarded_rev) && ! h2->is_fastq)
        fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  FILE * fp_fastaout = 0;
  FILE * fp_fastqout = 0;
  FILE * fp_fastaout_discarded = 0;
  FILE * fp_fastqout_discarded = 0;

  FILE * fp_fastaout_rev = 0;
  FILE * fp_fastqout_rev = 0;
  FILE * fp_fastaout_discarded_rev = 0;
  FILE * fp_fastqout_discarded_rev = 0;

  if (opt_fastaout)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (!fp_fastaout)
        fatal("Unable to open FASTA output file for writing");
    }

  if (opt_fastqout)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (!fp_fastqout)
        fatal("Unable to open FASTQ output file for writing");
    }

  if (opt_fastaout_discarded)
    {
      fp_fastaout_discarded = fopen_output(opt_fastaout_discarded);
      if (!fp_fastaout_discarded)
        fatal("Unable to open FASTA output file for writing");
    }

  if (opt_fastqout_discarded)
    {
      fp_fastqout_discarded = fopen_output(opt_fastqout_discarded);
      if (!fp_fastqout_discarded)
        fatal("Unable to open FASTQ output file for writing");
    }

  if (h2)
    {
      if (opt_fastaout_rev)
        {
          fp_fastaout_rev = fopen_output(opt_fastaout_rev);
          if (!fp_fastaout_rev)
            fatal("Unable to open FASTA output file for writing");
        }

      if (opt_fastqout_rev)
        {
          fp_fastqout_rev = fopen_output(opt_fastqout_rev);
          if (!fp_fastqout_rev)
            fatal("Unable to open FASTQ output file for writing");
        }

      if (opt_fastaout_discarded_rev)
        {
          fp_fastaout_discarded_rev = fopen_output(opt_fastaout_discarded_rev);
          if (!fp_fastaout_discarded_rev)
            fatal("Unable to open FASTA output file for writing");
        }

      if (opt_fastqout_discarded_rev)
        {
          fp_fastqout_discarded_rev = fopen_output(opt_fastqout_discarded_rev);
          if (!fp_fastqout_discarded_rev)
            fatal("Unable to open FASTQ output file for writing");
        }
    }

  progress_init("Reading input file", filesize);

  int64_t kept = 0;
  int64_t discarded = 0;
  int64_t truncated = 0;

  while(fastx_next(h1, 0, chrmap_no_change))
    {
      if (h2 && ! fastx_next(h2, 0, chrmap_no_change))
        fatal("More forward reads than reverse reads");

      struct analysis_res res1 = { false, false, 0, 0, 0.0 } ;
      struct analysis_res res2 = { false, false, 0, 0, -1.0 } ;

      res1 = analyse(h1);
      if (h2)
        res2 = analyse(h2);

      if (res1.discarded || res2.discarded)
        {
          /* discard the sequence(s) */

          discarded++;

          if (opt_fastaout_discarded)
            fasta_print_general(fp_fastaout_discarded,
                                0,
                                fastx_get_sequence(h1) + res1.start,
                                res1.length,
                                fastx_get_header(h1),
                                fastx_get_header_length(h1),
                                fastx_get_abundance(h1),
                                discarded,
                                res1.ee,
                                -1,
                                -1,
                                0,
                                0.0);

          if (opt_fastqout_discarded)
            fastq_print_general(fp_fastqout_discarded,
                                fastx_get_sequence(h1) + res1.start,
                                res1.length,
                                fastx_get_header(h1),
                                fastx_get_header_length(h1),
                                fastx_get_quality(h1) + res1.start,
                                fastx_get_abundance(h1),
                                discarded,
                                res1.ee);

          if (h2)
            {
              if (opt_fastaout_discarded_rev)
                fasta_print_general(fp_fastaout_discarded_rev,
                                    0,
                                    fastx_get_sequence(h2) + res2.start,
                                    res2.length,
                                    fastx_get_header(h2),
                                    fastx_get_header_length(h2),
                                    fastx_get_abundance(h2),
                                    discarded,
                                    res2.ee,
                                    -1,
                                    -1,
                                    0,
                                    0.0);

              if (opt_fastqout_discarded_rev)
                fastq_print_general(fp_fastqout_discarded_rev,
                                    fastx_get_sequence(h2) + res2.start,
                                    res2.length,
                                    fastx_get_header(h2),
                                    fastx_get_header_length(h2),
                                    fastx_get_quality(h2) + res2.start,
                                    fastx_get_abundance(h2),
                                    discarded,
                                    res2.ee);
            }
        }
      else
        {
          /* keep the sequence(s) */

          kept++;

          if (res1.truncated || res2.truncated)
            truncated++;

          if (opt_fastaout)
            fasta_print_general(fp_fastaout,
                                0,
                                fastx_get_sequence(h1) + res1.start,
                                res1.length,
                                fastx_get_header(h1),
                                fastx_get_header_length(h1),
                                fastx_get_abundance(h1),
                                kept,
                                res1.ee,
                                -1,
                                -1,
                                0,
                                0.0);

          if (opt_fastqout)
            fastq_print_general(fp_fastqout,
                                fastx_get_sequence(h1) + res1.start,
                                res1.length,
                                fastx_get_header(h1),
                                fastx_get_header_length(h1),
                                fastx_get_quality(h1) + res1.start,
                                fastx_get_abundance(h1),
                                kept,
                                res1.ee);

          if (h2)
            {
              if (opt_fastaout_rev)
                fasta_print_general(fp_fastaout_rev,
                                    0,
                                    fastx_get_sequence(h1) + res2.start,
                                    res2.length,
                                    fastx_get_header(h2),
                                    fastx_get_header_length(h2),
                                    fastx_get_abundance(h2),
                                    kept,
                                    res2.ee,
                                    -1,
                                    -1,
                                    0,
                                    0.0);

              if (opt_fastqout_rev)
                fastq_print_general(fp_fastqout_rev,
                                    fastx_get_sequence(h2) + res2.start,
                                    res2.length,
                                    fastx_get_header(h2),
                                    fastx_get_header_length(h2),
                                    fastx_get_quality(h2) + res2.start,
                                    fastx_get_abundance(h2),
                                    kept,
                                    res2.ee);
            }
        }

      progress_update(fastx_get_position(h1));
    }

  progress_done();

  if (h2 && fastx_next(h2, 0, chrmap_no_change))
    fatal("More reverse reads than forward reads");

  if (! opt_quiet)
    fprintf(stderr,
            "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
            kept,
            truncated,
            discarded);

  if (opt_log)
    fprintf(fp_log,
            "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
            kept,
            truncated,
            discarded);

  if (h2)
    {
      if (opt_fastaout_rev)
        fclose(fp_fastaout_rev);

      if (opt_fastqout_rev)
        fclose(fp_fastqout_rev);

      if (opt_fastaout_discarded_rev)
        fclose(fp_fastaout_discarded_rev);

      if (opt_fastqout_discarded_rev)
        fclose(fp_fastqout_discarded_rev);

      fastx_close(h2);
    }

  if (opt_fastaout)
    fclose(fp_fastaout);

  if (opt_fastqout)
    fclose(fp_fastqout);

  if (opt_fastaout_discarded)
    fclose(fp_fastaout_discarded);

  if (opt_fastqout_discarded)
    fclose(fp_fastqout_discarded);

  fastx_close(h1);
}

void fastq_filter()
{
  filter(1, opt_fastq_filter);
}

void fastx_filter()
{
  filter(0, opt_fastx_filter);
}
