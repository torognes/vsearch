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

void filter(bool fastq_only, char * filename)
{
  fastx_handle h = fastx_open(filename);

  if (!h)
    fatal("Unrecognized file type (not proper FASTA or FASTQ format)");

  if (fastq_only && ! h->is_fastq)
    fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");

  if ((opt_fastqout || opt_fastqout_discarded) && ! h->is_fastq)
    fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");

  uint64_t filesize = fastx_get_size(h);

  FILE * fp_fastaout = 0;
  FILE * fp_fastqout = 0;
  FILE * fp_fastaout_discarded = 0;
  FILE * fp_fastqout_discarded = 0;

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

  uint64_t header_alloc = 0;
  char * header = 0;
  if (opt_relabel)
    {
      header_alloc = strlen(opt_relabel) + 25;
      header = (char*) xmalloc(header_alloc);
    }

  progress_init("Reading input file", filesize);

  int64_t kept = 0;
  int64_t discarded = 0;
  int64_t truncated = 0;

  while(fastx_next(h, 0, chrmap_no_change))
    {
      int64_t length = fastx_get_sequence_length(h);
      char * d = fastx_get_header(h);
      char * p = fastx_get_sequence(h);
      char * q = fastx_get_quality(h);
      int64_t abundance = fastx_get_abundance(h);

      /* strip initial part */
      if (opt_fastq_stripleft > 0)
        {
          if (opt_fastq_stripleft < length)
            {
              p += opt_fastq_stripleft;
              q += opt_fastq_stripleft;
              length -= opt_fastq_stripleft;
            }
          else
            {
              p += length;
              q += length;
              length = 0;
            }
        }

      /* strip right end */
      if (opt_fastq_stripright > 0)
        {
          if (opt_fastq_stripright < length)
            length -= opt_fastq_stripright;
          else
            length = 0;
        }

      /* truncate trailing part */
      if (opt_fastq_trunclen >= 0)
        {
          if (length >= opt_fastq_trunclen)
            length = opt_fastq_trunclen;
          else
            length = 0;
        }

      /* truncate trailing part, but keep if short */
      if ((opt_fastq_trunclen_keep >= 0) && (length > opt_fastq_trunclen_keep))
        length = opt_fastq_trunclen_keep;

      /* quality and ee truncation */
      double ee = 0.0;
      if (h->is_fastq)
        {
          for (int64_t i = 0; i < length; i++)
            {
              int qual = fastq_get_qual(q[i]);
              ee += exp10(- qual / 10.0);

              if ((qual <= opt_fastq_truncqual) ||
                  (ee > opt_fastq_truncee))
                {
                  ee -= exp10(- qual / 10.0);
                  length = i;
                  break;
                }
            }
        }

      /* count n's */
      int64_t ncount = 0;
      for (int64_t i = 0; i < length; i++)
        {
          int pc = p[i];
          if ((pc == 'N') || (pc == 'n'))
            ncount++;
        }

      if ((length >= opt_fastq_minlen) &&
          (length <= opt_fastq_maxlen) &&
          ((opt_fastq_trunclen < 0) || (length >= opt_fastq_trunclen)) &&
          (ncount <= opt_fastq_maxns) &&
          (ee <= opt_fastq_maxee) &&
          ((length == 0) || (ee / length <= opt_fastq_maxee_rate)) &&
          ((opt_minsize == 0) || (abundance >= opt_minsize)) &&
          ((opt_maxsize == 0) || (abundance <= opt_maxsize)))
        {
          /* keep the sequence */

          kept++;

          if ((uint64_t)(length) < fastx_get_sequence_length(h))
            {
              truncated++;
              p[length] = 0;
              if (h->is_fastq)
                q[length] = 0;
            }

          if (opt_fastaout)
            fasta_print_general(fp_fastaout,
                                0,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                abundance,
                                kept,
                                -1,
                                -1,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);

          if (opt_fastqout)
            fastq_print_general(fp_fastqout,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                q,
                                abundance,
                                kept,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);
        }
      else
        {
          /* discard the sequence */

          discarded++;

          if (opt_fastaout_discarded)
            fasta_print_general(fp_fastaout_discarded,
                                0,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                abundance,
                                discarded,
                                -1,
                                -1,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);

          if (opt_fastqout_discarded)
            fastq_print_general(fp_fastqout_discarded,
                                p,
                                length,
                                d,
                                fastx_get_header_length(h),
                                q,
                                abundance,
                                discarded,
                                (opt_eeout || opt_fastq_eeout) ? "ee" : 0,
                                ee);
        }

      progress_update(fastx_get_position(h));
    }
  progress_done();

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

  if (header)
    xfree(header);

  if (opt_fastaout)
    fclose(fp_fastaout);

  if (opt_fastqout)
    fclose(fp_fastqout);

  if (opt_fastaout_discarded)
    fclose(fp_fastaout_discarded);

  if (opt_fastqout_discarded)
    fclose(fp_fastqout_discarded);

  fastx_close(h);
}

void fastq_filter()
{
  filter(1, opt_fastq_filter);
}

void fastx_filter()
{
  filter(0, opt_fastx_filter);
}
