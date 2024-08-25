/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "maps.h"
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstring>  // std::memset
#include <limits>
#include <vector>


auto q2p(double quality_value) -> double
{
  static constexpr auto base = 10.0;
  return std::pow(base, -quality_value / base);
}


auto fastq_stats() -> void
{
  auto input_handle = fastq_open(opt_fastq_stats);

  auto const filesize = fastq_get_size(input_handle);

  progress_init("Reading FASTQ file", filesize);

  uint64_t seq_count = 0;
  uint64_t symbols = 0;

  int64_t read_length_alloc = 512;

  auto * read_length_table = (uint64_t *) xmalloc(sizeof(uint64_t) * read_length_alloc);
  memset(read_length_table, 0, sizeof(uint64_t) * read_length_alloc);

  auto * qual_length_table = (uint64_t *) xmalloc(sizeof(uint64_t) * read_length_alloc * 256);
  memset(qual_length_table, 0, sizeof(uint64_t) * read_length_alloc * 256);

  auto * ee_length_table = (uint64_t *) xmalloc(sizeof(uint64_t) * read_length_alloc * 4);
  memset(ee_length_table, 0, sizeof(uint64_t) * read_length_alloc * 4);

  auto * q_length_table = (uint64_t *) xmalloc(sizeof(uint64_t) * read_length_alloc * 4);
  memset(q_length_table, 0, sizeof(uint64_t) * read_length_alloc * 4);

  auto * sumee_length_table = (double *) xmalloc(sizeof(double) * read_length_alloc);
  memset(sumee_length_table, 0, sizeof(double) * read_length_alloc);

  int64_t len_min = std::numeric_limits<long>::max();
  int64_t len_max = 0;

  auto qmin = std::numeric_limits<int>::max();
  auto qmax = std::numeric_limits<int>::min();

  std::vector<uint64_t> quality_chars(256);

  while (fastq_next(input_handle, false, chrmap_upcase))
    {
      ++seq_count;

      auto const len = static_cast<int64_t>(fastq_get_sequence_length(input_handle));
      auto * q = fastq_get_quality(input_handle);

      /* update length statistics */

      if (len + 1 > read_length_alloc)
        {
          read_length_table = (uint64_t *) xrealloc(read_length_table,
                                                   sizeof(uint64_t) * (len + 1));
          memset(read_length_table + read_length_alloc, 0,
                 sizeof(uint64_t) * (len + 1 - read_length_alloc));

          qual_length_table = (uint64_t *) xrealloc(qual_length_table,
                                                   sizeof(uint64_t) * (len + 1) * 256);
          memset(qual_length_table + (256 * read_length_alloc), 0,
                 sizeof(uint64_t) * (len + 1 - read_length_alloc) * 256);

          ee_length_table = (uint64_t *) xrealloc(ee_length_table,
                                                 sizeof(uint64_t) * (len + 1) * 4);
          memset(ee_length_table + (4 * read_length_alloc), 0,
                 sizeof(uint64_t) * (len + 1 - read_length_alloc) * 4);

          q_length_table = (uint64_t *) xrealloc(q_length_table,
                                                sizeof(uint64_t) * (len + 1) * 4);
          memset(q_length_table + (4 * read_length_alloc), 0,
                 sizeof(uint64_t) * (len + 1 - read_length_alloc) * 4);

          sumee_length_table = (double *) xrealloc(sumee_length_table,
                                                   sizeof(double) * (len + 1));
          memset(sumee_length_table + read_length_alloc, 0,
                 sizeof(double) * (len + 1 - read_length_alloc));

          read_length_alloc = len + 1;
        }

      ++read_length_table[len];

      len_min = std::min(len, len_min);
      len_max = std::max(len, len_max);

      /* update quality statistics */

      symbols += len;

      std::array<double, 4> const ee_limits = { 1.0, 0.5, 0.25, 0.1 };

      double ee = 0.0;
      int qmin_this = std::numeric_limits<int>::max();
      for (int64_t i = 0; i < len; i++)
        {
          int const qc = q[i];

          int const qual = qc - opt_fastq_ascii;
          if ((qual < opt_fastq_qmin) || (qual > opt_fastq_qmax))
            {
              char * msg = nullptr;
              if (xsprintf(& msg,
                           "FASTQ quality value (%d) out of range (%" PRId64 "-%" PRId64 ").\n"
                           "Please adjust the FASTQ quality base character or range with the\n"
                           "--fastq_ascii, --fastq_qmin or --fastq_qmax options. For a complete\n"
                           "diagnosis with suggested values, please run vsearch --fastq_chars file.",
                           qual, opt_fastq_qmin, opt_fastq_qmax) > 0)
                {
                  fatal(msg);
                }
              else
                {
                  fatal("Out of memory");
                }
              xfree(msg);
            }

          ++quality_chars[qc];
          qmin = std::min(qc, qmin);
          qmax = std::max(qc, qmax);

          ++qual_length_table[(256 * i) + qc];

          ee += q2p(qual);

          sumee_length_table[i] += ee;

          for (int z = 0; z < 4; z++)
            {
              if (ee <= ee_limits[z])
                {
                  ++ee_length_table[(4 * i) + z];
                }
              else
                {
                  break;
                }
            }

          qmin_this = std::min(qual, qmin_this);

          for (int z = 0; z < 4; z++)
            {
              if (qmin_this > 5 * (z + 1))
                {
                  ++q_length_table[(4 * i) + z];
                }
              else
                {
                  break;
                }
            }
        }

      progress_update(fastq_get_position(input_handle));
    }
  progress_done();

  /* compute various distributions */

  std::vector<uint64_t> length_dist(len_max + 1);
  std::vector<int64_t> symb_dist(len_max + 1);
  std::vector<double> rate_dist(len_max + 1);
  std::vector<double> avgq_dist(len_max + 1);
  std::vector<double> avgee_dist(len_max + 1);
  std::vector<double> avgp_dist(len_max + 1);

  int64_t length_accum = 0;
  int64_t symb_accum = 0;

  for (int64_t i = 0; i <= len_max; i++)
    {
      length_accum += read_length_table[i];
      length_dist[i] = length_accum;

      symb_accum += seq_count - length_accum;
      symb_dist[i] = symb_accum;

      int64_t q = 0;
      int64_t x = 0;
      double e_sum = 0.0;
      for (int c = qmin; c <= qmax; c++)
        {
          int const qual = c - opt_fastq_ascii;
          x += qual_length_table[(256 * i) + c];
          q += qual_length_table[(256 * i) + c] * qual;
          e_sum += qual_length_table[(256 * i) + c] * q2p(qual);
        }
      avgq_dist[i] = 1.0 * q / x;
      avgp_dist[i] = e_sum / x;
      avgee_dist[i] = sumee_length_table[i] / x;
      rate_dist[i] = avgee_dist[i] / (i + 1);
    }

  if (fp_log)
    {
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Read length distribution\n");
      fprintf(fp_log, "      L           N      Pct   AccPct\n");
      fprintf(fp_log, "-------  ----------  -------  -------\n");

      for (int64_t i = len_max; i >= len_min; i--)
        {
          if (read_length_table[i] > 0)
            {
              fprintf(fp_log, "%2s%5" PRId64 "  %10" PRIu64 "   %5.1lf%%   %5.1lf%%\n",
                      (i == len_max ? ">=" : "  "),
                      i,
                      read_length_table[i],
                      read_length_table[i] * 100.0 / seq_count,
                      100.0 * (seq_count - length_dist[i - 1]) / seq_count);
            }
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "Q score distribution\n");
      fprintf(fp_log, "ASCII    Q       Pe           N      Pct   AccPct\n");
      fprintf(fp_log, "-----  ---  -------  ----------  -------  -------\n");

      int64_t qual_accum = 0;
      for (int c = qmax ; c >= qmin ; c--)
        {
          if (quality_chars[c] > 0)
            {
              qual_accum += quality_chars[c];
              fprintf(fp_log,
                      "    %c  %3" PRId64 "  %7.5lf  %10" PRIu64 "  %6.1lf%%  %6.1lf%%\n",
                      c,
                      c - opt_fastq_ascii,
                      q2p(c - opt_fastq_ascii),
                      quality_chars[c],
                      100.0 * quality_chars[c] / symbols,
                      100.0 * qual_accum / symbols);
            }
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct\n");
      fprintf(fp_log, "-----  -------  ----  -------  --------  -----  ---------  --------\n");

      for (int64_t i = 2; i <= len_max; i++)
        {
          double const PctRecs = 100.0 * (seq_count - length_dist[i - 1]) / seq_count;
          double const AvgQ = avgq_dist[i - 1];
          double const AvgP = avgp_dist[i - 1];
          double const AvgEE = avgee_dist[i - 1];
          double const Rate = rate_dist[i - 1];

          fprintf(fp_log,
                  "%5" PRId64 "  %6.1lf%%  %4.1lf  %7.5lf  %8.6lf  %5.2lf  %9.6lf  %7.3lf%%\n",
                  i,
                  PctRecs,
                  AvgQ,
                  q2p(AvgQ),
                  AvgP,
                  AvgEE,
                  Rate,
                  100.0 * Rate);
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000\n");
      fprintf(fp_log, "-----  -------  -------  -------  -------  -------  -------  -------  -------\n");

      for (int64_t i = len_max; i >= 1; i--)
        {
          int64_t read_count[4];
          double read_percentage[4];

          for (int z = 0; z < 4; z++)
            {
              read_count[z] = ee_length_table[(4 * (i - 1)) + z];
              read_percentage[z] = 100.0 * read_count[z] / seq_count;
            }

          if (read_count[0] > 0)
            {
              fprintf(fp_log,
                      "%5" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  %7" PRId64 "  "
                      "%6.2lf%%  %6.2lf%%  %6.2lf%%  %6.2lf%%\n",
                      i,
                      read_count[0], read_count[1],
                      read_count[2], read_count[3],
                      read_percentage[0], read_percentage[1],
                      read_percentage[2], read_percentage[3]);
            }
        }


      fprintf(fp_log, "\n");
      fprintf(fp_log, "Truncate at first Q\n");
      fprintf(fp_log, "  Len     Q=5    Q=10    Q=15    Q=20\n");
      fprintf(fp_log, "-----  ------  ------  ------  ------\n");

      for (int64_t i = len_max; i >= MAX(1, len_max / 2); i--)
        {
          double read_percentage[4];

          for (int z = 0; z < 4; z++)
            {
              read_percentage[z] = 100.0 * q_length_table[(4 * (i - 1)) + z] / seq_count;
            }

          fprintf(fp_log, "%5" PRId64 "  %5.1lf%%  %5.1lf%%  %5.1lf%%  %5.1lf%%\n",
                  i,
                  read_percentage[0], read_percentage[1],
                  read_percentage[2], read_percentage[3]);
        }

      fprintf(fp_log, "\n");
      fprintf(fp_log, "%10" PRIu64 "  Recs (%.1lfM), 0 too long\n",
              seq_count, seq_count / 1.0e6);
      if (seq_count > 0)
        {
          fprintf(fp_log, "%10.1lf  Avg length\n", 1.0 * symbols / seq_count);
        }
      fprintf(fp_log, "%9.1lfM  Bases\n", symbols / 1.0e6);
    }

  xfree(read_length_table);
  xfree(qual_length_table);
  xfree(ee_length_table);
  xfree(q_length_table);
  xfree(sumee_length_table);

  fastq_close(input_handle);

  if (! opt_quiet)
    {
      fprintf(stderr, "Read %" PRIu64 " sequences.\n", seq_count);
    }
}


auto fastx_revcomp() -> void
{
  uint64_t buffer_alloc = 512;
  char * seq_buffer = (char*) xmalloc(buffer_alloc);
  char * qual_buffer = (char*) xmalloc(buffer_alloc);

  if ((! opt_fastaout) && (! opt_fastqout)) {
    fatal("No output files specified");
  }

  fastx_handle h = fastx_open(opt_fastx_revcomp);

  if (! h)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  if (opt_fastqout && ! (h->is_fastq || h->is_empty))
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  uint64_t const filesize = fastx_get_size(h);

  FILE * fp_fastaout = nullptr;
  FILE * fp_fastqout = nullptr;

  if (opt_fastaout)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (! fp_fastaout)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (! fp_fastqout)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (h->is_fastq)
    {
      progress_init("Reading FASTQ file", filesize);
    }
  else
    {
      progress_init("Reading FASTA file", filesize);
    }

  int count = 0;

  while (fastx_next(h, false, chrmap_no_change))
    {
      ++count;

      /* header */

      uint64_t const hlen = fastx_get_header_length(h);
      char * header = fastx_get_header(h);
      int64_t const abundance = fastx_get_abundance(h);


      /* sequence */

      uint64_t const length = fastx_get_sequence_length(h);

      if (length + 1 > buffer_alloc)
        {
          buffer_alloc = length + 1;
          seq_buffer = (char *) xrealloc(seq_buffer, buffer_alloc);
          qual_buffer = (char *) xrealloc(qual_buffer, buffer_alloc);
        }

      char * p = fastx_get_sequence(h);
      reverse_complement(seq_buffer, p, length);


      /* quality values */

      char * q = fastx_get_quality(h);

      if (fastx_is_fastq(h))
        {
          /* reverse quality values */
          for (uint64_t i = 0; i < length; i++)
            {
              qual_buffer[i] = q[length - 1 - i];
            }
          qual_buffer[length] = 0;
        }

      if (opt_fastaout)
        {
          fasta_print_general(fp_fastaout,
                              nullptr,
                              seq_buffer,
                              length,
                              header,
                              hlen,
                              abundance,
                              count,
                              -1.0,
                              -1, -1, nullptr, 0.0);
        }

      if (opt_fastqout)
        {
          fastq_print_general(fp_fastqout,
                              seq_buffer,
                              length,
                              header,
                              hlen,
                              qual_buffer,
                              abundance,
                              count,
                              -1.0);
        }

      progress_update(fastx_get_position(h));
    }
  progress_done();

  if (opt_fastaout)
    {
      fclose(fp_fastaout);
    }

  if (opt_fastqout)
    {
      fclose(fp_fastqout);
    }

  fastx_close(h);

  xfree(seq_buffer);
  xfree(qual_buffer);
}


auto fastq_convert() -> void
{
  if (! opt_fastqout) {
    fatal("No output file specified with --fastqout");
  }

  auto input_handle = fastq_open(opt_fastq_convert);

  if (! input_handle)
    {
      fatal("Unable to open FASTQ file");
    }

  auto const filesize = fastq_get_size(input_handle);

  FILE * fp_fastqout = nullptr;

  fp_fastqout = fopen_output(opt_fastqout);
  if (! fp_fastqout)
    {
      fatal("Unable to open FASTQ output file for writing");
    }

  progress_init("Reading FASTQ file", filesize);

  auto j = 1;
  while (fastq_next(input_handle, false, chrmap_no_change))
    {
      /* header */

      auto * header = fastq_get_header(input_handle);
      auto const abundance = fastq_get_abundance(input_handle);

      /* sequence */

      auto const length = fastq_get_sequence_length(input_handle);
      auto * sequence = fastq_get_sequence(input_handle);

      /* convert quality values */

      auto * quality = fastq_get_quality(input_handle);
      for (uint64_t i = 0; i < length; i++)
        {
          int q = quality[i] - opt_fastq_ascii;
          if (q < opt_fastq_qmin)
            {
              fprintf(stderr,
                      "\nFASTQ quality score (%d) below minimum (%" PRId64
                      ") in entry no %" PRIu64
                      " starting on line %" PRIu64 "\n",
                      q,
                      opt_fastq_qmin,
                      fastq_get_seqno(input_handle) + 1,
                      fastq_get_lineno(input_handle));
              fatal("FASTQ quality score too low");
            }
          if (q > opt_fastq_qmax)
            {
              fprintf(stderr,
                      "\nFASTQ quality score (%d) above maximum (%" PRId64
                      ") in entry no %" PRIu64
                      " starting on line %" PRIu64 "\n",
                      q,
                      opt_fastq_qmax,
                      fastq_get_seqno(input_handle) + 1,
                      fastq_get_lineno(input_handle));
              fatal("FASTQ quality score too high");
            }
          q = std::max<int64_t>(q, opt_fastq_qminout);
          q = std::min<int64_t>(q, opt_fastq_qmaxout);
          q += opt_fastq_asciiout;
          q = std::max(q, 33);
          q = std::min(q, 126);
          quality[i] = q;
        }
      quality[length] = 0;

      int const hlen = fastq_get_header_length(input_handle);
      fastq_print_general(fp_fastqout,
                          sequence,
                          length,
                          header,
                          hlen,
                          quality,
                          abundance,
                          j,
                          -1.0);

      ++j;
      progress_update(fastq_get_position(input_handle));
    }

  progress_done();

  fclose(fp_fastqout);
  fastq_close(input_handle);
}
