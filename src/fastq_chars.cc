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
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <vector>


constexpr unsigned int n_characters = 256;

struct statistics {
  std::vector<uint64_t> sequence_chars;
  std::vector<uint64_t> quality_chars;
  std::vector<uint64_t> tail_chars;
  std::vector<int> maxrun;
  uint64_t total_chars = 0;
  uint64_t seq_count = 0;
  int qmin_n = 255;
  int qmax_n = 0;
  char qmin = '\0';
  char qmax = '\0';
  char fastq_ascii = '\0';
  char fastq_qmin = '\0';
  char fastq_qmax = '\0';
};


auto fastq_chars(struct Parameters const & parameters) -> void
{
  struct statistics stats;
  stats.sequence_chars.resize(n_characters);
  stats.quality_chars.resize(n_characters);
  stats.tail_chars.resize(n_characters);
  stats.maxrun.resize(n_characters);

  fastx_handle h = fastq_open(parameters.opt_fastq_chars);

  auto const filesize = fastq_get_size(h);

  progress_init("Reading FASTQ file", filesize);

  while (fastq_next(h, false, chrmap_upcase))
    {
      int64_t const len = fastq_get_sequence_length(h);
      char * p = fastq_get_sequence(h);
      char * q = fastq_get_quality(h);

      stats.seq_count++;
      stats.total_chars += len;

      int run_char = -1;
      int run = 0;

      int64_t i = 0;
      while (i < len)
        {
          int const pc = *p++;
          int const qc = *q++;
          stats.sequence_chars[pc]++;
          stats.quality_chars[qc]++;

          if ((pc == 'N') || (pc == 'n'))
            {
              if (qc < stats.qmin_n)
                {
                  stats.qmin_n = qc;
                }
              if (qc > stats.qmax_n)
                {
                  stats.qmax_n = qc;
                }
            }

          if (pc == run_char)
            {
              run++;
              if (run > stats.maxrun[run_char])
                {
                  stats.maxrun[run_char] = run;
                }
            }
          else
            {
              run_char = pc;
              run = 0;
            }

          i++;
        }

      if (len >= parameters.opt_fastq_tail)
        {
          q = fastq_get_quality(h) + len - 1;
          int const tail_char = *q--;
          int tail_len = 1;
          while (*q-- == tail_char)
            {
              tail_len++;
              if (tail_len >= parameters.opt_fastq_tail)
                {
                  break;
                }
            }
          if (tail_len >= parameters.opt_fastq_tail)
            {
              stats.tail_chars[tail_char]++;
            }
        }

      progress_update(fastq_get_position(h));
    }
  progress_done();

  fastq_close(h);

  // refactor: find first non-null
  for (int c = 0; c <= 255; c++)
    {
      if (stats.quality_chars[c])
        {
          stats.qmin = c;
          break;
        }
    }

  // refactor: find last non-null
  for(int c = 255; c >= 0; c--)
    {
      if (stats.quality_chars[c])
        {
          stats.qmax = c;
          break;
        }
    }

  char fastq_qmax = '\0';

  if ((stats.qmin < 59) || (stats.qmax < 75))
    {
      stats.fastq_ascii = 33;
    }
  else
    {
      stats.fastq_ascii = 64;
    }

  fastq_qmax = stats.qmax - stats.fastq_ascii;
  stats.fastq_qmin = stats.qmin - stats.fastq_ascii;

  if (! parameters.opt_quiet)
    {
      fprintf(stderr, "Read %" PRIu64 " sequences.\n", stats.seq_count);

      if (stats.seq_count > 0)
        {
          fprintf(stderr, "Qmin %d, Qmax %d, Range %d\n",
                  stats.qmin, stats.qmax, stats.qmax - stats.qmin + 1);

          fprintf(stderr, "Guess: -fastq_qmin %d -fastq_qmax %d -fastq_ascii %d\n",
                  stats.fastq_qmin, fastq_qmax, stats.fastq_ascii);

          if (stats.fastq_ascii == 64)
            {
              if (stats.qmin < 64)
                {
                  fprintf(stderr, "Guess: Solexa format (phred+64)\n");
                }
              else if (stats.qmin < 66)
                {
                  fprintf(stderr, "Guess: Illumina 1.3+ format (phred+64)\n");
                }
              else
                {
                  fprintf(stderr, "Guess: Illumina 1.5+ format (phred+64)\n");
                }
            }
          else
            {
              if (stats.qmax > 73)
                {
                  fprintf(stderr, "Guess: Illumina 1.8+ format (phred+33)\n");
                }
              else
                {
                  fprintf(stderr, "Guess: Original Sanger format (phred+33)\n");
                }
            }

          fprintf(stderr, "\n");
          fprintf(stderr, "Letter          N   Freq MaxRun\n");
          fprintf(stderr, "------ ---------- ------ ------\n");

          for (int c = 0; c < 256; c++)
            {
              if (stats.sequence_chars[c] > 0)
                {
                  fprintf(stderr, "     %c %10" PRIu64 " %5.1f%% %6d",
                          c,
                          stats.sequence_chars[c],
                          100.0 * stats.sequence_chars[c] / stats.total_chars,
                          stats.maxrun[c]);
                  if ((c == 'N') || (c == 'n'))
                    {
                      if (stats.qmin_n < stats.qmax_n)
                        {
                          fprintf(stderr, "  Q=%c..%c", stats.qmin_n, stats.qmax_n);
                        }
                      else
                        {
                          fprintf(stderr, "  Q=%c", stats.qmin_n);
                        }
                    }
                  fprintf(stderr, "\n");
                }
            }

          fprintf(stderr, "\n");
          fprintf(stderr, "Char  ASCII    Freq       Tails\n");
          fprintf(stderr, "----  -----  ------  ----------\n");

          for (int c = stats.qmin; c <= stats.qmax; c++)
            {
              if (stats.quality_chars[c] > 0)
                {
                  fprintf(stderr, " '%c'  %5d  %5.1f%%  %10" PRIu64 "\n",
                          c,
                          c,
                          100.0 * stats.quality_chars[c] / stats.total_chars,
                          stats.tail_chars[c]);
                }
            }
        }
    }

  if (parameters.opt_log)
    {
      fprintf(fp_log, "Read %" PRIu64 " sequences.\n", stats.seq_count);

      if (stats.seq_count > 0)
        {
          fprintf(fp_log, "Qmin %d, Qmax %d, Range %d\n",
                  stats.qmin, stats.qmax, stats.qmax-stats.qmin + 1);

          fprintf(fp_log, "Guess: -fastq_qmin %d -fastq_qmax %d -fastq_ascii %d\n",
                  stats.fastq_qmin, fastq_qmax, stats.fastq_ascii);

          if (stats.fastq_ascii == 64)
            {
              if (stats.qmin < 64)
                {
                  fprintf(fp_log, "Guess: Solexa format (phred+64)\n");
                }
              else if (stats.qmin < 66)
                {
                  fprintf(fp_log, "Guess: Illumina 1.3+ format (phred+64)\n");
                }
              else
                {
                  fprintf(fp_log, "Guess: Illumina 1.5+ format (phred+64)\n");
                }
            }
          else
            {
              if (stats.qmax > 73)
                {
                  fprintf(fp_log, "Guess: Illumina 1.8+ format (phred+33)\n");
                }
              else
                {
                  fprintf(fp_log, "Guess: Original Sanger format (phred+33)\n");
                }
            }

          fprintf(fp_log, "\n");
          fprintf(fp_log, "Letter          N   Freq MaxRun\n");
          fprintf(fp_log, "------ ---------- ------ ------\n");

          for (int c = 0; c < 256; c++)
            {
              if (stats.sequence_chars[c] > 0)
                {
                  fprintf(fp_log, "     %c %10" PRIu64 " %5.1f%% %6d",
                          c,
                          stats.sequence_chars[c],
                          100.0 * stats.sequence_chars[c] / stats.total_chars,
                          stats.maxrun[c]);
                  if ((c == 'N') || (c == 'n'))
                    {
                      if (stats.qmin_n < stats.qmax_n)
                        {
                          fprintf(fp_log, "  Q=%c..%c", stats.qmin_n, stats.qmax_n);
                        }
                      else
                        {
                          fprintf(fp_log, "  Q=%c", stats.qmin_n);
                        }
                    }
                  fprintf(fp_log, "\n");
                }
            }

          fprintf(fp_log, "\n");
          fprintf(fp_log, "Char  ASCII    Freq       Tails\n");
          fprintf(fp_log, "----  -----  ------  ----------\n");

          for (int c = stats.qmin; c <= stats.qmax; c++)
            {
              if (stats.quality_chars[c] > 0)
                {
                  fprintf(fp_log, " '%c'  %5d  %5.1f%%  %10" PRIu64 "\n",
                          c,
                          c,
                          100.0 * stats.quality_chars[c] / stats.total_chars,
                          stats.tail_chars[c]);
                }
            }
        }
    }
}
