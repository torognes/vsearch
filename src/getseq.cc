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

/* Implement fastx_getseq, fastx_getseqs and fastx_getsubseq as described here:
   https://drive5.com/usearch/manual/cmd_fastx_getseqs.html                  */

#include "vsearch.h"
#include <cctype>  // isalnum
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::snprintf, std::fileno, std::fgets
#include <cstdlib>  // std::realloc, std::free
#include <cstring>  // std::strlen, std::memset, std::strcpy, std::strstr
#include <string.h>  // strdup


static int labels_alloc = 0;
static int labels_count = 0;
static int labels_longest = 0;
static char * * labels_data = nullptr;

auto read_labels_file(char * filename) -> void
{
  FILE * fp_labels = fopen_input(filename);
  if (! fp_labels)
    {
      fatal("Unable to open labels file (%s)", filename);
    }

  xstat_t fs;
  if (xfstat(fileno(fp_labels), & fs))
    {
      fatal("Unable to get status for labels file (%s)", filename);
    }

  bool const is_pipe = S_ISFIFO(fs.st_mode);
  uint64_t file_size = 0;
  if (! is_pipe)
    {
      file_size = fs.st_size;
    }

  progress_init("Reading labels", file_size);

  while (true)
    {
      const int buffer_size = 1024;
      char buffer[buffer_size];
      char * ret = fgets(buffer, buffer_size, fp_labels);
      if (ret)
        {
          int len = strlen(buffer);
          if ((len > 0) && (buffer[len - 1] == '\n'))
            {
              buffer[len - 1] = 0;
              len--;
            }

          if (len > labels_longest)
            {
              labels_longest = len;
            }

          if (labels_count + 1 > labels_alloc)
            {
              labels_alloc += 1024;
              labels_data = (char * *) realloc(labels_data,
                                               labels_alloc * sizeof (char*));
              if (! labels_data)
                {
                  fatal("Unable to allocate memory for labels");
                }
            }
          labels_data[labels_count++] = strdup(buffer);
        }
      else
        {
          break;
        }
    }

  fclose(fp_labels);
  progress_done();

  if (labels_longest >= 1023)
    {
      if (! opt_quiet)
        {
          fprintf(stderr, "WARNING: Labels longer than 1023 characters are not supported\n");
        }

      if (opt_log)
        {
          fprintf(fp_log, "WARNING: Labels longer than 1023 characters are not supported\n");
        }
    }
}

auto free_labels() -> void
{
  for(int i = 0; i < labels_count; i++)
    {
      free(labels_data[i]);
    }
  free(labels_data);
  labels_data = nullptr;
}

auto test_label_match(fastx_handle h) -> bool
{
  char * header = fastx_get_header(h);
  int const hlen = fastx_get_header_length(h);
  char * field_buffer = nullptr;
  int field_len = 0;
  if (opt_label_field)
    {
      field_len = strlen(opt_label_field);
      int field_buffer_size = field_len + 2;
      if (opt_label_word)
        {
          field_buffer_size += strlen(opt_label_word);
        }
      else
        {
          field_buffer_size += labels_longest;
        }
      field_buffer = (char *) xmalloc(field_buffer_size);
      snprintf(field_buffer, field_buffer_size, "%s=", opt_label_field);
    }

  if (opt_label)
    {
      char * needle = opt_label;
      int const wlen = strlen(needle);
      if (opt_label_substr_match)
        {
          return xstrcasestr(header, needle);
        }
      else
        {
          return (hlen == wlen) && ! strcasecmp(header, needle);
        }
    }
  else if (opt_labels)
    {
      if (opt_label_substr_match)
        {
          for (int i = 0; i < labels_count; i++)
            {
              if (xstrcasestr(header, labels_data[i]))
                {
                  return true;
                }
            }
        }
      else
        {
          for (int i = 0; i < labels_count; i++)
            {
              char * needle = labels_data[i];
              int const wlen = strlen(needle);
              if ((hlen == wlen) && ! strcasecmp(header, needle))
                {
                  return true;
                }
            }
        }
    }
  else if (opt_label_word)
    {
      char * needle = opt_label_word;
      if (opt_label_field)
        {
          strcpy(field_buffer + field_len + 1, needle);
          needle = field_buffer;
        }
      int const wlen = strlen(needle);
      char * hit = header;
      while (true)
        {
          hit = strstr(hit, needle);
          if (hit)
            {
              if (opt_label_field)
                {
                  /* check of field */
                  if (((hit == header) ||
                       (*(hit - 1) == ';')) &&
                      ((hit + wlen == header + hlen) ||
                       (*(hit + wlen) == ';')))
                    {
                      return true;
                    }
                }
              else
                {
                  /* check of full word */
                  if (((hit == header) ||
                       (! isalnum(*(hit - 1)))) &&
                      ((hit + wlen == header + hlen) ||
                       (! isalnum(*(hit + wlen)))))
                    {
                      return true;
                    }
                }
              hit++;
            }
          else
            {
              break;
            }
        }
    }
  else if (opt_label_words)
    {
      for (int i = 0; i < labels_count; i++)
        {
          char * needle = labels_data[i];
          if (opt_label_field)
            {
              strcpy(field_buffer + field_len + 1, needle);
              needle = field_buffer;
            }
          int const wlen = strlen(needle);
          char * hit = header;
          while (true)
            {
              hit = strstr(hit, needle);
              if (hit)
                {
                  if (opt_label_field)
                    {
                      /* check of field */
                      if (((hit == header) ||
                           (*(hit - 1) == ';')) &&
                          ((hit + wlen == header + hlen) ||
                           (*(hit + wlen) == ';')))
                        {
                          return true;
                        }
                    }
                  else
                    {
                      /* check of full word */
                      if (((hit == header) ||
                           (! isalnum(*(hit - 1)))) &&
                          ((hit + wlen == header + hlen) ||
                           (! isalnum(*(hit + wlen)))))
                        {
                          return true;
                        }
                    }
                  hit++;
                }
              else
                {
                  break;
                }
            }
        }
    }
  return false;
}

auto getseq(char * filename) -> void
{
  if ((! opt_fastqout) && (! opt_fastaout) &&
      (! opt_notmatched) && (! opt_notmatchedfq))
    {
      fatal("No output files specified");
    }

  if (opt_fastx_getseq)
    {
      if (! opt_label)
        {
          fatal("Missing label option");
        }
    }
  else if (opt_fastx_getsubseq)
    {
      if (! opt_label)
        {
          fatal("Missing label option");
        }

      if ((opt_subseq_start < 1) || (opt_subseq_end < 1))
        {
          fatal("The argument to options subseq_start and subseq_end must be at least 1");
        }

      if (opt_subseq_start > opt_subseq_end)
        {
          fatal("The argument to option subseq_start must be equal or less than to subseq_end");
        }
    }
  else if (opt_fastx_getseqs)
    {
      int label_options = 0;
      if (opt_label)
        {
          label_options++;
        }
      if (opt_labels)
        {
          label_options++;
        }
      if (opt_label_word)
        {
          label_options++;
        }
      if (opt_label_words)
        {
          label_options++;
        }

      if (label_options != 1)
        {
          fatal("Specify one label option (label, labels, label_word or label_words)");
        }

      if (opt_labels)
        {
          read_labels_file(opt_labels);
        }

      if (opt_label_words)
        {
          read_labels_file(opt_label_words);
        }
    }

  fastx_handle h1 = nullptr;

  h1 = fastx_open(filename);

  if (! h1)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  if ((opt_fastqout || opt_notmatchedfq) && ! (h1->is_fastq || h1->is_empty))
    {
      fatal("Cannot write FASTQ output from FASTA input");
    }

  uint64_t const filesize = fastx_get_size(h1);

  FILE * fp_fastaout = nullptr;
  FILE * fp_fastqout = nullptr;
  FILE * fp_notmatched = nullptr;
  FILE * fp_notmatchedfq = nullptr;

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

  if (opt_notmatched)
    {
      fp_notmatched = fopen_output(opt_notmatched);
      if (! fp_notmatched)
        {
          fatal("Unable to open FASTA output file (notmatched) for writing");
        }
    }

  if (opt_notmatchedfq)
    {
      fp_notmatchedfq = fopen_output(opt_notmatchedfq);
      if (! fp_notmatchedfq)
        {
          fatal("Unable to open FASTQ output file (notmatchedfq) for writing");
        }
    }

  progress_init("Extracting sequences", filesize);

  int64_t kept = 0;
  int64_t discarded = 0;

  while(fastx_next(h1, ! opt_notrunclabels, chrmap_no_change))
    {
      bool const match = test_label_match(h1);

      int64_t start = 1;
      int64_t end = fastx_get_sequence_length(h1);
      if (opt_fastx_getsubseq)
        {
          if (opt_subseq_start > start)
            {
              start = opt_subseq_start;
            }
          if (opt_subseq_end < end)
            {
              end = opt_subseq_end;
            }
        }
      int64_t const length = end - start + 1;

      if (match)
        {
          /* keep the sequence(s) */

          kept++;

          if (opt_fastaout)
            {
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  fastx_get_sequence(h1) + start - 1,
                                  length,
                                  fastx_get_header(h1),
                                  fastx_get_header_length(h1),
                                  fastx_get_abundance(h1),
                                  kept,
                                  -1.0,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout)
            {
              fastq_print_general(fp_fastqout,
                                  fastx_get_sequence(h1) + start - 1,
                                  length,
                                  fastx_get_header(h1),
                                  fastx_get_header_length(h1),
                                  fastx_get_quality(h1) + start - 1,
                                  fastx_get_abundance(h1),
                                  kept,
                                  -1.0);
            }
        }
      else
        {
          /* discard the sequence */

          discarded++;

          if (opt_notmatched)
            {
              fasta_print_general(fp_notmatched,
                                  nullptr,
                                  fastx_get_sequence(h1) + start - 1,
                                  length,
                                  fastx_get_header(h1),
                                  fastx_get_header_length(h1),
                                  fastx_get_abundance(h1),
                                  discarded,
                                  -1.0,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_notmatchedfq)
            {
              fastq_print_general(fp_notmatchedfq,
                                  fastx_get_sequence(h1) + start - 1,
                                  length,
                                  fastx_get_header(h1),
                                  fastx_get_header_length(h1),
                                  fastx_get_quality(h1) + start - 1,
                                  fastx_get_abundance(h1),
                                  discarded,
                                  -1.0);
            }
        }

      progress_update(fastx_get_position(h1));
    }

  progress_done();

  if (! opt_quiet)
    {
      fprintf(stderr,
              "%" PRId64 " of %" PRId64 " sequences extracted",
              kept,
              kept + discarded);
      if (kept + discarded > 0)
        {
          fprintf(stderr,
                  " (%.1lf%%)",
                  100.0 * kept / (kept + discarded));
        }
      fprintf(stderr, "\n");
    }

  if (opt_log)
    {
      fprintf(fp_log,
              "%" PRId64 " of %" PRId64 " sequences extracted",
              kept,
              kept + discarded);
      if (kept + discarded > 0)
        {
          fprintf(fp_log,
                  " (%.1lf%%)",
                  100.0 * kept / (kept + discarded));
        }
      fprintf(fp_log, "\n");
    }

  if (opt_fastaout)
    {
      fclose(fp_fastaout);
    }

  if (opt_fastqout)
    {
      fclose(fp_fastqout);
    }

  if (opt_notmatched)
    {
      fclose(fp_notmatched);
    }

  if (opt_notmatchedfq)
    {
      fclose(fp_notmatchedfq);
    }

  fastx_close(h1);

  if (opt_labels || opt_label_words)
    {
      free_labels();
    }
}

auto fastx_getseq() -> void
{
  getseq(opt_fastx_getseq);
}

auto fastx_getseqs() -> void
{
  getseq(opt_fastx_getseqs);
}

auto fastx_getsubseq() -> void
{
  getseq(opt_fastx_getsubseq);
}
