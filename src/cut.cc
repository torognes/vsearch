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
#include <cassert>
#include <cinttypes>  // macros PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstring>  // std::strlen


static uint64_t fragment_no = 0;
static uint64_t fragment_rev_no = 0;
static uint64_t fragment_discarded_no = 0;
static uint64_t fragment_discarded_rev_no = 0;

struct a_file {
  char * name = nullptr;
  std::FILE * handle = nullptr;
};

struct a_strand {
  a_file forward;
  a_file reverse;
};

struct file_purpose {
  a_strand kept;
  a_strand discarded;
};


auto cut_one(fastx_handle h,
            std::FILE * fp_fastaout,
            std::FILE * fp_fastaout_discarded,
            std::FILE * fp_fastaout_rev,
            std::FILE * fp_fastaout_discarded_rev,
            char * pattern,
            int pattern_length,
            int cut_fwd,
            int cut_rev) -> int64_t
{
  char * seq  = fasta_get_sequence(h);
  auto const seq_length = static_cast<int>(fasta_get_sequence_length(h));

  /* get reverse complement */
  char * rc = (char *) xmalloc(seq_length + 1);
  reverse_complement(rc, seq, seq_length);

  int frag_start = 0;
  int frag_length = seq_length;
  int64_t matches = 0;

  int rc_start = seq_length;
  int rc_length = 0;

  for (int i = 0; i < seq_length - pattern_length + 1; i++)
    {
      auto match = true;
      for (int j = 0; j < pattern_length; j++)
        {
          if ((chrmap_4bit[(unsigned char) (pattern[j])] &
               chrmap_4bit[(unsigned char) (seq[i + j])]) == 0)
            {
              match = false;
              break;
            }
        }

      if (match)
        {
          ++matches;

          frag_length = i + cut_fwd - frag_start;

          rc_length = rc_start - (seq_length - (i + cut_rev));
          rc_start -= rc_length;

          if (frag_length > 0)
            {
              if (opt_fastaout)
                {
                  fasta_print_general(fp_fastaout,
                                      nullptr,
                                      fasta_get_sequence(h) + frag_start,
                                      frag_length,
                                      fasta_get_header(h),
                                      fasta_get_header_length(h),
                                      fasta_get_abundance(h),
                                      ++fragment_no,
                                      -1.0,
                                      -1,
                                      -1,
                                      nullptr,
                                      0.0);
                }
            }

          if (rc_length > 0)
            {
              if (opt_fastaout_rev)
                {
                  fasta_print_general(fp_fastaout_rev,
                                      nullptr,
                                      rc + rc_start,
                                      rc_length,
                                      fasta_get_header(h),
                                      fasta_get_header_length(h),
                                      fasta_get_abundance(h),
                                      ++fragment_rev_no,
                                      -1.0,
                                      -1,
                                      -1,
                                      nullptr,
                                      0.0);
                }
            }

          frag_start += frag_length;
        }
    }

  if (matches > 0)
    {
      frag_length = seq_length - frag_start;

      if (frag_length > 0)
        {
          if (opt_fastaout)
            {
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  fasta_get_sequence(h) + frag_start,
                                  frag_length,
                                  fasta_get_header(h),
                                  fasta_get_header_length(h),
                                  fasta_get_abundance(h),
                                  ++fragment_no,
                                  -1.0,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }
        }

      rc_length = rc_start;
      rc_start = 0;

      if (rc_length > 0)
        {
          if (opt_fastaout_rev)
            {
              fasta_print_general(fp_fastaout_rev,
                                  nullptr,
                                  rc + rc_start,
                                  rc_length,
                                  fasta_get_header(h),
                                  fasta_get_header_length(h),
                                  fasta_get_abundance(h),
                                  ++fragment_rev_no,
                                  -1.0,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }
        }
    }
  else
    {
      if (opt_fastaout_discarded)
        {
          fasta_print_general(fp_fastaout_discarded,
                              nullptr,
                              fasta_get_sequence(h),
                              seq_length,
                              fasta_get_header(h),
                              fasta_get_header_length(h),
                              fasta_get_abundance(h),
                              ++fragment_discarded_no,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }

      if (opt_fastaout_discarded_rev)
        {
          fasta_print_general(fp_fastaout_discarded_rev,
                              nullptr,
                              rc,
                              seq_length,
                              fasta_get_header(h),
                              fasta_get_header_length(h),
                              fasta_get_abundance(h),
                              ++fragment_discarded_rev_no,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }
    }

  xfree(rc);

  return matches;
}


auto ckeck_if_output_is_set(struct Parameters const & parameters) -> void
{
  if ((parameters.opt_fastaout == nullptr) and
      (parameters.opt_fastaout_discarded == nullptr) and
      (parameters.opt_fastaout_rev == nullptr) and
      (parameters.opt_fastaout_discarded_rev == nullptr))
    {
      fatal("No output files specified");
    }
}


auto cut(struct Parameters const & parameters) -> void
{
  ckeck_if_output_is_set(parameters);

  struct file_purpose ouput;
  ouput.kept.forward.name = parameters.opt_fastaout;
  ouput.discarded.forward.name = parameters.opt_fastaout_discarded;
  ouput.kept.reverse.name = parameters.opt_fastaout_rev;
  ouput.discarded.reverse.name = parameters.opt_fastaout_discarded_rev;

  fastx_handle h = fasta_open(parameters.opt_cut);
  assert(h != nullptr);  // verified by fasta_open()

  auto const filesize = fasta_get_size(h);

  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastaout_discarded = nullptr;
  std::FILE * fp_fastaout_rev = nullptr;
  std::FILE * fp_fastaout_discarded_rev = nullptr;

  if (parameters.opt_fastaout)
    {
      fp_fastaout = fopen_output(parameters.opt_fastaout);
      if (not fp_fastaout)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (parameters.opt_fastaout_rev)
    {
      fp_fastaout_rev = fopen_output(parameters.opt_fastaout_rev);
      if (not fp_fastaout_rev)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (parameters.opt_fastaout_discarded)
    {
      fp_fastaout_discarded = fopen_output(parameters.opt_fastaout_discarded);
      if (not fp_fastaout_discarded)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (parameters.opt_fastaout_discarded_rev)
    {
      fp_fastaout_discarded_rev = fopen_output(parameters.opt_fastaout_discarded_rev);
      if (not fp_fastaout_discarded_rev)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  char * pattern = parameters.opt_cut_pattern;
  assert(pattern != nullptr);  // verified by <getopt.h>

  int const n = strlen(pattern);

  if (n == 0)
    {
      fatal("Empty cut pattern string");
    }

  int cut_fwd = -1;
  int cut_rev = -1;

  int j = 0;  // number of nucleotides (pattern minus cutting sites)
  for (int i = 0; i < n ; i++)
    {
      unsigned char const x = pattern[i];
      if (x == '^')
        {
          if (cut_fwd != -1)
            {
              fatal("Multiple cut sites not supported");

            }
          cut_fwd = j;
        }
      else if (x == '_')
        {
          if (cut_rev != -1)
            {
              fatal("Multiple cut sites not supported");

            }
          cut_rev = j;
        }
      else if (chrmap_4bit[(unsigned int) x])
        {
          pattern[j] = x;
          ++j;
        }
      else
        {
          fatal("Illegal character in cut pattern");
        }
    }

  if (cut_fwd < 0)
    {
      fatal("No forward sequence cut site (^) found in pattern");
    }

  if (cut_rev < 0)
    {
      fatal("No reverse sequence cut site (_) found in pattern");
    }

  progress_init("Cutting sequences", filesize);

  int64_t cut = 0;
  int64_t uncut = 0;
  int64_t matches = 0;

  while (fasta_next(h, false, chrmap_no_change))
    {
      auto const a_match = cut_one(h,
                                   fp_fastaout,
                                   fp_fastaout_discarded,
                                   fp_fastaout_rev,
                                   fp_fastaout_discarded_rev,
                                   pattern,
                                   n - 2,
                                   cut_fwd,
                                   cut_rev);
      matches += a_match;
      if (a_match > 0)
        {
          ++cut;
        }
      else
        {
          ++uncut;
        }

      progress_update(fasta_get_position(h));
    }

  progress_done();

  if (not parameters.opt_quiet)
    {
      std::fprintf(stderr,
              "%" PRId64 " sequence(s) cut %" PRId64 " times, %" PRId64 " sequence(s) never cut.\n",
              cut, matches, uncut);
    }

  if (parameters.opt_log)
    {
      std::fprintf(fp_log,
              "%" PRId64 " sequence(s) cut %" PRId64 " times, %" PRId64 " sequence(s) never cut.\n",
              cut, matches, uncut);
    }

  if (parameters.opt_fastaout)
    {
      fclose(fp_fastaout);
    }

  if (parameters.opt_fastaout_rev)
    {
      fclose(fp_fastaout_rev);
    }

  if (parameters.opt_fastaout_discarded)
    {
      fclose(fp_fastaout_discarded);
    }

  if (parameters.opt_fastaout_discarded_rev)
    {
      fclose(fp_fastaout_discarded_rev);
    }

  fasta_close(h);
}
