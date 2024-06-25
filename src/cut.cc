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
#include <vector>


struct statistics {
  uint64_t fragment_no = 0;
  uint64_t fragment_rev_no = 0;
  uint64_t fragment_discarded_no = 0;
  uint64_t fragment_discarded_rev_no = 0;
  int64_t cut = 0;
  int64_t uncut = 0;
};

struct a_file {
  char * name = nullptr;
  std::FILE * handle = nullptr;
};

struct a_strand {
  a_file forward;
  a_file reverse;
};

struct file_purpose {
  a_strand cut;
  a_strand discarded;
};


auto cut_one(fastx_handle input_handle,
             char * pattern,
             int pattern_length,
             int cut_fwd,
             int cut_rev,
             struct file_purpose const & fastaout,
             struct statistics & counters) -> int64_t
{
  char * seq  = fasta_get_sequence(input_handle);
  auto const seq_length = static_cast<int>(fasta_get_sequence_length(input_handle));

  /* get reverse complement */
  std::vector<char> rc_buffer_v(seq_length + 1);
  reverse_complement(rc_buffer_v.data(), seq, seq_length);

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
              if (fastaout.cut.forward.name != nullptr)
                {
                  fasta_print_general(fastaout.cut.forward.handle,
                                      nullptr,
                                      fasta_get_sequence(input_handle) + frag_start,
                                      frag_length,
                                      fasta_get_header(input_handle),
                                      fasta_get_header_length(input_handle),
                                      fasta_get_abundance(input_handle),
                                      ++counters.fragment_no,
                                      -1.0,
                                      -1,
                                      -1,
                                      nullptr,
                                      0.0);
                }
            }

          if (rc_length > 0)
            {
              if (fastaout.cut.reverse.name != nullptr)
                {
                  fasta_print_general(fastaout.cut.reverse.handle,
                                      nullptr,
                                      rc_buffer_v.data() + rc_start,
                                      rc_length,
                                      fasta_get_header(input_handle),
                                      fasta_get_header_length(input_handle),
                                      fasta_get_abundance(input_handle),
                                      ++counters.fragment_rev_no,
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
          if (fastaout.cut.forward.name != nullptr)
            {
              fasta_print_general(fastaout.cut.forward.handle,
                                  nullptr,
                                  fasta_get_sequence(input_handle) + frag_start,
                                  frag_length,
                                  fasta_get_header(input_handle),
                                  fasta_get_header_length(input_handle),
                                  fasta_get_abundance(input_handle),
                                  ++counters.fragment_no,
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
          if (fastaout.cut.reverse.name != nullptr)
            {
              fasta_print_general(fastaout.cut.reverse.handle,
                                  nullptr,
                                  rc_buffer_v.data() + rc_start,
                                  rc_length,
                                  fasta_get_header(input_handle),
                                  fasta_get_header_length(input_handle),
                                  fasta_get_abundance(input_handle),
                                  ++counters.fragment_rev_no,
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
      if (fastaout.discarded.forward.name != nullptr)
        {
          fasta_print_general(fastaout.discarded.forward.handle,
                              nullptr,
                              fasta_get_sequence(input_handle),
                              seq_length,
                              fasta_get_header(input_handle),
                              fasta_get_header_length(input_handle),
                              fasta_get_abundance(input_handle),
                              ++counters.fragment_discarded_no,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }

      if (fastaout.discarded.reverse.name != nullptr)
        {
          fasta_print_general(fastaout.discarded.reverse.handle,
                              nullptr,
                              rc_buffer_v.data(),
                              seq_length,
                              fasta_get_header(input_handle),
                              fasta_get_header_length(input_handle),
                              fasta_get_abundance(input_handle),
                              ++counters.fragment_discarded_rev_no,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }
    }

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


auto open_output_files(struct file_purpose & fastaout) -> void {
  if (fastaout.cut.forward.name != nullptr) {
    fastaout.cut.forward.handle = fopen_output(fastaout.cut.forward.name);
  }
  if (fastaout.discarded.forward.name != nullptr) {
    fastaout.discarded.forward.handle = fopen_output(fastaout.discarded.forward.name);
  }
  if (fastaout.cut.reverse.name != nullptr) {
    fastaout.cut.reverse.handle = fopen_output(fastaout.cut.reverse.name);
  }
  if (fastaout.discarded.reverse.name != nullptr) {
    fastaout.discarded.reverse.handle = fopen_output(fastaout.discarded.reverse.name);
  }
}


auto check_output_files(struct file_purpose const & fastaout) -> void {
  if (fastaout.cut.forward.name != nullptr) {
    if (fastaout.cut.forward.handle == nullptr) {
      fatal("Unable to open FASTA output file for writing");
    }
  }
  if (fastaout.discarded.forward.name != nullptr) {
    if (fastaout.discarded.forward.handle == nullptr) {
      fatal("Unable to open FASTA output file for writing");
    }
  }
  if (fastaout.cut.reverse.name != nullptr) {
    if (fastaout.cut.reverse.handle == nullptr) {
      fatal("Unable to open FASTQ output file for writing");
    }
  }
  if (fastaout.discarded.reverse.name != nullptr) {
    if (fastaout.discarded.reverse.handle == nullptr) {
      fatal("Unable to open FASTQ output file for writing");
    }
  }
}


auto close_output_files(struct file_purpose const & fastaout) -> void {
  for (auto * fp_outputfile : {
           fastaout.cut.forward.handle, fastaout.discarded.forward.handle,
           fastaout.cut.reverse.handle, fastaout.discarded.reverse.handle}) {
    if (fp_outputfile != nullptr) {
      static_cast<void>(std::fclose(fp_outputfile));
    }
  }
}


auto cut(struct Parameters const & parameters) -> void
{
  ckeck_if_output_is_set(parameters);

  struct statistics counters;
  struct file_purpose fastaout;
  fastaout.cut.forward.name = parameters.opt_fastaout;
  fastaout.discarded.forward.name = parameters.opt_fastaout_discarded;
  fastaout.cut.reverse.name = parameters.opt_fastaout_rev;
  fastaout.discarded.reverse.name = parameters.opt_fastaout_discarded_rev;

  fastx_handle input_handle = fasta_open(parameters.opt_cut);
  assert(input_handle != nullptr);  // verified by fasta_open()

  auto const filesize = fasta_get_size(input_handle);

  open_output_files(fastaout);
  check_output_files(fastaout);

  char * pattern = parameters.opt_cut_pattern;
  assert(pattern != nullptr);  // verified by <getopt.h>

  auto const pattern_length = static_cast<int>(strlen(pattern));

  if (pattern_length == 0)
    {
      fatal("Empty cut pattern string");
    }

  int cut_fwd = -1;
  int cut_rev = -1;

  int j = 0;  // number of nucleotides (pattern minus cutting sites)
  for (int i = 0; i < pattern_length ; i++)
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

  int64_t matches = 0;

  while (fasta_next(input_handle, false, chrmap_no_change))
    {
      auto const a_match = cut_one(input_handle,
                                   pattern,
                                   pattern_length - 2,
                                   cut_fwd,
                                   cut_rev,
                                   fastaout,
                                   counters);
      matches += a_match;
      if (a_match > 0)
        {
          ++counters.cut;
        }
      else
        {
          ++counters.uncut;
        }

      progress_update(fasta_get_position(input_handle));
    }

  progress_done();

  if (not parameters.opt_quiet)
    {
      static_cast<void>(std::fprintf(stderr,
              "%" PRId64 " sequence(s) cut %" PRId64 " times, %" PRId64 " sequence(s) never cut.\n",
                                     counters.cut, matches, counters.uncut));
    }

  if (parameters.opt_log != nullptr)
    {
      static_cast<void>(std::fprintf(fp_log,
              "%" PRId64 " sequence(s) cut %" PRId64 " times, %" PRId64 " sequence(s) never cut.\n",
                                     counters.cut, matches, counters.uncut));
    }

  close_output_files(fastaout);
  fasta_close(input_handle);
}
