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
#include "utils/maps.hpp"
#include <algorithm>  // std::count, std::for_each, std::equal
#include <cassert>
#include <cinttypes>  // macros PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <iterator>  // std::next
#include <string>
#include <utility>  // std::move
#include <vector>


struct statistics {
  int fragment_no = 0;
  int fragment_rev_no = 0;
  int fragment_discarded_no = 0;
  int fragment_discarded_rev_no = 0;
  int64_t cut = 0;
  int64_t uncut = 0;
  int64_t matches = 0;
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

struct restriction_pattern {
  std::string pattern;
  std::string coded_pattern;
  int cut_fwd;
  int cut_rev;
};


auto cut_a_sequence(fastx_handle input_handle,
             struct restriction_pattern const & restriction,
             struct file_purpose const & fastaout,
             struct statistics & counters,
             std::vector<char> & rc_buffer) -> void
{
  auto const pattern_length = static_cast<int>(restriction.pattern.size());
  char * seq = fasta_get_sequence(input_handle);
  auto const seq_length = static_cast<int>(fasta_get_sequence_length(input_handle));
  // failed refactoring: use transform to create a coded std::string
  // and find() to search for pattern occurrences, IUPAC chars make it
  // harder to compare sequences

  /* get reverse complement */
  rc_buffer.clear();
  rc_buffer.resize(seq_length + 1);
  reverse_complement(rc_buffer.data(), seq, seq_length);

  int64_t local_matches = 0;
  int frag_start = 0;
  int frag_length = seq_length;
  int rc_start = seq_length;
  int rc_length = 0;

  for (int i = 0; i < seq_length - pattern_length + 1; ++i)
    {
      auto const match = std::equal(restriction.coded_pattern.cbegin(),
                                    restriction.coded_pattern.cend(),
                                    std::next(seq, i),
                                    [](char const & lhs, char const & rhs) -> bool {
                                      auto const lhs_unsigned = static_cast<unsigned char>(lhs);
                                      auto const rhs_unsigned = chrmap_4bit_vector[static_cast<unsigned char>(rhs)];
                                      return ((lhs_unsigned & rhs_unsigned) != 0);  // explanation needed
                                    });

      if (not match) {
        continue;
      }

      ++local_matches;

      frag_length = i + restriction.cut_fwd - frag_start;

      rc_length = rc_start - (seq_length - (i + restriction.cut_rev));
      rc_start -= rc_length;

      if ((frag_length > 0) and (fastaout.cut.forward.name != nullptr))
        {
          fasta_print_general(fastaout.cut.forward.handle,
                              nullptr,
                              std::next(seq, frag_start),
                              frag_length,
                              fasta_get_header(input_handle),
                              static_cast<int>(fasta_get_header_length(input_handle)),
                              fasta_get_abundance(input_handle),
                              ++counters.fragment_no,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }

      if ((rc_length > 0) and (fastaout.cut.reverse.name != nullptr))
        {
          fasta_print_general(fastaout.cut.reverse.handle,
                              nullptr,
                              &rc_buffer[rc_start],
                              rc_length,
                              fasta_get_header(input_handle),
                              static_cast<int>(fasta_get_header_length(input_handle)),
                              fasta_get_abundance(input_handle),
                              ++counters.fragment_rev_no,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }

      frag_start += frag_length;
    }

  if (local_matches > 0)
    {
      ++counters.cut;
      frag_length = seq_length - frag_start;
      rc_length = rc_start;
      rc_start = 0;
    }

  if ((local_matches > 0) and (frag_length > 0) and (fastaout.cut.forward.name != nullptr))
    {
      fasta_print_general(fastaout.cut.forward.handle,
                          nullptr,
                          std::next(seq, frag_start),
                          frag_length,
                          fasta_get_header(input_handle),
                          static_cast<int>(fasta_get_header_length(input_handle)),
                          fasta_get_abundance(input_handle),
                          ++counters.fragment_no,
                          -1.0,
                          -1,
                          -1,
                          nullptr,
                          0.0);
    }

  if ((local_matches > 0) and (rc_length > 0) and (fastaout.cut.reverse.name != nullptr))
    {
      fasta_print_general(fastaout.cut.reverse.handle,
                          nullptr,
                          &rc_buffer[rc_start],
                          rc_length,
                          fasta_get_header(input_handle),
                          static_cast<int>(fasta_get_header_length(input_handle)),
                          fasta_get_abundance(input_handle),
                          ++counters.fragment_rev_no,
                          -1.0,
                          -1,
                          -1,
                          nullptr,
                          0.0);
    }

  if (local_matches == 0)
    {
      ++counters.uncut;
    }

  if ((local_matches == 0) and (fastaout.discarded.forward.name != nullptr))
    {
      fasta_print_general(fastaout.discarded.forward.handle,
                          nullptr,
                          seq,
                          seq_length,
                          fasta_get_header(input_handle),
                          static_cast<int>(fasta_get_header_length(input_handle)),
                          fasta_get_abundance(input_handle),
                          ++counters.fragment_discarded_no,
                          -1.0,
                          -1,
                          -1,
                          nullptr,
                          0.0);
    }

  if ((local_matches == 0) and (fastaout.discarded.reverse.name != nullptr))
    {
      fasta_print_general(fastaout.discarded.reverse.handle,
                          nullptr,
                          rc_buffer.data(),
                          seq_length,
                          fasta_get_header(input_handle),
                          static_cast<int>(fasta_get_header_length(input_handle)),
                          fasta_get_abundance(input_handle),
                          ++counters.fragment_discarded_rev_no,
                          -1.0,
                          -1,
                          -1,
                          nullptr,
                          0.0);
    }

  counters.matches += local_matches;
}


auto ckeck_if_output_is_set(struct Parameters const & parameters) -> void {
  if ((parameters.opt_fastaout == nullptr) and
      (parameters.opt_fastaout_discarded == nullptr) and
      (parameters.opt_fastaout_rev == nullptr) and
      (parameters.opt_fastaout_discarded_rev == nullptr)) {
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


auto check_if_contains_circumflex(std::string const & pattern) -> void {
  auto const occurrences = std::count(pattern.cbegin(), pattern.cend(), '^');
  if (occurrences == 0) {
    fatal("No forward sequence cut site (^) found in pattern");
  }
  if (occurrences > 1) {
    fatal("Multiple cut sites not supported");
  }
}


auto check_if_contains_underscore(std::string const & pattern) -> void {
  auto const occurrences = std::count(pattern.cbegin(), pattern.cend(), '_');
  if (occurrences == 0) {
    fatal("No reverse sequence cut site (_) found in pattern");
  }
  if (occurrences > 1) {
    fatal("Multiple cut sites not supported");
  }
}


auto locate_forward_restriction_site(std::string pattern) -> int {
  auto const underscore_position = pattern.find('_');
  pattern.erase(underscore_position, 1);
  return static_cast<int>(pattern.find('^'));
}


auto locate_reverse_restriction_site(std::string pattern) -> int {
  auto const circumflex_position = pattern.find('^');
  pattern.erase(circumflex_position, 1);
  return static_cast<int>(pattern.find('_'));
}


auto remove_restriction_sites(std::string pattern) -> std::string {
  auto const circumflex_position = pattern.find('^');
  pattern.erase(circumflex_position, 1);
  auto const underscore_position = pattern.find('_');
  return pattern.erase(underscore_position, 1);
}


auto reencode_restriction_pattern(std::string raw_pattern) -> std::string {
  auto pattern = remove_restriction_sites(std::move(raw_pattern));
  auto encode_characters = [](char const & character) -> char {
    auto const symbol_uchar = static_cast<unsigned char>(character);
    auto const coded_symbol_uchar = chrmap_4bit_vector[symbol_uchar];
    return static_cast<char>(coded_symbol_uchar);
  };
  std::transform(pattern.cbegin(), pattern.cend(),
                 pattern.begin(), encode_characters);
  return pattern;
}


auto search_illegal_characters(std::string const & pattern) -> void {
  auto character_is_illegal = [](char const & character) {
    auto const unsigned_character = static_cast<unsigned char>(character);
    if (chrmap_4bit_vector[unsigned_character] == 0) {
      fatal("Illegal character in cut pattern");
    }
  };
  std::for_each(pattern.cbegin(), pattern.cend(), character_is_illegal);
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

  auto const raw_pattern = parameters.opt_cut_pattern;
  // assert(pattern != nullptr);  // verified by <getopt.h>

  // check for the expected number of restriction sites
  check_if_contains_circumflex(raw_pattern);
  check_if_contains_underscore(raw_pattern);

  // locate restriction sites and trim pattern
  struct restriction_pattern const restriction = {
    remove_restriction_sites(raw_pattern),
    reencode_restriction_pattern(raw_pattern),
    locate_forward_restriction_site(raw_pattern),
    locate_reverse_restriction_site(raw_pattern)
  };

  search_illegal_characters(restriction.pattern);

  if (restriction.pattern.empty())
    {
      fatal("Empty cut pattern string");
    }

  progress_init("Cutting sequences", filesize);

  std::vector<char> rc_buffer;
  while (fasta_next(input_handle, false, chrmap_no_change_array.data()))
    {
      cut_a_sequence(input_handle,
              restriction,
              fastaout,
              counters,
              rc_buffer);

      progress_update(fasta_get_position(input_handle));
    }

  progress_done();

  if (not parameters.opt_quiet)
    {
      static_cast<void>(std::fprintf(stderr,
              "%" PRId64 " sequence(s) cut %" PRId64 " times, %" PRId64 " sequence(s) never cut.\n",
                                     counters.cut, counters.matches, counters.uncut));
    }

  if (parameters.opt_log != nullptr)
    {
      static_cast<void>(std::fprintf(fp_log,
              "%" PRId64 " sequence(s) cut %" PRId64 " times, %" PRId64 " sequence(s) never cut.\n",
                                     counters.cut, counters.matches, counters.uncut));
    }

  close_output_files(fastaout);
  fasta_close(input_handle);
}
