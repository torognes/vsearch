/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "maps.h"
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"
#include "utils/span.hpp"
#include <algorithm>  // std::copy, std::max, std::min, std::search, std::equal
#include <array>
#include <cassert>
#include <cctype>  // isalnum, std::toupper
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::snprintf, std::fileno, std::fgets, EOF, std::size_t
#include <cstring>  // std::strlen, std::strcpy, std::strstr
#include <vector>


std::vector<std::vector<char>> labels_data;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto compare_chars = [](char const lhs, char const rhs) -> bool {
    assert((lhs >= 0) or (lhs == EOF));
    auto const lhs_unsigned = static_cast<unsigned char>(lhs);
    auto const rhs_unsigned = static_cast<unsigned char>(rhs);
    return std::toupper(lhs_unsigned) == std::toupper(rhs_unsigned);
  };


  auto contains_substring(Span<char> const haystack, Span<char> const needle) -> bool {
    // case insensitive
    auto const hit = std::search(haystack.begin(), haystack.end(),
                                 needle.begin(), needle.end(),
                                 compare_chars);
    return (hit != haystack.end());
  }


  auto are_same_string(Span<char> const haystack, std::vector<char> const & needle) -> bool {
    // case insensitive
    if (haystack.size() != needle.size()) {
      return false;
    }
    return std::equal(haystack.begin(), haystack.end(),
                      needle.begin(), compare_chars);
  }


  auto are_same_string(Span<char> const haystack, Span<char> const needle) -> bool {
    // case insensitive
    if (haystack.size() != needle.size()) {
      return false;
    }
    return std::equal(haystack.begin(), haystack.end(),
                      needle.begin(), compare_chars);
  }


  auto find_length_longest_label(std::vector<std::vector<char>> const & labels)
    -> std::size_t {
    auto longest = std::size_t{0};
    for (auto const & label : labels) {
      longest = std::max(longest, label.size());
    }
    return longest;
  }

}  // end of anonymous namespace


auto read_labels_file(char * filename) -> void
{
  auto labels_alloc = 0U;
  auto labels_count = 0U;
  auto labels_longest = std::size_t{0};
  auto fp_labels = open_input_file(filename);
  if (fp_labels.get() == nullptr)
    {
      fatal("Unable to open labels file (%s)", filename);
    }

  xstat_t file_status;
  if (xfstat(fileno(fp_labels.get()), & file_status) != 0)
    {
      fatal("Unable to get status for labels file (%s)", filename);
    }

  auto const is_pipe = S_ISFIFO(file_status.st_mode);  // linuxism
  uint64_t const file_size = is_pipe ? 0: file_status.st_size;

  progress_init("Reading labels", file_size);

  static constexpr auto a_memory_chunck = 1024U;
  while (true)
    {
      static constexpr auto buffer_size = 1024U;
      std::array<char, buffer_size> buffer {{}};
      auto * return_value = std::fgets(buffer.data(), buffer_size, fp_labels.get());
      if (return_value == nullptr) { break; }

      auto length = std::strlen(buffer.data());
      if ((length != 0) and (buffer[length - 1] == '\n'))
        {
          buffer[length - 1] = '\0';
          --length;
        }

      labels_longest = std::max(length, labels_longest);

      if (labels_count + 1 > labels_alloc)
        {
          labels_alloc += a_memory_chunck;
          labels_data.resize(labels_alloc);
        }
      labels_data[labels_count].resize(length);
      std::copy(buffer.begin(), std::next(buffer.begin(), length),
                labels_data[labels_count].begin());
      ++labels_count;
    }

  progress_done();

  // definitive number of labels is known
  labels_data.resize(labels_count);

  static constexpr auto max_label_length = std::size_t{1023};
  if (labels_longest >= max_label_length)
    {
      if (not opt_quiet)
        {
          fprintf(stderr, "WARNING: Labels longer than 1023 characters are not supported\n");
        }

      if (opt_log != nullptr)
        {
          fprintf(fp_log, "WARNING: Labels longer than 1023 characters are not supported\n");
        }
    }
}


auto test_label_match(fastx_handle input_handle) -> bool
{
  char * header = fastx_get_header(input_handle);
  auto const header_length = fastx_get_header_length(input_handle);
  int const hlen = static_cast<int>(header_length);
  auto const header_view = Span<char>{header, header_length};
  auto const longest_label = find_length_longest_label(labels_data);
  std::vector<char> field_buffer;
  int field_len = 0;
  if (opt_label_field != nullptr)
    {
      field_len = std::strlen(opt_label_field);
      int field_buffer_size = field_len + 2;
      if (opt_label_word != nullptr)
        {
          field_buffer_size += std::strlen(opt_label_word);
        }
      else
        {
          field_buffer_size += static_cast<int>(longest_label);
        }
      field_buffer.resize(field_buffer_size);
      snprintf(field_buffer.data(), field_buffer_size, "%s=", opt_label_field);
    }

  if (opt_label != nullptr)
    {
      char * needle = opt_label;
      auto const needle_view = Span<char>{opt_label, std::strlen(needle)};
      if (opt_label_substr_match)
        {
          return (contains_substring(header_view, needle_view));
        }
      return (are_same_string(header_view, needle_view));
    }
  if (opt_labels != nullptr)
    {
      if (opt_label_substr_match)
        {
          for (auto & label: labels_data) {
            auto const label_view = Span<char>{label.data(), label.size()};
            if (contains_substring(header_view, label_view)) {
              return true;
            }
          }
        }
      else
        {
          for (auto const & label: labels_data) {
            if (are_same_string(header_view, label)) {
              return true;
            }
          }
        }
    }
  else if (opt_label_word != nullptr)
    {
      char * needle = opt_label_word;
      if (opt_label_field != nullptr)
        {
          std::strcpy(&field_buffer[field_len + 1], needle);
          needle = field_buffer.data();
        }
      int const wlen = std::strlen(needle);
      char * hit = header_view.data();
      while (true)
        {
          hit = std::strstr(hit, needle);
          if (hit == nullptr) {
            break;
          }
          if (opt_label_field != nullptr)
            {
              /* check of field */
              if (((hit == header) or
                   (*(hit - 1) == ';')) and
                  ((hit + wlen == header + hlen) or
                   (*(hit + wlen) == ';')))
                {
                  return true;
                }
            }
          else
            {
              /* check of full word */
              if (((hit == header) or
                   (isalnum(*(hit - 1)) == 0)) and
                  ((hit + wlen == header + hlen) or
                   (isalnum(*(hit + wlen)) == 0)))
                {
                  return true;
                }
            }
          ++hit;
        }
    }
  else if (opt_label_words != nullptr)
    {
      for (auto & label: labels_data) {
        char * needle = label.data();
        if (opt_label_field != nullptr)
          {
            std::copy(label.begin(), label.end(), &field_buffer[field_len + 1]);
            needle = field_buffer.data();
          }
        int const wlen = std::strlen(needle);
        char * hit = header;
        while (true)
          {
            hit = std::strstr(hit, needle);
            if (hit == nullptr) {
              break;
            }
            if (opt_label_field != nullptr)
              {
                /* check of field */
                if (((hit == header) or
                     (*(hit - 1) == ';')) and
                    ((hit + wlen == header + hlen) or
                     (*(hit + wlen) == ';')))
                  {
                    return true;
                  }
              }
            else
              {
                /* check of full word */
                if (((hit == header) or
                     (isalnum(*(hit - 1)) == 0)) and
                    ((hit + wlen == header + hlen) or
                     (isalnum(*(hit + wlen)) == 0)))
                  {
                    return true;
                  }
              }
            ++hit;
          }
      }  // end of range-for loop
    }
  return false;
}


auto getseq(struct Parameters const & parameters, char * filename) -> void
{
  if ((opt_fastqout == nullptr) and (opt_fastaout == nullptr) and
      (opt_notmatched == nullptr) and (opt_notmatchedfq == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_fastx_getseq != nullptr)
    {
      if (opt_label == nullptr)
        {
          fatal("Missing label option");
        }
    }
  else if (parameters.opt_fastx_getsubseq != nullptr)
    {
      if (opt_label == nullptr)
        {
          fatal("Missing label option");
        }

      if ((opt_subseq_start < 1) or (opt_subseq_end < 1))
        {
          fatal("The argument to options subseq_start and subseq_end must be at least 1");
        }

      if (opt_subseq_start > opt_subseq_end)
        {
          fatal("The argument to option subseq_start must be equal or less than to subseq_end");
        }
    }
  else if (parameters.opt_fastx_getseqs != nullptr)
    {
      int label_options = 0;
      if (opt_label != nullptr)
        {
          ++label_options;
        }
      if (opt_labels != nullptr)
        {
          ++label_options;
        }
      if (opt_label_word != nullptr)
        {
          ++label_options;
        }
      if (opt_label_words != nullptr)
        {
          ++label_options;
        }

      if (label_options != 1)
        {
          fatal("Specify one label option (label, labels, label_word or label_words)");
        }

      if (opt_labels != nullptr)
        {
          read_labels_file(opt_labels);
        }

      if (opt_label_words != nullptr)
        {
          read_labels_file(opt_label_words);
        }
    }

  fastx_handle h1 = nullptr;

  h1 = fastx_open(filename);

  if (h1 == nullptr)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  if (((opt_fastqout != nullptr) or (opt_notmatchedfq != nullptr)) and not (h1->is_fastq or h1->is_empty))
    {
      fatal("Cannot write FASTQ output from FASTA input");
    }

  uint64_t const filesize = fastx_get_size(h1);

  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_notmatched = nullptr;
  std::FILE * fp_notmatchedfq = nullptr;

  if (opt_fastaout != nullptr)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (fp_fastaout == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout != nullptr)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (fp_fastqout == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (opt_notmatched != nullptr)
    {
      fp_notmatched = fopen_output(opt_notmatched);
      if (fp_notmatched == nullptr)
        {
          fatal("Unable to open FASTA output file (notmatched) for writing");
        }
    }

  if (opt_notmatchedfq != nullptr)
    {
      fp_notmatchedfq = fopen_output(opt_notmatchedfq);
      if (fp_notmatchedfq == nullptr)
        {
          fatal("Unable to open FASTQ output file (notmatchedfq) for writing");
        }
    }

  progress_init("Extracting sequences", filesize);

  int64_t kept = 0;
  int64_t discarded = 0;

  while (fastx_next(h1, opt_notrunclabels == 0, chrmap_no_change))
    {
      bool const match = test_label_match(h1);

      int64_t start = 1;
      int64_t end = fastx_get_sequence_length(h1);
      if (parameters.opt_fastx_getsubseq != nullptr)
        {
          start = std::max(opt_subseq_start, start);
          end = std::min(opt_subseq_end, end);
        }
      int64_t const length = end - start + 1;

      if (match)
        {
          /* keep the sequence(s) */

          ++kept;

          if (opt_fastaout != nullptr)
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

          if (opt_fastqout != nullptr)
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

          ++discarded;

          if (opt_notmatched != nullptr)
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

          if (opt_notmatchedfq != nullptr)
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

  if (not opt_quiet)
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

  if (opt_log != nullptr)
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

  if (opt_fastaout != nullptr)
    {
      fclose(fp_fastaout);
    }

  if (opt_fastqout != nullptr)
    {
      fclose(fp_fastqout);
    }

  if (opt_notmatched != nullptr)
    {
      fclose(fp_notmatched);
    }

  if (opt_notmatchedfq != nullptr)
    {
      fclose(fp_notmatchedfq);
    }

  fastx_close(h1);
}


auto fastx_getseq(struct Parameters const & parameters) -> void
{
  getseq(parameters, parameters.opt_fastx_getseq);
}


auto fastx_getseqs(struct Parameters const & parameters) -> void
{
  getseq(parameters, parameters.opt_fastx_getseqs);
}


auto fastx_getsubseq(struct Parameters const & parameters) -> void
{
  getseq(parameters, parameters.opt_fastx_getsubseq);
}
