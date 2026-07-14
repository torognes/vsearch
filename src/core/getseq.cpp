/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "core/getseq.hpp"
#include "vsearch.h"
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "os/system.hpp"  // xstat_t, xfstat, S_ISFIFO
#include "utils/progress.hpp"
#include "utils/compare_strings_nocase.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/span.hpp"
#include <algorithm>  // std::copy, std::max, std::min, std::search, std::equal
#include <array>
#include <cassert>
#include <cctype>  // isalnum
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::snprintf, std::fileno, std::fgets, EOF, std::size_t
#include <cstring>  // std::strlen, std::strcpy, std::strstr
#include <vector>


static std::vector<std::vector<char>> labels_data;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto find_length_longest_label(std::vector<std::vector<char>> const & labels)
    -> std::size_t {
    auto longest = std::size_t{0};
    for (auto const & label : labels) {
      longest = std::max(longest, label.size());
    }
    return longest;
  }

}  // end of anonymous namespace


auto read_labels_file(char const * filename, struct Parameters const & parameters) -> void
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
  uint64_t const file_size = is_pipe ? 0: static_cast<uint64_t>(file_status.st_size);

  {
    Progress const progress("Reading labels", file_size, parameters);

    static constexpr auto a_memory_chunck = 1024U;
    while (true)
      {
        static constexpr auto buffer_size = 1024U;
        std::array<char, buffer_size> buffer {{}};
        auto const * return_value = std::fgets(buffer.data(), buffer_size, fp_labels.get());
        if (return_value == nullptr) { break; }

        auto length = std::strlen(buffer.data());
        if ((length != 0) and (buffer[length - 1] == '\n'))
          {
            buffer[length - 1] = '\0';
            --length;
          }

        // silently skip empty lines: an empty label would match every
        // header, and storing an empty std::vector<char> would make
        // label.data() return nullptr, crashing downstream scanners
        if (length == 0) { continue; }

        labels_longest = std::max(length, labels_longest);

        if (labels_count + 1 > labels_alloc)
          {
            labels_alloc += a_memory_chunck;
            labels_data.resize(labels_alloc);
          }
        labels_data[labels_count].resize(length);
        std::copy(buffer.begin(), std::next(buffer.begin(), static_cast<std::ptrdiff_t>(length)),
                  labels_data[labels_count].begin());
        ++labels_count;
      }
  }

  // definitive number of labels is known
  labels_data.resize(labels_count);

  static constexpr auto max_label_length = std::size_t{1023};
  if (labels_longest >= max_label_length)
    {
      if (not parameters.opt_quiet)
        {
          std::fprintf(stderr, "WARNING: Labels longer than 1023 characters are not supported\n");
        }

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log, "WARNING: Labels longer than 1023 characters are not supported\n");
        }
    }
}


auto test_label_match(fastx_handle input_handle, struct Parameters const & parameters) -> bool
{
  char const * header = fastx_get_header(input_handle);
  auto const header_length = fastx_get_header_length(input_handle);
  auto const hlen = static_cast<std::size_t>(header_length);
  auto const header_view = Span<char>{header, header_length};
  auto const longest_label = find_length_longest_label(labels_data);
  std::vector<char> field_buffer;
  std::size_t field_len = 0;
  if (parameters.opt_label_field != nullptr)
    {
      field_len = std::strlen(parameters.opt_label_field);
      std::size_t field_buffer_size = field_len + 2;
      if (parameters.opt_label_word != nullptr)
        {
          field_buffer_size += std::strlen(parameters.opt_label_word);
        }
      else
        {
          field_buffer_size += longest_label;
        }
      field_buffer.resize(field_buffer_size);
      std::snprintf(field_buffer.data(), field_buffer_size, "%s=", parameters.opt_label_field);
    }

  if (parameters.opt_label != nullptr)
    {
      auto const needle_view = Span<char>{parameters.opt_label, std::strlen(parameters.opt_label)};
      if (parameters.opt_label_substr_match)
        {
          return (contains_substring(header_view, needle_view));
        }
      return (are_same_string(header_view, needle_view));
    }
  if (parameters.opt_labels != nullptr)
    {
      if (parameters.opt_label_substr_match)
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
  else if (parameters.opt_label_word != nullptr)
    {
      char const * needle = parameters.opt_label_word;
      if (parameters.opt_label_field != nullptr)
        {
          std::strcpy(&field_buffer[field_len + 1], needle);
          needle = field_buffer.data();
        }
      auto const wlen = std::strlen(needle);
      char const * hit = header_view.data();
      while (true)
        {
          hit = std::strstr(hit, needle);
          if (hit == nullptr) {
            break;
          }
          if (parameters.opt_label_field != nullptr)
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
                   (std::isalnum(*(hit - 1)) == 0)) and
                  ((hit + wlen == header + hlen) or
                   (std::isalnum(*(hit + wlen)) == 0)))
                {
                  return true;
                }
            }
          ++hit;
        }
    }
  else if (parameters.opt_label_words != nullptr)
    {
      char const * const header_end = header + hlen;
      for (auto const & label: labels_data) {
        // labels read from a file are stored as std::vector<char>
        // without a trailing '\0', so the needle length must come from
        // label.size() and the search must be range-based; strlen and
        // strstr would read past the vector's storage
        char const * needle = label.data();
        std::size_t wlen = label.size();
        if (parameters.opt_label_field != nullptr)
          {
            std::copy(label.begin(), label.end(), &field_buffer[field_len + 1]);
            needle = field_buffer.data();
            wlen = field_len + 1 + label.size();
          }
        char const * const needle_end = needle + wlen;
        char const * hit = header;
        while (true)
          {
            hit = std::search(hit, header_end, needle, needle_end);
            if (hit == header_end) {
              break;
            }
            if (parameters.opt_label_field != nullptr)
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
                     (std::isalnum(*(hit - 1)) == 0)) and
                    ((hit + wlen == header + hlen) or
                     (std::isalnum(*(hit + wlen)) == 0)))
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


auto getseq(struct Parameters const & parameters, char const * filename) -> void
{
  if ((parameters.opt_fastqout == nullptr) and (parameters.opt_fastaout == nullptr) and
      (parameters.opt_notmatched == nullptr) and (parameters.opt_notmatchedfq == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_fastx_getseq != nullptr)
    {
      if (parameters.opt_label == nullptr)
        {
          fatal("Missing label option");
        }
    }
  else if (parameters.opt_fastx_getsubseq != nullptr)
    {
      if (parameters.opt_label == nullptr)
        {
          fatal("Missing label option");
        }

      if ((parameters.opt_subseq_start < 1) or (parameters.opt_subseq_end < 1))
        {
          fatal("The argument to options subseq_start and subseq_end must be at least 1");
        }

      if (parameters.opt_subseq_start > parameters.opt_subseq_end)
        {
          fatal("The argument to option subseq_start must be equal or less than to subseq_end");
        }
    }
  else if (parameters.opt_fastx_getseqs != nullptr)
    {
      int label_options = 0;
      if (parameters.opt_label != nullptr)
        {
          ++label_options;
        }
      if (parameters.opt_labels != nullptr)
        {
          ++label_options;
        }
      if (parameters.opt_label_word != nullptr)
        {
          ++label_options;
        }
      if (parameters.opt_label_words != nullptr)
        {
          ++label_options;
        }

      if (label_options != 1)
        {
          fatal("Specify one label option (label, labels, label_word or label_words)");
        }

      if (parameters.opt_labels != nullptr)
        {
          read_labels_file(parameters.opt_labels, parameters);
        }

      if (parameters.opt_label_words != nullptr)
        {
          read_labels_file(parameters.opt_label_words, parameters);
        }
    }

  fastx_handle h1 = nullptr;

  h1 = fastx_open(filename, parameters);

  if (((parameters.opt_fastqout != nullptr) or (parameters.opt_notmatchedfq != nullptr)) and not (h1->is_fastq or h1->is_empty))
    {
      fatal("Cannot write FASTQ output from FASTA input");
    }

  uint64_t const filesize = fastx_get_size(h1);

  auto fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  auto fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
  auto notmatched_handle = open_optional_output_file(parameters.opt_notmatched, OutputOption{"--notmatched"});
  auto notmatchedfq_handle = open_optional_output_file(parameters.opt_notmatchedfq, OutputOption{"--notmatchedfq"});
  std::FILE * const fp_fastaout = fastaout_handle.get();
  std::FILE * const fp_fastqout = fastqout_handle.get();
  std::FILE * const fp_notmatched = notmatched_handle.get();
  std::FILE * const fp_notmatchedfq = notmatchedfq_handle.get();

  int64_t kept = 0;
  int64_t discarded = 0;

  {
    Progress progress("Extracting sequences", filesize, parameters);
    while (fastx_next(h1, not parameters.opt_notrunclabels, chrmap_no_change()))
      {
        bool const match = test_label_match(h1, parameters);

        if (match)
          {
            /* keep the sequence(s) */

            ++kept;

            int64_t start = 1;
            int64_t end = static_cast<int64_t>(fastx_get_sequence_length(h1));
            if (parameters.opt_fastx_getsubseq != nullptr)
              {
                start = std::max(parameters.opt_subseq_start, start);
                end = std::min(parameters.opt_subseq_end, end);
              }
            /* When --subseq_start is past this sequence's length (trivially hit on
               a mixed-length file), end < start and the subsequence is empty. Emit
               it empty with in-bounds pointers rather than offsetting past the end
               with a negative length (S4). The guard must precede the offset of
               BOTH the sequence and the quality pointer. */
            int64_t length = end - start + 1;
            int64_t offset = start - 1;
            if (length <= 0)
              {
                length = 0;
                offset = 0;
              }

            if (parameters.opt_fastaout != nullptr)
              {
                fasta_print_general(fp_fastaout,
                                    nullptr,
                                    fastx_get_sequence(h1) + offset,
                                    static_cast<int>(length),
                                    fastx_get_header(h1),
                                    static_cast<int>(fastx_get_header_length(h1)),
                                    static_cast<uint64_t>(fastx_get_abundance(h1)),
                                    kept,
                                    -1.0,
                                    -1,
                                    -1,
                                    nullptr,
                                    0.0,
                                    0,
                                    parameters);
              }

            if (parameters.opt_fastqout != nullptr)
              {
                fastq_print_general(fp_fastqout,
                                    fastx_get_sequence(h1) + offset,
                                    static_cast<int>(length),
                                    fastx_get_header(h1),
                                    static_cast<int>(fastx_get_header_length(h1)),
                                    fastx_get_quality(h1) + offset,
                                    static_cast<uint64_t>(fastx_get_abundance(h1)),
                                    kept,
                                    -1.0,
                                    parameters);
              }
          }
        else
          {
            /* discard the sequence: non-matching sequences are always
               written in full, even when --subseq_start/--subseq_end
               are set (see vsearch-fastx_getsubseq(1)) */

            ++discarded;

            int64_t const length = static_cast<int64_t>(fastx_get_sequence_length(h1));

            if (parameters.opt_notmatched != nullptr)
              {
                fasta_print_general(fp_notmatched,
                                    nullptr,
                                    fastx_get_sequence(h1),
                                    static_cast<int>(length),
                                    fastx_get_header(h1),
                                    static_cast<int>(fastx_get_header_length(h1)),
                                    static_cast<uint64_t>(fastx_get_abundance(h1)),
                                    discarded,
                                    -1.0,
                                    -1,
                                    -1,
                                    nullptr,
                                    0.0,
                                    0,
                                    parameters);
              }

            if (parameters.opt_notmatchedfq != nullptr)
              {
                fastq_print_general(fp_notmatchedfq,
                                    fastx_get_sequence(h1),
                                    static_cast<int>(length),
                                    fastx_get_header(h1),
                                    static_cast<int>(fastx_get_header_length(h1)),
                                    fastx_get_quality(h1),
                                    static_cast<uint64_t>(fastx_get_abundance(h1)),
                                    discarded,
                                    -1.0,
                                    parameters);
              }
          }

        progress.update(fastx_get_position(h1));
      }
  }

  if (not parameters.opt_quiet)
    {
      std::fprintf(stderr,
              "%" PRId64 " of %" PRId64 " sequences extracted",
              kept,
              kept + discarded);
      if (kept + discarded > 0)
        {
          std::fprintf(stderr,
                  " (%.1lf%%)",
                  100.0 * static_cast<double>(kept) / static_cast<double>(kept + discarded));
        }
      std::fprintf(stderr, "\n");
    }

  if (parameters.opt_log != nullptr)
    {
      std::fprintf(parameters.fp_log,
              "%" PRId64 " of %" PRId64 " sequences extracted",
              kept,
              kept + discarded);
      if (kept + discarded > 0)
        {
          std::fprintf(parameters.fp_log,
                  " (%.1lf%%)",
                  100.0 * static_cast<double>(kept) / static_cast<double>(kept + discarded));
        }
      std::fprintf(parameters.fp_log, "\n");
    }

  if (parameters.opt_fastaout != nullptr)
    {
      fastaout_handle.reset();
    }

  if (parameters.opt_fastqout != nullptr)
    {
      fastqout_handle.reset();
    }

  if (parameters.opt_notmatched != nullptr)
    {
      notmatched_handle.reset();
    }

  if (parameters.opt_notmatchedfq != nullptr)
    {
      notmatchedfq_handle.reset();
    }

  fastx_close(h1, parameters);
}
