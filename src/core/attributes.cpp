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

#include "utils/view.hpp"
#include "core/attributes.hpp"  // View<char>
#include "vsearch.hpp"  // struct Parameters
#include "utils/fatal.hpp"
#include "utils/header_attribute.hpp"  // Attribute, header_find_attribute
#include "utils/print_view.hpp"  // fprint
#include "utils/sequence_digest.hpp"  // fprint_seq_digest_sha1, fprint_seq_digest_md5
#include <algorithm>  // std::find_if, std::search, std::sort
#include <array>
#include <cerrno>  // errno
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fputs
#include <cstdlib>  // std::strtoll
#include <cstring>  // std::strlen
#include <iterator>  // std::next, std::distance


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  using vsearch::Attribute;
  using vsearch::Attribute_span;
  using vsearch::Value_chars;
  using vsearch::header_find_attribute;


  /* The three numeric annotations vsearch reads back from a header. All three
     reject an empty value -- "size=" with nothing after it has never been an
     abundance -- which is what tells them apart from the four names
     core/otutable.cpp and core/tax.cpp read (see
     utils/header_attribute.hpp). */
  struct Attributes {
    Attribute ee {"ee=", 3, Value_chars::digits_and_dot, false};
    Attribute length {"length=", 7, Value_chars::digits, false};
    Attribute size {"size=", 5, Value_chars::digits, false};
  };

  constexpr Attributes attributes;


  constexpr auto n_expected_attributes = std::size_t{3};  // 3 attributes: size, ee, length


}  // end of anonymous namespace


auto header_get_size(View<char> const header) -> int64_t {
  /* read size/abundance annotation */
  static constexpr auto decimal_base = 10;
  auto const annotation = header_find_attribute(header, attributes.size);
  if (not annotation.present) {
    return 0;  // refactoring: return 1 by default?
  }

  char * next_character = nullptr;
  // C++17 refactoring: replace strtoll with std::from_chars
  /* attribute_value() rather than std::next(header.data(), ...): the offset
     comes from a match position inside the header, and the subspan it builds
     is where an out-of-range one would be caught. std::strtoll still needs the
     bare pointer, and stops at the first non-digit -- the ';' or the
     terminator that follows the annotation in the reader's and the Database's
     buffers alike, which is also why the view's length is not passed on. */
  auto const * const value =
    vsearch::attribute_value(header, annotation, attributes.size).data();
  auto const abundance = std::strtoll(value, &next_character, decimal_base);
  auto const range_error = (errno == ERANGE);

  if (range_error) {
    fatal("Invalid (range error) abundance annotation in FASTA file header");
  }

  if (abundance == 0) {
    fatal("Invalid (zero) abundance annotation in FASTA file header");
  }

  return abundance;
}


auto annotation_separator(bool & trailing_separator) -> char const * {
  /*
    Return the separator to prepend to the next annotation. When the text
    printed so far already ends with the separator ';' (e.g. a header or a
    label suffix ending with ';'), reuse it rather than emit a second one,
    so that annotations are merged with a single ';' (see issue #271).
  */
  if (trailing_separator) {
    trailing_separator = false;
    return "";
  }
  return ";";
}


auto header_fprint_strip(std::FILE * output_handle,
                         View<char> const header_view,
                         bool const strip_size,
                         bool const strip_ee,
                         bool const strip_length) -> bool
{
  /* the attributes found, ordered by position in the header by the sort below;
     one array of spans, where two parallel start/end arrays used to be kept
     side by side */
  std::array<Attribute_span, n_expected_attributes> found {{}};
  auto nth_attribute = std::size_t{0};

  auto collect = [&](bool const wanted, Attribute const & attribute) -> void {
    if (not wanted) { return; }
    auto const span = header_find_attribute(header_view, attribute);
    if (span.present) {
      found[nth_attribute] = span;
      ++nth_attribute;
    }
  };

  /* look for size attribute */
  collect(strip_size, attributes.size);

  /* look for ee attribute */
  collect(strip_ee, attributes.ee);

  /* look for length attribute */
  collect(strip_length, attributes.length);

  /* sort */

  /* by position in the header; the attribute names differ, so no two spans
     share a start and the order among equals never comes up */

  /* the whole array is sorted, not just its [0, nth_attribute) prefix, and the
     absent slots are ordered last so that the prefix is unchanged. Sorting a
     runtime-length sub-range of a 3-element array makes GCC inline the
     >_S_threshold (16) arm of std::__final_insertion_sort and then warn that
     found[16] is out of bounds (-Warray-bounds), even though nth_attribute is
     at most 3 by construction. A compile-time constant length folds that dead
     arm away. */
  std::sort(found.begin(), found.end(),
            [](Attribute_span const & lhs, Attribute_span const & rhs) -> bool
            {
              if (lhs.present != rhs.present) { return lhs.present; }
              return lhs.start < rhs.start;
            });

  /* print */

  auto last_emitted = View<char>{};  // the last chunk printed, empty until one is

  auto emit = [output_handle, &last_emitted](View<char> const chunk) -> void {
    if (chunk.empty()) { return; }
    fprint(output_handle, chunk);
    last_emitted = chunk;
  };

  if (nth_attribute == 0)
    {
      emit(header_view);
    }
  else
    {
      auto prev_end = std::size_t{0};
      for (std::size_t i = 0; i < nth_attribute; ++i)
        {
          /* print part of header in front of this attribute */
          if (found[i].start > prev_end + 1)
            {
              emit(header_view.subspan(prev_end, found[i].start - prev_end - 1));
            }
          prev_end = found[i].end;
        }

      /* print the rest, if any */
      if (header_view.size() > prev_end + 1)
        {
          emit(header_view.drop(prev_end));
        }
    }

  /* report whether the last emitted character is the annotation separator */
  return (not last_emitted.empty()) and (last_emitted.back() == ';');
}


auto fprint_header_annotations(std::FILE * const output_handle,
                               View<char> const sequence,
                               View<char> const header,
                               OutputAnnotations const & annotations,
                               struct Parameters const & parameters) -> void
{
  // track whether the text printed so far ends with the annotation
  // separator ';', so that appended annotations are merged with a single
  // separator instead of producing ";;" (see issue #271)
  auto trailing_separator = false;

  if (parameters.opt_relabel_self)
    {
      /* normalize first? */
      fprint(output_handle, sequence);
    }
  else if (parameters.opt_relabel_sha1)
    {
      fprint_seq_digest_sha1(output_handle, sequence);
    }
  else if (parameters.opt_relabel_md5)
    {
      fprint_seq_digest_md5(output_handle, sequence);
    }
  else if ((parameters.opt_relabel != nullptr) and (annotations.ordinal > 0))
    {
      std::fputs(parameters.opt_relabel, output_handle);
      fprint_integer(output_handle, annotations.ordinal);
    }
  else
    {
      bool const strip_size = parameters.opt_xsize or (parameters.opt_sizeout and (annotations.abundance > 0));
      bool const strip_ee = parameters.opt_xee or ((parameters.opt_eeout or parameters.opt_fastq_eeout) and (annotations.expected_error >= 0.0));
      bool const strip_length = parameters.opt_xlength or parameters.opt_lengthout;
      trailing_separator = header_fprint_strip(output_handle,
                                               header,
                                               strip_size,
                                               strip_ee,
                                               strip_length);
    }

  if (parameters.opt_label_suffix != nullptr)
    {
      std::fputs(parameters.opt_label_suffix, output_handle);
      if (*parameters.opt_label_suffix != '\0')
        {
          trailing_separator = (parameters.opt_label_suffix[std::strlen(parameters.opt_label_suffix) - 1] == ';');
        }
    }

  if (parameters.opt_sample != nullptr)
    {
      std::fputs(annotation_separator(trailing_separator), output_handle);
      fprint(output_handle, "sample=");
      std::fputs(parameters.opt_sample, output_handle);
    }

  if (annotations.clustersize > 0)
    {
      std::fputs(annotation_separator(trailing_separator), output_handle);
      fprint(output_handle, "seqs=");
      fprint_integer(output_handle, annotations.clustersize);
    }

  if (annotations.clusterid >= 0)
    {
      std::fputs(annotation_separator(trailing_separator), output_handle);
      fprint(output_handle, "clusterid=");
      fprint_integer(output_handle, annotations.clusterid);
    }

  if (parameters.opt_sizeout and (annotations.abundance > 0))
    {
      std::fputs(annotation_separator(trailing_separator), output_handle);
      fprint(output_handle, "size=");
      fprint_integer(output_handle, annotations.abundance);
    }

  if (parameters.opt_centroid_sizeout and (annotations.centroid_size > 0))
    {
      std::fputs(annotation_separator(trailing_separator), output_handle);
      fprint(output_handle, "centroid_size=");
      fprint_integer(output_handle, annotations.centroid_size);
    }

  if ((parameters.opt_eeout or parameters.opt_fastq_eeout) and (annotations.expected_error >= 0.0))
    {
      auto const expected_error = annotations.expected_error;
      auto const * separator = annotation_separator(trailing_separator);
      if (expected_error < 0.000000001) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.13lf", expected_error);
      } else if (expected_error < 0.00000001) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.12lf", expected_error);
      } else if (expected_error < 0.0000001) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.11lf", expected_error);
      } else if (expected_error < 0.000001) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.10lf", expected_error);
      } else if (expected_error < 0.00001) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.9lf", expected_error);
      } else if (expected_error < 0.0001) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.8lf", expected_error);
      } else if (expected_error < 0.001) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.7lf", expected_error);
      } else if (expected_error < 0.01) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.6lf", expected_error);
      } else if (expected_error < 0.1) {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.5lf", expected_error);
      } else {
        std::fputs(separator, output_handle);
        fprint(output_handle, "ee=");
        std::fprintf(output_handle, "%.4lf", expected_error);
      }
    }

  if (parameters.opt_lengthout)
    {
      /* widened by assignment, not by a cast: std::size_t already is
         uint64_t on a 64-bit target, where a cast would be flagged useless */
      uint64_t const sequence_length = sequence.size();
      std::fputs(annotation_separator(trailing_separator), output_handle);
      fprint(output_handle, "length=");
      fprint_integer(output_handle, sequence_length);
    }

  if (annotations.score_name != nullptr)
    {
      std::fputs(annotation_separator(trailing_separator), output_handle);
      std::fputs(annotations.score_name, output_handle);
      fprint(output_handle, '=');
      std::fprintf(output_handle, "%.4lf", annotations.score);
    }

  if (parameters.opt_relabel_keep and
      (((parameters.opt_relabel != nullptr) and (annotations.ordinal > 0)) or parameters.opt_relabel_sha1 or parameters.opt_relabel_md5 or parameters.opt_relabel_self))
    {
      fprint(output_handle, ' ');
      fprint(output_handle, header);
    }
}
