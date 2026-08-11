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

#pragma once

#include "utils/view.hpp"  // View<char>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE


auto header_get_size(View<char> header) -> int64_t;

auto annotation_separator(bool & trailing_separator) -> char const *;

auto header_fprint_strip(std::FILE * output_handle,
                         View<char> header,
                         bool strip_size,
                         bool strip_ee,
                         bool strip_length) -> bool;


/* The values a print helper may append to a record's header, as opposed to the
   record's own stored fields. Deliberately separate from SeqRecord for the
   reason SeqRecord's own comment gives about abundance: these are all output
   decisions (a match count, a summed cluster size, a per-run ordinal, a
   detection score), not properties of the record that was read. The defaults
   are the "not applicable" sentinels, so a caller names only what it means. */
struct OutputAnnotations
{
  uint64_t abundance = 0;
  int64_t  ordinal = 0;
  double   expected_error = -1.0;
  int64_t  clustersize = -1;
  int      clusterid = -1;
  char const * score_name = nullptr;
  double   score = 0.0;
  uint64_t centroid_size = 0;

  /* The default member initializers above cost the C++11 aggregate-init that
     SeqRecord keeps, so the two-value form most of the call sites want is a
     constructor instead. Everything else default-constructs and assigns the
     field it means by name. */
  OutputAnnotations() = default;
  OutputAnnotations(uint64_t const new_abundance, int64_t const new_ordinal)
    : abundance {new_abundance}, ordinal {new_ordinal} {}
};


/* Emit the annotated header for one record: everything between the '>' or '@'
   the caller has already written and the newline before the payload. Returns
   nothing -- the two trailers differ and stay with their own printer.

   'sequence' is a parameter even though the header does not print it: the
   --relabel_self label and the --relabel_sha1/--relabel_md5 digests are
   computed from it, and --lengthout reports its length. */
auto fprint_header_annotations(std::FILE * output_handle,
                               View<char> sequence,
                               View<char> header,
                               OutputAnnotations const & annotations,
                               struct Parameters const & parameters) -> void;
