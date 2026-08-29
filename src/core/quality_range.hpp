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
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
  OF THE POSSIBILITY OF SUCH DAMAGE.

*/



#pragma once

#include "utils/quality_encoding.hpp"  // QualitySymbolRange
#include <cstdint>  // int64_t, uint64_t
#include <string>  // std::string

struct Parameters;


namespace vsearch {

/* One wording for every quality score outside --fastq_qmin/--fastq_qmax, in
   one place.

   Six sites test that window while reading FASTQ, and four of them --
   core/filter.cpp, core/eestats.cpp, core/mergepairs.cpp and core/derep.cpp
   -- carried character-for-character identical text in four copies. The two
   that differ, commands/fastq_stats.cpp and commands/fastq_convert.cpp, are
   left alone here: migrating them changes what a user sees, which is a
   separate decision (see DONE_20260825_quality_range.md).

   Detection is separable from reporting because core/mergepairs.cpp tests
   the window on a worker thread and cannot fatal() in place: it records a
   verdict and the main thread turns it into text after the join (see the
   deferred-error note in fastx.hpp). classify_quality() is that test,
   quality_out_of_range_message() is that text, and check_quality_score() is
   the two together for the callers that can stop where they stand. */

enum struct QualityBound { in_range, below_qmin, above_qmax };


/* Where in the input the offending score was found, when the caller can say.

   Four of the six sites hold the reader handle and can (filter, derep,
   fastq_stats, fastq_convert); core/mergepairs.cpp cannot, because it tests
   the window on a worker thread where the only ordinal to hand is a pair
   number spanning two files, with no line number at all. So the location is
   optional, and a default-constructed QualityLocation means "not known" --
   the message then reads exactly as it did before. record is 1-based, the
   way commands/fastq_convert.cpp has always printed it.

   A plain aggregate, deliberately: default member initializers would make it
   a non-aggregate before C++14, and QualityLocation{} is how a caller says
   "not known". Build it with fastx_s::quality_location() rather than by hand,
   so the two uint64_t members cannot be filled in the wrong order. */
struct QualityLocation {
  uint64_t record;
  uint64_t line;

  auto known() const noexcept -> bool { return record != 0; }
};


/* Pure test, no output: safe to call from a worker thread. */
auto classify_quality(int64_t quality_score,
                      struct Parameters const & parameters) noexcept -> QualityBound;

/* The message a verdict deserves. Never called with QualityBound::in_range;
   the caller has already decided there is something to report. */
auto quality_out_of_range_message(QualityBound bound,
                                  int64_t quality_score,
                                  struct Parameters const & parameters,
                                  QualityLocation location = QualityLocation{}) -> std::string;

/* Test, and fatal() with that message if the score is outside the window.
   Not noexcept: fatal() throws VsearchError in a library session. */
auto check_quality_score(int64_t quality_score,
                         struct Parameters const & parameters,
                         QualityLocation location = QualityLocation{}) -> void;

/* The same check for a caller that has a whole record's symbol range rather
   than one score: decodes both extremes with --fastq_ascii and checks them,
   lowest first, which is the order the per-symbol callers would have reached
   them in. An empty range (FASTA record, empty quality line) checks nothing. */
auto check_quality_range(QualitySymbolRange range,
                         struct Parameters const & parameters) -> void;

}  // namespace vsearch
