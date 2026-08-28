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

#include "utils/print_view.hpp"  // fprint, fprint_integer
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fputs


namespace vsearch {

/* One wording for "this option discarded that many sequences", in one place.

   Four options report it -- --minseqlength, --maxseqlength and --minsize from
   Database::read() (core/db.cpp), and --minseqlength/--maxseqlength again from
   the derep commands, which use their own reader (core/derep_stats.cpp). The
   two files had grown an identical writer each, character for character apart
   from the parameter names and the width of the count.

   The sentence is:

     <option> <threshold>: <n> sequence(s) discarded.

   and it is a warning, so callers write it to stderr with no --quiet gate:
   --quiet suppresses "messages to stdout and stderr, except for warnings and
   error messages" (man/commands/fragments/option_quiet.md). The caller keeps
   the choice of destinations and, where one is wanted, the log copy's
   trailing separator.

   Named for the sentence rather than for the criterion: the derep copy was
   called print_length_filtered(), but --minsize compares an abundance, not a
   length.

   Distinct parameter types throughout, so no two arguments can be swapped
   silently: the option name is a pointer, the threshold is a signed int64_t
   (it is an option value, and --minseqlength defaults to -1), and the count
   is an unsigned uint64_t. */
inline auto print_discarded(std::FILE * output_stream,
                            char const * const option_name,
                            int64_t const threshold,
                            uint64_t const discarded) -> void
{
  std::fputs(option_name, output_stream);
  fprint(output_stream, ' ');
  fprint_integer(output_stream, threshold);
  fprint(output_stream, ": ");
  fprint_integer(output_stream, discarded);
  fprint(output_stream, ' ');
  std::fputs((discarded == 1 ? "sequence" : "sequences"), output_stream);
  fprint(output_stream, " discarded.\n");
}

}  // namespace vsearch
