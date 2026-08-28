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
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf


namespace vsearch {

/* One writer for the "how many queries matched" report, in one place.

   commands/usearch_global.cpp and commands/search_exact.cpp print the same
   two lines and had grown a copy each -- character for character identical,
   differing only in the type of the state object they read the counters from.
   Neither file owns a header the other includes, hence this one.

   The report is:

     Matching unique query sequences: <n> of <n> (<pct>%)
     Matching total query sequences: <n> of <n> (<pct>%)

   the second line only under --sizein, and each percentage only when its
   denominator is non-zero.

   The counters arrive in a struct rather than as four adjacent integers,
   which no signature could keep un-swappable, and the two commands build it
   from their own state at the call site. That is deliberately more explicit
   than templating on the state type: a template would compile today only
   because both structs happen to name their members the same way, and would
   say so nowhere. */
struct MatchCounts
{
  int unique_matched;
  int unique_total;
  uint64_t abundance_matched;
  uint64_t abundance_total;
};


/* The percentage of the first line divides an int by an int through 100.0,
   the second divides two uint64_t after an explicit cast to double. Both are
   copied from the call sites they replace rather than unified, because they
   are what the two commands have always printed. */
inline auto print_match_counts(std::FILE * output_stream,
                               MatchCounts const & counts,
                               bool const report_abundances) -> void
{
  fprint(output_stream, "Matching unique query sequences: ");
  fprint_integer(output_stream, counts.unique_matched);
  fprint(output_stream, " of ");
  fprint_integer(output_stream, counts.unique_total);
  if (counts.unique_total > 0)
    {
      fprint(output_stream, " (");
      std::fprintf(output_stream, "%.2f", 100.0 * counts.unique_matched / counts.unique_total);
      fprint(output_stream, "%)");
    }
  fprint(output_stream, '\n');
  if (report_abundances)
    {
      fprint(output_stream, "Matching total query sequences: ");
      fprint_integer(output_stream, counts.abundance_matched);
      fprint(output_stream, " of ");
      fprint_integer(output_stream, counts.abundance_total);
      if (counts.abundance_total > 0)
        {
          fprint(output_stream, " (");
          std::fprintf(output_stream, "%.2f", 100.0 * static_cast<double>(counts.abundance_matched) / static_cast<double>(counts.abundance_total));
          fprint(output_stream, "%)");
        }
      fprint(output_stream, '\n');
    }
}

}  // namespace vsearch
