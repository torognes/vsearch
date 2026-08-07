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

#include <cstdint>  // int64_t, uint64_t
#include <limits>


/* The run statistics the dereplication engines accumulate, and the reporting
   of them. Shared because all three engines -- core/derep.cpp,
   commands/derep_prefix.cpp and commands/derep_smallmem.cpp -- print the same
   messages, in the same wording, and had each spelled them out separately.
   The engines themselves are standalone and stay that way; this is only their
   common vocabulary for saying what they did. */

struct Derep_stats
{
  uint64_t sequencecount = 0;
  uint64_t nucleotidecount = 0;
  int64_t shortest = std::numeric_limits<int64_t>::max();
  int64_t longest = 0;
  uint64_t discarded_short = 0;
  uint64_t discarded_long = 0;
  uint64_t clusters = 0;
  int64_t sumsize = 0;
  uint64_t maxsize = 0;
};


// statistics / summary reporting: each helper folds the former
// stderr-then-log duplicate blocks. The message text is written verbatim to
// stderr and, when --log is in effect, to the log; the log copy is followed
// by the extra blank line the original code emitted (the "\n\n" endings).

auto report_input_stats(Derep_stats const & stats,
                        struct Parameters const & parameters) -> void;

auto report_length_filtered(struct Parameters const & parameters,
                            char const * option_name,
                            int64_t const length_limit,
                            uint64_t const discarded) -> void;

auto report_unique_summary(Derep_stats const & stats,
                           double const average,
                           double const median,
                           struct Parameters const & parameters) -> void;

auto report_selected(uint64_t const selected,
                     Derep_stats const & stats,
                     struct Parameters const & parameters) -> void;
