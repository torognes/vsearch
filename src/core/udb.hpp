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

#include "utils/span.hpp"  // Span
#include <cstdint>  // uint64_t

/* What the caller wants out of a UDB file. This replaced two adjacent bool
   parameters, create_bitmaps and parse_abundances, that every one of the six
   callers set together: true for the four commands that go on to search the
   index, false for the two that only report on the file. Two same-typed
   neighbours that must never disagree are exactly the shape a caller can swap
   by mistake, so the one fact they carried between them is now one value. */
/* What the caller intends to do with the file, which decides how much of the
   k-mer index is kept. The k-mer sections are 72 % to 99 % of a UDB file, so
   the difference is not a detail: a word length 13 database costs 1.07 GB of
   dense tables to load for search and 21 MB to convert to fasta. Every value
   reads and validates the whole file either way -- what varies is what
   survives the read. */
enum struct UdbUse : unsigned char {
  search,     /* the whole index, the bitmaps, and the ;size= annotations */
  word_stats, /* the word counts, for reporting on them; no bitmaps, no
                 abundances, and no per-k-mer offset table */
  sequences,  /* no k-mer data at all: the sequences and their headers */
};

auto udb_detect_isudb(char const * filename) -> bool;
auto udb_read(char const * filename,
              UdbUse usage,
              struct Dbindex & dbindex,
              struct Database & db,
              struct Parameters const & parameters) -> void;

/* Fill `entries` from the word list, starting at entry number `first`.

   For a UdbUse::word_stats caller, which does not hold that section: the whole
   word list is 4 bytes per index entry and can be a gigabyte, while the report
   shows at most eight entries for each of eleven k-mers. So it is fetched on
   demand rather than kept, by reopening the file and seeking -- the same way
   udb_detect_isudb() and --udbinfo already read a UDB independently of
   udb_read().

   Every value is checked against `seqcount` here, not only when the file was
   loaded, so a file rewritten in between cannot put an out-of-range sequence
   number into the report. The section's position comes from the index
   udb_read() just filled, not from the session's configuration -- so the word
   length used to find it is the file's own. */
auto udb_read_word_entries(char const * filename,
                           struct Dbindex const & dbindex,
                           unsigned int seqcount,
                           uint64_t first,
                           Span<unsigned int> entries) -> void;
