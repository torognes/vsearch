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

#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"  // fasta_print_general
#include "core/fastq.hpp"  // fastq_print_general
#include "utils/open_file.hpp"  // OutputFileHandle
#include "utils/view.hpp"  // View<char>


struct Parameters;
struct SeqRecord;


namespace vsearch {

/* The two destinations one record may be written to. Nine commands offer a
   --fasta*out and a --fastq*out for the same category of records (kept and
   discarded, merged and not merged, matched and not matched), and each such
   pair is one OutputPair. A command with several categories holds several.

   Either handle may be null: the user asked for neither format, or only one.
   The handles are what the emit helpers gate on, rather than the matching
   parameters.opt_* name, so a destination stays a single value that carries
   its own "is this wanted" answer. That also covers the pairs whose two
   option names share no stem, such as getseq's --notmatched/--notmatchedfq. */
struct OutputPair
{
  OutputFileHandle fasta;
  OutputFileHandle fastq;

  /* Will write_record() write anything at all? It tests both handles itself,
     so this is not needed to be correct -- it is for a caller whose argument
     views cost something to build, because as arguments they are built before
     write_record() can decline. One test at the call site instead of the two
     the migrated if/if pair used. Measured on fastq_mergepairs: building the
     six views of a not-merged pair that nobody asked for cost 42 instructions
     a record.

     No user-provided constructor and no default member initializers, so
     OutputPair stays a C++11 aggregate and callers can still brace-init it. */
  auto wanted() const -> bool {
    return (fasta != nullptr) or (fastq != nullptr);
  }
};


/* Emit one record to whichever of the pair's two destinations is open.

   FASTA first, then FASTQ. The order is observable and not incidental: two
   output options naming the same file share one std::FILE (see
   utils/open_file.hpp), so `--fastaout out.txt --fastqout out.txt` writes the
   two renderings of each record interleaved in this order. It is the order
   every migrated call site already used.

   inline, and in the header: core/filter.cpp and commands/fastq_mergepairs.cpp
   call this once per record on paths measured at 1.9-3.0x, where an
   out-of-line call would cost about 21 instructions a record. */
inline auto write_record(OutputPair const & destination,
                         View<char> const sequence,
                         View<char> const header,
                         View<char> const quality,
                         OutputAnnotations const & annotations,
                         struct Parameters const & parameters) -> void
{
  if (destination.fasta != nullptr)
    {
      fasta_print_general(destination.fasta.get(), sequence, header,
                          annotations, parameters);
    }

  if (destination.fastq != nullptr)
    {
      fastq_print_general(destination.fastq.get(), sequence, header, quality,
                          annotations, parameters);
    }
}


/* As above, for a caller that holds the whole record (a reader's record() or
   Database::record()) rather than the separate views.

   The quality view of a FASTA record is empty, but no caller can reach the
   FASTQ arm with one: every command offering a --fastq*out rejects FASTA input
   up front, either with its own fatal() or by opening the input with
   fastq_open() instead of fastx_open(). Audited 2026-09-11 across all nine. */
inline auto write_record(OutputPair const & destination,
                         SeqRecord const & record,
                         OutputAnnotations const & annotations,
                         struct Parameters const & parameters) -> void
{
  if (destination.fasta != nullptr)
    {
      fasta_print_general(destination.fasta.get(), record, annotations, parameters);
    }

  if (destination.fastq != nullptr)
    {
      fastq_print_general(destination.fastq.get(), record, annotations, parameters);
    }
}

}  // namespace vsearch
