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

// Declarations for the internal functions defined in parameters.cpp:
// thread-count validation and Parameters sentinel/range resolution. The
// library session lifecycle (the VsearchSession ctor/dtor, which calls
// vsearch_apply_defaults_fixups below) is public API defined in vsearch_api.cpp
// and declared in the public header vsearch_api.h instead. struct Parameters is
// (for now) still defined in vsearch.hpp; the functions that take it do so by
// reference, so a forward declaration suffices.

#include <cstdint>  // int64_t (validate_thread_count)

struct Parameters;

// Fatal unless the requested thread count is within the accepted range
// (see the upper bound local to validate_thread_count()).
auto validate_thread_count(int64_t threads) -> void;

/* Resolve the values derived purely from other options, with no command-
   specific variant: the maxhits/minwordmatches sentinels and the once-only
   gap-open adjustment (plus its folded infinite-penalty flag). Idempotent per
   struct. Both the CLI and the library call this. */
auto parameters_resolve_derived(struct Parameters & parameters) -> void;

/* Where a written quality symbol comes from. Two independent facts hang off
   it, which is why one enum carries both: the ceiling that applies to the
   symbol, and the ASCII offset that encodes it.

     merged          --fastq_mergepairs. The Edgar & Flyvbjerg posterior of a
                     merged base is computed rather than read, so it keeps the
                     pre-3.0 ceiling of 41; and it is written with
                     --fastq_ascii, because the command does not accept
                     --fastq_asciiout.
     generated       --fasta2fastq. Also a produced score, so also 41, but
                     written with --fastq_asciiout like everything else.
     passed_through  a symbol read from the input and re-encoded
                     (--fastq_convert, --fastx_uniques, --sff_convert). Capped
                     only by what the output offset can carry.

   The two axes really are independent -- 'generated' is a produced score
   written with asciiout, 'merged' a produced score written with ascii -- which
   is why neither can be derived from the other. --fastq_mergepairs is the only
   command on the ascii side: the others that reject --fastq_asciiout
   (--fastx_filter, --cut, --orient, --fastx_revcomp, --fastx_subsample) copy
   quality through verbatim and never clamp it, and --derep_fulllength has no
   --fastqout at all. */
enum struct QualityOrigin { merged, generated, passed_through };

/* The ASCII offset a written quality symbol carries. Never derive this from
   the opt_<command> pointers: a library session sets none of them, so a merge
   driven through MergePairs would read as a pass-through and be encoded with
   the wrong offset. The caller states its origin instead. */
auto fastq_output_offset(struct Parameters const & parameters,
                         QualityOrigin origin) -> int64_t;

/* The resolved --fastq_qmaxout ceiling: the caller's value whenever one was
   given, otherwise the default for this origin. The CLI resolves once at parse
   time and writes the answer back into the struct, so by the time a command
   runs the member holds a real value and a direct read is correct. The library
   cannot -- vsearch_apply_defaults_fixups() is command-agnostic by
   construction and does not know which entry point is coming -- so it leaves
   the sentinel in place and every library-reachable consumer calls this. */
auto resolve_fastq_qmaxout(struct Parameters const & parameters,
                           QualityOrigin origin) -> int64_t;

/* Fatal unless the thread count and the
   maxaccepts/maxrejects/wordlength/iddef/chimeras_parents_max values are in
   range. Shared range validation for both the CLI and the library, so a
   library caller is held to the same bounds the CLI enforces. */
auto parameters_validate(struct Parameters const & parameters) -> void;

/* Apply every default fix-up to a Parameters in the order the compute engines
   expect: parameters_resolve_derived(), then the command-agnostic sentinel
   defaults (weak_id clamp, threads, maxrejects, wordlength), then
   parameters_validate(). The library's single entry point (called by the
   VsearchSession constructor); the CLI supplies its own command-aware defaults
   and calls the two helpers above directly. Idempotent per struct. */
auto vsearch_apply_defaults_fixups(struct Parameters & parameters) -> void;
