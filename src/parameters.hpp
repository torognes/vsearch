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
// library-session entry points (vsearch_session_begin/vsearch_session_end),
// also defined in parameters.cpp, are public API and are declared in the
// public header vsearch_api.h instead. struct Parameters is (for now) still
// defined in vsearch.hpp; the functions that take it do so by reference, so a
// forward declaration suffices.

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

/* Fatal unless the thread count and the maxaccepts/maxrejects/wordlength
   values are in range. Shared range validation for both the CLI and the
   library; command-specific validation (e.g. --chimeras_parents_max) lives in
   the CLI's own validators. */
auto parameters_validate(struct Parameters const & parameters) -> void;

/* Apply every default fix-up to a Parameters in the order the compute engines
   expect: parameters_resolve_derived(), then the command-agnostic sentinel
   defaults (weak_id clamp, threads, maxrejects, wordlength), then
   parameters_validate(). The library's single entry point (called by
   vsearch_session_begin()); the CLI supplies its own command-aware defaults
   and calls the two helpers above directly. Idempotent per struct. */
auto vsearch_apply_defaults_fixups(struct Parameters & parameters) -> void;
