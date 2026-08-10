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

#include <string>  // std::string


/* The process-wide warning reporter, the non-fatal counterpart of fatal()
   (utils/fatal.hpp): it emits the message and returns, where fatal() exits or
   throws. Both write the same two places, so a warning reaches the --log file
   whenever one is open -- which is what the fifteen hand-written warning sites
   this replaced only did at six of them.

   In namespace vsearch, unlike fatal(): "warn" is a far more common identifier
   than "fatal", the header is reachable from library consumers, and an
   unqualified call in consumer code could otherwise be captured silently by
   ADL. A new header costs nothing to be born namespaced (no call sites to
   convert later), so it is the first file to adopt the decision recorded in
   TBD_20260725_vsearch_namespace.md. */
namespace vsearch
{

  /* Emits "\nWARNING: " + message + "\n" to stderr, and the same to the --log
     file when one is open (log_file::handle(), utils/logfile.hpp).

     The leading newline is unconditional: it separates the warning from
     whatever preceded it, typically a progress line still missing its own
     terminator, which is what two of the original sites spelled out by hand and
     the other thirteen wanted too.

     Never consults --quiet: warnings survive it by documented contract
     ("Suppress messages to stdout and stderr, except for warnings and error
     messages", man/commands/fragments/option_quiet.md), which is also why this
     needs no Parameters. fatal() has never consulted it either.

     The caller keeps any interior newlines: exactly one "WARNING: " prefix and
     exactly one trailing newline are added, so a multi-line message such as
     "... does not support multithreading.\nOnly 1 thread used." arrives
     unchanged, and a companion line that is *not* a warning (fastx.cpp's
     "REMINDER: ...") stays a separate emission rather than gaining the prefix.

     A string literal selects the char const * overload, not the std::string
     one: array-to-pointer is an exact match and the std::string conversion is
     user-defined. The pair mirrors fatal.hpp's, for the same reason -- the
     literal sites stay literal and only the composed ones allocate. */
  auto warn(char const * message) -> void;
  auto warn(std::string const & message) -> void;

}  // namespace vsearch
