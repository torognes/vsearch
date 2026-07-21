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

// Implementation of the non-compute parts of the library API declared in
// vsearch_api.h: the API-version accessors and the session lifecycle (the
// VsearchSession constructor/destructor). The compute entry points declared in
// vsearch_api.h (search, cluster, chimera, ...) are implemented under core/.
// Parameter resolution lives in parameters.cpp; this TU reaches it only through
// vsearch_apply_defaults_fixups() when a session opens.

#include "vsearch_api.h"  // VSEARCH_API_VERSION*, vsearch_api_version*, VsearchSession
#include "parameters.hpp"  // vsearch_apply_defaults_fixups
#include "utils/fatal.hpp"  // fatal_detail::throw_on_fatal


auto vsearch_api_version() -> int
{
  return VSEARCH_API_VERSION;
}

auto vsearch_api_version_string() -> const char *
{
  return VSEARCH_API_VERSION_STRING;
}


/* Reference to the thread_local flag that puts fatal() into throwing mode.
   Defined here, next to its only mutators (the VsearchSession ctor/dtor below),
   rather than in fatal_throw.cpp: fatal_throw.cpp merely reads it, and keeping
   the storage out of that translation unit avoids a cppcheck false positive (it
   would otherwise inline the accessor, see the default-false initial value with
   no local write, and flag the throw branch as dead). Function-local
   thread_local so the default (false) is guaranteed on every thread — worker
   threads keep it false and keep using cooperative abort; only a thread that
   opened a session flips it to true. */
namespace fatal_detail {
  auto throw_on_fatal() -> bool &
  {
    thread_local bool enabled = false;
    return enabled;
  }
}


/* A library session is now just a caller-owned object: no process-wide lock and
   no begin/end pair, because vsearch keeps no shared mutable state to serialize
   (so independent sessions can run concurrently in different threads). The
   constructor puts fatal() into throwing mode on this thread, then resolves the
   struct's sentinels/ranges; enabling throwing first means a configuration
   error surfacing during the fixups is itself a catchable VsearchError, the
   same way in-session fatals are. The previous mode is saved and restored
   (rather than forced back to false) so nested sessions on one thread compose,
   and worker threads — which never construct a session — keep the default,
   non-throwing behaviour (an exception must not escape a std::thread). */
VsearchSession::VsearchSession(struct Parameters & parameters)
  : previous_throw_mode(fatal_detail::throw_on_fatal())
{
  fatal_detail::throw_on_fatal() = true;
  vsearch_apply_defaults_fixups(parameters);
}


VsearchSession::~VsearchSession()
{
  fatal_detail::throw_on_fatal() = previous_throw_mode;
}
