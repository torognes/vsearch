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

#include <cstdint>  // std::uint8_t

/* Internal seam between the CLI chimera-detection engine (core/chimera.cpp)
   and its five thin command wrappers (commands/uchime_denovo.cpp,
   commands/uchime2_denovo.cpp, commands/uchime3_denovo.cpp,
   commands/uchime_ref.cpp, commands/chimeras_denovo.cpp). Not part of the
   public library API, so it lives here rather than in core/chimera.hpp.

   The five commands all run this one engine; it selects reference-vs-de-novo
   detection and the algorithm variant from the mode. */

/* Which of the five chimera commands the shared engine runs as. They differ in
   where the parent candidates come from (a separate --db for uchime_ref, the
   query set itself for the four de-novo modes), in the scoring rule, in the
   output options they honour, and in what the log records. chimera() used to
   tell them apart by asking which opt_<command> pointer was non-null; a
   library caller sets none of them, so the caller states its mode instead --
   the same reasoning as QualityOrigin (parameters.hpp) and Derep_mode
   (derep_internal.hpp). */
enum struct ChimeraMode : std::uint8_t {
  uchime_ref, uchime_denovo, uchime2_denovo, uchime3_denovo, chimeras_denovo
};

/* De-novo detection: the four modes that draw their candidate parents from the
   sequences already processed, rather than from the separate --db reference
   that --uchime_ref searches. */
inline auto chimera_is_denovo(ChimeraMode const mode) -> bool {
  return mode != ChimeraMode::uchime_ref;
}

auto chimera(ChimeraMode mode, struct Parameters const & parameters) -> void;
