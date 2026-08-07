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


#include "utils/taxonomic_fields.h"  // tax_levels
#include <array>


struct Database;

/* Where one taxonomic rank's name sits inside a database header: the offset of
   its first byte and its length, so that db.header_view(seqno).subspan(start,
   length) is the name. A zero length means the header does not carry that rank.

   Offsets rather than a View<char> into the header, deliberately.  Every
   consumer does build that view, and commands/sintax.cpp builds it immediately
   -- but results_show_lcaout keeps its candidate ranks across iterations and
   re-derives them from a different seqno each round, so a view stored in round
   one would be compared against a header fetched in round three. The offsets
   are the safe currency inside tax_split; the views are the right thing for the
   caller to build right afterwards. */
struct TaxLevel {
  int start = 0;
  int length = 0;
};

/* Fills `levels` with the ranks the header of sequence `seqno` carries, leaving
   the ranks it does not carry as they were -- so callers pass a freshly
   default-constructed array. One array of a two-field record rather than the
   two parallel int arrays this used to take as separate out-parameters: a
   wrong-size buffer now fails to compile instead of running off the end, and
   the two halves can no longer be filled to different depths. */
auto tax_split(int seqno, std::array<TaxLevel, tax_levels> & levels,
               struct Database const & db) -> void;
