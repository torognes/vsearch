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

#include "core/bitmap.hpp"
#include "utils/view.hpp"
#include <cstdint>  // uint64_t
#include <vector>


struct Database;
struct Parameters;

struct dbhash_bucket_s
{
  uint64_t hash = 0;
  uint64_t seqno = 0;
};

struct dbhash_search_info_s
{
  View<char> seq;
  uint64_t hash = 0;
  uint64_t index = 0;
};


/* The exact-match dedup hash index over a Database (used by --search_exact).
   Owns its bucket table and open-addressing bitmap (RAII: released by clear()
   and the destructor) and is non-copyable/non-movable. The search API
   (search_first/search_next) is const, so the worker threads can query one
   shared index concurrently; add()/add_all() build it single-threaded first. */
struct Dbhash
{
private:
  Bitmap bitmap_;  /* one occupancy bit per bucket slot */
  uint64_t mask_ = 0;  /* table size is a power of two; index = hash & mask_ */
  std::vector<struct dbhash_bucket_s> table_;

public:
  Dbhash() = default;
  ~Dbhash();
  Dbhash(Dbhash const &) = delete;
  auto operator=(Dbhash const &) -> Dbhash & = delete;
  Dbhash(Dbhash &&) = delete;
  auto operator=(Dbhash &&) -> Dbhash & = delete;

  auto open(uint64_t maxelements) -> void;
  auto clear() -> void;

  auto add(View<char> seq, uint64_t seqno, struct Database const & db) -> void;
  auto add_all(struct Database const & db, struct Parameters const & parameters) -> void;

  auto search_first(View<char> seq,
                    struct dbhash_search_info_s & info,
                    struct Database const & db) const -> int64_t;
  auto search_next(struct dbhash_search_info_s & info, struct Database const & db) const -> int64_t;

private:
  /* Shared probe of the two entry points above, which held a character-identical
     copy of it each. Walk forward from 'index' to the first slot that is either
     free or holds info.seq, record where the scan stopped, and report the
     sequence found there (-1 for a free slot). Both callers have already put the
     query's hash and sequence in info -- search_first computes them,
     search_next carries them over from the previous call -- so the only thing
     that differs is where the walk starts, and that is the argument: the hash's
     own slot, or one past the last stop. */
  auto probe_from(uint64_t index,
                  struct dbhash_search_info_s & info,
                  struct Database const & db) const -> int64_t;
};
