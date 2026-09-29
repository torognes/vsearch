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

#include <map>  // std::map
#include <utility>  // std::move


/* Writing chunks of work in the order they were claimed.

   A command that parses its input inside the claim (under the input lock)
   gives each chunk a rank, its claim order, which is also the input order.
   Chunks are then processed in parallel and finish in any order. Submitted
   with the output lock held, a chunk is written at once if its turn has come,
   followed by every waiting chunk whose turn follows; otherwise it waits here,
   and the worker moves on to claim another. The output is therefore
   identical to a single-threaded run.

   First written in --search_exact (commit 6f8f4e4d), where the claim order is
   the input order; shared since with --search_oligodb.

   Not thread-safe by itself: every call is made under the caller's output
   lock. submit() cannot be noexcept: std::map allocates, and the writer is
   the caller's. */
template <typename Chunk>
class ChunkReorder
{
public:
  /* A chunk that has to wait is moved out of `chunk`, which the caller then
     holds empty (default-constructed) and can reuse for its next claim. */
  template <typename Write>
  auto submit(unsigned long const rank, Chunk & chunk, Write write) -> void
  {
    if (rank != next_rank_)
      {
        waiting_.emplace(rank, std::move(chunk));
        chunk = Chunk{};
        return;
      }
    write(chunk);
    ++next_rank_;
    auto next = waiting_.find(next_rank_);
    while (next != waiting_.end())
      {
        write(next->second);
        waiting_.erase(next);
        ++next_rank_;
        next = waiting_.find(next_rank_);
      }
  }

  /* true once every submitted chunk was written */
  auto empty() const noexcept -> bool { return waiting_.empty(); }

private:
  unsigned long next_rank_ = 0;
  std::map<unsigned long, Chunk> waiting_;
};
