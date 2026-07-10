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

#include <mutex>


/*
  run_worker_loop drives one worker thread's self-scheduling loop, the
  pattern shared by every ThreadRunner-based command in vsearch: repeatedly
  claim the next unit of work under a shared input mutex, then process that
  unit with the mutex released so the other workers can claim in parallel.

  claim is called while input_mutex is held. It should claim the next unit
  of work (for a counter-driven source: read and advance the counter; for a
  reader-driven source: pull the next record and copy it into thread-local
  storage before another worker overwrites the shared reader buffers) and
  return true, or return false when the input is exhausted. Returning false
  ends the loop.

  work is called with no lock held. It performs the computation for the unit
  the preceding claim call selected, communicating through state the two
  callables share (typically captured by reference). Any output
  synchronisation is work's own responsibility -- several callers need none
  because they write to disjoint per-unit slots.

  Centralising the loop makes "release the input mutex before doing the
  work" a property of one place instead of an invariant each caller
  re-establishes by hand-placing input_lock.unlock().

  This cannot be marked noexcept: locking the mutex may throw
  std::system_error and the supplied callables are unconstrained. This
  matches ThreadRunner, whose run() is likewise not noexcept.
*/
template <typename ClaimFn, typename WorkFn>
auto run_worker_loop(std::mutex & input_mutex,
                     ClaimFn claim,
                     WorkFn work) -> void {
  while (true) {
    {
      std::lock_guard<std::mutex> const lock(input_mutex);
      if (not claim()) {
        break;
      }
    }
    work();
  }
}
