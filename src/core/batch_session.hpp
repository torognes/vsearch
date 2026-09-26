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

#include "utils/round_pool.hpp"  // RoundPool
#include <cstddef>  // std::size_t
#include <functional>  // std::function
#include <memory>  // std::unique_ptr


/* What a library session keeps between calls to the batch functions
   (search_batch, chimera_detect_batch, cluster_assign_batch): one pool of
   worker threads, shared by the three, and the per-thread working state of
   search_batch and chimera_detect_batch. Both are built at the first batch
   call that needs them and kept until the VsearchSession is destroyed, so a
   caller that submits many small batches no longer pays for thread creation
   and per-thread initialization on every call (0.1 to 2.1 ms per call,
   measured at 8 to 24 threads on 2026-09-25). The per-thread state of
   cluster_assign_batch belongs to its cluster_session_s, which only borrows
   the threads.

   The batch functions find the session through a thread_local pointer to
   the innermost VsearchSession open on the calling thread. A batch call
   from a thread without a session creates its threads and state for that
   call only, as all batch calls did before. So two threads cannot share
   one session's pool: batch calls on one session come from the thread that
   opened it, one at a time. */

/* Base of the per-function state kept by a batch_session_s: the session
   destroys it without knowing which function built it. */
struct batch_state_s
{
  batch_state_s() = default;
  virtual ~batch_state_s() = default;
  batch_state_s(batch_state_s const &) = delete;
  batch_state_s(batch_state_s &&) = delete;
  auto operator=(batch_state_s const &) -> batch_state_s & = delete;
  auto operator=(batch_state_s &&) -> batch_state_s & = delete;
};


struct batch_session_s
{
  /* each is built, and replaced when its inputs change, by its function;
     the state never reads the database or the index when destroyed, so
     the session may outlive them */
  std::unique_ptr<batch_state_s> search;   // core/search.cpp
  std::unique_ptr<batch_state_s> chimera;  // core/chimera.cpp

  /* the session's worker threads; created at the first call, and again
     when a call asks for a different number (opt_threads changed). Not
     noexcept: creating threads may throw (see RoundPool). */
  auto pool(std::size_t worker_count) -> RoundPool &;

private:
  /* last: destroyed first, joining idle threads before the states go */
  std::unique_ptr<RoundPool> pool_;
};


/* The batch session of the innermost VsearchSession open on this thread,
   nullptr when none is. */
auto current_batch_session() noexcept -> batch_session_s *;

/* One round of worker_count tasks, task(0) to task(worker_count - 1), each on
   its own thread, the calling thread waiting: on the session's pool, or on
   threads created for this call outside a session. Not noexcept: creating
   threads may throw (see RoundPool). */
auto run_batch_round(std::size_t worker_count,
                     std::function<void(std::size_t)> const & task) -> void;


namespace batch_session_detail {
  /* the pointer behind current_batch_session(), set and restored by the
     VsearchSession constructor and destructor (vsearch_api.cpp) */
  auto current() noexcept -> batch_session_s *&;
}
