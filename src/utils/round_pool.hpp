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

#include "utils/round_gate.hpp"  // RoundGate
#include <cassert>  // assert
#include <cstddef>  // std::size_t
#include <cstdint>  // std::uint8_t
#include <functional>  // std::function
#include <thread>  // std::thread
#include <vector>


/* Persistent worker threads that run rounds of work for one driving thread,
   through a RoundGate (utils/round_gate.hpp). Each call to run() is one
   round: every worker calls the task once, with its own index, and run()
   returns when all of them are done. The threads are created once and
   joined by the destructor, so a caller that runs many short rounds (the
   library batch functions, called once per batch, or the clustering pool,
   once per round) pays for thread creation only once. */

/* Whether the driving thread takes a share of a round's work. The library
   batch functions keep it out: in a library session fatal() throws on the
   thread that built the session (and only there), and an exception leaving
   the task in the middle of a round would leave the workers running a round
   whose data is being unwound. Their workers call std::exit on a fatal(),
   as they did before the pool. */
enum struct Participation : std::uint8_t { caller_waits, caller_joins };


class RoundPool
{
public:
  /* not noexcept: std::thread's constructor throws std::system_error if a
     thread cannot be created, and the vector may throw std::bad_alloc; in
     the CLI (no exception handling) either terminates the program */
  explicit RoundPool(std::size_t const worker_count) : gate_(worker_count)
  {
    threads_.reserve(worker_count);
    for (std::size_t worker = 0; worker < worker_count; ++worker)
      {
        threads_.emplace_back([this, worker]() -> void {
          unsigned long seen = 0;
          while (gate_.wait_open(seen))
            {
              (*task_)(worker);
              gate_.finish();
            }
        });
      }
  }

  ~RoundPool()
  {
    gate_.close();
    for (auto & thread : threads_)
      {
        thread.join();
      }
  }

  RoundPool(RoundPool const &) = delete;
  RoundPool(RoundPool &&) = delete;
  auto operator=(RoundPool const &) -> RoundPool & = delete;
  auto operator=(RoundPool &&) -> RoundPool & = delete;

  auto worker_count() const noexcept -> std::size_t { return threads_.size(); }

  /* one round: each worker runs task(its index, 0 .. worker_count() - 1);
     with Participation::caller_joins the calling thread also runs
     task(worker_count()). Not noexcept: with caller_joins, the task runs on
     the calling thread and may throw (in the CLI, fatal() exits instead). */
  auto run(std::function<void(std::size_t)> const & task,
           Participation const participation) -> void
  {
    /* with nobody to run it, a round would end before it started */
    assert(participation == Participation::caller_joins or worker_count() > 0);
    /* rounds do not nest: run() is called by one thread, and not from a task */
    assert(task_ == nullptr);
    task_ = &task;  // published to the workers by open_round()
    gate_.open_round();
    if (participation == Participation::caller_joins)
      {
        task(worker_count());
      }
    gate_.wait_round();
    task_ = nullptr;
  }

private:
  RoundGate gate_;
  std::function<void(std::size_t)> const * task_ = nullptr;
  std::vector<std::thread> threads_;  // last: their lambda reads the members above
};
