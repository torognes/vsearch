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

#include <atomic>  // std::atomic
#include <chrono>  // std::chrono::steady_clock, std::chrono::microseconds
#include <condition_variable>  // std::condition_variable
#include <cstddef>  // std::size_t
#include <mutex>  // std::mutex, std::lock_guard, std::unique_lock
#include <thread>  // std::this_thread::yield


/* Rounds of parallel work driven by one thread: the main thread opens a
   round, the workers (and usually the main thread itself) do its work,
   and the main thread waits until every worker is done before it goes on
   alone (committing results, extending an index, ...). Used by the denovo
   chimera batches (core/chimera.cpp) and the clustering pool
   (core/cluster.cpp), whose rounds range from one to eight queries per
   thread, sized adaptively by BatchSizer (utils/batch_sizer.hpp).

   Both sides spin, yielding, for a while before they sleep on a condition
   variable. A thread woken from sleep starts a few hundred microseconds
   late (measured on an 8 P-core + 16 E-core machine: 229 us on average
   after a chimera commit phase, and the last clustering worker of a round
   420 us late at 24 threads), against 100 to 700 us for a whole query. A
   round driven by condition variables alone therefore lasted two to five
   queries instead of one. The spin is bounded, so that workers do not burn
   cores through a long serial phase.
   // C++20 refactoring: std::barrier and std::atomic::wait could replace this */


/* how long RoundGate waits by spinning before it sleeps (C++11: a
   namespace-scope constant, as a static member bound to a reference would
   need an out-of-class definition) */
constexpr std::chrono::microseconds round_gate_spin {2000};


class RoundGate
{
public:
  explicit RoundGate(std::size_t const worker_count) noexcept : worker_count_(worker_count) {}

  /* main thread: start a round; everything written before is visible to
     the workers that see it */
  auto open_round() -> void
  {
    running_.store(worker_count_);
    {
      std::lock_guard<std::mutex> const lock(mutex_);
      ++generation_;
    }
    workers_cv_.notify_all();
  }

  /* main thread: wait until every worker has finished the round */
  auto wait_round() -> void
  {
    if (spin_until([this]() -> bool { return running_.load() == 0; }))
      {
        return;
      }
    std::unique_lock<std::mutex> lock(mutex_);
    main_cv_.wait(lock, [this]() -> bool { return running_.load() == 0; });
  }

  /* main thread: let the workers leave */
  auto close() -> void
  {
    {
      std::lock_guard<std::mutex> const lock(mutex_);
      closed_ = true;
      ++generation_;
    }
    workers_cv_.notify_all();
  }

  /* worker: wait for the round after `seen`; false once the gate is closed */
  auto wait_open(unsigned long & seen) -> bool
  {
    auto const opened = [this, seen]() -> bool { return generation_.load() != seen; };
    if (not spin_until(opened))
      {
        std::unique_lock<std::mutex> lock(mutex_);
        workers_cv_.wait(lock, opened);
      }
    seen = generation_.load();
    std::lock_guard<std::mutex> const lock(mutex_);
    return not closed_;
  }

  /* worker: done with this round */
  auto finish() -> void
  {
    if (running_.fetch_sub(1) == 1)
      {
        std::lock_guard<std::mutex> const lock(mutex_);
        main_cv_.notify_one();
      }
  }

private:
  template <typename Condition>
  static auto spin_until(Condition const & condition) -> bool
  {
    auto const deadline = std::chrono::steady_clock::now() + round_gate_spin;
    while (not condition())
      {
        if (std::chrono::steady_clock::now() > deadline)
          {
            return false;
          }
        std::this_thread::yield();
      }
    return true;
  }

  std::size_t const worker_count_;
  std::atomic<std::size_t> running_ {0};
  std::atomic<unsigned long> generation_ {0};
  bool closed_ = false;  /* guarded by mutex_ */
  std::mutex mutex_;
  std::condition_variable workers_cv_;
  std::condition_variable main_cv_;
};
