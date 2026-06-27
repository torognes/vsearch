/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>


/*
  ThreadRunner manages a fixed pool of persistent worker threads. The
  threads are created once at construction and parked in stand-by mode.
  Each call to run() wakes all workers, lets each execute the supplied
  function once (receiving its own thread index), and blocks until all
  have finished. The threads are joined and destroyed at destruction.

  This replaces the hand-rolled pthread worker pools previously found
  throughout vsearch (a pthread_t plus a pthread_mutex_t and a
  pthread_cond_t per thread). RAII removes all the manual init/destroy
  bookkeeping: std::mutex, std::condition_variable and std::thread clean
  up after themselves.
*/

class ThreadRunner {
private:

  enum struct Work_state : std::int8_t { wait = 0, work = 1, quit = -1 };

  struct thread_s {
    uint64_t thread_id {0};
    std::function<void(uint64_t)> fun;
    std::thread thread;
    std::mutex workmutex;
    std::condition_variable workcond;
    Work_state work {Work_state::wait};
  };

  std::vector<struct thread_s> thread_array;

  static auto worker(struct thread_s & tip) -> void {
    std::unique_lock<std::mutex> lock(tip.workmutex);

    /* loop until signalled to quit */
    while (tip.work != Work_state::quit) {
      /* wait for work available */
      if (tip.work == Work_state::wait) {
        tip.workcond.wait(lock);
      }

      if (tip.work == Work_state::work) {
        tip.fun(tip.thread_id);
        tip.work = Work_state::wait;
        tip.workcond.notify_one();
      }
    }
  }


public:

  // The constructor cannot be marked noexcept: std::thread's
  // constructor throws std::system_error if a thread cannot be created,
  // and std::vector allocation may throw std::bad_alloc. In practice
  // thread creation almost never fails; if it does, the exception
  // propagates and terminates the (non-exception-handling) program,
  // which matches the previous fatal() behaviour on pthread_create
  // failure.
  ThreadRunner(std::size_t const thread_count,
               std::function<void(uint64_t nth_thread)> const & function) :
      thread_array(thread_count) {
    /* init and create worker threads */
    // std::ref is required: std::thread decays its arguments by
    // default, so passing `tip` directly would copy thread_s — which
    // contains non-copyable members (std::mutex, std::condition_variable,
    // std::thread) and would fail to compile. std::ref preserves the
    // reference semantics so worker() receives the live thread_s.
    uint64_t counter {0};
    for (auto & tip : thread_array) {
      tip.thread_id = counter;
      tip.fun = function;
      tip.thread = std::thread(worker, std::ref(tip));
      ++counter;
    }
  }


  ~ThreadRunner() {
    /* ask threads to quit and wait for them to join */
    for (auto & tip : thread_array) {
      /* tell worker to quit */
      {
        std::lock_guard<std::mutex> const lock(tip.workmutex);
        tip.work = Work_state::quit;
        tip.workcond.notify_one();
      }
      /* wait for worker to quit */
      tip.thread.join();
    }
  }


  ThreadRunner(ThreadRunner const &) = delete; // copy constructor
  ThreadRunner(ThreadRunner &&) = delete; // move constructor
  auto operator=(ThreadRunner const &) -> ThreadRunner & = delete; // copy assignment
  auto operator=(ThreadRunner &&) -> ThreadRunner & = delete; // move assignment


  auto run() -> void {
    /* wake up threads */
    for (auto & tip : thread_array) {
      std::lock_guard<std::mutex> const lock(tip.workmutex);
      tip.work = Work_state::work;
      tip.workcond.notify_one();
    }

    /* wait for threads to finish their work */
    for (auto & tip : thread_array) {
      std::unique_lock<std::mutex> lock(tip.workmutex);
      tip.workcond.wait(lock, [&tip]() -> bool { return tip.work != Work_state::work; });
    }
  }
};
