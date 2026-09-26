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

#include "core/batch_session.hpp"
#include "utils/make_unique.hpp"  // make_unique
#include "utils/fatal.hpp"  // fatal_detail::ExitOnFatal
#include "utils/round_pool.hpp"  // RoundPool
#include <cassert>  // assert
#include <cstddef>  // std::size_t
#include <functional>  // std::function


auto batch_session_s::pool(std::size_t const thread_count) -> RoundPool &
{
  assert(thread_count >= 1);
  auto const worker_count = thread_count - 1;
  if (pool_ == nullptr or pool_->worker_count() != worker_count)
    {
      pool_.reset();  // join the old threads before creating new ones
      pool_ = make_unique<RoundPool>(worker_count);
    }
  return *pool_;
}


namespace batch_session_detail {
  auto current() noexcept -> batch_session_s *&
  {
    thread_local batch_session_s * session = nullptr;
    return session;
  }
}


auto current_batch_session() noexcept -> batch_session_s *
{
  return batch_session_detail::current();
}


auto run_batch_round(std::size_t const thread_count,
                     std::function<void(std::size_t)> const & task) -> void
{
  assert(thread_count >= 1);
  auto const share = [&task](std::size_t const thread) -> void {
    fatal_detail::ExitOnFatal const exit_on_fatal;
    task(thread);
  };
  auto * const session = current_batch_session();
  if (session == nullptr)
    {
      RoundPool pool(thread_count - 1);
      pool.run(share);
      return;
    }
  session->pool(thread_count).run(share);
}
