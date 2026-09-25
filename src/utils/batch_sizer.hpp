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

#include <algorithm>  // std::max
#include <cstddef>  // std::size_t


/* Chooses the batch size of chimera_denovo_batches while it runs.

   The batch size trades two losses against each other:
   - idle threads: a round lasts as long as its slowest query, so the
     threads that finish early wait. That loss comes from the round's tail,
     about the same for every round, so it shrinks as 1/size;
   - serial checks and recomputations: each query is checked against the
     earlier queries of its batch, and a query found stale is detected
     again, serially, while the other threads wait. Both the number of
     checks and the stale share grow about linearly with the batch size,
     and so does that loss.
   A sum of the form a/size + b x size is smallest where its two terms are
   equal, so the size grows while the measured idle time exceeds the
   measured checking and recomputation loss, and shrinks in the opposite
   case (with 10 %
   hysteresis). Both are measured on the same batches, so the comparison
   does not suffer from the per-query cost growing along the run.

   Measured with 8, 16 and 24 threads (2026-09-25), the best fixed size was
   one query per thread on 61k 16S V4 amplicons, and two to four per thread
   on 219k 18S V9 ones. Sizes move by half the thread count, between one
   and eight queries per thread. The output does not depend on the size. */

/* BatchSizer's bounds and window (namespace-scope: C++11 static members
   bound to a reference, as by std::max, would need an out-of-class
   definition) */
constexpr unsigned int batch_sizer_largest_per_thread = 8;
constexpr unsigned int batch_sizer_window = 16;
constexpr double batch_sizer_hysteresis = 1.1;
/* Only half of the idle time a larger batch recovers turns into
   throughput: the rest goes to contention, as more queries running at
   once slow each other down (measured: 13 % slower per query with 8
   cores busy, 24 % with 16). */
constexpr double batch_sizer_idle_weight = 0.5;


/* The measurements of one batch that BatchSizer weighs, in seconds. An
   aggregate (no default member initializers, which C++11 does not allow
   in one), built from its three values at once. */
struct batch_times_s {
  double round_wall;  /* wall time of the parallel round */
  double busy;  /* time the threads spent detecting in it, summed */
  double serial_wall;  /* wall time of the serial checks and recomputations */
};


class BatchSizer
{
public:
  explicit BatchSizer(std::size_t const threads) noexcept
    : threads_(static_cast<double>(threads)),
      minimum_(static_cast<unsigned int>(threads)),
      maximum_(static_cast<unsigned int>(batch_sizer_largest_per_thread * threads)),
      step_(std::max(1U, static_cast<unsigned int>(threads / 2))),
      size_(static_cast<unsigned int>(threads)) {}

  auto size() const noexcept -> unsigned int { return size_; }
  auto largest() const noexcept -> unsigned int { return maximum_; }

  /* after each batch */
  auto record(struct batch_times_s const & times) noexcept -> void
  {
    idle_ += batch_sizer_idle_weight * std::max(0.0, (threads_ * times.round_wall) - times.busy);
    recompute_ += threads_ * times.serial_wall;
    ++batches_;
    if (batches_ < batch_sizer_window)
      {
        return;
      }
    if ((idle_ > batch_sizer_hysteresis * recompute_) and (size_ + step_ <= maximum_))
      {
        size_ += step_;
      }
    else if ((recompute_ > batch_sizer_hysteresis * idle_) and (size_ >= minimum_ + step_))
      {
        size_ -= step_;
      }
    batches_ = 0;
    idle_ = 0.0;
    recompute_ = 0.0;
  }

private:
  double const threads_;
  unsigned int const minimum_;
  unsigned int const maximum_;
  unsigned int const step_;
  unsigned int size_;
  unsigned int batches_ = 0;
  double idle_ = 0.0;  /* thread-seconds */
  double recompute_ = 0.0;  /* thread-seconds */
};
