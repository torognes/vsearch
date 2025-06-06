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

#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf


constexpr auto one_hundred_percent = 100UL;


class Progress {
public:
  explicit Progress(char const * prompt, std::uint64_t const max_size,
                    struct Parameters const & parameters)
    : prompt_{prompt}, max_size_{max_size}, parameters_{parameters} {
    assert(prompt != nullptr);
    is_visible_ = check_if_visible();
    if (parameters_.opt_quiet) { return; }
    static_cast<void>(std::fprintf(stderr, "%s", prompt));
    if (not is_visible_) { return; }
    static_cast<void>(std::fprintf(stderr, " %d%%", 0));
    if (max_size_ == 0) {
      static_cast<void>(std::fprintf(stderr, "  \r%s 0%%", prompt_));
      return;
    }
    current_percentage_ = calculate_percentage();
    next_threshold_ = calculate_next_threshold();
  }

  Progress(Progress const &) = delete;  // copy constructor: no copies
  auto operator=(Progress const &) -> Progress & = delete;  // assignment operator: no self-assignment
  Progress(Progress&&) = delete;  // move constructor
  auto operator=(Progress&&) -> Progress & = delete; // move assignment operator
  
  auto update(std::uint64_t const counter) -> void {
    counter_ = counter;
    if ((not is_visible_) or (counter_ < next_threshold_)) { return; }
    current_percentage_ = calculate_percentage();
    static_cast<void>(std::fprintf(stderr,
                                   "  \r%s %" PRIu64 "%%",
                                   prompt_,
                                   current_percentage_));
    next_threshold_ = calculate_next_threshold();    
  };

  auto update() -> void {
    ++counter_;
    update(counter_);
  }

  ~Progress() {
    done();
  }
  
private:
  char const * prompt_ {};
  std::uint64_t max_size_ {};
  struct Parameters const & parameters_ {};

  // Internal
  std::uint64_t counter_ {};
  std::uint64_t current_percentage_ {};  // integer, and not a double
  std::uint64_t next_threshold_ {};
  bool is_visible_ {};
  
  // Helpers
  auto check_if_visible() const -> bool {
    return (parameters_.opt_stderr_is_tty)
      and (not parameters_.opt_quiet)
      and (not parameters_.opt_no_progress);
  };

  auto calculate_percentage() const -> std::uint64_t {
    assert(max_size_ != 0);
    return counter_ * one_hundred_percent / max_size_;
  };

  auto calculate_next_threshold() const -> std::uint64_t {
    static constexpr auto nighty_nine_percent = 99UL;
    return ((current_percentage_ + 1) * max_size_ + nighty_nine_percent) /
      one_hundred_percent;
  };

  auto done() const -> void {
    if (parameters_.opt_quiet) { return; }
    if (is_visible_) {
      static_cast<void>(std::fprintf(stderr, "  \r%s", prompt_));
    }
    static_cast<void>(std::fprintf(stderr, " %lu%%\n", one_hundred_percent));
  }
};
