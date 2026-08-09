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

#include "print_view.hpp"  // fprint
#include "vsearch.hpp"  // struct Parameters

#include "view.hpp"  // View, make_view
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf
#include <string>  // std::string
#include <utility>  // std::move


class Progress {
public:
  /* The prompt is owned, not borrowed: it is printed again on every update and
     once more from the destructor, so a Progress built from a std::string used
     to require that string to outlive it -- a requirement no signature stated
     and four callers met only by accident of scope. A Progress is built once
     per phase, never per record, so the copy is unmeasurable; the ~60 callers
     passing a string literal are unchanged in source. */
  explicit Progress(std::string prompt, std::uint64_t const max_size,
                    struct Parameters const & parameters)
      : prompt_{std::move(prompt)},
        max_size_{max_size},
        stderr_is_tty_(parameters.opt_stderr_is_tty),
        is_quiet_(parameters.opt_quiet),
        no_progress_(parameters.opt_no_progress) {
    is_visible_ = check_if_visible();
    if (is_quiet_) { return; }
    fprint(stderr, make_view(prompt_));
    if (not is_visible_) { return; }
    fprint(stderr, " 0%");
    if (max_size_ == 0) {
      fprint(stderr, "  \r");
      fprint(stderr, make_view(prompt_));
      fprint(stderr, " 0%");
      return;
    }
    current_percentage_ = calculate_percentage();
    next_threshold_ = calculate_next_threshold();
  }

  Progress(Progress const &) = delete;  // copy constructor: no copies
  auto operator=(Progress const &) -> Progress & = delete;  // assignment operator: no self-assignment
  Progress(Progress&&) = delete;  // no move constructor
  auto operator=(Progress&&) -> Progress & = delete; // no move assignment operator
  
  auto update(std::uint64_t const counter) -> void {
    counter_ = counter;
    if ((not is_visible_) or (counter_ < next_threshold_)) { return; }
    current_percentage_ = calculate_percentage();
    fprint(stderr, "  \r");
    fprint(stderr, make_view(prompt_));
    fprint(stderr, ' ');
    fprint_integer(stderr, current_percentage_);
    fprint(stderr, '%');
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
  static constexpr auto one_hundred_percent = 100UL;

  // Construction-time parameters
  std::string prompt_;
  std::uint64_t max_size_ {};
  bool stderr_is_tty_ {};
  bool is_quiet_ {};
  bool no_progress_ {};

  // Internal parameters
  std::uint64_t counter_ {};
  std::uint64_t current_percentage_ {};  // integer, and not a double
  std::uint64_t next_threshold_ {};
  bool is_visible_ {};
  
  // Helper functions
  auto check_if_visible() const -> bool {
    return stderr_is_tty_
      and (not is_quiet_)
      and (not no_progress_);
  };

  auto calculate_percentage() const -> std::uint64_t {
    if (max_size_ == 0) { return 0; }  // when reading from a pipe
    return counter_ * one_hundred_percent / max_size_;
  };

  auto calculate_next_threshold() const -> std::uint64_t {
    static constexpr auto nighty_nine_percent = 99UL;
    return (((current_percentage_ + 1) * max_size_) + nighty_nine_percent) /
      one_hundred_percent;
  };

  auto done() const -> void {
    if (is_quiet_) { return; }
    if (is_visible_) {
      fprint(stderr, "  \r");
      fprint(stderr, make_view(prompt_));
    }
    fprint(stderr, ' ');
    fprint_integer(stderr, one_hundred_percent);
    fprint(stderr, "%\n");
  }
};
