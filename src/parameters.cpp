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

// Shared parameter normalization and validation: resolves the sentinel and
// range values of a Parameters struct (vsearch_apply_defaults_fixups). This
// logic is common to both the command-line front-end (cli.cc) and the library
// API (the VsearchSession constructor calls vsearch_apply_defaults_fixups), so
// it lives in its own translation unit rather than in the argument parser or in
// main(). validate_thread_count() lives here too, next to its fixups caller.

#include "parameters.hpp"  // validate_thread_count, parameters_resolve_derived, parameters_validate, vsearch_apply_defaults_fixups
#include "vsearch.hpp"  // struct Parameters
#include "core/chimera.hpp"  // maxparents
#include "core/searchcore.hpp"  // minwordmatches_defaults
#include "os/system.hpp"  // system_get_cores
#include "utils/fatal.hpp"  // fatal
#include <array>  // std::array::size
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t
#include <limits>  // std::numeric_limits
#include <string>  // std::to_string, std::string


auto validate_thread_count(int64_t const threads) -> void
{
  // Upper bound for the --threads option. Kept local to this validator
  // (rather than as a Parameters field) so the option's range limit lives
  // with the check that enforces it instead of in the configuration struct.
  constexpr int64_t n_threads_max = 1024;
  if ((threads < 0) or (threads > n_threads_max))
    {
      std::string const message =
        "The argument to --threads must be in the range 0 (default) to "
        + std::to_string(n_threads_max);
      fatal(message.c_str());
    }
}


/* Sentinel/range resolution operating on a Parameters struct: resolves the
   values the parser (or a library caller) left as sentinels, applies the
   default thread count, and range-checks a few options. Shared by the CLI
   (apply_command_defaults) and the library session (the VsearchSession
   constructor). The gap-open adjustment is guarded by the struct's own
   gap_penalties_adjusted so a repeated call stays idempotent. */
auto parameters_resolve_derived(struct Parameters & parameters) -> void
{
  if (parameters.opt_maxhits == 0)
    {
      parameters.opt_maxhits = std::numeric_limits<int64_t>::max();
    }

  if (parameters.opt_minwordmatches < 0)
    {
      if (parameters.opt_wordlength >= 0 and
          parameters.opt_wordlength < static_cast<int64_t>(minwordmatches_defaults.size()))
        {
          parameters.opt_minwordmatches = minwordmatches_defaults[static_cast<size_t>(parameters.opt_wordlength)];
        }
      else
        {
          parameters.opt_minwordmatches = 0;
        }
    }

  if (not parameters.gap_penalties_adjusted)
    {
      parameters.opt_gap_open_query_left -= parameters.opt_gap_extension_query_left;
      parameters.opt_gap_open_target_left -= parameters.opt_gap_extension_target_left;
      parameters.opt_gap_open_query_interior -= parameters.opt_gap_extension_query_interior;
      parameters.opt_gap_open_target_interior -= parameters.opt_gap_extension_target_interior;
      parameters.opt_gap_open_query_right -= parameters.opt_gap_extension_query_right;
      parameters.opt_gap_open_target_right -= parameters.opt_gap_extension_target_right;
      parameters.gap_penalties_adjusted = true;
    }

  /* Fold the twelve '*' gap-penalty sentinels into one flag so the accept path
     can skip the cigar scan on the common finite-penalty runs. Idempotent. */
  parameters.opt_gap_penalty_has_infinite =
    parameters.opt_gap_open_query_left_infinite or
    parameters.opt_gap_open_query_interior_infinite or
    parameters.opt_gap_open_query_right_infinite or
    parameters.opt_gap_open_target_left_infinite or
    parameters.opt_gap_open_target_interior_infinite or
    parameters.opt_gap_open_target_right_infinite or
    parameters.opt_gap_extension_query_left_infinite or
    parameters.opt_gap_extension_query_interior_infinite or
    parameters.opt_gap_extension_query_right_infinite or
    parameters.opt_gap_extension_target_left_infinite or
    parameters.opt_gap_extension_target_interior_infinite or
    parameters.opt_gap_extension_target_right_infinite;
}


auto parameters_validate(struct Parameters const & parameters) -> void
{
  validate_thread_count(parameters.opt_threads);

  if (parameters.opt_maxaccepts < 0)
    {
      fatal("The argument to --maxaccepts must not be negative");
    }
  if (parameters.opt_maxrejects < 0)
    {
      fatal("The argument to --maxrejects must not be negative");
    }

  if ((parameters.opt_wordlength < 3) or (parameters.opt_wordlength > 15))
    {
      fatal("The argument to --wordlength must be in the range 3 to 15");
    }

  if ((parameters.opt_iddef < 0) or (parameters.opt_iddef > 4))
    {
      fatal("The argument to --iddef must be in the range 0 to 4");
    }

  if ((parameters.opt_chimeras_parents_max < 2) or
      (parameters.opt_chimeras_parents_max > maxparents))
    {
      std::string const message =
        "The argument to --chimeras_parents_max must be in the range 2 to "
        + std::to_string(maxparents);
      fatal(message.c_str());
    }
}


auto vsearch_apply_defaults_fixups(struct Parameters & parameters) -> void
{
  parameters_resolve_derived(parameters);

  /* Command-agnostic sentinel defaults. The CLI overrides these with
     command-aware values (validate_option_values / configure_threads) and then
     calls parameters_resolve_derived() and parameters_validate() directly, so
     it never runs this block; the library reaches it through this umbrella. */
  if (parameters.opt_id >= 0.0 and parameters.opt_weak_id > parameters.opt_id)
    {
      parameters.opt_weak_id = parameters.opt_id;
    }
  if (parameters.opt_threads == 0)
    {
      parameters.opt_threads = system_get_cores();
    }
  if (parameters.opt_maxrejects == -1)
    {
      parameters.opt_maxrejects = 32;
    }
  if (parameters.opt_wordlength == 0)
    {
      parameters.opt_wordlength = 8;
    }

  parameters_validate(parameters);
}
