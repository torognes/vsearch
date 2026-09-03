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

#include "parameters.hpp"  // validate_thread_count, parameters_resolve_derived, parameters_validate, vsearch_apply_defaults_fixups, QualityOrigin, fastq_output_offset, resolve_fastq_qmaxout
#include "vsearch.hpp"  // struct Parameters
#include "core/chimera.hpp"  // maxparents
#include "core/searchcore.hpp"  // minwordmatches_defaults
#include "os/system.hpp"  // system_get_available_cores
#include "utils/fatal.hpp"  // fatal
#include "utils/quality_encoding.hpp"  // lowest_printable_ascii, highest_printable_ascii, sanger_ascii_offset, legacy_max_quality, fastq_qmaxout_unset
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
      fatal(message);
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

  if (not parameters.runtime.gap_penalties_adjusted)
    {
      parameters.opt_gap_open_query_left -= parameters.opt_gap_extension_query_left;
      parameters.opt_gap_open_target_left -= parameters.opt_gap_extension_target_left;
      parameters.opt_gap_open_query_interior -= parameters.opt_gap_extension_query_interior;
      parameters.opt_gap_open_target_interior -= parameters.opt_gap_extension_target_interior;
      parameters.opt_gap_open_query_right -= parameters.opt_gap_extension_query_right;
      parameters.opt_gap_open_target_right -= parameters.opt_gap_extension_target_right;
      parameters.runtime.gap_penalties_adjusted = true;
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


auto fastq_output_offset(struct Parameters const & parameters,
                         QualityOrigin const origin) -> int64_t
{
  /* --fastq_mergepairs writes with the input offset because it has no
     --fastq_asciiout to write with; every other producer of quality uses the
     output offset. See QualityOrigin in parameters.hpp for why this cannot be
     read off the opt_<command> pointers. */
  return (origin == QualityOrigin::merged) ? parameters.opt_fastq_ascii
                                           : parameters.opt_fastq_asciiout;
}


auto resolve_fastq_qmaxout(struct Parameters const & parameters,
                           QualityOrigin const origin) -> int64_t
{
  if (parameters.opt_fastq_qmaxout != fastq_qmaxout_unset)
    {
      return parameters.opt_fastq_qmaxout;
    }
  /* a score vsearch produced keeps the pre-3.0 ceiling; one that passed
     through is capped only by what the output offset can carry */
  if (origin == QualityOrigin::passed_through)
    {
      return highest_printable_ascii - fastq_output_offset(parameters, origin);
    }
  return legacy_max_quality;
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
      fatal(message);
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
      /* The cores this run may use, not the machine's: a job confined to one
         core by a cgroup would otherwise default to one thread per core the
         host has (issue #584). See os/system.hpp. */
      parameters.opt_threads = system_get_available_cores();
    }
  if (parameters.opt_maxrejects == -1)
    {
      parameters.opt_maxrejects = 32;
    }
  if (parameters.opt_wordlength == 0)
    {
      parameters.opt_wordlength = 8;
    }

  /* The input quality ceiling follows the ASCII offset, the way cli.cc's
     resolve_quality_bound_defaults() does it for the command line. Without
     this, a library caller who sets opt_fastq_ascii = 64 keeps the Sanger
     ceiling of 93, and 64 + 93 = 157 is the very combination the CLI's sum
     rule refuses (torognes/vsearch#564). The struct's own comment says
     "cli.cc lowers the unset default accordingly"; the library path simply
     never got the same treatment.

     "Still the member-initializer value" is the test for "not set by the
     caller", as it already is for opt_maxrejects and opt_wordlength above --
     the library has no options_selected vector to consult. A caller who
     deliberately wrote 93 next to offset 64 asked for the impossible
     combination anyway, so resolving it is a fix rather than a loss.

     No symbol changes verdict either way: 126 - 64 = 62 is exactly the
     highest value a printable byte can carry at that offset, so nothing was
     being wrongly accepted before. What this corrects is the value the
     caller reads back, and the divergence from the CLI.

     opt_fastq_qmaxout is still NOT resolved here, for the same reason as
     before: its default is entry-point-specific (see QualityOrigin), and this
     function is command-agnostic by construction. It now carries a sentinel
     rather than a plausible-looking 93, so the consumer that does know its
     origin -- MergePairs -- resolves it, and the two cannot be confused. */
  if (parameters.opt_fastq_qmax == highest_printable_ascii - sanger_ascii_offset)
    {
      parameters.opt_fastq_qmax =
        highest_printable_ascii - parameters.opt_fastq_ascii;
    }

  /* The three output-quality rules the CLI applies in
     validate_option_values() but the library had nowhere. Without the ceiling
     one, a caller who sets opt_fastq_ascii = 64 and leaves the ceiling at 93
     gets merged quality bytes up to 64 + 93 = 157 -- an invalid FASTQ file,
     and a SIGSEGV in vsearch 2.31.0.

     They are checked here rather than in parameters_validate(), which the CLI
     also calls, because the offset a symbol is written with depends on the
     command: a legal CLI run (--fastq_convert --fastq_ascii 64
     --fastq_asciiout 33 --fastq_qmaxout 93) writes with asciiout and would be
     rejected by any command-agnostic reading of the rule. A library session
     has no command, and MergePairs is its only consumer of these bounds, so
     the merged origin is the one that applies -- which also means the offset
     is opt_fastq_ascii, and the messages name it.

     The ceiling has to be the RESOLVED one: the member still holds the
     sentinel unless the caller set it, and comparing qminout against
     INT64_MIN would refuse every default configuration. */
  auto const merged_ceiling =
    resolve_fastq_qmaxout(parameters, QualityOrigin::merged);
  auto const merged_offset =
    fastq_output_offset(parameters, QualityOrigin::merged);

  if (parameters.opt_fastq_qminout > merged_ceiling)
    {
      fatal("The argument to --fastq_qminout cannot be larger than --fastq_qmaxout");
    }

  if (merged_offset + parameters.opt_fastq_qminout < lowest_printable_ascii)
    {
      fatal("Sum of arguments to --fastq_ascii and --fastq_qminout must be no less than 33");
    }

  if (merged_offset + merged_ceiling > highest_printable_ascii)
    {
      fatal("Sum of arguments to --fastq_ascii and --fastq_qmaxout must be no more than 126");
    }

  parameters_validate(parameters);
}
