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

#include "vsearch.h"
#include "vsearch_api.h"
#include "commands/allpairs_global.hpp"
#include "core/chimera.hpp"  // maxparents
#include "core/searchcore.hpp"  // minwordmatches_defaults
#include "commands/uchime_denovo.hpp"
#include "commands/uchime2_denovo.hpp"
#include "commands/uchime3_denovo.hpp"
#include "commands/uchime_ref.hpp"
#include "commands/chimeras_denovo.hpp"
#include "cli.h"
#include "commands/cluster_fast.hpp"
#include "commands/cluster_smallmem.hpp"
#include "commands/cluster_size.hpp"
#include "commands/cluster_unoise.hpp"
#include "commands/cut.hpp"
#include "commands/derep_fulllength.hpp"
#include "commands/derep_id.hpp"
#include "commands/fastx_uniques.hpp"
#include "commands/derep_prefix.hpp"
#include "commands/derep_smallmem.hpp"
#include "os/dynlibs.hpp"
#include "commands/fastq_eestats.hpp"
#include "commands/fastq_eestats2.hpp"
#include "commands/fasta2fastq.hpp"
#include "commands/fastq_chars.hpp"
#include "commands/fastq_join.hpp"
#include "commands/fastq_mergepairs.hpp"
#include "commands/fastq_convert.hpp"
#include "commands/fastq_stats.hpp"
#include "commands/fastx_revcomp.hpp"
#include "commands/fastx_subsample.hpp"
#include "commands/fastx_syncpairs.hpp"
#include "commands/fastq_filter.hpp"
#include "commands/fastx_filter.hpp"
#include "commands/fastx_getseq.hpp"
#include "commands/fastx_getseqs.hpp"
#include "commands/fastx_getsubseq.hpp"
#include "commands/help.hpp"
#include "commands/fastx_mask.hpp"
#include "commands/maskfasta.hpp"
#include "commands/orient.hpp"
#include "commands/rereplicate.hpp"
#include "commands/usearch_global.hpp"
#include "commands/search_exact.hpp"
#include "commands/sff_convert.hpp"
#include "commands/shuffle.hpp"
#include "commands/sintax.hpp"
#include "commands/sortbylength.hpp"
#include "commands/sortbysize.hpp"
#include "commands/makeudb_usearch.hpp"
#include "commands/udb2fasta.hpp"
#include "commands/udbinfo.hpp"
#include "commands/udbstats.hpp"
#include "commands/version.hpp"
#include "utils/compare_strings_nocase.hpp"
#ifdef __x86_64__
#include "arch/x86_64/cpu_features.hpp"
#endif
#include "utils/fatal.hpp"
#include "utils/logfile.hpp"  // LogFile
#include "utils/prog_id.hpp"  // PROG_NAME, PROG_VERSION, PROG_ARCH
#include "utils/random.hpp"
#include <algorithm>  // std::count, std::any_of
#include <array>
#include <cerrno>  // errno, ERANGE
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t, std::snprintf, std::printf
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <getopt.h>  // getopt_long_only, optarg, optind, opterr, struct
                     // option (no_argument, required_argument)
#include <limits>
#include <mutex>  // std::mutex
#include <new>  // std::set_new_handler
#include <string>
#include <unistd.h>  // write, _exit, STDERR_FILENO
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // performance: lambda compiles to a single x86-64 instruction (shr al, 7),
  // equivalent to (c & ~0x7F) == 0
  auto is_not_ASCII(std::string const & user_string) -> bool {
    static constexpr auto ascii_max = static_cast<unsigned char>(std::numeric_limits<signed char>::max());
    auto is_not_in_range = [](char const user_char) -> bool {
      auto const unsigned_user_char = static_cast<unsigned char>(user_char);
      return (unsigned_user_char > ascii_max);
    };
    return std::any_of(user_string.begin(), user_string.end(), is_not_in_range);
  }

}  // end of anonymous namespace


/* Serializes library sessions: acquired by vsearch_session_begin() and released
   by vsearch_session_end(), so only one session runs at a time. */
static std::mutex session_mutex;


auto vsearch_api_version() -> int
{
  return VSEARCH_API_VERSION;
}

auto vsearch_api_version_string() -> const char *
{
  return VSEARCH_API_VERSION_STRING;
}



auto vsearch_session_end() -> void
{
  session_mutex.unlock();
}




/* Parameters-based sentinel/range resolution: an exact mirror of the global
   overload above, operating on the struct. Introduced for the E1 migration to
   Parameters-primary configuration (F2); the two are kept in lockstep until
   the globals are removed. The gap-open adjustment is guarded by the struct's
   own gap_penalties_adjusted so a repeated call stays idempotent. */
auto vsearch_apply_defaults_fixups(struct Parameters & parameters) -> void
{
  if (parameters.opt_maxhits == 0)
    {
      parameters.opt_maxhits = Parameters::int64_max;
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

  if (parameters.opt_id >= 0.0 and parameters.opt_weak_id > parameters.opt_id)
    {
      parameters.opt_weak_id = parameters.opt_id;
    }

  if ((parameters.opt_threads < 0) or (parameters.opt_threads > Parameters::n_threads_max))
    {
      fatal("The argument to --threads must be in the range 0 (default) to 1024");
    }
  if (parameters.opt_threads == 0)
    {
      parameters.opt_threads = system_get_cores();
    }

  if (parameters.opt_maxrejects == -1)
    {
      parameters.opt_maxrejects = 32;
    }
  if (parameters.opt_maxaccepts < 0)
    {
      fatal("The argument to --maxaccepts must not be negative");
    }
  if (parameters.opt_maxrejects < 0)
    {
      fatal("The argument to --maxrejects must not be negative");
    }

  if (parameters.opt_wordlength == 0)
    {
      parameters.opt_wordlength = 8;
    }
  if ((parameters.opt_wordlength < 3) or (parameters.opt_wordlength > 15))
    {
      fatal("The argument to --wordlength must be in the range 3 to 15");
    }

  if ((parameters.opt_chimeras_parents_max < 2) or (parameters.opt_chimeras_parents_max > maxparents))
    {
      std::string const message =
        "The argument to --chimeras_parents_max must be in the range 2 to "
        + std::to_string(maxparents);
      fatal(message.c_str());
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


/* Begin a library session from a Parameters (E1/F2, shape A). Acquires the
   session lock (same non-blocking semantics as the retired
   vsearch_init_defaults) and resolves the struct's sentinels/ranges. Pair with
   vsearch_session_end(). */
auto vsearch_session_begin(struct Parameters & parameters) -> void
{
  if (not session_mutex.try_lock())
    {
      fatal("A vsearch library session is already active: a previous "
            "vsearch_session_begin() was not paired with vsearch_session_end(). "
            "Call vsearch_session_end() before starting a new session.");
    }
  vsearch_apply_defaults_fixups(parameters);
}


auto cmd_allpairs_global(struct Parameters const & parameters) -> void
{
  /* check options */

  if ((parameters.opt_alnout == nullptr) and (parameters.opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (parameters.opt_blast6out == nullptr) and
      (parameters.opt_matched == nullptr) and (parameters.opt_notmatched == nullptr) and
      (parameters.opt_samout == nullptr) and (parameters.opt_fastapairs == nullptr))
    {
      fatal("No output files specified");
    }

  if (not ((parameters.opt_acceptall != 0) or ((parameters.opt_id >= 0.0) and (parameters.opt_id <= 1.0))))
    {
      fatal("Specify either --acceptall or --id with an identity from 0.0 to 1.0");
    }

  allpairs_global(parameters);
}


auto cmd_usearch_global(struct Parameters const & parameters) -> void
{
  /* check options */

  if ((parameters.opt_alnout == nullptr) and (parameters.opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (parameters.opt_blast6out == nullptr) and
      (parameters.opt_matched == nullptr) and (parameters.opt_notmatched == nullptr) and
      (parameters.opt_dbmatched == nullptr) and (parameters.opt_dbnotmatched == nullptr) and
      (parameters.opt_samout == nullptr) and (parameters.opt_otutabout == nullptr) and
      (parameters.opt_biomout == nullptr) and (parameters.opt_mothur_shared_out == nullptr) and
      (parameters.opt_fastapairs == nullptr) and (parameters.opt_lcaout == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_db == nullptr)
    {
      fatal("Database filename not specified with --db");
    }

  if ((parameters.opt_id < 0.0) or (parameters.opt_id > 1.0))
    {
      fatal("Identity between 0.0 and 1.0 must be specified with --id");
    }

  usearch_global(parameters);
}


auto cmd_search_exact(struct Parameters const & parameters) -> void
{
  /* check options */

  if ((parameters.opt_alnout == nullptr) and (parameters.opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (parameters.opt_blast6out == nullptr) and
      (parameters.opt_matched == nullptr) and (parameters.opt_notmatched == nullptr) and
      (parameters.opt_dbmatched == nullptr) and (parameters.opt_dbnotmatched == nullptr) and
      (parameters.opt_samout == nullptr) and (parameters.opt_otutabout == nullptr) and
      (parameters.opt_biomout == nullptr) and (parameters.opt_mothur_shared_out == nullptr) and
      (parameters.opt_fastapairs == nullptr) and (parameters.opt_lcaout == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_db == nullptr)
    {
      fatal("Database filename not specified with --db");
    }

  search_exact(parameters);
}


auto cmd_subsample(struct Parameters const & parameters) -> void
{
  if ((parameters.opt_fastaout == nullptr) and (parameters.opt_fastqout == nullptr))
    {
      fatal("Specify output files for subsampling with --fastaout and/or --fastqout");
    }

  if ((parameters.opt_sample_pct > 0) == (parameters.opt_sample_size > 0))
    {
      fatal("Specify either --sample_pct or --sample_size");
    }

  subsample(parameters);
}


auto cmd_none(struct Parameters const & parameters) -> void {
  if (parameters.opt_quiet) { return ; }
  std::fprintf(stderr,
          "For more help, please enter: %s --help\n"
          "For further details, please consult the manual by entering: man vsearch\n"
          "\n"
          "Selected command examples:\n"
          "\n"
          "vsearch --allpairs_global FILENAME --id 0.5 --alnout FILENAME\n"
          "vsearch --cluster_size FILENAME --id 0.97 --centroids FILENAME\n"
          "vsearch --cut FILENAME --cut_pattern G^AATT_C --fastaout FILENAME\n"
          "vsearch --fastq_chars FILENAME\n"
          "vsearch --fastq_convert FILENAME --fastqout FILENAME --fastq_ascii 64\n"
          "vsearch --fastq_eestats FILENAME --output FILENAME\n"
          "vsearch --fastq_eestats2 FILENAME --output FILENAME\n"
          "vsearch --fastq_mergepairs FILENAME --reverse FILENAME --fastqout FILENAME\n"
          "vsearch --fastq_stats FILENAME --log FILENAME\n"
          "vsearch --fastx_filter FILENAME --fastaout FILENAME --fastq_trunclen 100\n"
          "vsearch --fastx_getseq FILENAME --label LABEL --fastaout FILENAME\n"
          "vsearch --fastx_mask FILENAME --fastaout FILENAME\n"
          "vsearch --fastx_revcomp FILENAME --fastqout FILENAME\n"
          "vsearch --fastx_subsample FILENAME --fastaout FILENAME --sample_pct 1\n"
          "vsearch --fastx_uniques FILENAME --fastaout FILENAME\n"
          "vsearch --makeudb_usearch FILENAME --output FILENAME\n"
          "vsearch --search_exact FILENAME --db FILENAME --alnout FILENAME\n"
          "vsearch --sff_convert FILENAME --output FILENAME --sff_clip\n"
          "vsearch --shuffle FILENAME --output FILENAME\n"
          "vsearch --sintax FILENAME --db FILENAME --tabbedout FILENAME\n"
          "vsearch --sortbylength FILENAME --output FILENAME\n"
          "vsearch --sortbysize FILENAME --output FILENAME\n"
          "vsearch --uchime_denovo FILENAME --nonchimeras FILENAME\n"
          "vsearch --uchime_ref FILENAME --db FILENAME --nonchimeras FILENAME\n"
          "vsearch --usearch_global FILENAME --db FILENAME --id 0.97 --alnout FILENAME\n"
          "\n"
          "Other commands: cluster_fast, cluster_smallmem, cluster_unoise, cut,\n"
          "                derep_id, derep_fulllength, derep_prefix, derep_smallmem,\n"
          "                fasta2fastq, fastq_filter, fastq_join, fastx_getseqs,\n"
          "                fastx_getsubseq, fastx_syncpairs, maskfasta, orient, rereplicate,\n"
          "                uchime2_denovo, uchime3_denovo, udb2fasta, udbinfo, udbstats,\n"
          "                version\n"
          "\n",
          parameters.progname);
}


auto cmd_cluster(struct Parameters const & parameters) -> void
{
  if ((parameters.opt_alnout == nullptr) and (parameters.opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (parameters.opt_blast6out == nullptr) and
      (parameters.opt_matched == nullptr) and (parameters.opt_notmatched == nullptr) and
      (parameters.opt_centroids == nullptr) and (parameters.opt_clusters == nullptr) and
      (parameters.opt_consout == nullptr) and (parameters.opt_msaout == nullptr) and
      (parameters.opt_samout == nullptr) and (parameters.opt_profile == nullptr) and
      (parameters.opt_otutabout == nullptr) and (parameters.opt_biomout == nullptr) and
      (parameters.opt_mothur_shared_out == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_cluster_unoise == nullptr)
    {
      if ((parameters.opt_id < 0.0) or (parameters.opt_id > 1.0))
        {
          fatal("Identity between 0.0 and 1.0 must be specified with --id");
        }
    }

  if (parameters.opt_cluster_fast != nullptr)
    {
      cluster_fast(parameters);
    }
  else if (parameters.opt_cluster_smallmem != nullptr)
    {
      cluster_smallmem(parameters);
    }
  else if (parameters.opt_cluster_size != nullptr)
    {
      cluster_size(parameters);
    }
  else if (parameters.opt_cluster_unoise != nullptr)
    {
      cluster_unoise(parameters);
    }
}


auto cmd_fastq_join(struct Parameters & parameters) -> void
{
  if (is_not_ASCII(parameters.opt_join_padgap)) {
    fatal("Option --join_padgap contains non-ASCII characters");
  }
  if (is_not_ASCII(parameters.opt_join_padgapq)) {
    fatal("Option --join_padgapq contains non-ASCII characters");
  }
  if ((not parameters.opt_join_padgapq_set_by_user) and
      (parameters.opt_fastq_ascii != Parameters::default_ascii_offset)) {
    std::string const alternative_quality_padding = "hhhhhhhh";  // Q40 with an offset of 64
    parameters.opt_join_padgapq = alternative_quality_padding;
  }
  fastq_join(parameters);
}


auto cmd_fastq_mergepairs(struct Parameters const & parameters) -> void
{
  if (parameters.opt_reverse == nullptr)
    {
      fatal("No reverse reads file specified with --reverse");
    }
  if ((parameters.opt_fastqout == nullptr) and
      (parameters.opt_fastaout == nullptr) and
      (parameters.opt_fastqout_notmerged_fwd == nullptr) and
      (parameters.opt_fastqout_notmerged_rev == nullptr) and
      (parameters.opt_fastaout_notmerged_fwd == nullptr) and
      (parameters.opt_fastaout_notmerged_rev == nullptr) and
      (parameters.opt_eetabbedout == nullptr))
    {
      fatal("No output files specified");
    }
  if (parameters.opt_fastq_maxdiffs < 0) {
    fatal("Argument to --fastq_maxdiffs must be positive");
  }
  if (parameters.opt_fastq_maxee <= 0.0) {
    /* expected error is the sum of per-base error probabilities;
       probabilities are strictly positive (min quality score is 93,
       corresponding to probability ~1e-9.3), so the sum cannot be
       zero or negative. A null or negative threshold would always
       reject every read and is almost certainly a user mistake. */
    fatal("Argument to --fastq_maxee must be a strictly positive number");
  }
  if (parameters.opt_fastq_maxlen < 1) {
    fatal("Argument to --fastq_maxlen must be a positive integer");
  }
  if (parameters.opt_fastq_minlen < 1) {
    fatal("Argument to --fastq_minlen must be a positive integer");
  }
  if (parameters.opt_fastq_maxns < 0) {
    fatal("Argument to --fastq_maxns must be a non-negative integer");
  }
  if (parameters.opt_fastq_maxmergelen < 1) {
    fatal("Argument to --fastq_maxmergelen must be a positive integer");
  }
  if (parameters.opt_fastq_minmergelen < 0) {
    fatal("Argument to --fastq_minmergelen must be a non-negative integer");
  }
  {
    /* Quality score range: 0..93 (Phred scores encoded in fastq).
       The default value is std::numeric_limits<long>::min(), meaning
       "no truncation"; skip the range check in that case so the default
       is preserved. */
    static constexpr auto long_min = std::numeric_limits<long>::min();
    if ((parameters.opt_fastq_truncqual != long_min) and
        ((parameters.opt_fastq_truncqual < 0) or (parameters.opt_fastq_truncqual > 93))) {
      fatal("Argument to --fastq_truncqual must be in range 0..93");
    }
  }
  fastq_mergepairs(parameters);
}


auto fill_prog_header(struct Parameters & parameters) -> void
{
  static constexpr auto max_line_length = std::size_t{80};
  static constexpr auto one_gigabyte = double{1024 * 1024 * 1024};
  auto const * const format = "%s v%s_%s, %.1fGB RAM, %ld cores";
  std::array<char, max_line_length> buffer {{}};
  static_cast<void>(std::snprintf(
      buffer.data(), max_line_length, format, PROG_NAME, PROG_VERSION,
      PROG_ARCH, static_cast<double>(system_get_memtotal()) / one_gigabyte,
      system_get_cores()));
  parameters.prog_header = buffer.data();
}


auto getentirecommandline(int argc, char** argv) -> std::string
{
  std::string command_line;
  for (int i = 0; i < argc; i++)
    {
      if (i > 0)
        {
          command_line += ' ';
        }
      command_line += argv[i];
    }
  return command_line;
}


auto show_header(struct Parameters const & parameters) -> void {
  if (parameters.opt_quiet) { return ; }
  std::fprintf(stderr, "%s\n", parameters.prog_header.c_str());
  std::fprintf(stderr, "https://github.com/torognes/vsearch\n");
  std::fprintf(stderr, "\n");
}


/* When building as a library (VSEARCH_NO_MAIN), exclude main() to avoid
   symbol conflicts with the embedding application's entry point.
   All global variable definitions and helper functions above remain
   available — only the CLI entry point is excluded. */
#ifndef VSEARCH_NO_MAIN
/* Installed via std::set_new_handler so that a failed C++ allocation
   (std::vector, std::string, ...) reports the same friendly message as the
   xmalloc path instead of throwing std::bad_alloc, which under -fno-exceptions
   would call std::terminate()/abort(). The handler runs while memory is
   exhausted, so it must not allocate: it writes a fixed string with write(2)
   and leaves via _exit() without flushing stdio. This catches a refused
   allocation (operator new gets nullptr); it cannot catch a kernel OOM kill
   of an overcommitted allocation. */
static auto vsearch_new_handler() -> void
{
  static char const message[] = "\n\nFatal error: Unable to allocate enough memory.\n";
  ssize_t const written = write(STDERR_FILENO, message, sizeof(message) - 1);
  (void) written;
  _exit(EXIT_FAILURE);
}


/* Run the command selected on the command line. Exactly one command option is
   set (enforced during parsing); the chain below picks it out and calls the
   matching handler. A few commands adjust parameters just before running, and
   an unrecognised (or empty) command falls through to cmd_none(). */
auto dispatch_command(struct Parameters & parameters) -> void
{
  if (parameters.opt_help)
    {
      help(parameters);
    }
  else if (parameters.opt_allpairs_global != nullptr)
    {
      parameters.opt_strand = false;
      parameters.opt_uc_allhits = true;
      cmd_allpairs_global(parameters);
    }
  else if (parameters.opt_usearch_global != nullptr)
    {
      cmd_usearch_global(parameters);
    }
  else if (parameters.opt_sortbysize != nullptr)
    {
      sortbysize(parameters);
    }
  else if (parameters.opt_sortbylength != nullptr)
    {
      sortbylength(parameters);
    }
  else if (parameters.opt_derep_fulllength != nullptr)
    {
      derep_fulllength(parameters);
    }
  else if (parameters.opt_derep_prefix != nullptr)
    {
      derep_prefix(parameters);
    }
  else if (parameters.opt_derep_smallmem != nullptr)
    {
      derep_smallmem(parameters);
    }
  else if (parameters.opt_derep_id != nullptr)
    {
      derep_id(parameters);
    }
  else if (parameters.opt_shuffle != nullptr)
    {
      shuffle(parameters);
    }
  else if (parameters.opt_fastx_subsample != nullptr)
    {
      cmd_subsample(parameters);
    }
  else if (parameters.opt_maskfasta != nullptr)
    {
      maskfasta(parameters);
    }
  else if ((parameters.opt_cluster_smallmem != nullptr) or (parameters.opt_cluster_fast != nullptr) or (parameters.opt_cluster_size != nullptr) or (parameters.opt_cluster_unoise != nullptr))
    {
      cmd_cluster(parameters);
    }
  else if (parameters.opt_uchime_denovo != nullptr)
    {
      uchime_denovo(parameters);
    }
  else if (parameters.opt_uchime_ref != nullptr)
    {
      uchime_ref(parameters);
    }
  else if (parameters.opt_uchime2_denovo != nullptr)
    {
      uchime2_denovo(parameters);
    }
  else if (parameters.opt_uchime3_denovo != nullptr)
    {
      uchime3_denovo(parameters);
    }
  else if (parameters.opt_chimeras_denovo != nullptr)
    {
      chimeras_denovo(parameters);
    }
  else if (parameters.opt_fastq_chars != nullptr)
    {
      fastq_chars(parameters);
    }
  else if (parameters.opt_fastq_stats != nullptr)
    {
      fastq_stats(parameters);
    }
  else if (parameters.opt_fastq_filter != nullptr)
    {
      fastq_filter(parameters);
    }
  else if (parameters.opt_fastx_filter != nullptr)
    {
      fastx_filter(parameters);
    }
  else if (parameters.opt_fastx_syncpairs != nullptr)
    {
      fastx_syncpairs(parameters);
    }
  else if (parameters.opt_fastx_revcomp != nullptr)
    {
      fastx_revcomp(parameters);
    }
  else if (parameters.opt_search_exact != nullptr)
    {
      cmd_search_exact(parameters);
    }
  else if (parameters.opt_fastx_mask != nullptr)
    {
      fastx_mask(parameters);
    }
  else if (parameters.opt_fastq_convert != nullptr)
    {
      fastq_convert(parameters);
    }
  else if (parameters.opt_fastq_mergepairs != nullptr)
    {
      cmd_fastq_mergepairs(parameters);
    }
  else if (parameters.opt_fastq_eestats != nullptr)
    {
      fastq_eestats(parameters);
    }
  else if (parameters.opt_fastq_eestats2 != nullptr)
    {
      fastq_eestats2(parameters);
    }
  else if (parameters.opt_fastq_join != nullptr)
    {
      cmd_fastq_join(parameters);
    }
  else if (parameters.opt_rereplicate != nullptr)
    {
      parameters.opt_xsize = true;
      rereplicate(parameters);
    }
  else if (parameters.opt_version)
    {
      version(parameters);
    }
  else if (parameters.opt_makeudb_usearch != nullptr)
    {
      makeudb_usearch(parameters);
    }
  else if (parameters.opt_udb2fasta != nullptr)
    {
      udb2fasta(parameters);
    }
  else if (parameters.opt_udbinfo != nullptr)
    {
      udbinfo(parameters);
    }
  else if (parameters.opt_udbstats != nullptr)
    {
      udbstats(parameters);
    }
  else if (parameters.opt_sintax != nullptr)
    {
      sintax(parameters);
    }
  else if (parameters.opt_sff_convert != nullptr)
    {
      sff_convert(parameters);
    }
  else if (parameters.opt_fastx_getseq != nullptr)
    {
      fastx_getseq(parameters);
    }
  else if (parameters.opt_fastx_getseqs != nullptr)
    {
      fastx_getseqs(parameters);
    }
  else if (parameters.opt_fastx_getsubseq != nullptr)
    {
      fastx_getsubseq(parameters);
    }
  else if (parameters.opt_cut != nullptr)
    {
      cut(parameters);
    }
  else if (parameters.opt_orient != nullptr)
    {
      orient(parameters);
    }
  else if (parameters.opt_fasta2fastq != nullptr)
    {
      fasta2fastq(parameters);
    }
  else if (parameters.opt_fastx_uniques != nullptr)
    {
      fastx_uniques(parameters);
    }
  else
    {
      cmd_none(parameters);
    }
}


auto main(int argc, char** argv) -> int
{
  std::set_new_handler(vsearch_new_handler);

  struct Parameters parameters;

  fill_prog_header(parameters);

  parameters.command_line = getentirecommandline(argc, argv);

#ifdef __x86_64__
  cpu_features_detect(parameters);
#endif

  args_init(argc, argv, parameters);

  /* The log file's Started banner and Finished/elapsed/max-memory footer are
     written by the LogFile constructor and destructor; the explicit scope
     makes the destructor (and hence the footer) run before the stdout flush
     check below, matching the original ordering. */
  {
    LogFile const log_file(parameters);

    random_init(parameters);

    show_header(parameters);

    /* RAII: loads the optional compression libraries here and closes them
       when this scope ends (replaces the former dynlibs_open/close pair,
       whose close ran outside this block). */
    DynamicLibraries const dynamic_libraries;
    parameters.dyn_libs = &dynamic_libraries;

#ifdef __x86_64__
    cpu_features_test(parameters);
#endif

    dispatch_command(parameters);
  }

  /* Output written directly to stdout is not closed through an OutputFileHandle
     (whose CheckedCloseOutputHandle deleter surfaces deferred write errors), so
     surface any deferred write error here (full disk, quota, or a broken pipe
     such as `vsearch ... | head`) rather than exiting 0 with truncated output
     (I1). */
  if ((std::fflush(stdout) != 0) or (std::ferror(stdout) != 0))
    {
      fatal("Unable to write to standard output (disk full, quota exceeded, or broken pipe?)");
    }
}
#endif /* VSEARCH_NO_MAIN */
