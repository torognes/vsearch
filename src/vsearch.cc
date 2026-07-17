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

#include "vsearch.hpp"
#include "commands/allpairs_global.hpp"
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
#include "os/system.hpp"  // system_get_cores, system_get_memtotal
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
#include "arch/cpu_features.hpp"
#include "utils/fatal.hpp"
#include "utils/logfile.hpp"  // LogFile
#include "utils/prog_id.hpp"  // PROG_NAME, PROG_VERSION, PROG_ARCH
#include "utils/random.hpp"
#include <array>
#include <cerrno>  // errno, ERANGE
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t, std::snprintf, std::printf
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::strlen
#include <getopt.h>  // getopt_long_only, optarg, optind, opterr, struct
                     // option (no_argument, required_argument)
#include <new>  // std::set_new_handler
#include <string>
#include <unistd.h>  // write, _exit, STDERR_FILENO
#include <vector>


/* Standalone CLI entry point: main(), the command dispatcher, and their
   internal-linkage helpers. Listed only in the vsearch executable's sources
   (never libvsearch_core), so the library defines no main() and cannot clash
   with an embedding application's entry point. */
namespace {
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
  if (argc <= 0)
    {
      return command_line;
    }

  // Size the result up front so it is built in a single allocation: the sum
  // of every argument's length plus one separating space between adjacent
  // arguments.
  auto total_length = static_cast<std::string::size_type>(argc - 1);
  for (int i = 0; i < argc; ++i)
    {
      total_length += std::strlen(argv[i]);
    }
  command_line.reserve(total_length);

  for (int i = 0; i < argc; ++i)
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
}  // end of anonymous namespace


namespace {

/* Installed via std::set_new_handler so that a failed C++ allocation
   (std::vector, std::string, ...) reports the same friendly message as the
   xmalloc path instead of throwing std::bad_alloc, which under -fno-exceptions
   would call std::terminate()/abort(). The handler runs while memory is
   exhausted, so it must not allocate: it writes a fixed string with write(2)
   and leaves via _exit() without flushing stdio. This catches a refused
   allocation (operator new gets nullptr); it cannot catch a kernel OOM kill
   of an overcommitted allocation. */
auto vsearch_new_handler() -> void
{
  static char const message[] = "\n\nFatal error: Unable to allocate enough memory.\n";
  ssize_t const written = write(STDERR_FILENO, message, sizeof(message) - 1);
  (void) written;
  _exit(EXIT_FAILURE);
}


/* Run the command selected on the command line. The command was resolved once
   during parsing (args_init) and is passed in, so dispatch is a plain table
   lookup rather than a chain of opt_* pointer tests. Configuration is fully
   resolved by args_init(), so this reads `parameters` without modifying it; an
   unrecognised (or empty) command falls through to cmd_none() via default. */
auto dispatch_command(struct Parameters const & parameters, Command const command) -> void
{
  switch (command)
    {
    case Command::help:             help(parameters);             break;
    case Command::version:          version(parameters);          break;
    case Command::allpairs_global:  allpairs_global(parameters);  break;
    case Command::usearch_global:   usearch_global(parameters);   break;
    case Command::search_exact:     search_exact(parameters);     break;
    case Command::sintax:           sintax(parameters);           break;
    case Command::orient:           orient(parameters);           break;
    case Command::cluster_fast:     cluster_fast(parameters);     break;
    case Command::cluster_smallmem: cluster_smallmem(parameters); break;
    case Command::cluster_size:     cluster_size(parameters);     break;
    case Command::cluster_unoise:   cluster_unoise(parameters);   break;
    case Command::uchime_denovo:    uchime_denovo(parameters);    break;
    case Command::uchime2_denovo:   uchime2_denovo(parameters);   break;
    case Command::uchime3_denovo:   uchime3_denovo(parameters);   break;
    case Command::uchime_ref:       uchime_ref(parameters);       break;
    case Command::chimeras_denovo:  chimeras_denovo(parameters);  break;
    case Command::derep_fulllength: derep_fulllength(parameters); break;
    case Command::derep_prefix:     derep_prefix(parameters);     break;
    case Command::derep_id:         derep_id(parameters);         break;
    case Command::derep_smallmem:   derep_smallmem(parameters);   break;
    case Command::fastq_chars:      fastq_chars(parameters);      break;
    case Command::fastq_stats:      fastq_stats(parameters);      break;
    case Command::fastq_filter:     fastq_filter(parameters);     break;
    case Command::fastx_filter:     fastx_filter(parameters);     break;
    case Command::fastq_convert:    fastq_convert(parameters);    break;
    case Command::fastq_eestats:    fastq_eestats(parameters);    break;
    case Command::fastq_eestats2:   fastq_eestats2(parameters);   break;
    case Command::fastq_join:       fastq_join(parameters);       break;
    case Command::fastq_mergepairs: fastq_mergepairs(parameters); break;
    case Command::fastx_uniques:    fastx_uniques(parameters);    break;
    case Command::fastx_mask:       fastx_mask(parameters);       break;
    case Command::fastx_revcomp:    fastx_revcomp(parameters);    break;
    case Command::fastx_syncpairs:  fastx_syncpairs(parameters);  break;
    case Command::fastx_getseq:     fastx_getseq(parameters);     break;
    case Command::fastx_getseqs:    fastx_getseqs(parameters);    break;
    case Command::fastx_getsubseq:  fastx_getsubseq(parameters);  break;
    case Command::fastx_subsample:  subsample(parameters);        break;
    case Command::fasta2fastq:      fasta2fastq(parameters);      break;
    case Command::cut:              cut(parameters);              break;
    case Command::shuffle:          shuffle(parameters);          break;
    case Command::sortbylength:     sortbylength(parameters);     break;
    case Command::sortbysize:       sortbysize(parameters);       break;
    case Command::rereplicate:      rereplicate(parameters);      break;
    case Command::maskfasta:        maskfasta(parameters);        break;
    case Command::sff_convert:      sff_convert(parameters);      break;
    case Command::makeudb_usearch:  makeudb_usearch(parameters);  break;
    case Command::udb2fasta:        udb2fasta(parameters);        break;
    case Command::udbinfo:          udbinfo(parameters);          break;
    case Command::udbstats:         udbstats(parameters);         break;
    default:                        cmd_none(parameters);         break;
    }
}


/* Output written directly to stdout is not closed through an OutputFileHandle
   (whose CheckedCloseOutputHandle deleter surfaces deferred write errors), so
   surface any deferred write error here (full disk, quota, or a broken pipe
   such as `vsearch ... | head`) rather than exiting 0 with truncated output
   (I1). */
auto flush_stdout() -> void
{
  if ((std::fflush(stdout) != 0) or (std::ferror(stdout) != 0))
    {
      fatal("Unable to write to standard output (disk full, quota exceeded, or broken pipe?)");
    }
}

}  // end of anonymous namespace


auto main(int argc, char** argv) -> int
{
  std::set_new_handler(vsearch_new_handler);

  struct Parameters parameters;

  fill_prog_header(parameters);

  parameters.command_line = getentirecommandline(argc, argv);

  cpu_features_detect(parameters);

  Command const command = args_init(argc, argv, parameters);

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

    cpu_features_test(parameters);

    dispatch_command(parameters, command);
  }

  flush_stdout();
}
