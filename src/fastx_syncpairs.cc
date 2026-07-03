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
#include "fasta.h"  // fasta_print_general
#include "fastq.h"  // fastq_print_general
#include "fastx.h"  // fastx_handle
#include "util.h"  // progress_init, progress_update, progress_done, fopen_output
#include "utils/fatal.hpp"
#include "utils/maps.hpp"  // chrmap_no_change_vector
#include <cinttypes>  // macros PRIu64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <string>
#include <unordered_map>
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // default set of characters that may introduce the mate number 1 or 2
  // (used when the user does not provide --read_separators). A whitespace
  // is always treated as a separator, regardless of this set.
  char const * const default_read_separators = "/";

  struct output_file {
    char * name = nullptr;
    std::FILE * handle = nullptr;
  };

  // a destination for a category of reads (synced or orphaned), which may
  // be written in fasta and/or fastq format depending on user options
  struct output_pair {
    output_file fasta;
    output_file fastq;
  };

  struct output_files {
    output_pair synced_fwd;
    output_pair synced_rev;
    output_pair orphans_fwd;
    output_pair orphans_rev;
  };

  // a single read kept in memory while indexing the reverse file
  struct read_record {
    std::string header;
    std::string sequence;
    std::string quality;  // empty when the input is in fasta format
    int64_t abundance = 1;
  };

  using read_index = std::unordered_map<std::string, std::size_t>;


  auto check_parameters(struct Parameters const & parameters) -> void {
    if (parameters.opt_reverse == nullptr) {
      fatal("No reverse reads file specified with --reverse");
    }

    if ((parameters.opt_fastaout == nullptr) and
        (parameters.opt_fastqout == nullptr) and
        (parameters.opt_fastaout_rev == nullptr) and
        (parameters.opt_fastqout_rev == nullptr) and
        (parameters.opt_fastaout_orphans == nullptr) and
        (parameters.opt_fastqout_orphans == nullptr) and
        (parameters.opt_fastaout_orphans_rev == nullptr) and
        (parameters.opt_fastqout_orphans_rev == nullptr)) {
      fatal("No output files specified");
    }
  }


  auto requests_fastq_output(struct Parameters const & parameters) -> bool {
    return (parameters.opt_fastqout != nullptr) or
      (parameters.opt_fastqout_rev != nullptr) or
      (parameters.opt_fastqout_orphans != nullptr) or
      (parameters.opt_fastqout_orphans_rev != nullptr);
  }


  auto open_output_file(char * name) -> output_file {
    output_file outfile;
    outfile.name = name;
    if (outfile.name != nullptr) {
      outfile.handle = fopen_output(outfile.name);
      if (outfile.handle == nullptr) {
        fatal("Unable to open file for writing (%s)", outfile.name);
      }
    }
    return outfile;
  }


  auto open_output_files(struct Parameters const & parameters) -> output_files {
    output_files outfiles;
    outfiles.synced_fwd.fasta = open_output_file(parameters.opt_fastaout);
    outfiles.synced_fwd.fastq = open_output_file(parameters.opt_fastqout);
    outfiles.synced_rev.fasta = open_output_file(parameters.opt_fastaout_rev);
    outfiles.synced_rev.fastq = open_output_file(parameters.opt_fastqout_rev);
    outfiles.orphans_fwd.fasta = open_output_file(parameters.opt_fastaout_orphans);
    outfiles.orphans_fwd.fastq = open_output_file(parameters.opt_fastqout_orphans);
    outfiles.orphans_rev.fasta = open_output_file(parameters.opt_fastaout_orphans_rev);
    outfiles.orphans_rev.fastq = open_output_file(parameters.opt_fastqout_orphans_rev);
    return outfiles;
  }


  auto close_output_files(output_files const & outfiles) -> void {
    for (auto const * pair : {& outfiles.synced_fwd, & outfiles.synced_rev,
                              & outfiles.orphans_fwd, & outfiles.orphans_rev}) {
      for (auto * handle : {pair->fasta.handle, pair->fastq.handle}) {
        if (handle != nullptr) {
          static_cast<void>(fclose_output(handle));
        }
      }
    }
  }


  // Derive the key shared by the two mates of a pair. Casava 1.8+ headers
  // ("instrument... 1:N:0:" and "... 2:N:0:") already share the substring
  // before the first whitespace, so truncating there is enough. Older
  // headers carry the mate number as a "/1" or "/2" suffix, removed here
  // when its separator belongs to the configured set.
  auto matching_key(char const * header,
                    uint64_t const header_length,
                    std::string const & separators) -> std::string {
    std::string key {header, header_length};

    auto const blank = key.find_first_of(" \t");
    if (blank != std::string::npos) {
      key.resize(blank);
    }

    if (key.size() >= 2) {
      auto const last = key.back();
      auto const separator = key[key.size() - 2];
      if (((last == '1') or (last == '2')) and
          (separators.find(separator) != std::string::npos)) {
        key.resize(key.size() - 2);
      }
    }

    return key;
  }


  auto store_record(fastx_handle handle, bool const is_fastq) -> read_record {
    read_record record;
    record.header.assign(fastx_get_header(handle), fastx_get_header_length(handle));
    record.sequence.assign(fastx_get_sequence(handle), fastx_get_sequence_length(handle));
    if (is_fastq) {
      record.quality.assign(fastx_get_quality(handle), fastx_get_sequence_length(handle));
    }
    record.abundance = fastx_get_abundance(handle);
    return record;
  }


  auto write_record(output_pair const & destination,
                    read_record const & record,
                    int64_t const ordinal) -> void {
    auto const length = static_cast<int>(record.sequence.length());
    auto const header_length = static_cast<int>(record.header.length());
    if (destination.fastq.handle != nullptr) {
      fastq_print_general(destination.fastq.handle,
                          record.sequence.c_str(),
                          length,
                          record.header.c_str(),
                          header_length,
                          record.quality.c_str(),
                          static_cast<uint64_t>(record.abundance),
                          ordinal,
                          -1.0);
    }
    if (destination.fasta.handle != nullptr) {
      fasta_print_general(destination.fasta.handle,
                          nullptr,
                          record.sequence.c_str(),
                          length,
                          record.header.c_str(),
                          header_length,
                          static_cast<uint64_t>(record.abundance),
                          ordinal,
                          -1.0,
                          -1,
                          -1,
                          nullptr,
                          0,
                          0);
    }
  }


  // Read the reverse file once, keeping every record in memory (in file
  // order, for orphan output) and mapping its matching key to its position.
  auto index_reverse(fastx_handle reverse_handle,
                     bool const is_fastq,
                     std::string const & separators,
                     std::vector<read_record> & records,
                     read_index & index) -> void {
    progress_init("Indexing reverse reads", fastx_get_size(reverse_handle));
    while (fastx_next(reverse_handle, false, chrmap_no_change_vector.data())) {
      auto key = matching_key(fastx_get_header(reverse_handle),
                              fastx_get_header_length(reverse_handle),
                              separators);
      auto const position = records.size();
      auto const inserted = index.emplace(std::move(key), position);
      if (not inserted.second) {
        fatal("Duplicate read label in reverse file");
      }
      records.push_back(store_record(reverse_handle, is_fastq));
      progress_update(fastx_get_position(reverse_handle));
    }
    progress_done();
  }


  auto stats_message(std::FILE * output_stream,
                     uint64_t const pairs,
                     uint64_t const orphans_fwd,
                     uint64_t const orphans_rev) -> void {
    static_cast<void>(std::fprintf(output_stream,
                                   "%" PRIu64 " pairs synchronized, "
                                   "%" PRIu64 " forward and "
                                   "%" PRIu64 " reverse orphan reads\n",
                                   pairs, orphans_fwd, orphans_rev));
  }

}  // end of anonymous namespace


auto fastx_syncpairs(struct Parameters const & parameters) -> void
{
  /* check parameters */

  check_parameters(parameters);

  /* open and check input files */

  auto * forward_handle = fastx_open(parameters.opt_fastx_syncpairs);
  auto * reverse_handle = fastx_open(parameters.opt_reverse);

  auto const forward_empty = fastx_is_empty(forward_handle);
  auto const reverse_empty = fastx_is_empty(reverse_handle);
  auto const forward_is_fastq = fastx_is_fastq(forward_handle);
  auto const reverse_is_fastq = fastx_is_fastq(reverse_handle);

  if ((not forward_empty) and (not reverse_empty) and
      (forward_is_fastq != reverse_is_fastq)) {
    fatal("Forward and reverse files must both be FASTA or both FASTQ");
  }

  /* the effective format follows the non-empty input file */
  auto const is_fastq = forward_empty ? reverse_is_fastq : forward_is_fastq;

  if ((not is_fastq) and requests_fastq_output(parameters)) {
    fatal("Cannot write FASTQ output from FASTA input (no quality scores)");
  }

  /* open output files */

  auto const outfiles = open_output_files(parameters);

  std::string const separators =
    (parameters.opt_read_separators != nullptr) ?
    parameters.opt_read_separators : default_read_separators;

  /* index the reverse file (read once, kept in memory) */

  std::vector<read_record> reverse_records;
  read_index reverse_index;
  index_reverse(reverse_handle, is_fastq, separators, reverse_records, reverse_index);

  /* stream the forward file, emitting synced pairs in forward order */

  std::vector<bool> reverse_used(reverse_records.size(), false);
  uint64_t pairs = 0;
  uint64_t orphans_fwd = 0;

  progress_init("Synchronizing reads", fastx_get_size(forward_handle));
  while (fastx_next(forward_handle, false, chrmap_no_change_vector.data())) {
    auto const key = matching_key(fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  separators);
    auto const match = reverse_index.find(key);
    if (match == reverse_index.end()) {
      write_record(outfiles.orphans_fwd, store_record(forward_handle, is_fastq),
                   static_cast<int64_t>(orphans_fwd + 1));
      ++orphans_fwd;
    }
    else {
      auto const position = match->second;
      // a reverse read already claimed by an earlier forward read means
      // two forward reads share the same matching key: the pairing is
      // ambiguous. Forward orphans that share a key are harmless and are
      // not detected here (they are simply written out twice).
      if (reverse_used[position]) {
        fatal("Duplicate read label in forward file");
      }
      reverse_used[position] = true;
      ++pairs;
      write_record(outfiles.synced_fwd, store_record(forward_handle, is_fastq),
                   static_cast<int64_t>(pairs));
      write_record(outfiles.synced_rev, reverse_records[position],
                   static_cast<int64_t>(pairs));
    }
    progress_update(fastx_get_position(forward_handle));
  }
  progress_done();

  /* write the reverse reads that had no forward mate, in reverse order */

  uint64_t orphans_rev = 0;
  for (std::size_t position = 0; position < reverse_records.size(); ++position) {
    if (not reverse_used[position]) {
      write_record(outfiles.orphans_rev, reverse_records[position],
                   static_cast<int64_t>(orphans_rev + 1));
      ++orphans_rev;
    }
  }

  /* report */

  if (not parameters.opt_quiet) {
    stats_message(stderr, pairs, orphans_fwd, orphans_rev);
  }
  if (parameters.opt_log != nullptr) {
    stats_message(parameters.fp_log, pairs, orphans_fwd, orphans_rev);
  }

  /* clean up */

  close_output_files(outfiles);
  fastx_close(forward_handle);
  fastx_close(reverse_handle);
}
