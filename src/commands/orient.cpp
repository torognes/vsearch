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

#include "utils/span.hpp"
#include "utils/view.hpp"
#include "vsearch.hpp"
#include <memory>  // std::unique_ptr
#include "utils/progress.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/db.hpp"
#include "core/dbindex.hpp"
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "core/mask.hpp"
#include "core/udb.hpp"
#include "core/unique.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/print_view.hpp"  // fprint
#include <algorithm>  // std::reverse_copy
#include <array>
#include <cassert>
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fputs, std::size_t
#include <vector>


// documentation: assuming opt_wordlength = 3 (6 bits)
// input       output
// 0b000000 -> 0b111111
// 0b111111 -> 0b000000
// 0b111100 -> 0b110000
// 0b110000 -> 0b111100
// 0b001100 -> 0b110011
// 0b000011 -> 0b001111
// 0b001111 -> 0b000011
// 0b100001 -> 0b101101
// 0b011110 -> 0b010010
// 0b101010 -> 0b010101
// 0b010101 -> 0b101010
namespace {
/* reverse complement of the four bases packed in one byte, tabulated once;
   each entry is produced by rc_kmer's original per-base loop applied to a
   4-base word, so a lookup is exactly that loop's result */
auto rc_byte_table() -> std::array<unsigned char, 256> const &
{
  static std::array<unsigned char, 256> const table = []() -> std::array<unsigned char, 256> {
    std::array<unsigned char, 256> values {{}};
    for (auto byte = 0U; byte < 256U; ++byte)
      {
        auto fwd = byte;
        auto rev = 0U;
        for (auto rank = 0U; rank < 4U; ++rank)
          {
            // compute complement of the last two bits
            auto const complement_bits = (fwd & 3U) ^ 3U;
            fwd = fwd >> 2U;
            rev = rev << 2U;
            rev |= complement_bits;
          }
        values[byte] = static_cast<unsigned char>(rev);
      }
    return values;
  }();
  return table;
}


auto rc_kmer(unsigned int const kmer, unsigned int const wordlength) -> unsigned int
{
  /* reverse complement a kmer where k = wordlength */

  /* wordlength is the effective index word length (Dbindex::wordlength, set by
     Dbindex::prepare for a FASTA db, or by udb_read for a UDB db whose stored
     width overrides the configured one). Query kmers must be extracted at this
     width to match the index; reading parameters.opt_wordlength here would use
     the wrong width against a UDB index (mismatch, out-of-bounds when wider). */
  assert(wordlength != 0);
  assert(wordlength * 2 <= 32);
  auto const & table = rc_byte_table();
  /* reverse complement all sixteen 2-bit slots of the 32-bit word, one byte
     at a time; the empty slots above the kmer come out as low-order 3s and
     the final shift drops exactly those */
  auto const rc_full_word =
      (static_cast<unsigned int>(table[kmer & 0xFFU]) << 24U)
    | (static_cast<unsigned int>(table[(kmer >> 8U) & 0xFFU]) << 16U)
    | (static_cast<unsigned int>(table[(kmer >> 16U) & 0xFFU]) << 8U)
    |  static_cast<unsigned int>(table[kmer >> 24U]);
  return rc_full_word >> (32U - (2U * wordlength));
}


/* the orientation assigned to one query, printed as-is in --tabbedout */
enum struct Strand : char
{
  positive = '+',
  negative = '-',
  undetermined = '?',
};


/* what the orienting loop counted over the whole input (no default member
   initializers: they would make this a non-aggregate before C++14) */
struct OrientationTotals
{
  int queries;
  int qmatches;
  int matches_fwd;
  int matches_rev;
  int notmatched;
};


/* the same summary goes to stderr (unless --quiet) and to --log (when set) */
auto report_orient(std::FILE * output_stream,
                   struct OrientationTotals const & totals) -> void
{
  auto const queries = totals.queries;
  auto const print_line = [output_stream, queries](char const * label,
                                                   int const count) -> void {
    std::fputs(label, output_stream);
    fprint_integer(output_stream, count);
    if (queries > 0)
      {
        fprint(output_stream, " (");
        std::fprintf(output_stream, "%.2f", 100.0 * count / queries);
        fprint(output_stream, "%)");
      }
    fprint(output_stream, '\n');
  };
  print_line("Forward oriented sequences: ", totals.matches_fwd);
  print_line("Reverse oriented sequences: ", totals.matches_rev);
  print_line("All oriented sequences:     ", totals.qmatches);
  print_line("Not oriented sequences:     ", totals.notmatched);
  fprint(output_stream, "Total number of sequences:  ");
  fprint_integer(output_stream, queries);
  fprint(output_stream, '\n');
}
}  // anonymous namespace


auto orient(struct Parameters const & parameters) -> void
{
  std::unique_ptr<fastx_s> query_h;
  // refactoring: use struct, like in subsample
  OutputFileHandle fastaout_handle;
  OutputFileHandle fastqout_handle;
  OutputFileHandle tabbedout_handle;
  OutputFileHandle notmatched_handle;
  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_tabbedout = nullptr;
  std::FILE * fp_notmatched = nullptr;

  int queries = 0;
  int qmatches = 0;
  int matches_fwd = 0;
  int matches_rev = 0;
  int notmatched = 0;

  /* check arguments */

  if (parameters.opt_db == nullptr)
    {
      fatal("Database not specified with --db");
    }

  if (not ((parameters.opt_fastaout != nullptr) or (parameters.opt_fastqout != nullptr) or (parameters.opt_notmatched != nullptr) or (parameters.opt_tabbedout != nullptr)))
    {
      fatal("Output file not specified with --fastaout, --fastqout, --notmatched or --tabbedout");
    }

  /* prepare reading of queries */

  query_h = fastx_open(parameters.opt_orient, parameters);

  /* open output files */

  fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  fp_fastaout = fastaout_handle.get();

  if (parameters.opt_fastqout != nullptr)
    {
      if (not query_h->is_fastq_input())
        {
          fatal("Cannot write FASTQ output with FASTA input");
        }

      fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
      fp_fastqout = fastqout_handle.get();
    }

  notmatched_handle = open_optional_output_file(parameters.opt_notmatched, OutputOption{"--notmatched"});
  fp_notmatched = notmatched_handle.get();
  tabbedout_handle = open_optional_output_file(parameters.opt_tabbedout, OutputOption{"--tabbedout"});
  fp_tabbedout = tabbedout_handle.get();

  /* the sequence database this run owns (RAII) */
  Database db;
  /* the k-mer index this run owns (RAII) */
  Dbindex dbindex;

  /* check if it may be an UDB file */

  auto const is_udb = udb_detect_isudb(parameters.opt_db);

  if (is_udb)
    {
      udb_read(parameters.opt_db, UdbUse::search, dbindex, db, parameters);
    }
  else
    {
      db.read(parameters.opt_db, 0, parameters);
    }

  if (not is_udb)
    {
      if (parameters.opt_dbmask == Masking::dust)
        {
          dust_all(db, parameters);
        }
      else if ((parameters.opt_dbmask == Masking::soft) and parameters.opt_hardmask)
        {
          hardmask_all(db);
        }
      dbindex.prepare(parameters.opt_dbmask, db, parameters);
      dbindex.add_all_sequences(parameters.opt_dbmask, db, parameters);
    }

  Uniquer uh_fwd;

  std::size_t alloc = 0;
  std::vector<char> qseq_rev;
  std::vector<char> query_qual_rev;

  {
    Progress progress_bar("Orienting sequences", query_h->get_size(), parameters);

    while (query_h->next(
                      (not parameters.opt_notrunclabels),
                      chrmap_no_change()))
      {
        auto const query_head = query_h->header_view();
        auto const query_sequence = query_h->sequence_view();
        auto const qseqlen = static_cast<int>(query_sequence.size());
        int64_t const qsize = query_h->get_abundance();
        auto const query_qual_fwd = query_h->quality_view();

        /* find kmers in query sequence */

        /* dbindex.wordlength: the effective index width (see rc_kmer) */
        auto const kmer_list_fwd = uh_fwd.count(static_cast<int>(dbindex.wordlength),
                                                query_sequence, parameters.opt_qmask);

        /* count kmers matching on each strand */

        unsigned int count_fwd = 0;
        unsigned int count_rev = 0;
        constexpr auto hits_factor = 8U;

        for (auto const kmer_fwd : kmer_list_fwd)
          {
            unsigned int const kmer_rev = rc_kmer(kmer_fwd, dbindex.wordlength);

            unsigned int const hits_fwd = dbindex.getmatchcount(kmer_fwd);
            unsigned int const hits_rev = dbindex.getmatchcount(kmer_rev);

            /* require 8 times as many matches on one stand than the other */

            if (hits_fwd > hits_factor * hits_rev)
              {
                ++count_fwd;
              }
            else if (hits_rev > hits_factor * hits_fwd)
              {
                ++count_rev;
              }
          }

        /* get progress as amount of input file read */

        auto const progress = query_h->get_position();

        /* update stats */

        ++queries;

        auto strand = Strand::undetermined;
        auto const min_count = 1U;
        auto const min_factor = 4U;

        if ((count_fwd >= min_count) and (count_fwd >= min_factor * count_rev))
          {
            /* fwd */

            strand = Strand::positive;
            ++matches_fwd;
            ++qmatches;

            if (parameters.opt_fastaout != nullptr)
              {
                fasta_print_general(fp_fastaout,
                                    nullptr,
                                    query_h->record(),
                                    OutputAnnotations{static_cast<uint64_t>(qsize), qmatches},
                                    parameters);
              }

            if (parameters.opt_fastqout != nullptr)
              {
                fastq_print_general(fp_fastqout,
                                    query_h->record(),
                                    OutputAnnotations{static_cast<uint64_t>(qsize), qmatches},
                                    parameters);
              }
          }
        else if ((count_rev >= min_count) and (count_rev >= min_factor * count_fwd))
          {
            /* rev */

            strand = Strand::negative;
            ++matches_rev;
            ++qmatches;

            /* alloc more mem if necessary to keep reverse sequence and qual */
            assert(qseqlen > 0);
            static_assert(sizeof(std::size_t) >= sizeof(int), "size_t is too small");
            const std::size_t requirements = static_cast<std::size_t>(qseqlen) + 1;
            // refactoring: unsigned int qseqlen
            if (requirements > alloc)
              {
                alloc = requirements;
                qseq_rev.resize(alloc);
                if (query_h->is_fastq_input())
                  {
                    query_qual_rev.resize(alloc);
                  }
              }

            /* get reverse complementary sequence */

            reverse_complement(make_span(qseq_rev).first(query_sequence.size() + 1), query_sequence);
            auto const rc_sequence = make_view(qseq_rev).first(query_sequence.size());

            if (parameters.opt_fastaout != nullptr)
              {
                fasta_print_general(fp_fastaout,
                                    nullptr,
                                    rc_sequence,
                                    query_head,
                                    OutputAnnotations{static_cast<uint64_t>(qsize), qmatches},
                                    parameters);
              }

            if (parameters.opt_fastqout != nullptr)
              {
                /* reverse quality scores */

                if (query_h->is_fastq_input())
                  {
                    // copy query string in reverse order
                    std::reverse_copy(query_qual_fwd.cbegin(), query_qual_fwd.cend(),
                                      query_qual_rev.begin());
                    query_qual_rev[query_qual_fwd.size()] = '\0';
                  }

                fastq_print_general(fp_fastqout,
                                    rc_sequence,
                                    query_head,
                                    make_view(query_qual_rev).first(query_sequence.size()),
                                    OutputAnnotations{static_cast<uint64_t>(qsize), qmatches},
                                    parameters);
              }
          }
        else
          {
            /* undecided */

            strand = Strand::undetermined;
            ++notmatched;

            if (parameters.opt_notmatched != nullptr)
              {
                if (query_h->is_fastq_input())
                  {
                    fastq_print_general(fp_notmatched,
                                        query_h->record(),
                                        OutputAnnotations{static_cast<uint64_t>(qsize), notmatched},
                                        parameters);
                  }
                else
                  {
                    fasta_print_general(fp_notmatched,
                                        nullptr,
                                        query_h->record(),
                                        OutputAnnotations{static_cast<uint64_t>(qsize), notmatched},
                                        parameters);
                  }
              }
          }

        if (parameters.opt_tabbedout != nullptr)
          {
            fprint(fp_tabbedout, query_head);
            fprint(fp_tabbedout, '\t');
            fprint(fp_tabbedout, static_cast<char>(strand));
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, count_fwd);
            fprint(fp_tabbedout, '\t');
            fprint_integer(fp_tabbedout, count_rev);
            fprint(fp_tabbedout, '\n');
          }

        /* show progress */

        progress_bar.update(progress);
      }
  }


  /* clean up */

  /* uh_fwd (a Uniquer) frees its buffers when it goes out of scope (RAII) */

  dbindex.clear();
  db.clear();

  if (parameters.opt_tabbedout != nullptr)
    {
      tabbedout_handle.reset();
    }
  if (parameters.opt_notmatched != nullptr)
    {
      notmatched_handle.reset();
    }
  if (parameters.opt_fastqout != nullptr)
    {
      fastqout_handle.reset();
    }
  if (parameters.opt_fastaout != nullptr)
    {
      fastaout_handle.reset();
    }

  query_h->report_stripped_warning(parameters);

  OrientationTotals const totals = {queries, qmatches, matches_fwd, matches_rev, notmatched};

  if (not parameters.opt_quiet)
    {
      report_orient(stderr, totals);
    }

  if (parameters.fp_log != nullptr)
    {
      report_orient(parameters.fp_log, totals);
    }
}
