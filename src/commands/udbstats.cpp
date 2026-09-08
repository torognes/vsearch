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
#include "commands/udbstats.hpp"
#include "core/db.hpp"
#include "core/udb.hpp"
#include "core/dbindex.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/span.hpp"  // make_span
#include "utils/view.hpp"
#include <algorithm>  // std::max, std::min, std::sort
#include <cassert>  // assert
#include <cmath>  // std::lround
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::fprintf
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct wordfreq
  {
    unsigned int kmer;
    unsigned int count;
  };

  using wordfreq_t = struct wordfreq;


  /* The order the report's word rows are in. It used to be implied rather than
     stated: wc_compare() sorted the whole table by count ascending and by
     k-mer descending among equal counts, and the printing loop then walked the
     result backwards, which is this. */
  auto ranks_above(wordfreq_t const & lhs, wordfreq_t const & rhs) -> bool
  {
    if (lhs.count != rhs.count)
      {
        return lhs.count > rhs.count;
      }
    return lhs.kmer < rhs.kmer;
  }

  /* How many word rows the report prints (the row loop breaks after i == 10),
     and how many of each row's matching sequence numbers it shows before the
     ellipsis. Named because the word list is no longer held in memory: these
     two numbers are exactly how much of it has to be fetched. */
  constexpr auto reported_rows = std::size_t{11};
  constexpr auto reported_entries_per_row = std::size_t{8};

}  // end of anonymous namespace


auto udbstats(struct Parameters const & parameters) -> void
{
  /* show word statistics for an UDB file */

  Database db;  /* the sequence database this run owns (RAII) */
  Dbindex dbindex;  /* the k-mer index this run owns (RAII) */

  /* read UDB file */

  udb_read(parameters.input_filename, UdbUse::word_stats, dbindex, db, parameters);

  /* Every line this command reports goes to the log file (documented in
     man/commands/vsearch-udbstats.1.md), so without --log there is nothing to
     produce: the analysis below builds a table of 4^wordlength entries and
     sorts it, and then prints none of it. Returning here leaves the loading
     summary udb_read() already wrote on stderr as the whole output, which is
     what such a run produces today. */

  if (parameters.fp_log == nullptr)
    {
      dbindex.clear();
      db.clear();
      return;
    }

  /* dbindex.wordlength below is the effective index width that udb_read() just
     published from this UDB file's header (which may differ from the configured
     parameters.opt_wordlength); read it, not the config (E1). */

  auto const seqcount = static_cast<unsigned int>(db.getsequencecount());
  auto const nt = db.getnucleotidecount();

  /* analyze word counts

     Everything reported below is a property of one distribution -- how many
     k-mer slots hold each count -- plus the highest-ranking slots: Max size is
     its largest value, Median size its two central ones, the Size lo / Size hi
     table its totals over doubling buckets, and the word rows the eleven slots
     that rank above all others. So the counts are summarised into that
     distribution in a single pass, rather than materialised as 4^wordlength
     {k-mer, count} pairs and sorted. A count is a number of database sequences
     (udb_read() rejects a larger one), so the distribution has seqcount + 1
     cells however wide the index is, and every figure taken from it is exact.

     What that removes, at word length 13: an 8-bytes-per-slot table (537 MB)
     and the sort of it, which was 74 % of the run's instructions at word
     length 10 and 37 % at the default 8. */

  std::vector<uint64_t> histogram(static_cast<std::size_t>(seqcount) + 1, 0);

  std::vector<wordfreq_t> top;  /* the reported rows, kept in report order */
  top.reserve(reported_rows);

  for (auto kmer = 0U; kmer < dbindex.hashsize; ++kmer)
    {
      auto const count = dbindex.kmercount[kmer];
      assert(count <= seqcount);  /* checked when the UDB was read */
      ++histogram[count];

      /* Strictly greater, so among equal counts the k-mer seen first -- the
         smallest -- keeps its place, which is the order the printing loop used
         to get by reading the sorted table from its far end. */
      if ((top.size() < reported_rows) or (count > top.back().count))
        {
          wordfreq_t const candidate {kmer, count};
          auto const where = std::upper_bound(top.begin(), top.end(), candidate, ranks_above);
          top.insert(where, candidate);
          if (top.size() > reported_rows)
            {
              top.pop_back();
            }
        }
    }

  assert(top.size() == std::min<std::size_t>(reported_rows, dbindex.hashsize));

  /* Max size, and the k-mer the report names beside it: the largest count, and
     the smallest k-mer holding it -- the first reported row, by construction. */
  auto const wcmax = top.front().count;
  auto const wcmax_kmer = top.front().kmer;

  /* the two central counts, as the sorted table gave them: the (hashsize/2)-th
     and (hashsize/2 + 1)-th smallest */
  auto const count_at_rank = [&histogram](uint64_t const rank) -> unsigned int
  {
    uint64_t seen = 0;
    for (std::size_t count = 0; count < histogram.size(); ++count)
      {
        seen += histogram[count];
        if (seen >= rank)
          {
            return static_cast<unsigned int>(count);
          }
      }
    assert(false);  /* the ranks asked for are below hashsize */
    return 0U;
  };
  auto const wcmedian = (count_at_rank(dbindex.hashsize / 2)
                         + count_at_rank((dbindex.hashsize / 2) + 1)) / 2;

  /* show stats */

  if (parameters.fp_log != nullptr)
    {
      fprint(parameters.fp_log, "      Alphabet  nt\n");
      fprint(parameters.fp_log, "    Word width  ");
      fprint_integer(parameters.fp_log, dbindex.wordlength);
      fprint(parameters.fp_log, '\n');
      fprint(parameters.fp_log, "     Word ones  ");
      fprint_integer(parameters.fp_log, dbindex.wordlength);
      fprint(parameters.fp_log, '\n');
      fprint(parameters.fp_log, "        Spaced  No\n");
      fprint(parameters.fp_log, "        Hashed  No\n");
      fprint(parameters.fp_log, "         Coded  No\n");
      fprint(parameters.fp_log, "       Stepped  No\n");
      fprint(parameters.fp_log, "         Slots  ");
      fprint_integer(parameters.fp_log, dbindex.hashsize);
      fprint(parameters.fp_log, " (");
      std::fprintf(parameters.fp_log, "%.1f", 1.0 * dbindex.hashsize / 1000.0);
      fprint(parameters.fp_log, "k)\n");
      fprint(parameters.fp_log, "       DBAccel  ");
      fprint_integer(parameters.fp_log, dbindex.dbaccel);
      fprint(parameters.fp_log, "%\n");
      fprint(parameters.fp_log, '\n');

      fprint_integer(parameters.fp_log, nt, 10);
      fprint(parameters.fp_log, "  DB size (");
      std::fprintf(parameters.fp_log, "%.1f", 1.0 * static_cast<double>(nt) / 1000.0);
      fprint(parameters.fp_log, "k)\n");
      fprint_integer(parameters.fp_log, dbindex.indexsize, 10);
      fprint(parameters.fp_log, "  Words\n");
      fprint_integer(parameters.fp_log, wcmedian, 10);
      fprint(parameters.fp_log, "  Median size\n");
      std::fprintf(parameters.fp_log, "%10.1f", 1.0 * static_cast<double>(dbindex.indexsize) / dbindex.hashsize);
      fprint(parameters.fp_log, "  Mean size\n");
      fprint(parameters.fp_log, '\n');

      fprint(parameters.fp_log, "     iWord         sWord         Cap        Size  Row\n");
      fprint(parameters.fp_log, "----------  ------------  ----------  ----------  ---\n");

      /* The matching sequence numbers the rows below show. udb_read() did not
         keep the word list (UdbUse::word_stats) nor the per-k-mer offset
         table, so both are recovered here for these rows alone: the offsets by
         one pass over the counts, and the entries by reading the file at those
         offsets. The rows are visited in report order, but the counts must be
         walked in k-mer order, so the requests are sorted by k-mer first and
         the results land back in row order. */

      std::vector<unsigned int> row_entries(top.size() * reported_entries_per_row, 0U);
      std::vector<std::size_t> row_shown(top.size(), 0);
      {
        std::vector<std::size_t> by_kmer;
        by_kmer.reserve(reported_rows);
        for (std::size_t row = 0; row < top.size(); ++row)
          {
            by_kmer.push_back(row);
          }
        auto const kmer_of_row = [&top](std::size_t const row) -> unsigned int
        { return top[row].kmer; };
        std::sort(by_kmer.begin(), by_kmer.end(),
                  [&kmer_of_row](std::size_t const lhs, std::size_t const rhs) -> bool
                  { return kmer_of_row(lhs) < kmer_of_row(rhs); });

        uint64_t running = 0;
        std::size_t next = 0;
        for (auto kmer = 0U; (kmer < dbindex.hashsize) and (next < by_kmer.size()); ++kmer)
          {
            while ((next < by_kmer.size()) and (kmer_of_row(by_kmer[next]) == kmer))
              {
                auto const row = by_kmer[next];
                auto const count = static_cast<std::size_t>(dbindex.kmercount[kmer]);
                row_shown[row] = std::min(count, reported_entries_per_row);
                udb_read_word_entries(parameters.input_filename,
                                      dbindex.wordlength,
                                      seqcount,
                                      running,
                                      make_span(row_entries)
                                        .subspan(row * reported_entries_per_row,
                                                 row_shown[row]));
                ++next;
              }
            running += dbindex.kmercount[kmer];
          }
        assert(next == by_kmer.size());
      }

      for (std::size_t i = 0; i < top.size(); ++i)
        {
          fprint_integer(parameters.fp_log, top[i].kmer, 10);
          fprint(parameters.fp_log, "  ");

          /* pad the k-mer column out to 12 characters */
          static constexpr char twelve_spaces[] = "            ";
          auto const padding = static_cast<std::size_t>(
            std::max(12 - static_cast<int>(dbindex.wordlength), 0));
          fprint(parameters.fp_log, View<char>{twelve_spaces, padding});

          fprint_kmer(parameters.fp_log, dbindex.wordlength, top[i].kmer);

          fprint(parameters.fp_log, "  ");
          fprint_integer(parameters.fp_log, 0U, 10);
          fprint(parameters.fp_log, "  ");
          fprint_integer(parameters.fp_log, top[i].count, 10);

          fprint(parameters.fp_log, ' ');

          for (std::size_t j = 0; j < row_shown[i]; ++j)
            {
              fprint(parameters.fp_log, ' ');
              fprint_integer(parameters.fp_log, row_entries[(i * reported_entries_per_row) + j]);
            }


          if (top[i].count > reported_entries_per_row)
            {
              fprint(parameters.fp_log, "...");
            }

          fprint(parameters.fp_log, '\n');
        }

      fprint(parameters.fp_log, "\n\n");

      fprint(parameters.fp_log, "Word width  ");
      fprint_integer(parameters.fp_log, dbindex.wordlength);
      fprint(parameters.fp_log, '\n');
      fprint(parameters.fp_log, "Slots       ");
      fprint_integer(parameters.fp_log, dbindex.hashsize);
      fprint(parameters.fp_log, '\n');
      fprint(parameters.fp_log, "Words       ");
      fprint_integer(parameters.fp_log, dbindex.indexsize);
      fprint(parameters.fp_log, '\n');
      fprint(parameters.fp_log, "Max size    ");
      fprint_integer(parameters.fp_log, wcmax);
      fprint(parameters.fp_log, " (");
      fprint_kmer(parameters.fp_log, dbindex.wordlength, wcmax_kmer);
      fprint(parameters.fp_log, ")\n\n");

      fprint(parameters.fp_log, "   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct\n");
      fprint(parameters.fp_log, "----------  ----------  ----------  ----------  ------  ------\n");


      auto size_lo = 0U;
      auto size_hi = 0U;
      auto x = 0U;
      auto totpct = 0.0;

      while (size_lo < seqcount)
        {

          /* the slots whose count falls in this bucket, read straight off the
             distribution: x is the lowest count not yet reported */
          uint64_t count = 0;
          uint64_t size = 0;
          for (auto value = x; value <= size_hi; ++value)
            {
              count += histogram[value];
              size += value * histogram[value];
            }
          x = size_hi + 1;

          auto const pct = 100.0 * static_cast<double>(count) / dbindex.hashsize;
          totpct += pct;

          if (size_lo < size_hi)
            {
              fprint_integer(parameters.fp_log, size_lo, 10);
            }
          else
            {
              fprint(parameters.fp_log, "          ");
            }

          fprint(parameters.fp_log, "  ");
          fprint_integer(parameters.fp_log, size_hi, 10);

          if (size >= 10000)
            {
              fprint(parameters.fp_log, "  ");
              std::fprintf(parameters.fp_log, "%9.1f", static_cast<double>(size) * 0.001);
              fprint(parameters.fp_log, 'k');
            }
          else
            {
              fprint(parameters.fp_log, "  ");
              std::fprintf(parameters.fp_log, "%10.1f", static_cast<double>(size));
            }

          if (count >= 10000)
            {
              fprint(parameters.fp_log, "  ");
              std::fprintf(parameters.fp_log, "%9.1f", static_cast<double>(count) * 0.001);
              fprint(parameters.fp_log, 'k');
            }
          else
            {
              fprint(parameters.fp_log, "  ");
              std::fprintf(parameters.fp_log, "%10.1f", static_cast<double>(count));
            }

          fprint(parameters.fp_log, "  ");
          std::fprintf(parameters.fp_log, "%5.1f", pct);
          fprint(parameters.fp_log, "%  ");
          std::fprintf(parameters.fp_log, "%5.1f", totpct);
          fprint(parameters.fp_log, '%');

          static constexpr auto divider = 3.0;
          const auto dots = std::lround(pct / divider);

          if (dots > 0)
            {
              fprint(parameters.fp_log, "  ");
            }

          for (auto i = 0L; i < dots ; i++)
            {
              fprint(parameters.fp_log, '*');
            }

          fprint(parameters.fp_log, '\n');

          size_lo = size_hi + 1;
          if (size_hi > 0)
            {
              size_hi *= 2;
            }
          else
            {
              size_hi = 1;
            }
          size_hi = std::min(size_hi, seqcount);
        }

      fprint(parameters.fp_log, "----------  ----------  ----------  ----------\n");
      fprint(parameters.fp_log, "                      ");

      if (dbindex.indexsize >= 10000)
        {
          fprint(parameters.fp_log, "  ");
          std::fprintf(parameters.fp_log, "%9.1f", static_cast<double>(dbindex.indexsize) * 0.001);
          fprint(parameters.fp_log, 'k');
        }
      else
        {
          fprint(parameters.fp_log, "  ");
          std::fprintf(parameters.fp_log, "%10.1f", static_cast<double>(dbindex.indexsize) * 1.0);
        }

      if (dbindex.hashsize >= 10000)
        {
          fprint(parameters.fp_log, "  ");
          std::fprintf(parameters.fp_log, "%9.1f", dbindex.hashsize * 0.001);
          fprint(parameters.fp_log, 'k');
        }
      else
        {
          fprint(parameters.fp_log, "  ");
          std::fprintf(parameters.fp_log, "%10.1f", dbindex.hashsize * 1.0);
        }

      fprint(parameters.fp_log, "\n\n");

      fprint_integer(parameters.fp_log, nt, 10);
      fprint(parameters.fp_log, "  Upper\n");
      fprint_integer(parameters.fp_log, 0U, 10);
      fprint(parameters.fp_log, "  Lower (");
      std::fprintf(parameters.fp_log, "%.1f", 0.0);
      fprint(parameters.fp_log, "%)\n");
      fprint_integer(parameters.fp_log, nt, 10);
      fprint(parameters.fp_log, "  Total\n");
      fprint_integer(parameters.fp_log, dbindex.indexsize, 10);
      fprint(parameters.fp_log, "  Indexed words\n");
    }

  dbindex.clear();
  db.clear();
}
