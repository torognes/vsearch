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


  auto wc_compare(wordfreq_t const & lhs, wordfreq_t const & rhs) -> bool
  {
    if (lhs.count != rhs.count)
      {
        return lhs.count < rhs.count;
      }
    return lhs.kmer > rhs.kmer;
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

  /* analyze word counts */

  std::vector<wordfreq_t> freqtable(dbindex.hashsize);

  for (auto i = 0U; i < dbindex.hashsize; i++)
    {
      freqtable[i].kmer = i;
      freqtable[i].count = dbindex.kmercount[i];
    }

  std::sort(freqtable.begin(), freqtable.end(), wc_compare);

  auto const wcmax = freqtable[dbindex.hashsize-1].count;
  auto const wcmedian = ( freqtable[(dbindex.hashsize / 2) - 1].count +
                            freqtable[dbindex.hashsize / 2].count ) / 2;

  auto const seqcount = static_cast<unsigned int>(db.getsequencecount());
  auto const nt = db.getnucleotidecount();

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

      std::vector<unsigned int> row_entries(reported_rows * reported_entries_per_row, 0U);
      std::vector<std::size_t> row_shown(reported_rows, 0);
      {
        std::vector<std::size_t> by_kmer;
        by_kmer.reserve(reported_rows);
        for (std::size_t row = 0; row < reported_rows; ++row)
          {
            by_kmer.push_back(row);
          }
        auto const kmer_of_row = [&freqtable, &dbindex](std::size_t const row) -> unsigned int
        { return freqtable[dbindex.hashsize - 1 - row].kmer; };
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

      for (auto i = 0U; i < dbindex.hashsize; i++)
        {
          fprint_integer(parameters.fp_log, freqtable[dbindex.hashsize - 1 - i].kmer, 10);
          fprint(parameters.fp_log, "  ");

          /* pad the k-mer column out to 12 characters */
          static constexpr char twelve_spaces[] = "            ";
          auto const padding = static_cast<std::size_t>(
            std::max(12 - static_cast<int>(dbindex.wordlength), 0));
          fprint(parameters.fp_log, View<char>{twelve_spaces, padding});

          fprint_kmer(parameters.fp_log, dbindex.wordlength, freqtable[dbindex.hashsize - 1 - i].kmer);

          fprint(parameters.fp_log, "  ");
          fprint_integer(parameters.fp_log, 0U, 10);
          fprint(parameters.fp_log, "  ");
          fprint_integer(parameters.fp_log, freqtable[dbindex.hashsize - 1 - i].count, 10);

          fprint(parameters.fp_log, ' ');

          for (std::size_t j = 0; j < row_shown[i]; ++j)
            {
              fprint(parameters.fp_log, ' ');
              fprint_integer(parameters.fp_log, row_entries[(i * reported_entries_per_row) + j]);
            }


          if (freqtable[dbindex.hashsize-1-i].count > reported_entries_per_row)
            {
              fprint(parameters.fp_log, "...");
            }

          fprint(parameters.fp_log, '\n');

          if (i + 1 == reported_rows)
            {
              break;
            }
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
      fprint_kmer(parameters.fp_log, dbindex.wordlength, freqtable[dbindex.hashsize - 1].kmer);
      fprint(parameters.fp_log, ")\n\n");

      fprint(parameters.fp_log, "   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct\n");
      fprint(parameters.fp_log, "----------  ----------  ----------  ----------  ------  ------\n");


      auto size_lo = 0U;
      auto size_hi = 0U;
      auto x = 0U;
      auto totpct = 0.0;

      while (size_lo < seqcount)
        {

          auto count = 0;
          auto size = 0U;
          while ((x < dbindex.hashsize) and (freqtable[x].count <= size_hi))
            {
              count++;
              size += freqtable[x].count;
              x++;
            }

          auto const pct = 100.0 * count / dbindex.hashsize;
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
              std::fprintf(parameters.fp_log, "%9.1f", size * 0.001);
              fprint(parameters.fp_log, 'k');
            }
          else
            {
              fprint(parameters.fp_log, "  ");
              std::fprintf(parameters.fp_log, "%10.1f", size * 1.0);
            }

          if (count >= 10000)
            {
              fprint(parameters.fp_log, "  ");
              std::fprintf(parameters.fp_log, "%9.1f", count * 0.001);
              fprint(parameters.fp_log, 'k');
            }
          else
            {
              fprint(parameters.fp_log, "  ");
              std::fprintf(parameters.fp_log, "%10.1f", count * 1.0);
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
