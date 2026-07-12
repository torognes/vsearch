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
#include "commands/udbstats.hpp"
#include "core/udb.hpp"
#include "dbindex.h"
#include <algorithm>  // std::max, std::min
#include <cinttypes>  // macro PRIu64
#include <cmath>  // std::lround
#include <cstdint>  // uint64_t
#include <cstdio>  // std::fprintf
#include <cstdlib>  // std::qsort
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  struct wordfreq
  {
    unsigned int kmer;
    unsigned int count;
  };

  using wordfreq_t = struct wordfreq;


  auto wc_compare(const void * a, const void * b) -> int
  {
    auto const * lhs = static_cast<wordfreq_t const *>(a);
    auto const * rhs = static_cast<wordfreq_t const *>(b);
    if (lhs->count < rhs->count)
      {
        return -1;
      }
    if (lhs->count > rhs->count)
      {
        return +1;
      }

    if (lhs->kmer < rhs->kmer)
      {
        return +1;
      }
    if (lhs->kmer > rhs->kmer)
      {
        return -1;
      }
    return 0;
  }

}  // end of anonymous namespace


auto udb_stats(struct Parameters const & parameters) -> void
{
  /* show word statistics for an UDB file */

  Dbindex dbindex;  /* the k-mer index this run owns (RAII) */

  /* read UDB file */

  udb_read(parameters.opt_udbstats, false, false, dbindex, parameters);

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

  std::qsort(freqtable.data(), dbindex.hashsize, sizeof(wordfreq_t), wc_compare);

  auto const wcmax = freqtable[dbindex.hashsize-1].count;
  auto const wcmedian = ( freqtable[(dbindex.hashsize / 2) - 1].count +
                            freqtable[dbindex.hashsize / 2].count ) / 2;

  unsigned int const seqcount = static_cast<unsigned int>(db_getsequencecount());
  auto const nt = db_getnucleotidecount();

  /* show stats */

  if (parameters.opt_log != nullptr)
    {
      std::fprintf(parameters.fp_log, "      Alphabet  nt\n");
      std::fprintf(parameters.fp_log, "    Word width  %u\n", dbindex.wordlength);
      std::fprintf(parameters.fp_log, "     Word ones  %u\n", dbindex.wordlength);
      std::fprintf(parameters.fp_log, "        Spaced  No\n");
      std::fprintf(parameters.fp_log, "        Hashed  No\n");
      std::fprintf(parameters.fp_log, "         Coded  No\n");
      std::fprintf(parameters.fp_log, "       Stepped  No\n");
      std::fprintf(parameters.fp_log,
              "         Slots  %u (%.1fk)\n",
              dbindex.hashsize,
              1.0 * dbindex.hashsize / 1000.0);
      std::fprintf(parameters.fp_log, "       DBAccel  %u%%\n", dbindex.dbaccel);
      std::fprintf(parameters.fp_log, "\n");

      std::fprintf(parameters.fp_log,
              "%10" PRIu64 "  DB size (%.1fk)\n",
              nt,
              1.0 * static_cast<double>(nt) / 1000.0);
      std::fprintf(parameters.fp_log, "%10" PRIu64 "  Words\n", dbindex.indexsize);
      std::fprintf(parameters.fp_log, "%10u  Median size\n", wcmedian);
      std::fprintf(parameters.fp_log,
              "%10.1f  Mean size\n",
              1.0 * static_cast<double>(dbindex.indexsize) / dbindex.hashsize);
      std::fprintf(parameters.fp_log, "\n");

      std::fprintf(parameters.fp_log,
              "     iWord         sWord         Cap        Size  Row\n");
      std::fprintf(parameters.fp_log,
              "----------  ------------  ----------  ----------  ---\n");

      for (auto i = 0U; i < dbindex.hashsize; i++)
        {
          std::fprintf(parameters.fp_log,
                  "%10u  ",
                  freqtable[dbindex.hashsize - 1 - i].kmer);

          std::fprintf(parameters.fp_log,
                  "%.*s", std::max(12 - static_cast<int>(dbindex.wordlength), 0), "            ");

          fprint_kmer(parameters.fp_log, dbindex.wordlength, freqtable[dbindex.hashsize - 1 - i].kmer);

          std::fprintf(parameters.fp_log,
                  "  %10u  %10u",
                  0U,
                  freqtable[dbindex.hashsize - 1 - i].count);

          std::fprintf(parameters.fp_log, " ");

          for (auto j = 0U; j < freqtable[dbindex.hashsize - 1 - i].count; j++)
            {
              std::fprintf(parameters.fp_log,
                      " %u", dbindex.kmerindex[dbindex.kmerhash[freqtable[dbindex.hashsize - 1 - i].kmer] + j]);

              if (j == 7)
                {
                  break;
                }
            }


          if (freqtable[dbindex.hashsize-1-i].count > 8)
            {
              std::fprintf(parameters.fp_log, "...");
            }

          std::fprintf(parameters.fp_log, "\n");

          if (i == 10)
            {
              break;
            }
        }

      std::fprintf(parameters.fp_log, "\n\n");

      std::fprintf(parameters.fp_log, "Word width  %u\n", dbindex.wordlength);
      std::fprintf(parameters.fp_log, "Slots       %u\n", dbindex.hashsize);
      std::fprintf(parameters.fp_log, "Words       %" PRIu64 "\n", dbindex.indexsize);
      std::fprintf(parameters.fp_log, "Max size    %u (", wcmax);
      fprint_kmer(parameters.fp_log, dbindex.wordlength, freqtable[dbindex.hashsize - 1].kmer);
      std::fprintf(parameters.fp_log, ")\n\n");

      std::fprintf(parameters.fp_log, "   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct\n");
      std::fprintf(parameters.fp_log, "----------  ----------  ----------  ----------  ------  ------\n");


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
              std::fprintf(parameters.fp_log, "%10u", size_lo);
            }
          else
            {
              std::fprintf(parameters.fp_log, "          ");
            }

          std::fprintf(parameters.fp_log, "  %10u", size_hi);

          if (size >= 10000)
            {
              std::fprintf(parameters.fp_log, "  %9.1fk", size * 0.001);
            }
          else
            {
              std::fprintf(parameters.fp_log, "  %10.1f", size * 1.0);
            }

          if (count >= 10000)
            {
              std::fprintf(parameters.fp_log, "  %9.1fk", count * 0.001);
            }
          else
            {
              std::fprintf(parameters.fp_log, "  %10.1f", count * 1.0);
            }

          std::fprintf(parameters.fp_log, "  %5.1f%%  %5.1f%%", pct, totpct);

          static constexpr auto divider = 3.0;
          const auto dots = std::lround(pct / divider);

          if (dots > 0)
            {
              std::fprintf(parameters.fp_log, "  ");
            }

          for (auto i = 0L; i < dots ; i++)
            {
              std::fprintf(parameters.fp_log, "*");
            }

          std::fprintf(parameters.fp_log, "\n");

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

      std::fprintf(parameters.fp_log, "----------  ----------  ----------  ----------\n");
      std::fprintf(parameters.fp_log, "                      ");

      if (dbindex.indexsize >= 10000)
        {
          std::fprintf(parameters.fp_log, "  %9.1fk", static_cast<double>(dbindex.indexsize) * 0.001);
        }
      else
        {
          std::fprintf(parameters.fp_log, "  %10.1f", static_cast<double>(dbindex.indexsize) * 1.0);
        }

      if (dbindex.hashsize >= 10000)
        {
          std::fprintf(parameters.fp_log, "  %9.1fk", dbindex.hashsize * 0.001);
        }
      else
        {
          std::fprintf(parameters.fp_log, "  %10.1f", dbindex.hashsize * 1.0);
        }

      std::fprintf(parameters.fp_log, "\n\n");

      std::fprintf(parameters.fp_log, "%10" PRIu64 "  Upper\n", nt);
      std::fprintf(parameters.fp_log, "%10u  Lower (%.1f%%)\n", 0U, 0.0);
      std::fprintf(parameters.fp_log, "%10" PRIu64 "  Total\n", nt);
      std::fprintf(parameters.fp_log, "%10" PRIu64 "  Indexed words\n", dbindex.indexsize);
    }

  dbindex.clear();
  db_free();
}
