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
#include "commands/fastx_mask.hpp"
#include "core/db.hpp"
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/mask.hpp"
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"
#include "utils/progress.hpp"
#include <cctype>  // std::isupper
#include <cstdint>  // uint64_t
#include <cstdio>  // std::fprintf


auto fastx_mask(struct Parameters const & parameters) -> void
{
  if ((parameters.opt_fastaout == nullptr) && (parameters.opt_fastqout == nullptr)) {
    fatal("Specify output files for masking with --fastaout and/or --fastqout");
  }

  auto fp_fastaout = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  auto fp_fastqout = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});

  Database db;
  db.read(parameters.opt_fastx_mask, 0, parameters);
  // memory-intensive: the entire database is now held in memory

  if ((fp_fastqout != nullptr) && ! db.is_fastq())
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  uint64_t const seqcount = db.getsequencecount();

  if (parameters.opt_qmask == Masking::dust)
    {
      dust_all(db, parameters);
    }
  else if ((parameters.opt_qmask == Masking::soft) && parameters.opt_hardmask)
    {
      hardmask_all(db);
    }

  auto kept = 0;
  auto discarded_less = 0;
  auto discarded_more = 0;
  {
    Progress progress("Writing output", seqcount, parameters);
    for (uint64_t i = 0; i < seqcount; i++)
      {
        auto unmasked = 0;
        auto const * seq = db.getsequence(i);
        const int len = static_cast<int>(db.getsequencelen(i));
        if (parameters.opt_qmask == Masking::none)
          {
            unmasked = len;
          }
        else if (parameters.opt_hardmask)
          {
            for (auto j = 0; j < len; j++)
              {
                if (seq[j] != 'N')
                  {
                    ++unmasked;
                  }
              }
          }
        else
          {
            for (auto j = 0; j < len; j++)
              {
                if (std::isupper(seq[j]) != 0)
                  {
                    ++unmasked;
                  }
              }
          }
        auto const unmasked_pct = 100.0 * unmasked / len;

        if (unmasked_pct < parameters.opt_min_unmasked_pct)
          {
            ++discarded_less;
          }
        else if (unmasked_pct >  parameters.opt_max_unmasked_pct)
          {
            ++discarded_more;
          }
        else
          {
            ++kept;

            if (parameters.opt_fastaout != nullptr)
              {
                fasta_print_general(fp_fastaout.get(),
                                    nullptr,
                                    seq,
                                    len,
                                    db.getheader(i),
                                    static_cast<int>(db.getheaderlen(i)),
                                    db.getabundance(i),
                                    kept,
                                    -1.0,
                                    -1, -1, nullptr, 0.0,
                                    0,
                                    parameters);
              }

            if (parameters.opt_fastqout != nullptr)
              {
                fastq_print_general(fp_fastqout.get(),
                                    seq,
                                    len,
                                    db.getheader(i),
                                    static_cast<int>(db.getheaderlen(i)),
                                    db.getquality(i),
                                    db.getabundance(i),
                                    kept,
                                    -1.0,
                                    parameters);
              }
          }

        progress.update(i);
      }
  }

  if (! parameters.opt_quiet)
    {
      if (parameters.opt_min_unmasked_pct > 0.0)
        {
          std::fprintf(stderr, "%d sequences with less than %.1lf%% unmasked residues discarded\n", discarded_less, parameters.opt_min_unmasked_pct);
        }
      if (parameters.opt_max_unmasked_pct < 100.0)
        {
          std::fprintf(stderr, "%d sequences with more than %.1lf%% unmasked residues discarded\n", discarded_more, parameters.opt_max_unmasked_pct);
        }
      std::fprintf(stderr, "%d sequences kept\n", kept);
    }

  if (parameters.opt_log != nullptr)
    {
      if (parameters.opt_min_unmasked_pct > 0.0)
        {
          std::fprintf(parameters.fp_log, "%d sequences with less than %.1lf%% unmasked residues discarded\n", discarded_less, parameters.opt_min_unmasked_pct);
        }
      if (parameters.opt_max_unmasked_pct < 100.0)
        {
          std::fprintf(parameters.fp_log, "%d sequences with more than %.1lf%% unmasked residues discarded\n", discarded_more, parameters.opt_max_unmasked_pct);
        }
      std::fprintf(parameters.fp_log, "%d sequences kept\n", kept);
    }

  db.clear();

  /* reset() in the original fclose order (fastaout before fastqout); a no-op on
     an unopened handle, and scope-exit would reverse the flush order for the
     degenerate case of both outputs sharing stdout. */
  fp_fastaout.reset();
  fp_fastqout.reset();
}
