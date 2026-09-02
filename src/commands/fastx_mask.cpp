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
#include "commands/fastx_mask.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/db.hpp"
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/mask.hpp"
#include "utils/ascii_case.hpp"  // is_upper
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include <algorithm>  // std::count_if
#include <cstdint>  // uint64_t
#include <cstdio>  // std::fprintf


namespace {

  /* The three counts travel together and are all the same type, so they are
     named in a struct rather than passed as three adjacent swappable
     arguments. No default member initializers: they would make this a
     non-aggregate before C++14, and the call site brace-initializes it. */
  struct MaskCounts
  {
    int kept;
    int discarded_less;
    int discarded_more;
  };


  /* Both option-driven branches (--min_unmasked_pct, --max_unmasked_pct) stay
     inside the writer, which is what made this report 18 lines twice over. */
  auto print_mask_summary(std::FILE * output_stream,
                          struct Parameters const & parameters,
                          MaskCounts const & counts) -> void
  {
    if (parameters.opt_min_unmasked_pct > 0.0)
      {
        fprint_integer(output_stream, counts.discarded_less);
        fprint(output_stream, " sequences with less than ");
        std::fprintf(output_stream, "%.1lf", parameters.opt_min_unmasked_pct);
        fprint(output_stream, "% unmasked residues discarded\n");
      }
    if (parameters.opt_max_unmasked_pct < 100.0)
      {
        fprint_integer(output_stream, counts.discarded_more);
        fprint(output_stream, " sequences with more than ");
        std::fprintf(output_stream, "%.1lf", parameters.opt_max_unmasked_pct);
        fprint(output_stream, "% unmasked residues discarded\n");
      }
    fprint_integer(output_stream, counts.kept);
    fprint(output_stream, " sequences kept\n");
  }

}  // anonymous namespace


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
        auto const seq = db.sequence_view(i);
        const int len = static_cast<int>(seq.size());
        if (parameters.opt_qmask == Masking::none)
          {
            unmasked = len;
          }
        else if (parameters.opt_hardmask)
          {
            unmasked = static_cast<int>(std::count_if(seq.cbegin(), seq.cend(),
                                                      [](char const nucleotide) -> bool
                                                      { return nucleotide != 'N'; }));
          }
        else
          {
            unmasked = static_cast<int>(std::count_if(seq.cbegin(), seq.cend(),
                                                      [](char const nucleotide) -> bool
                                                      { return is_upper(nucleotide); }));
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
                                    db.record(i),
                                    OutputAnnotations{db.getabundance(i), kept},
                                    parameters);
              }

            if (parameters.opt_fastqout != nullptr)
              {
                fastq_print_general(fp_fastqout.get(),
                                    db.record(i),
                                    OutputAnnotations{db.getabundance(i), kept},
                                    parameters);
              }
          }

        progress.update(i);
      }
  }

  auto const summary = MaskCounts{kept, discarded_less, discarded_more};

  if (! parameters.opt_quiet)
    {
      print_mask_summary(stderr, parameters, summary);
    }

  if (parameters.fp_log != nullptr)
    {
      print_mask_summary(parameters.fp_log, parameters, summary);
    }

  db.clear();

}
