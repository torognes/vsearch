/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <algorithm>  // std::count_if
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::floor
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <functional>  // std::minus
#include <numeric>  // std::fill
#include <vector>

// refactoring:
// - accept sample_size = 0 and sample_pct = 0.0?
// - fastaout should be empty, all reads should be in fastaout_discarded

// refactoring:
// - store parameters in a struct, pass struct to reduce list of arguments
// - use uint64_t for abundance values?
// - create discarded vector only if needed


namespace {
  // anonymous namespace to avoid linker error (multiple definitions
  // of function with identical names and parameters)
  auto create_deck(bool const opt_sizein) -> std::vector<int> {
    auto const dbsequencecount = db_getsequencecount();
    std::vector<int> deck(dbsequencecount, 1);
    if (opt_sizein) {
      auto counter = std::size_t{0};
      for (auto & abundance : deck) {
        abundance = db_getabundance(counter);
        ++counter;
      }
    }
    return deck;
  }
}


auto number_of_reads_to_sample(int64_t const opt_sample_size,
                               double const opt_sample_pct,
                               uint64_t const mass_total) -> uint64_t {
  // assert(mass_total < max_uint64 / opt_sample_pct)
  if (opt_sample_size != 0) {
    return static_cast<uint64_t>(opt_sample_size);
  }
  return static_cast<uint64_t>(std::floor(mass_total * opt_sample_pct / 100.0));
}


auto random_subsampling(std::vector<int> & deck, uint64_t const mass_total,
                        uint64_t const n_reads) -> void {
  auto n_reads_left = n_reads;
  auto amplicon_number = 0;
  uint64_t n_read_being_checked = 0;
  uint64_t accumulated_mass = 0;
  auto amplicon_mass = opt_sizein ? db_getabundance(0) : 1;

  // refactoring C++17: std::sample()
  progress_init("Subsampling", mass_total);
  while (n_reads_left > 0)
    {
      auto const random = random_ulong(mass_total - n_read_being_checked);

      if (random < n_reads_left)
        {
          /* selected read r from amplicon a */
          ++deck[amplicon_number];
          --n_reads_left;
        }

      ++n_read_being_checked;
      ++accumulated_mass;
      if (accumulated_mass >= amplicon_mass)
        {
          /* next amplicon */
          ++amplicon_number;
          amplicon_mass = opt_sizein ? db_getabundance(amplicon_number) : 1;
          accumulated_mass = 0;
        }
      progress_update(n_read_being_checked);
    }
  progress_done();
}


auto writing_fasta_output(std::vector<int> const & deck,
                          char * ptr_fasta_file_name,
                          std::FILE * ptr_fasta_file) -> void {
  if (ptr_fasta_file_name == nullptr) {
    return;
  }
  int amplicons_printed = 0;
  progress_init("Writing fasta output", deck.size());
  auto counter = 0U;
  for (auto const abundance_value : deck) {
    int64_t const new_abundance = abundance_value;
      if (new_abundance == 0) {
        ++counter;
        continue;
      }
      ++amplicons_printed;
      fasta_print_general(ptr_fasta_file,
                          nullptr,
                          db_getsequence(counter),
                          db_getsequencelen(counter),
                          db_getheader(counter),
                          db_getheaderlen(counter),
                          new_abundance,
                          amplicons_printed,
                          -1.0,
                          -1, -1, nullptr, 0.0);
      progress_update(counter);
      ++counter;
    }
  progress_done();
}


auto writing_fastq_output(std::vector<int> const & deck,
                          char * ptr_fastq_file_name,
                          std::FILE * ptr_fastq_file) -> void {
  if (ptr_fastq_file_name == nullptr) {
    return;
  }
  int amplicons_printed = 0;
  progress_init("Writing fastq output", deck.size());
  auto counter = 0U;
  for (auto const abundance_value : deck) {
    int64_t const new_abundance = abundance_value;
      if (new_abundance == 0) {
        ++counter;
        continue;
      }
      ++amplicons_printed;
      fastq_print_general(ptr_fastq_file,
                          db_getsequence(counter),
                          db_getsequencelen(counter),
                          db_getheader(counter),
                          db_getheaderlen(counter),
                          db_getquality(counter),
                          new_abundance,
                          amplicons_printed,
                          -1.0);
      progress_update(counter);
      ++counter;
  }
  progress_done();
}


auto subsample() -> void
{
  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastaout_discarded = nullptr;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_fastqout_discarded = nullptr;

  if (opt_fastaout != nullptr)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (fp_fastaout == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastaout_discarded != nullptr)
    {
      fp_fastaout_discarded = fopen_output(opt_fastaout_discarded);
      if (fp_fastaout_discarded == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout != nullptr)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (fp_fastqout == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (opt_fastqout_discarded != nullptr)
    {
      fp_fastqout_discarded = fopen_output(opt_fastqout_discarded);
      if (fp_fastqout_discarded == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  db_read(opt_fastx_subsample, 0);
  show_rusage();

  if ((fp_fastqout != nullptr or fp_fastqout_discarded != nullptr) and not db_is_fastq())
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  auto original_abundances = create_deck(opt_sizein);
  auto const mass_total = std::accumulate(original_abundances.cbegin(), original_abundances.cend(), uint64_t{0});
  auto subsampled_abundances = original_abundances;
  std::fill(subsampled_abundances.begin(), subsampled_abundances.end(), 0);  // temporary fix: reset vector to zero

  if (not opt_quiet)
    {
      std::fprintf(stderr, "Got %" PRIu64 " reads from %d amplicons\n",
                   mass_total, static_cast<int>(original_abundances.size()));
    }

  if (opt_log != nullptr)
    {
      std::fprintf(fp_log, "Got %" PRIu64 " reads from %d amplicons\n",
                   mass_total, static_cast<int>(original_abundances.size()));
    }

  auto const n_reads = number_of_reads_to_sample(opt_sample_size, opt_sample_pct, mass_total);

  if (n_reads > mass_total)
    {
      fatal("Cannot subsample more reads than in the original sample");
    }

  random_subsampling(subsampled_abundances, mass_total, n_reads);

  writing_fasta_output(subsampled_abundances, opt_fastaout, fp_fastaout);
  writing_fastq_output(subsampled_abundances, opt_fastqout, fp_fastqout);

  // refactoring: extract to a function, make discarded_abundances const
  std::vector<int> discarded_abundances(original_abundances.size());
  std::transform(original_abundances.cbegin(), original_abundances.cend(),
                 subsampled_abundances.cbegin(), discarded_abundances.begin(),
                 std::minus<int>());

  writing_fasta_output(discarded_abundances, opt_fastaout_discarded, fp_fastaout_discarded);
  writing_fastq_output(discarded_abundances, opt_fastqout_discarded, fp_fastqout_discarded);

  int const samples = std::count_if(subsampled_abundances.cbegin(),
                                    subsampled_abundances.cend(), [](int abundance) { return abundance != 0; });
  if (not opt_quiet)
    {
      std::fprintf(stderr, "Subsampled %" PRIu64 " reads from %d amplicons\n", n_reads, samples);
    }
  if (opt_log != nullptr)
    {
      std::fprintf(fp_log, "Subsampled %" PRIu64 " reads from %d amplicons\n", n_reads, samples);
    }

  db_free();

  if ((opt_fastaout != nullptr) and (fp_fastaout != nullptr))
    {
      static_cast<void>(std::fclose(fp_fastaout));
    }

  if ((opt_fastqout != nullptr) and (fp_fastqout != nullptr))
    {
      static_cast<void>(std::fclose(fp_fastqout));
    }

  if ((opt_fastaout_discarded != nullptr) and (fp_fastaout_discarded != nullptr))
    {
      static_cast<void>(std::fclose(fp_fastaout_discarded));
    }

  if ((opt_fastqout_discarded != nullptr) and (fp_fastqout_discarded != nullptr))
    {
      static_cast<void>(std::fclose(fp_fastqout_discarded));
    }
}
