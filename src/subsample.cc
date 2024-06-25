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
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::floor
#include <cstdint>  // int64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <functional>  // std::minus
#include <numeric>  // std::fill
#include <vector>


#ifndef NDEBUG
// all contiguous integers from 0 to 2^53 can be represented in the
// mantissa of a double
constexpr uint64_t contiguous_mantissa = 9007199254740992;  // 9 x 10^15 reads
#endif

struct a_file {
  char * name = nullptr;
  std::FILE * handle = nullptr;
};


struct file_purposes {
  a_file kept;
  a_file lost;
};


struct file_types {
  file_purposes fasta;
  file_purposes fastq;
};

// refactoring:
// - accept sample_size = 0 and sample_pct = 0.0?
// - fastaout should be empty, all reads should be in fastaout_discarded

// refactoring:
// - store parameters in a struct, pass struct to reduce list of arguments
// - use uint64_t for abundance values?
// - create discarded vector only if needed


// // discrete_distribution
// #include <iostream>
// #include <random>
// #include <vector>

// int main()
// {
//   const auto nreads = 500ULL; // target number of reads

//   std::random_device r;
//   std::default_random_engine generator(r());
//   std::vector<int> v1 = {1000, 100, 10, 1};
//   std::discrete_distribution<int> distribution(v1.begin(), v1.end());

//   std::vector<int> v2(v1.size());

//   for (auto i=0ULL; i<nreads; ++i) {
//     int const number = distribution(generator);
//     ++v2[number];
//   }

//   for (auto i=0UL; i<v2.size(); ++i)
//     std::cout << i << ": " << v1[i] << " " << v2[i] << "\n";

//   return 0;
// }


auto open_output_files(struct file_types & ouput_files) -> void {
  if (ouput_files.fasta.kept.name != nullptr) {
    ouput_files.fasta.kept.handle = fopen_output(ouput_files.fasta.kept.name);
  }
  if (ouput_files.fasta.lost.name != nullptr) {
    ouput_files.fasta.lost.handle = fopen_output(ouput_files.fasta.lost.name);
  }
  if (ouput_files.fastq.kept.name != nullptr) {
    ouput_files.fastq.kept.handle = fopen_output(ouput_files.fastq.kept.name);
  }
  if (ouput_files.fastq.lost.name != nullptr) {
    ouput_files.fastq.lost.handle = fopen_output(ouput_files.fastq.lost.name);
  }
}


auto abort_if_fastq_out_of_fasta(struct file_types const & ouput_files) -> void {
  auto const output_is_fastq = (ouput_files.fastq.kept.handle != nullptr
                                or ouput_files.fastq.lost.handle != nullptr);
  auto const input_is_fasta = not db_is_fastq();
  if (input_is_fasta and output_is_fastq) {
    fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
  }
}


auto check_output_files(struct file_types const & ouput_files) -> void {
  if (ouput_files.fasta.kept.name != nullptr) {
    if (ouput_files.fasta.kept.handle == nullptr) {
      fatal("Unable to open FASTA output file for writing");
    }
  }
  if (ouput_files.fasta.lost.name != nullptr) {
    if (ouput_files.fasta.lost.handle == nullptr) {
      fatal("Unable to open FASTA output file for writing");
    }
  }
  if (ouput_files.fastq.kept.name != nullptr) {
    if (ouput_files.fastq.kept.handle == nullptr) {
      fatal("Unable to open FASTQ output file for writing");
    }
  }
  if (ouput_files.fastq.lost.name != nullptr) {
    if (ouput_files.fastq.lost.handle == nullptr) {
      fatal("Unable to open FASTQ output file for writing");
    }
  }
}


namespace {
  // anonymous namespace to avoid linker error (multiple definitions
  // of function with identical names and parameters)
  auto create_deck(bool const sizein_requested) -> std::vector<int> {
    auto const dbsequencecount = db_getsequencecount();
    std::vector<int> deck(dbsequencecount, 1);
    if (sizein_requested) {
      auto counter = std::size_t{0};
      for (auto & abundance : deck) {
        abundance = db_getabundance(counter);
        ++counter;
      }
    }
    return deck;
  }
}


auto write_original_stats(std::vector<int> const & deck,
                          uint64_t const mass_total,
                          struct Parameters const & parameters) -> void {
  if (not parameters.opt_quiet) {
    std::fprintf(stderr, "Got %" PRIu64 " reads from %d amplicons\n",
                 mass_total, static_cast<int>(deck.size()));
  }
  if (parameters.opt_log != nullptr) {
    std::fprintf(parameters.fp_log, "Got %" PRIu64 " reads from %d amplicons\n",
                 mass_total, static_cast<int>(deck.size()));
  }
}


auto number_of_reads_to_sample(struct Parameters const & parameters,
                               uint64_t const mass_total) -> uint64_t {
  assert(mass_total <= contiguous_mantissa);
  if (parameters.opt_sample_size != 0) {
    return static_cast<uint64_t>(parameters.opt_sample_size);
  }
  return static_cast<uint64_t>(std::floor(static_cast<double>(mass_total) * parameters.opt_sample_pct / 100.0));
}


auto write_subsampling_stats(std::vector<int> const &deck,
                             uint64_t const n_reads,
                             struct Parameters const & parameters) -> void {
  int const samples = std::count_if(deck.begin(),
                                    deck.end(), [](int abundance) -> bool { return abundance != 0; });
  if (not parameters.opt_quiet) {
    std::fprintf(stderr, "Subsampled %" PRIu64 " reads from %d amplicons\n", n_reads, samples);
  }
  if (parameters.opt_log != nullptr) {
    std::fprintf(parameters.fp_log, "Subsampled %" PRIu64 " reads from %d amplicons\n", n_reads, samples);
  }
}


auto random_subsampling(std::vector<int> & deck, uint64_t const mass_total,
                        uint64_t const n_reads, bool const sizein_requested) -> void {
  auto n_reads_left = n_reads;
  auto amplicon_number = 0;
  uint64_t n_read_being_checked = 0;
  uint64_t accumulated_mass = 0;
  auto amplicon_mass = sizein_requested ? db_getabundance(0) : 1;

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
          amplicon_mass = sizein_requested ? db_getabundance(amplicon_number) : 1;
          accumulated_mass = 0;
        }
      progress_update(n_read_being_checked);
    }
  progress_done();
}


auto substract_two_decks(std::vector<int> const & original_deck,
                         std::vector<int> const & subsampled_deck) -> std::vector<int> {
  std::vector<int> difference_deck(original_deck.size());
  std::transform(original_deck.cbegin(), original_deck.cend(),
                 subsampled_deck.cbegin(), difference_deck.begin(),
                 std::minus<int>());
  return difference_deck;
}


auto writing_fasta_output(std::vector<int> const & deck,
                          struct a_file const & fasta_file) -> void {
  if (fasta_file.name == nullptr) {
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
      fasta_print_general(fasta_file.handle,
                          nullptr,
                          db_getsequence(counter),
                          static_cast<int>(db_getsequencelen(counter)),
                          db_getheader(counter),
                          static_cast<int>(db_getheaderlen(counter)),
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
                          struct a_file const & fastq_file) -> void {
  if (fastq_file.name == nullptr) {
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
      fastq_print_general(fastq_file.handle,
                          db_getsequence(counter),
                          static_cast<int>(db_getsequencelen(counter)),
                          db_getheader(counter),
                          static_cast<int>(db_getheaderlen(counter)),
                          db_getquality(counter),
                          static_cast<int>(new_abundance),
                          amplicons_printed,
                          -1.0);
      progress_update(counter);
      ++counter;
  }
  progress_done();
}


auto close_output_files(struct file_types const & ouput_files) -> void {
  for (auto * fp_outputfile : {
           ouput_files.fasta.kept.handle, ouput_files.fastq.kept.handle,
           ouput_files.fasta.lost.handle, ouput_files.fastq.lost.handle}) {
    if (fp_outputfile != nullptr) {
      static_cast<void>(std::fclose(fp_outputfile));
    }
  }
}


auto subsample(struct Parameters const & parameters) -> void {
  struct file_types ouput_files = {};
  ouput_files.fasta.kept.name = parameters.opt_fastaout;
  ouput_files.fasta.lost.name = parameters.opt_fastaout_discarded;
  ouput_files.fastq.kept.name = parameters.opt_fastqout;
  ouput_files.fastq.lost.name = parameters.opt_fastqout_discarded;
  open_output_files(ouput_files);
  check_output_files(ouput_files);

  db_read(parameters.opt_fastx_subsample, 0);
  show_rusage();

  abort_if_fastq_out_of_fasta(ouput_files);

  // subsampling
  auto const original_abundances = create_deck(parameters.opt_sizein);
  auto const mass_total = std::accumulate(original_abundances.cbegin(), original_abundances.cend(), uint64_t{0});
  auto subsampled_abundances = original_abundances;
  std::fill(subsampled_abundances.begin(), subsampled_abundances.end(), 0);  // temporary fix: reset vector to zero

  write_original_stats(original_abundances, mass_total, parameters);  // refactoring: move up?

  auto const n_reads = number_of_reads_to_sample(parameters, mass_total);

  if (n_reads > mass_total)
    {
      fatal("Cannot subsample more reads than in the original sample");
    }

  random_subsampling(subsampled_abundances, mass_total, n_reads, parameters.opt_sizein);  // refactoring: pass & original, copy, subsample, return new (const) vector

  // write output files
  writing_fasta_output(subsampled_abundances, ouput_files.fasta.kept);
  writing_fastq_output(subsampled_abundances, ouput_files.fastq.kept);
  auto const discarded_output_is_requested = (ouput_files.fasta.lost.handle != nullptr) or (ouput_files.fastq.lost.handle != nullptr);
  if (discarded_output_is_requested) {
    auto const discarded_abundances = substract_two_decks(original_abundances,
                                                          subsampled_abundances);
    writing_fasta_output(discarded_abundances, ouput_files.fasta.lost);
    writing_fastq_output(discarded_abundances, ouput_files.fastq.lost);
  }

  write_subsampling_stats(subsampled_abundances, n_reads, parameters);

  // clean up
  db_free();
  close_output_files(ouput_files);
}
