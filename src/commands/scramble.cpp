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


#include "commands/scramble.hpp"
#include "utils/span.hpp"
#include "utils/view.hpp"
#include "vsearch.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/random.hpp"
#include <algorithm>  // std::copy
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::size_t
#include <iterator>  // std::next
#include <vector>


constexpr auto initial_memory_allocation = 512;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* One Scrambler per run owns the scramble operation. With
     --scramble_kmer fixed at 1 (the only accepted value for now) it
     forwards to the portable Fisher-Yates in utils/random.hpp; this
     class is the seam where the k >= 2 Eulerian-path sampler (a
     private de Bruijn graph component with buffers reused across
     records) plugs in later. */
  class Scrambler {
  public:
    template <typename URBG>
    auto scramble(Span<char> const sequence, URBG & generator) -> void {
      random_shuffle(sequence, generator);
    }
  };

}  // end of anonymous namespace


auto scramble(struct Parameters const & parameters) -> void
{
  uint64_t buffer_alloc = initial_memory_allocation;
  std::vector<char> seq_buffer(buffer_alloc);
  std::vector<char> qual_buffer(buffer_alloc);

  if ((parameters.opt_fastaout == nullptr) && (parameters.opt_fastqout == nullptr)) {
    fatal("No output files specified");
  }

  /* K >= 2 (preserving all j-mer counts for j <= K) awaits the
     Eulerian-path sampler; only mononucleotide scrambling is available */
  if (parameters.opt_scramble_kmer > 1) {
    fatal("--scramble_kmer values greater than 1 are not supported (yet)");
  }

  auto input_handle = fastx_open(parameters.opt_scramble, parameters);

  if ((parameters.opt_fastqout != nullptr) && ! input_handle->is_fastq_input())
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  auto const filesize = input_handle->get_size();

  auto fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
  auto fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
  std::FILE * const fp_fastaout = fastaout_handle.get();
  std::FILE * const fp_fastqout = fastqout_handle.get();

  /* the RandomSeed carries the full 64-bit --randseed (or an OS value
     when 0); every draw goes through random_bounded(), so a given seed
     yields the same output on any platform (see utils/random.hpp).
     RandomSeed's constructor may throw (std::random_device on the
     --randseed 0 path), so scramble() cannot be noexcept. */
  RandomSeed const seed(parameters);
  Scrambler scrambler;

  {
    int64_t count = 0;  // the ordinal fed to --relabel; int would wrap at 2^31 records
    Progress progress("Scrambling", filesize, parameters);
    while (input_handle->next(false, chrmap_no_change()))
      {
        ++count;

        /* header */

        auto const hlen = input_handle->get_header_length();
        auto const * header = input_handle->get_header();
        auto const abundance = input_handle->get_abundance();


        /* sequence */

        auto const length = input_handle->get_sequence_length();

        if (length + 1 > buffer_alloc)
          {
            buffer_alloc = length + 1;
            seq_buffer.resize(buffer_alloc);
            qual_buffer.resize(buffer_alloc);
          }

        auto const * seq = input_handle->get_sequence();
        std::copy(seq, std::next(seq, static_cast<std::ptrdiff_t>(length)),
                  seq_buffer.begin());
        seq_buffer[length] = '\0';

        /* per-record substream (sintax-style): each record's scramble
           depends only on the base seed, its ordinal, and its length,
           never on how many draws earlier records consumed -- so fasta
           and fastq runs over the same records scramble identically */
        SplitMix64 generator(seed.substream(static_cast<uint64_t>(count)));
        scrambler.scramble(make_span(seq_buffer).first(static_cast<std::size_t>(length)),
                           generator);


        /* quality values */

        auto const * qual = input_handle->get_quality();

        if (input_handle->is_fastq_input())
          {
            /* the quality string is copied through untouched: quality is
               never scrambled, so the positional quality profile of each
               record is preserved exactly while the base<->quality pairing
               is deliberately broken */
            std::copy(qual, std::next(qual, static_cast<std::ptrdiff_t>(length)),
                      qual_buffer.begin());
            qual_buffer[length] = '\0';
          }

        if (parameters.opt_fastaout != nullptr)
          {
            fasta_print_general(fp_fastaout,
                                nullptr,
                                make_view(seq_buffer).first(static_cast<std::size_t>(length)),
                                View<char>{header, static_cast<std::size_t>(hlen)},
                                OutputAnnotations{static_cast<uint64_t>(abundance), count},
                                parameters);
          }

        if (parameters.opt_fastqout != nullptr)
          {
            fastq_print_general(fp_fastqout,
                                make_view(seq_buffer).first(static_cast<std::size_t>(length)),
                                View<char>{header, static_cast<std::size_t>(hlen)},
                                make_view(qual_buffer).first(static_cast<std::size_t>(length)),
                                OutputAnnotations{static_cast<uint64_t>(abundance), count},
                                parameters);
          }

        progress.update(input_handle->get_position());
      }
  }

  if (parameters.opt_fastaout != nullptr)
    {
      fastaout_handle.reset();
    }

  if (parameters.opt_fastqout != nullptr)
    {
      fastqout_handle.reset();
    }

  input_handle->report_stripped_warning(parameters);
}
