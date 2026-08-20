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

#pragma once

/* Phred quality symbol -> error probability, precomputed once per run.

   Every FASTQ-aware command needs 10^-(q/10) for each quality symbol it
   reads. Evaluated per base that is a std::pow() call into libm, about 5 ns
   and not inlinable; but the domain is only the 128 low ASCII ordinals, so
   the entire function fits in a 1 KiB array that stays resident in L1 and
   costs a single load. Two commands had already arrived at this table
   independently -- core/mergepairs.cpp (QualityTables::q2p) and
   commands/fastq_stats.cpp (precompute_probability_values) -- with two
   different shapes and two different clamping rules; this is the shared
   one.

   Each entry is produced by the very expression it replaces, so a lookup is
   bit-identical to the call, not an approximation. */

#include <algorithm>  // std::min
#include <array>
#include <cassert>
#include <cmath>  // std::pow
#include <cstddef>  // std::size_t
#include <limits>  // std::numeric_limits


namespace vsearch {

/* Every ASCII ordinal a quality symbol can take. The FASTQ readers reject
   DEL and every byte above it (quality_policy in core/fastq.cpp), and the
   option checks in cli.cc require 33 <= fastq_ascii + fastq_qmin and
   fastq_ascii + fastq_qmax <= 126, so a symbol reaching a conversion always
   lies in [33, 126]. The table spans the whole low ASCII range anyway, so
   indexing needs neither a bias nor a bounds test. */
constexpr auto quality_symbol_count = std::size_t{128};

/* The largest error probability an entry may hold. The three values are the
   three rules already in the code, not a new policy:

   - none: 10^-(q/10) as written. A negative quality value yields more than
     1.0, which --fastq_ascii 64 --fastq_qmin -10 can reach (cli.cc only
     requires fastq_ascii + fastq_qmin >= 33). core/filter.cpp behaves this
     way today.
   - certain_error: 1.0, i.e. a negative quality value is read as quality
     zero. What core/eestats.cpp's callers get from their std::max(qual, 0).
   - random_guess: 0.75, the chance of being wrong when picking one of four
     bases blind, which is all a quality value below 2 says. Capping the
     probability here is exactly the "if (quality_value < 2) return 0.75"
     branch in core/derep.cpp and core/mergepairs.cpp: 10^0 is 1.0 and
     10^-0.1 is 0.794, both above the cap, while 10^-0.2 is 0.631, below
     it. USEARCH does the same. */
enum struct ProbabilityCap { none, certain_error, random_guess };

constexpr double certain_error_probability = 1.0;
constexpr double blind_guess_over_four_bases = 0.75;

/* Written as a ternary chain rather than a switch because a C++11 constexpr
   function may hold only one return statement; the same shape as
   sequence_policy and quality_policy in core/fastq.cpp. */
constexpr auto largest_probability(ProbabilityCap const cap) noexcept -> double
{
  return cap == ProbabilityCap::certain_error ? certain_error_probability
    : cap == ProbabilityCap::random_guess     ? blind_guess_over_four_bases
    : std::numeric_limits<double>::infinity();
}


class QualityTable {
public:
  QualityTable() = default;

  /* noexcept: std::pow neither allocates nor throws (C++11 [res.on.exception.handling]
     forbids the math functions from throwing), and the table is a member array,
     so there is nothing here that can fail. */
  QualityTable(int const ascii_offset, ProbabilityCap const cap) noexcept
  {
    static constexpr auto base = 10.0;
    auto const largest = largest_probability(cap);
    for (std::size_t ordinal = 0; ordinal < quality_symbol_count; ++ordinal)
      {
        auto const quality_value = static_cast<int>(ordinal) - ascii_offset;
        probabilities_[ordinal] = std::min(std::pow(base, -quality_value / base),
                                           largest);
      }
  }

  /* The error probability of one quality symbol. Takes the symbol, not the
     decoded quality value: the offset is already baked into the table. */
  auto operator[](char const quality_symbol) const noexcept -> double
  {
    auto const ordinal = static_cast<unsigned char>(quality_symbol);
    assert(ordinal < quality_symbol_count);
    return probabilities_[ordinal];
  }

private:
  std::array<double, quality_symbol_count> probabilities_ {{}};
};

}  // namespace vsearch
