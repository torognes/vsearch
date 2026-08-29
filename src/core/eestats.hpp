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

#include <array>
#include <cassert>
#include <cstddef>  // std::size_t

/* Declared at global scope on purpose: QualityScoreTable's constructor below
   names `struct Parameters`, and an elaborated-type-specifier that finds
   nothing *declares* the type in its own scope -- which would silently make
   it an incomplete `vsearch::Parameters`, distinct from the real one. This
   used to be provided as a side effect of fastq_get_qual_eestats' signature
   sitting here; it is stated outright now that the function is gone. */
struct Parameters;

/* Quality helper shared by the fastq_eestats and fastq_eestats2 commands:
   q2p converts a Phred quality value to its error probability.

   fastq_get_qual_eestats() used to sit here too, decoding one symbol and
   range-checking it. Both callers reached it only to die -- they discarded
   its int and took the decoded value from QualityScoreTable::score() -- so
   they now call vsearch::check_quality_score() themselves, which is the same
   check and lets them name the record from the handle they already hold. */
auto q2p(int quality_value) -> double;


namespace vsearch {

/* Quality symbol -> decoded quality score, precomputed once per run.

   Decoding and range-checking each base through a function call into
   another translation unit cost a third of fastq_eestats2's runtime, by
   profiling; but the inputs are only the 128 low ASCII ordinals (see
   quality_table.hpp for why a symbol reaching a conversion always lies in
   [33, 126]), so both the decoded score and the range verdict fit in small
   per-symbol arrays filled up front.

   accepts() is vsearch::classify_quality() evaluated once per symbol rather
   than once per base; a caller that gets false back hands the symbol to
   vsearch::check_quality_score(), which is the same test again and produces
   the message. score() folds in the std::max(qual, 0) both eestats commands
   applied to the decoded value. */
class QualityScoreTable {
public:
  static constexpr auto n_symbols = std::size_t{128};

  explicit QualityScoreTable(struct Parameters const & parameters) noexcept;

  auto accepts(char const quality_symbol) const noexcept -> bool
  {
    return accepted_[ordinal_of(quality_symbol)];
  }

  /* the decoded quality score, clamped to zero at the low end */
  auto score(char const quality_symbol) const noexcept -> int
  {
    return scores_[ordinal_of(quality_symbol)];
  }

private:
  static auto ordinal_of(char const quality_symbol) noexcept -> std::size_t
  {
    auto const ordinal = static_cast<unsigned char>(quality_symbol);
    assert(ordinal < n_symbols);
    return ordinal;
  }

  std::array<int, n_symbols> scores_ {{}};
  std::array<bool, n_symbols> accepted_ {{}};
};

}  // namespace vsearch
