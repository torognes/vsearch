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

/* FASTQ quality-score ASCII offsets: the byte value subtracted from a
   quality character to recover its Phred score. Two encodings are in use;
   sanger_ascii_offset is the default (and the default for --fastq_ascii). */

constexpr int sanger_ascii_offset = 33;  // Phred+33 (Sanger / Illumina 1.8+)
constexpr int solexa_ascii_offset = 64;  // Phred+64 (Solexa / Illumina 1.3+)

/* The printable ASCII range a quality symbol may occupy. Every quality bound
   is constrained by these two: an offset plus a score must land inside
   [33, 126], which caps the representable score at 126 - offset (93 with the
   Sanger offset, 62 with the Solexa one) and floors it at 33 - offset. The
   sum rules in cli.cc and the default of --fastq_qmax are both stated in
   terms of these. */
constexpr int lowest_printable_ascii = 33;   // '!'
constexpr int highest_printable_ascii = 126; // '~'

/* The quality ceiling vsearch used before 3.0, when 41 was the highest score
   an Illumina 1.8+ file was expected to carry. It survives as the
   --fastq_qmaxout default of --fasta2fastq alone, which fabricates quality
   rather than clamping it (see resolve_quality_bound_defaults() in cli.cc). */
constexpr int legacy_max_quality = 41;
