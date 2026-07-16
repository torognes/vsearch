/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

// Command-line interface: parse and validate the user-supplied options,
// populate the global opt_* variables and the Parameters struct, and report
// usage errors. Extracted from vsearch.cc to keep the argument-parsing
// machinery separate from the command dispatch and the main program.

#pragma once


/* The single command a run performs, resolved by the CLI parser (from the
   requested command option) and returned to main() for the command dispatcher.
   One enumerator per dispatch handler: the --h/--help and --v/--version option
   aliases each collapse to a single command. Command::none means no (or no
   valid) command was requested. The underlying type is fixed to int for a
   stable, non-narrowing representation. This is CLI dispatch state, so it is
   deliberately kept out of the public Parameters/library surface. */
enum struct Command : int
  {
    none,
    help,
    version,
    allpairs_global,
    usearch_global,
    search_exact,
    sintax,
    orient,
    cluster_fast,
    cluster_smallmem,
    cluster_size,
    cluster_unoise,
    uchime_denovo,
    uchime2_denovo,
    uchime3_denovo,
    uchime_ref,
    chimeras_denovo,
    derep_fulllength,
    derep_prefix,
    derep_id,
    derep_smallmem,
    fastq_chars,
    fastq_stats,
    fastq_filter,
    fastx_filter,
    fastq_convert,
    fastq_eestats,
    fastq_eestats2,
    fastq_join,
    fastq_mergepairs,
    fastx_uniques,
    fastx_mask,
    fastx_revcomp,
    fastx_syncpairs,
    fastx_getseq,
    fastx_getseqs,
    fastx_getsubseq,
    fastx_subsample,
    fasta2fastq,
    cut,
    shuffle,
    sortbylength,
    sortbysize,
    rereplicate,
    maskfasta,
    sff_convert,
    makeudb_usearch,
    udb2fasta,
    udbinfo,
    udbstats
  };

// Parse the command line, set the matching fields in parameters, validate the
// requested command and its options, and return the resolved command.
auto args_init(int argc, char ** argv, struct Parameters & parameters) -> Command;
