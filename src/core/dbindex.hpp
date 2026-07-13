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

#include "core/bitmap.hpp"
#include "core/mask.hpp"  // Masking
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t


struct uhandle_s;
struct Parameters;
struct Database;


/* The database k-mer index. Owns its buffers (RAII: released by clear() and the
   destructor) and is non-copyable. The read API (the getX members) is const so
   the search worker threads can query one shared index concurrently. Built
   either by prepare() + add_all_sequences() from the in-memory FASTA database,
   or by udb_read() straight from a UDB file. */
struct Dbindex
{
  unsigned int * kmercount = nullptr; /* number of matching seqnos for each kmer */
  uint64_t * kmerhash = nullptr;  /* index into the list below for each kmer */
  unsigned int * kmerindex = nullptr; /* the list of matching seqnos for kmers */
  struct bitmap_s * * kmerbitmap = nullptr;
  unsigned int * map = nullptr;  /* mapping from index element number to seqno */
  uhandle_s * uhandle = nullptr;  /* unique-kmer finder, used while building */
  unsigned int count = 0;  /* number of sequences added to the index */
  unsigned int hashsize = 0;  /* number of kmer slots, i.e. 4^wordlength */
  uint64_t indexsize = 0;  /* total number of entries in kmerindex */

  /* effective word length of the built k-mer index (derived index state, not
     config): set by prepare (from parameters.opt_wordlength) for a FASTA
     database, or by udb_read for a UDB database whose stored word length differs.
     Read by everything that extracts query k-mers, so they match the index width
     without consulting the (immutable) opt_wordlength config. */
  unsigned int wordlength = 0;

  /* DBAccel percentage stored in a UDB header (buffer[6]): index metadata
     carried for reporting only, never consulted during search. Set by udb_read
     from the file, reported by udbstats; stays 0 for a FASTA-built index. */
  unsigned int dbaccel = 0;

  Dbindex() = default;
  ~Dbindex();
  Dbindex(Dbindex const &) = delete;
  auto operator=(Dbindex const &) -> Dbindex & = delete;
  Dbindex(Dbindex &&) = delete;
  auto operator=(Dbindex &&) -> Dbindex & = delete;

  auto prepare(int use_bitmap, Masking seqmask, struct Database const & db, struct Parameters const & parameters) -> void;
  auto add_sequence(unsigned int seqno, Masking seqmask, struct Database const & db) -> void;
  auto add_all_sequences(Masking seqmask, struct Database const & db, struct Parameters const & parameters) -> void;
  auto clear() -> void;

  auto getbitmap(unsigned int kmer) const -> unsigned char *;
  auto getmatchcount(unsigned int kmer) const -> unsigned int;
  auto getmatchlist(unsigned int kmer) const -> unsigned int *;
  auto getmapping(unsigned int index) const -> unsigned int;
  auto getcount() const -> unsigned int;
};


auto fprint_kmer(std::FILE * output_handle, unsigned int kmer_length, uint64_t kmer) -> void;
