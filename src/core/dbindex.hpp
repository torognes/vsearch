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
#include "core/unique.hpp"  // Uniquer
#include "utils/fatal_allocator.hpp"  // FatalAllocator
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t
#include <vector>  // std::vector


struct Parameters;
struct Database;


/* One database sequence that is in the k-mer index but holds fewer distinct
   k-mers than the candidate threshold asks for: its index element number (the
   numbering of Dbindex::map, not the sequence number) and how many distinct
   k-mers it contributed. Named fields rather than two adjacent unsigned ints a
   call site would be free to swap. */
struct LowKmerTarget
{
  unsigned int index;
  unsigned int kmers;
};


/* The database k-mer index. Owns its buffers (RAII: released by clear() and the
   destructor) and is non-copyable. The read API (the getX members) is const so
   the search worker threads can query one shared index concurrently. Built
   either by prepare() + add_all_sequences() from the in-memory FASTA database,
   or by udb_read() straight from a UDB file. */
struct Dbindex
{
  std::vector<unsigned int, FatalAllocator<unsigned int>> kmercount; /* number of matching seqnos for each kmer */
  std::vector<uint64_t, FatalAllocator<uint64_t>> kmerhash;  /* index into the list below for each kmer */
  std::vector<unsigned int, FatalAllocator<unsigned int>> kmerindex; /* the list of matching seqnos for kmers */
  /* Which k-mers use the bitmap form, and their bit-sets. kmerbitmap_slot has
     one entry per k-mer slot holding a 1-based index into bitmap_pool, or zero
     when the k-mer uses the list form (kmerindex) instead; bitmap_pool holds
     only the k-mers that earn a bitmap. Storing a 4-byte index per slot rather
     than a whole Bitmap object (24 bytes, or 56 under _GLIBCXX_DEBUG) keeps
     the dense array cheap: at word length 12 it is 64 MB instead of 384 MB
     (896 MB debug), allocated and zeroed on every index build, however few
     k-mers the database really contains. Indices rather than pointers, so
     growing the pool cannot leave a dangling reference behind. */
  std::vector<unsigned int, FatalAllocator<unsigned int>> kmerbitmap_slot;
  std::vector<Bitmap> bitmap_pool;
  unsigned int bitmap_width = 0;  /* bits in each pooled bitmap; see set_bitmap_width */
  std::vector<unsigned int, FatalAllocator<unsigned int>> map;  /* mapping from index element number to seqno */
  /* The index elements holding fewer distinct k-mers than minwordmatches, with
     that count. A target cannot share more k-mers than it holds, so asking one
     of these for minwordmatches shared k-mers asks for more than it can supply
     and hides it from every query, an exact match included
     (torognes/vsearch#328); search_topscores() gives each the threshold it can
     meet instead. A list rather than one threshold per index element because
     search reads it once per query on top of the per-sequence counter loop,
     and that loop is hot enough for one extra load per indexed sequence to
     show (see search_topscores in core/searchcore.cpp) -- while the list is
     empty for a database of ordinary-length, lightly masked sequences, and
     then costs nothing at all. Sequences that yield no k-mer are left out:
     they are absent from the index, can share nothing, and stay unreachable
     (prepare() warns about them). */
  std::vector<LowKmerTarget, FatalAllocator<LowKmerTarget>> low_kmer_targets;
  Uniquer uhandle {};  /* unique-kmer finder, used while building */
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

  /* The candidate threshold low_kmer_targets was built against, i.e. the
     opt_minwordmatches of the session that built the index: index state, like
     wordlength above, because it decides which targets are on that list. Set
     by prepare() for a FASTA database and by udb_read for a UDB one, and
     asserted against the searching session's own parameters in
     search_topscores. */
  unsigned int minwordmatches = 0;

  Dbindex() = default;
  ~Dbindex();
  Dbindex(Dbindex const &) = delete;
  auto operator=(Dbindex const &) -> Dbindex & = delete;
  Dbindex(Dbindex &&) = delete;
  auto operator=(Dbindex &&) -> Dbindex & = delete;

  auto prepare(Masking seqmask, struct Database const & db, struct Parameters const & parameters) -> void;
  auto add_sequence(unsigned int seqno, Masking seqmask, struct Database const & db) -> void;
  auto add_all_sequences(Masking seqmask, struct Database const & db, struct Parameters const & parameters) -> void;
  auto clear() -> void;

  /* Size the k-mer -> bitmap lookup to slots entries, all on the list form, and
     drop any bitmaps from a previous build. */
  auto bitmap_slots_reset(unsigned int slots) -> void;

  /* Width every bitmap in this index will have, as a number of database
     sequences: one bit each, plus the padding the SIMD counter routines need to
     read whole xmm registers past the last sequence. Call before the first
     bitmap_create. Recording it once, rather than passing a size to each
     bitmap_create, is what makes the widths uniform -- search scans every
     bitmap to the same indexed_count bound, so a short one would be read out of
     bounds. */
  auto set_bitmap_width(unsigned int sequences) -> void;

  /* Give kmer a bitmap of the recorded width and return it, for the caller to
     fill. The reference is invalidated by the next bitmap_create, which may
     grow the pool, so fill it before creating the next one. */
  auto bitmap_create(unsigned int kmer) -> Bitmap &;

  auto has_bitmap(unsigned int kmer) const -> bool;

  /* The bitmap of a k-mer that has one (assert): for the index builders and the
     UDB writer. Search reads bitmaps through getbitmap() instead. */
  auto bitmap_of(unsigned int kmer) -> Bitmap &;
  auto bitmap_of(unsigned int kmer) const -> Bitmap const &;

  auto getbitmap(unsigned int kmer) const -> unsigned char const *;
  auto getmatchcount(unsigned int kmer) const -> unsigned int;
  auto getmatchlist(unsigned int kmer) const -> unsigned int const *;
  auto getmapping(unsigned int index) const -> unsigned int;
  auto getcount() const -> unsigned int;
};


auto fprint_kmer(std::FILE * output_handle, unsigned int kmer_length, uint64_t kmer) -> void;


/* Smallest number of matching sequences that earns a k-mer a bitmap instead of
   the list form (kmerindex). Shared by the two index builders -- prepare() for
   a FASTA database and udb_read() for a UDB one, which fills the same bitmaps
   from the stored index -- so the rule cannot drift between them.
   noexcept: integer arithmetic on a by-value argument. */
auto bitmap_min_matches(unsigned int seqcount) noexcept -> unsigned int;
