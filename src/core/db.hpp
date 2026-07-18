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

#include "utils/fatal_allocator.hpp"  // FatalAllocator
#include "utils/span.hpp"  // Span<char>
#include "utils/view.hpp"  // View<char>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::size_t
#include <vector>


struct seqinfo_s
{
  std::size_t header_p {};
  std::size_t seq_p {};
  std::size_t qual_p {};
  // headerlen and seqlen are intentionally kept 32-bit: --maxseqlength is
  // capped at INT_MAX - buffer_headroom (see cli.cc), well under the 32-bit field
  // maximum, so no accepted record can overflow them, and the two fields share
  // a single 8-byte slot here. Widening either to 64-bit grows seqinfo_s from
  // 40 to 48 bytes — +8 bytes per sequence in the seqindex array (alignment
  // padding makes widening just one cost the same as widening both). If the
  // --maxseqlength cap is ever raised above UINT32_MAX, widen both to uint64_t
  // (and drop the narrowing casts in Database::add()); until then the cap is the guard
  // for the N1(c) storage path.
  unsigned int headerlen {};
  unsigned int seqlen {};
  uint64_t size {};  // abundance; 64-bit so a per-sequence ;size= above UINT_MAX
                     // does not wrap (fits the existing struct tail padding)
};

using seqinfo_t = struct seqinfo_s;


/* A read-only bundle of one database record's fields, returned by
   Database::record(). Groups the header, sequence and (FASTQ-only) quality as
   non-owning View<char> windows into the database's shared buffer, plus the
   parsed abundance. It suits the record-emit paths (e.g. FASTA/FASTQ printing)
   that consume every field of a record together; the per-field view accessors
   remain for callers that need only one. The quality view is empty for a FASTA
   database. */
struct DbRecord
{
  View<char> header;
  View<char> sequence;
  View<char> quality;
  // no default member initializer on purpose: it would make DbRecord a
  // non-aggregate under C++11, breaking the aggregate initialization in
  // Database::record(). record() always sets every field, so leaving abundance
  // uninitialized here is safe (and View has no default constructor, so a
  // fieldless DbRecord cannot be built by accident anyway).
  uint64_t   abundance;
};


/* The in-memory sequence database. Owns its two heap buffers (RAII: released by
   clear() and the destructor) and is non-copyable/non-movable. The read API (the
   getX members) is const so search worker threads can query one shared database
   concurrently. It is populated in one of three ways: init() + add() when a
   library caller assembles a database programmatically, read() from a
   FASTA/FASTQ file, or udb_read() straight from a UDB file. Each command (and
   each library caller) owns its own instance and threads a reference through
   the code. */
struct Database
{
private:
  // packed headers, sequences and qualities; per-sequence offsets/lengths/
  // abundance. FatalAllocator keeps the out-of-memory behaviour of the former
  // raw xmalloc buffers (fatal(), not std::terminate).
  std::vector<char, FatalAllocator<char>>           data_;
  std::vector<seqinfo_t, FatalAllocator<seqinfo_t>> seqindex_;

  bool     fastq_format = false;  // read through the is_fastq() accessor
  uint64_t sequences = 0;
  uint64_t nucleotides = 0;
  uint64_t longest = 0;
  uint64_t shortest = 0;
  uint64_t longestheader = 0;

  /* udb_read is a second database loader that fills data_/seqindex_ in place
     (it bypasses add()); grant it access to the otherwise-private buffers. */
  friend auto udb_read(const char * filename,
                       bool create_bitmaps,
                       bool parse_abundances,
                       struct Dbindex & dbindex,
                       struct Database & db,
                       struct Parameters const & parameters) -> void;

public:
  Database() = default;
  ~Database();
  Database(Database const &) = delete;
  auto operator=(Database const &) -> Database & = delete;
  Database(Database &&) = delete;
  auto operator=(Database &&) -> Database & = delete;

  auto init() -> void;

  auto add(bool is_fastq_record,
           char const * header,
           char const * sequence,
           char const * quality,
           std::size_t headerlength,
           std::size_t sequencelength,
           int64_t abundance) -> void;

  auto read(char const * filename, int upcase, struct Parameters const & parameters) -> void;
  auto clear() -> void;

  /* UDB bulk-load seam, used only by udb_read (the second database loader, which
     bypasses add() and fills the reserved buffers in place, mirroring how it
     fills Dbindex's buffers). udb_reserve allocates the two buffers up front;
     udb_finalize runs the terminator-insertion memmove pass and records the
     summary statistics. */
  auto udb_reserve(uint64_t count, uint64_t datap_bytes) -> void;
  auto udb_finalize(uint64_t count,
                    uint64_t nucleotide_count,
                    uint64_t longest_sequence,
                    uint64_t shortest_sequence,
                    uint64_t longest_header,
                    struct Parameters const & parameters) -> void;

  auto getquality(uint64_t seqno) const -> char const *;

  /* Note: the sorting members below must be called after read(),
     but before Dbindex::prepare */
  auto sortbylength(struct Parameters const & parameters) -> void;
  auto sortbylength_shortest_first(struct Parameters const & parameters) -> void;
  auto sortbyabundance(struct Parameters const & parameters) -> void;

  auto getheader(uint64_t seqno) const -> char const *
  {
    return data_.data() + seqindex_[seqno].header_p;
  }

  auto getsequence(uint64_t seqno) const -> char const *
  {
    return data_.data() + seqindex_[seqno].seq_p;
  }

  /* Non-const companion to getsequence(): hands out a writable pointer to the
     stored sequence. It backs the in-place masking passes (dust/hardmask) and is
     borrowed by allpairs_global into its writable qsequence buffer. Only a
     non-const Database exposes it, so read-only holders (e.g. the search worker
     threads) cannot mutate the database. */
  auto mutatesequence(uint64_t seqno) -> char *
  {
    return data_.data() + seqindex_[seqno].seq_p;
  }

  auto getabundance(uint64_t seqno) const -> uint64_t
  {
    return seqindex_[seqno].size;
  }

  auto getsequencelen(uint64_t seqno) const -> uint64_t
  {
    return seqindex_[seqno].seqlen;
  }

  auto getheaderlen(uint64_t seqno) const -> uint64_t
  {
    return seqindex_[seqno].headerlen;
  }

  auto getsequencecount() const -> uint64_t { return sequences; }
  auto getnucleotidecount() const -> uint64_t { return nucleotides; }
  auto getlongestheader() const -> uint64_t { return longestheader; }
  auto getlongestsequence() const -> uint64_t { return longest; }
  auto getshortestsequence() const -> uint64_t { return shortest; }
  auto is_fastq() const -> bool { return fastq_format; }

  /* View/Span companions to the raw pointer + length accessors above: each
     bundles a record's interior pointer with its length into a single
     non-owning window into the shared buffer, so callers no longer thread a
     (pointer, length) pair. They are header-inline like the raw accessors, so
     returning the small (16-byte) value costs no more than the bare pointer. */
  auto sequence_view(uint64_t seqno) const -> View<char>
  {
    return View<char>{getsequence(seqno), getsequencelen(seqno)};
  }

  auto header_view(uint64_t seqno) const -> View<char>
  {
    return View<char>{getheader(seqno), getheaderlen(seqno)};
  }

  /* A FASTA database has no quality, so the view is empty; for FASTQ the
     quality string has the same length as the sequence. */
  auto quality_view(uint64_t seqno) const -> View<char>
  {
    return View<char>{getquality(seqno),
                      fastq_format ? getsequencelen(seqno) : uint64_t{0}};
  }

  /* Writable companion to mutatesequence(): only a non-const Database exposes
     it, so read-only holders (e.g. search worker threads) cannot obtain a
     mutable window into the stored sequence. */
  auto mutable_sequence(uint64_t seqno) -> Span<char>
  {
    return Span<char>{mutatesequence(seqno), getsequencelen(seqno)};
  }

  auto record(uint64_t seqno) const -> DbRecord
  {
    return DbRecord{header_view(seqno), sequence_view(seqno),
                    quality_view(seqno), getabundance(seqno)};
  }
};


