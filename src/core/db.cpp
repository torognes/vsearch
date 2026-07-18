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

#include "vsearch.hpp"
#include "core/db.hpp"
#include "core/fastx.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/progress.hpp"
#include "utils/span.hpp"  // Span<char> (for mutable_sequence)
#include <algorithm>  // std::min, std::max, std::sort
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <cstring>  // std::memcpy, std::strcmp
#include <limits>
#include <memory>  // std::unique_ptr
#include <string>  // std::string


constexpr uint64_t memchunk = 16777216;  // 2^24


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* Grow a vector's capacity in fixed memchunk-byte steps, reproducing the raw
     buffers' former xrealloc growth policy so that a large database's peak
     memory does not balloon under std::vector's geometric growth. */
  template <class Vector>
  auto reserve_in_chunks(Vector & vec, std::size_t const needed_items) -> void
  {
    if (vec.capacity() >= needed_items) { return; }
    auto const item_size = sizeof(typename Vector::value_type);
    std::size_t chunked_bytes = vec.capacity() * item_size;
    std::size_t const needed_bytes = needed_items * item_size;
    while (chunked_bytes < needed_bytes) { chunked_bytes += memchunk; }
    vec.reserve(chunked_bytes / item_size);
  }

}  // end of anonymous namespace


/* Reset database state for library use.
   Must be called before the first add() when not using read().
   Frees any previously allocated data, then resets all counters.
   Sets shortest to UINT64_MAX so min-tracking works correctly. */
auto Database::init() -> void
{
  clear();
  fastq_format = false;
  sequences = 0;
  nucleotides = 0;
  longest = 0;
  shortest = std::numeric_limits<uint64_t>::max();
  longestheader = 0;
}


auto Database::udb_reserve(uint64_t const count, uint64_t const datap_bytes) -> void
{
  /* udb_read fills these buffers in place (it bypasses add()); size them up
     front to the exact extents derived from the UDB header. resize() (not
     reserve()) so data()/[] are valid across the whole buffer during the fill;
     the value-initialised bytes are overwritten by the loader. */
  seqindex_.resize(static_cast<std::size_t>(count));
  data_.resize(static_cast<std::size_t>(datap_bytes));
}


auto Database::udb_finalize(uint64_t const count,
                            uint64_t const nucleotide_count,
                            uint64_t const longest_sequence,
                            uint64_t const shortest_sequence,
                            uint64_t const longest_header,
                            struct Parameters const & parameters) -> void
{
  /* move sequences and insert zero at end of each sequence */

  {
    Progress progress("Reorganizing data in memory", count, parameters);
    char * const buffer = data_.data();
    for (auto i = count - 1; i > 0; i--)
      {
        auto const old_p = seqindex_[i].seq_p;
        auto const new_p = seqindex_[i].seq_p + i;
        auto const len   = seqindex_[i].seqlen;
        std::memmove(buffer + new_p, buffer + old_p, len);
        *(buffer + new_p + len) = 0;
        seqindex_[i].seq_p = new_p;
        progress.update(count - i);
      }
    *(buffer + seqindex_[0].seq_p + seqindex_[0].seqlen) = 0;
  }

  /* record the summary statistics (a UDB database is never FASTQ) */

  fastq_format = false;
  sequences = count;
  nucleotides = nucleotide_count;
  longest = longest_sequence;
  shortest = shortest_sequence;
  longestheader = longest_header;
}


auto Database::getquality(uint64_t seqno) const -> char const *
{
  if (fastq_format)
    {
      return data_.data() + seqindex_[seqno].qual_p;
    }
  return nullptr;
}


auto Database::mutable_sequence(uint64_t seqno) -> Span<char>
{
  return Span<char>{mutatesequence(seqno), getsequencelen(seqno)};
}


auto Database::add(bool const is_fastq_record,
                   char const * header,
                   char const * sequence,
                   char const * quality,
                   size_t const headerlength,
                   size_t const sequencelength,
                   int64_t const abundance) -> void
{
  /* Add a sequence to the database. Assumes that the database has been initialized. */

  /* grow space for data, if necessary (memchunk-stepped, as the raw buffer was) */
  size_t needed = data_.size() + headerlength + 1 + sequencelength + 1;
  if (is_fastq_record)
    {
      needed += sequencelength + 1;
    }
  reserve_in_chunks(data_, needed);

  /* store the header */
  size_t const header_p = data_.size();
  data_.insert(data_.end(), header, header + headerlength + 1);

  /* store sequence */
  size_t const sequence_p = data_.size();
  data_.insert(data_.end(), sequence, sequence + sequencelength + 1);

  size_t const quality_p = data_.size();
  if (is_fastq_record)
    {
      /* store quality */
      data_.insert(data_.end(), quality, quality + sequencelength + 1);

      /* A FASTQ record makes this a FASTQ database. read() sets the is_fastq
         flag itself before adding records, but callers that build a database
         directly with add() rely on this assignment so that the is_fastq flag
         and getquality() can reach the stored quality. */
      fastq_format = true;
    }

  /* grow space for the index, if necessary, then append the entry */
  reserve_in_chunks(seqindex_, seqindex_.size() + 1);
  seqinfo_t record;
  record.headerlen = static_cast<unsigned int>(headerlength);
  record.seqlen = static_cast<unsigned int>(sequencelength);
  record.header_p = header_p;
  record.seq_p = sequence_p;
  record.qual_p = quality_p;
  record.size = static_cast<uint64_t>(abundance);
  seqindex_.push_back(record);

  /* update statistics */
  ++sequences;
  nucleotides += sequencelength;
  longest = std::max(static_cast<uint64_t>(sequencelength), longest);
  shortest = std::min(static_cast<uint64_t>(sequencelength), shortest);
  longestheader = std::max(static_cast<uint64_t>(headerlength), longestheader);
}


auto Database::read(const char * filename, int upcase, struct Parameters const & parameters) -> void
{
  fastx_handle h = fastx_open(filename, parameters);
  /* Own the handle for the duration of the read loop: fastx_next() can fatal()
     on a malformed record, which unwinds in a library session. The guard frees
     the handle on that path; on the normal path it is released to the explicit
     fastx_close() below (which also emits the stripped-character warning). */
  std::unique_ptr<fastx_s> handle_guard(h);

  fastq_format = fastx_is_fastq(h);

  int64_t const filesize = static_cast<int64_t>(fastx_get_size(h));

  std::string const prompt = std::string("Reading file ") + filename;

  longest = 0;
  shortest = std::numeric_limits<uint64_t>::max();  // refactoring: direct initialization
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;

  int64_t discarded_short = 0;
  int64_t discarded_long = 0;
  int64_t discarded_unoise = 0;

  /* start from empty buffers; the records are appended by add() below */
  data_.clear();
  seqindex_.clear();

  {
    Progress progress(prompt.c_str(), static_cast<uint64_t>(filesize), parameters);
    while (fastx_next(h,
                     not parameters.opt_notrunclabels,
                      (upcase != 0) ? chrmap_upcase() : chrmap_no_change()))
      {
        size_t const sequencelength = fastx_get_sequence_length(h);
        int64_t const abundance = fastx_get_abundance(h);

        /* opt_minseqlength defaults to the -1 "unset" sentinel, which the CLI
           resolves to a command-specific value (1 or 32) before read() runs.
           A library caller may leave it unset, so guard the cast: a non-positive
           minimum means "no lower bound". Without this, static_cast<size_t>(-1)
           is SIZE_MAX and every sequence is discarded as too short. */
        if ((parameters.opt_minseqlength > 0) and
            (sequencelength < static_cast<size_t>(parameters.opt_minseqlength)))
          {
            ++discarded_short;
          }
        else if (sequencelength > static_cast<size_t>(parameters.opt_maxseqlength))
          {
            ++discarded_long;
          }
        else if ((parameters.opt_cluster_unoise != nullptr) && (abundance < parameters.opt_minsize))
          {
            ++discarded_unoise;
          }
        else
          {
            add(fastq_format,
                fastx_get_header(h),
                fastx_get_sequence(h),
                fastq_format ? fastx_get_quality(h) : nullptr,
                fastx_get_header_length(h),
                sequencelength,
                abundance);
          }
        progress.update(fastx_get_position(h));
      }
  }
  static_cast<void>(handle_guard.release());  // normal path: hand the handle to fastx_close()
  fastx_close(h, parameters);

  if (not parameters.opt_quiet)
    {
      if (sequences > 0)
        {
          std::fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, "
                  "min %" PRIu64 ", max %" PRIu64 ", avg %.0f\n",
                  getnucleotidecount(),
                  getsequencecount(),
                  getshortestsequence(),
                  getlongestsequence(),
                  static_cast<double>(getnucleotidecount()) / static_cast<double>(getsequencecount()));
        }
      else
        {
          std::fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  getnucleotidecount(),
                  getsequencecount());
        }
    }

  if (parameters.opt_log != nullptr)
    {
      if (sequences > 0)
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, "
                  "min %" PRIu64 ", max %" PRIu64 ", avg %.0f\n\n",
                  getnucleotidecount(),
                  getsequencecount(),
                  getshortestsequence(),
                  getlongestsequence(),
                  static_cast<double>(getnucleotidecount()) / static_cast<double>(getsequencecount()));
        }
      else
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n\n",
                  getnucleotidecount(),
                  getsequencecount());
        }
    }

  /* Warn about discarded sequences */

  if (discarded_short != 0)
    {
      std::fprintf(stderr,
              "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              parameters.opt_minseqlength,
              discarded_short,
              (discarded_short == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  parameters.opt_minseqlength,
                  discarded_short,
                  (discarded_short == 1 ? "sequence" : "sequences"));
        }
    }

  if (discarded_long != 0)
    {
      std::fprintf(stderr,
              "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              parameters.opt_maxseqlength,
              discarded_long,
              (discarded_long == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  parameters.opt_maxseqlength,
                  discarded_long,
                  (discarded_long == 1 ? "sequence" : "sequences"));
        }
    }

    if (discarded_unoise != 0)
    {
      std::fprintf(stderr,
              "minsize %" PRId64 ": %" PRId64 " %s discarded.\n",
              parameters.opt_minsize,
              discarded_unoise,
              (discarded_unoise == 1 ? "sequence" : "sequences"));

      if (parameters.opt_log != nullptr)
        {
          std::fprintf(parameters.fp_log,
                  "minsize %" PRId64 ": %" PRId64 " %s discarded.\n",
                  parameters.opt_minsize,
                  discarded_unoise,
                  (discarded_unoise == 1 ? "sequence" : "sequences"));
        }
    }
}


auto Database::clear() -> void
{
  /* release the heap storage (clear() alone keeps capacity); the vectors' own
     destructors would do this too, but callers use clear() to free a database
     early. */
  data_.clear();
  data_.shrink_to_fit();
  seqindex_.clear();
  seqindex_.shrink_to_fit();
  sequences = 0;
  nucleotides = 0;
  longest = 0;
  shortest = std::numeric_limits<uint64_t>::max();
  longestheader = 0;
}


Database::~Database()
{
  clear();
}


/* The sort comparators need the data buffer to compare header strings. Each is
   a lambda capturing the Database (this), so no file-scope pointer is needed
   (the old std::qsort path used one because a C comparator cannot capture, and
   qsort_r is non-portable). header_p increases with input order, so comparing
   it is a deterministic, stable tie-break that keeps equal records in their
   original order (it replaces the element-address comparison the qsort version
   relied on). */

auto Database::sortbylength(struct Parameters const & parameters) -> void
{
  Progress const progress("Sorting by length", 100, parameters);

  /* longest first, then by abundance, then by label, otherwise keep order */
  auto const * const buffer = data_.data();
  auto const by_length = [buffer](seqinfo_t const & lhs, seqinfo_t const & rhs) -> bool
  {
    if (lhs.seqlen != rhs.seqlen) { return lhs.seqlen > rhs.seqlen; }
    if (lhs.size != rhs.size) { return lhs.size > rhs.size; }
    auto const order = std::strcmp(buffer + lhs.header_p, buffer + rhs.header_p);
    if (order != 0) { return order < 0; }
    return lhs.header_p < rhs.header_p;
  };

  std::sort(seqindex_.begin(), seqindex_.end(), by_length);
}


auto Database::sortbylength_shortest_first(struct Parameters const & parameters) -> void
{
  Progress const progress("Sorting by length", 100, parameters);

  /* shortest first, then by abundance, then by label, otherwise keep order */
  auto const * const buffer = data_.data();
  auto const by_length_shortest = [buffer](seqinfo_t const & lhs, seqinfo_t const & rhs) -> bool
  {
    if (lhs.seqlen != rhs.seqlen) { return lhs.seqlen < rhs.seqlen; }
    if (lhs.size != rhs.size) { return lhs.size > rhs.size; }
    auto const order = std::strcmp(buffer + lhs.header_p, buffer + rhs.header_p);
    if (order != 0) { return order < 0; }
    return lhs.header_p < rhs.header_p;
  };

  std::sort(seqindex_.begin(), seqindex_.end(), by_length_shortest);
}


auto Database::sortbyabundance(struct Parameters const & parameters) -> void
{
  Progress const progress("Sorting by abundance", 100, parameters);

  /* most abundant first, then by label, otherwise keep order */
  auto const * const buffer = data_.data();
  auto const by_abundance = [buffer](seqinfo_t const & lhs, seqinfo_t const & rhs) -> bool
  {
    if (lhs.size != rhs.size) { return lhs.size > rhs.size; }
    auto const order = std::strcmp(buffer + lhs.header_p, buffer + rhs.header_p);
    if (order != 0) { return order < 0; }
    return lhs.header_p < rhs.header_p;
  };

  std::sort(seqindex_.begin(), seqindex_.end(), by_abundance);
}
