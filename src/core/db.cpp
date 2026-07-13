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

#include "vsearch.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/progress.hpp"
#include "utils/string_alloc.hpp"
#include <algorithm>  // std::min, std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <cstdlib>  // std::qsort
#include <cstring>  // std::memcpy, std::strcmp
#include <limits>


constexpr uint64_t memchunk = 16777216;  // 2^24


/* Reset database state for library use.
   Must be called before the first add() when not using read().
   Frees any previously allocated data, then resets all counters.
   Sets shortest to UINT64_MAX so min-tracking works correctly. */
auto Database::init() -> void
{
  clear();
  is_fastq = false;
  sequences = 0;
  nucleotides = 0;
  longest = 0;
  shortest = std::numeric_limits<uint64_t>::max();
  longestheader = 0;
  dataalloc = 0;
  datalen = 0;
  seqindex_alloc = 0;
}


auto Database::udb_reserve(uint64_t const count, uint64_t const datap_bytes) -> void
{
  /* udb_read fills these buffers in place (it bypasses add()); allocate them
     up front at the exact sizes derived from the UDB header. */
  seqindex = static_cast<seqinfo_t *>(xmalloc(count * sizeof(seqinfo_t)));
  datap = static_cast<char *>(xmalloc(datap_bytes));
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
    for (auto i = count - 1; i > 0; i--)
      {
        auto const old_p = seqindex[i].seq_p;
        auto const new_p = seqindex[i].seq_p + i;
        auto const len   = seqindex[i].seqlen;
        std::memmove(datap + new_p, datap + old_p, len);
        *(datap + new_p + len) = 0;
        seqindex[i].seq_p = new_p;
        progress.update(count - i);
      }
    *(datap + seqindex[0].seq_p + seqindex[0].seqlen) = 0;
  }

  /* record the summary statistics (a UDB database is never FASTQ) */

  is_fastq = false;
  sequences = count;
  nucleotides = nucleotide_count;
  longest = longest_sequence;
  shortest = shortest_sequence;
  longestheader = longest_header;
}


auto Database::getquality(uint64_t seqno) const -> char *
{
  if (is_fastq)
    {
      return datap + seqindex[seqno].qual_p;
    }
  return nullptr;
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

  /* grow space for data, if necessary */

  size_t const dataalloc_old = dataalloc;

  size_t needed = datalen + headerlength + 1 + sequencelength + 1;
  if (is_fastq_record)
    {
      needed += sequencelength + 1;
    }
  while (dataalloc < needed)
    {
      dataalloc += memchunk;
    }
  if (dataalloc > dataalloc_old)
    {
      datap = static_cast<char *>(xrealloc(datap, dataalloc));
    }

  /* store the header */
  size_t const header_p = datalen;
  std::memcpy(datap + header_p,
              header,
              headerlength + 1);
  datalen += headerlength + 1;

  /* store sequence */
  size_t const sequence_p = datalen;
  std::memcpy(datap + sequence_p,
              sequence,
              sequencelength + 1);
  datalen += sequencelength + 1;

  size_t const quality_p = datalen;
  if (is_fastq_record)
    {
      /* store quality */
      std::memcpy(datap + quality_p,
                  quality,
                  sequencelength + 1);
      datalen += sequencelength + 1;

      /* A FASTQ record makes this a FASTQ database. read() sets the is_fastq
         flag itself before adding records, but callers that build a database
         directly with add() rely on this assignment so that the is_fastq flag
         and getquality() can reach the stored quality. */
      is_fastq = true;
    }

  /* grow space for index, if necessary */
  size_t const seqindex_alloc_old = seqindex_alloc;
  while ((sequences + 1) * sizeof(seqinfo_t) > seqindex_alloc)
    {
      seqindex_alloc += memchunk;
    }
  if (seqindex_alloc > seqindex_alloc_old)
    {
      seqindex = static_cast<seqinfo_t *>(xrealloc(seqindex, seqindex_alloc));
    }

  /* update index */
  seqinfo_t * seqindex_p = seqindex + sequences;
  seqindex_p->headerlen = static_cast<unsigned int>(headerlength);
  seqindex_p->seqlen = static_cast<unsigned int>(sequencelength);
  seqindex_p->header_p = header_p;
  seqindex_p->seq_p = sequence_p;
  seqindex_p->qual_p = quality_p;
  seqindex_p->size = static_cast<uint64_t>(abundance);

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

  is_fastq = fastx_is_fastq(h);

  int64_t const filesize = static_cast<int64_t>(fastx_get_size(h));

  char * prompt = nullptr;
  if (xsprintf(&prompt, "Reading file %s", filename) == -1)
    {
      fatal("Out of memory");
    }

  longest = 0;
  shortest = std::numeric_limits<uint64_t>::max();  // refactoring: direct initialization
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;

  int64_t discarded_short = 0;
  int64_t discarded_long = 0;
  int64_t discarded_unoise = 0;

  /* allocate space for data */
  dataalloc = 0;
  datap = nullptr;
  datalen = 0;

  /* allocate space for index */
  seqindex_alloc = 0;
  seqindex = nullptr;

  {
    Progress progress(prompt, static_cast<uint64_t>(filesize), parameters);
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
            add(is_fastq,
                fastx_get_header(h),
                fastx_get_sequence(h),
                is_fastq ? fastx_get_quality(h) : nullptr,
                fastx_get_header_length(h),
                sequencelength,
                abundance);
          }
        progress.update(fastx_get_position(h));
      }
  }
  xfree(prompt);
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
  if (datap != nullptr)
    {
      xfree(datap);
      datap = nullptr;
    }
  if (seqindex != nullptr)
    {
      xfree(seqindex);
      seqindex = nullptr;
    }
  sequences = 0;
  nucleotides = 0;
  longest = 0;
  shortest = std::numeric_limits<uint64_t>::max();
  longestheader = 0;
  dataalloc = 0;
  datalen = 0;
  seqindex_alloc = 0;
}


Database::~Database()
{
  clear();
}


/* The std::qsort comparators below need the data buffer to compare header
   strings, but a C comparator is a plain function pointer that cannot capture
   the Database being sorted (qsort_r, which passes context, is non-portable).
   The sort members set this file-scope pointer to their own datap immediately
   before calling qsort; sorting is single-threaded, so the transient sharing is
   safe. Switching to std::sort with a capturing comparator would remove it (see
   TBD_20260713_Database_polish.md). */
namespace {
  char const * sort_datap = nullptr;
}


auto compare_bylength(const void * a, const void * b) -> int
{
  auto const * lhs = static_cast<seqinfo_t const *>(a);
  auto const * rhs = static_cast<seqinfo_t const *>(b);

  /* longest first, then by abundance, then by label, otherwise keep order */

  if (lhs->seqlen < rhs->seqlen)
    {
      return +1;
    }
  if (lhs->seqlen > rhs->seqlen)
    {
      return -1;
    }

  if (lhs->size < rhs->size)
    {
      return +1;
    }
  if (lhs->size > rhs->size)
    {
      return -1;
    }

  auto const result = std::strcmp(sort_datap + lhs->header_p, sort_datap + rhs->header_p);
  if (result != 0)
    {
      return result;
    }

  if (lhs < rhs)
    {
      return -1;
    }
  if (lhs > rhs)
    {
      return +1;
    }
  return 0;
}


auto compare_bylength_shortest_first(const void * a, const void * b) -> int
{
  auto const * lhs = static_cast<seqinfo_t const *>(a);
  auto const * rhs = static_cast<seqinfo_t const *>(b);

  /* shortest first, then by abundance, then by label, otherwise keep order */

  if (lhs->seqlen < rhs->seqlen)
    {
      return -1;
    }
  if (lhs->seqlen > rhs->seqlen)
    {
      return +1;
    }

  if (lhs->size < rhs->size)
    {
      return +1;
    }
  if (lhs->size > rhs->size)
    {
      return -1;
    }

  auto const result = std::strcmp(sort_datap + lhs->header_p, sort_datap + rhs->header_p);
  if (result != 0)
    {
      return result;
    }

  if (lhs < rhs)
    {
      return -1;
    }
  if (lhs > rhs)
    {
      return +1;
    }
  return 0;
}


inline auto compare_byabundance(const void * a, const void * b) -> int
{
  auto const * lhs = static_cast<seqinfo_t const *>(a);
  auto const * rhs = static_cast<seqinfo_t const *>(b);

  /* most abundant first, then by label, otherwise keep order */

  if (lhs->size > rhs->size)
    {
      return -1;
    }
  if (lhs->size < rhs->size)
    {
      return +1;
    }

  auto const result = std::strcmp(sort_datap + lhs->header_p, sort_datap + rhs->header_p);
  if (result != 0)
    {
      return result;
    }

  if (lhs < rhs)
    {
      return -1;
    }
  if (lhs > rhs)
    {
      return +1;
    }
  return 0;
}


auto Database::sortbylength(struct Parameters const & parameters) -> void
{
  Progress const progress("Sorting by length", 100, parameters);
  sort_datap = datap;
  if (sequences > 0)  // qsort requires a non-null pointer even for zero elements
    {
      std::qsort(seqindex,
            sequences,
            sizeof(seqinfo_t),
            compare_bylength);
    }
}


auto Database::sortbylength_shortest_first(struct Parameters const & parameters) -> void
{
  Progress const progress("Sorting by length", 100, parameters);
  sort_datap = datap;
  if (sequences > 0)  // qsort requires a non-null pointer even for zero elements
    {
      std::qsort(seqindex,
            sequences,
            sizeof(seqinfo_t),
            compare_bylength_shortest_first);
    }
}


auto Database::sortbyabundance(struct Parameters const & parameters) -> void
{
  Progress const progress("Sorting by abundance", 100, parameters);
  sort_datap = datap;
  if (sequences > 0)  // qsort requires a non-null pointer even for zero elements
    {
      std::qsort(seqindex,
            sequences,
            sizeof(seqinfo_t),
            compare_byabundance);
    }
}
