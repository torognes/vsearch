/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "maps.h"
#include <cinttypes>  // macros PRIu64 and PRId64
#include <climits>  // LONG_MAX
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <cstdlib>  // std::qsort
#include <cstring>  // std::memcpy, std::strcmp


constexpr uint64_t memchunk = 16777216;  // 2^24

static fastx_handle h = nullptr;
static bool is_fastq = false;
static uint64_t sequences = 0;
static uint64_t nucleotides = 0;
static uint64_t longest = 0;
static uint64_t shortest = 0;
static uint64_t longestheader = 0;

static uint64_t dataalloc = 0;
static uint64_t datalen = 0;
static size_t seqindex_alloc = 0;

seqinfo_t * seqindex = nullptr;
char * datap = nullptr;


auto db_setinfo(bool new_is_fastq,
                uint64_t new_sequences,
                uint64_t new_nucleotides,
                uint64_t new_longest,
                uint64_t new_shortest,
                uint64_t new_longestheader) -> void
{
  is_fastq = new_is_fastq;
  sequences = new_sequences;
  nucleotides = new_nucleotides;
  longest = new_longest;
  shortest = new_shortest;
  longestheader = new_longestheader;
}


auto db_is_fastq() -> bool
{
  return is_fastq;
}


auto db_getquality(uint64_t seqno) -> char *
{
  if (is_fastq)
    {
      return datap + seqindex[seqno].qual_p;
    }
  return nullptr;
}


auto db_add(bool is_fastq,
            char * header,
            char * sequence,
            char * quality,
            size_t headerlength,
            size_t sequencelength,
            int64_t abundance) -> void
{
  /* Add a sequence to the database. Assumes that the database has been initialized. */

  /* grow space for data, if necessary */

  size_t const dataalloc_old = dataalloc;

  size_t needed = datalen + headerlength + 1 + sequencelength + 1;
  if (is_fastq)
    {
      needed += sequencelength + 1;
    }
  while (dataalloc < needed)
    {
      dataalloc += memchunk;
    }
  if (dataalloc > dataalloc_old)
    {
      datap = (char *) xrealloc(datap, dataalloc);
    }

  /* store the header */
  size_t const header_p = datalen;
  memcpy(datap + header_p,
         header,
         headerlength + 1);
  datalen += headerlength + 1;

  /* store sequence */
  size_t const sequence_p = datalen;
  memcpy(datap + sequence_p,
         sequence,
         sequencelength + 1);
  datalen += sequencelength + 1;

  size_t const quality_p = datalen;
  if (is_fastq)
    {
      /* store quality */
      memcpy(datap + quality_p,
             quality,
             sequencelength + 1);
      datalen += sequencelength + 1;
    }

  /* grow space for index, if necessary */
  size_t const seqindex_alloc_old = seqindex_alloc;
  while ((sequences + 1) * sizeof(seqinfo_t) > seqindex_alloc)
    {
      seqindex_alloc += memchunk;
    }
  if (seqindex_alloc > seqindex_alloc_old)
    {
      seqindex = (seqinfo_t *) xrealloc(seqindex, seqindex_alloc);
    }

  /* update index */
  seqinfo_t * seqindex_p = seqindex + sequences;
  seqindex_p->headerlen = headerlength;
  seqindex_p->seqlen = sequencelength;
  seqindex_p->header_p = header_p;
  seqindex_p->seq_p = sequence_p;
  seqindex_p->qual_p = quality_p;
  seqindex_p->size = abundance;

  /* update statistics */
  ++sequences;
  nucleotides += sequencelength;
  if (sequencelength > longest)
    {
      longest = sequencelength;
    }
  if (sequencelength < shortest)
    {
      shortest = sequencelength;
    }
  if (headerlength > longestheader)
    {
      longestheader = headerlength;
    }
}


auto db_read(const char * filename, int upcase) -> void
{
  h = fastx_open(filename);

  if (not h)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  is_fastq = fastx_is_fastq(h);

  int64_t const filesize = fastx_get_size(h);

  char * prompt = nullptr;
  if (xsprintf(& prompt, "Reading file %s", filename) == -1)
    {
      fatal("Out of memory");
    }

  progress_init(prompt, filesize);

  longest = 0;
  shortest = LONG_MAX;
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

  while(fastx_next(h,
                   not opt_notrunclabels,
                   upcase ? chrmap_upcase : chrmap_no_change))
    {
      size_t const sequencelength = fastx_get_sequence_length(h);
      int64_t const abundance = fastx_get_abundance(h);

      if (sequencelength < (size_t) opt_minseqlength)
        {
          ++discarded_short;
        }
      else if (sequencelength > (size_t) opt_maxseqlength)
        {
          ++discarded_long;
        }
      else if (opt_cluster_unoise && (abundance < opt_minsize))
        {
          ++discarded_unoise;
        }
      else
        {
          db_add(is_fastq,
                 fastx_get_header(h),
                 fastx_get_sequence(h),
                 is_fastq ? fastx_get_quality(h) : nullptr,
                 fastx_get_header_length(h),
                 sequencelength,
                 abundance);
        }
      progress_update(fastx_get_position(h));
    }

  progress_done();
  xfree(prompt);
  fastx_close(h);

  if (not opt_quiet)
    {
      if (sequences > 0)
        {
          fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, "
                  "min %" PRIu64 ", max %" PRIu64 ", avg %.0f\n",
                  db_getnucleotidecount(),
                  db_getsequencecount(),
                  db_getshortestsequence(),
                  db_getlongestsequence(),
                  db_getnucleotidecount() * 1.0 / db_getsequencecount());
        }
      else
        {
          fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  db_getnucleotidecount(),
                  db_getsequencecount());
        }
    }

  if (opt_log)
    {
      if (sequences > 0)
        {
          fprintf(fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, "
                  "min %" PRIu64 ", max %" PRIu64 ", avg %.0f\n\n",
                  db_getnucleotidecount(),
                  db_getsequencecount(),
                  db_getshortestsequence(),
                  db_getlongestsequence(),
                  db_getnucleotidecount() * 1.0 / db_getsequencecount());
        }
      else
        {
          fprintf(fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n\n",
                  db_getnucleotidecount(),
                  db_getsequencecount());
        }
    }

  /* Warn about discarded sequences */

  if (discarded_short)
    {
      fprintf(stderr,
              "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_minseqlength,
              discarded_short,
              (discarded_short == 1 ? "sequence" : "sequences"));

      if (opt_log)
        {
          fprintf(fp_log,
                  "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  opt_minseqlength,
                  discarded_short,
                  (discarded_short == 1 ? "sequence" : "sequences"));
        }
    }

  if (discarded_long)
    {
      fprintf(stderr,
              "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_maxseqlength,
              discarded_long,
              (discarded_long == 1 ? "sequence" : "sequences"));

      if (opt_log)
        {
          fprintf(fp_log,
                  "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
                  opt_maxseqlength,
                  discarded_long,
                  (discarded_long == 1 ? "sequence" : "sequences"));
        }
    }

    if (discarded_unoise)
    {
      fprintf(stderr,
              "minsize %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_minsize,
              discarded_unoise,
              (discarded_unoise == 1 ? "sequence" : "sequences"));

      if (opt_log)
        {
          fprintf(fp_log,
                  "minsize %" PRId64 ": %" PRId64 " %s discarded.\n",
                  opt_minsize,
                  discarded_unoise,
                  (discarded_unoise == 1 ? "sequence" : "sequences"));
        }
    }

  show_rusage();
}


auto db_getsequencecount() -> uint64_t
{
  return sequences;
}


auto db_getnucleotidecount() -> uint64_t
{
  return nucleotides;
}


auto db_getlongestheader() -> uint64_t
{
  return longestheader;
}


auto db_getlongestsequence() -> uint64_t
{
  return longest;
}


auto db_getshortestsequence() -> uint64_t
{
  return shortest;
}


auto db_free() -> void
{
  if (datap)
    {
      xfree(datap);
    }
  if (seqindex)
    {
      xfree(seqindex);
    }
}


auto compare_bylength(const void * a, const void * b) -> int
{
  auto * lhs = (seqinfo_t *) a;
  auto * rhs = (seqinfo_t *) b;

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

  auto const result = std::strcmp(datap + lhs->header_p, datap + rhs->header_p);
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
  auto * lhs = (seqinfo_t *) a;
  auto * rhs = (seqinfo_t *) b;

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

  auto const result = std::strcmp(datap + lhs->header_p, datap + rhs->header_p);
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
  auto * lhs = (seqinfo_t *) a;
  auto * rhs = (seqinfo_t *) b;

  /* most abundant first, then by label, otherwise keep order */

  if (lhs->size > rhs->size)
    {
      return -1;
    }
  if (lhs->size < rhs->size)
    {
      return +1;
    }

  int const r = std::strcmp(datap + lhs->header_p, datap + rhs->header_p);
  if (r != 0)
    {
      return r;
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


auto db_sortbylength() -> void
{
  progress_init("Sorting by length", 100);
  qsort(seqindex,
        sequences,
        sizeof(seqinfo_t),
        compare_bylength);
  progress_done();
}


auto db_sortbylength_shortest_first() -> void
{
  progress_init("Sorting by length", 100);
  qsort(seqindex,
        sequences,
        sizeof(seqinfo_t),
        compare_bylength_shortest_first);
  progress_done();
}


auto db_sortbyabundance() -> void
{
  progress_init("Sorting by abundance", 100);
  qsort(seqindex,
        sequences,
        sizeof(seqinfo_t),
        compare_byabundance);
  progress_done();
}
