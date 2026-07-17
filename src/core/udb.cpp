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
#include "core/buffer_headroom.hpp"
#include "utils/progress.hpp"
#include "core/attributes.hpp"
#include "core/bitmap.hpp"
#include "core/db.hpp"  // Database, seqinfo_t
#include "core/dbindex.hpp"
#include "core/unique.hpp"
#include "os/system.hpp"  // xstat_t, xstat, xfstat, xmalloc, S_ISREG, S_ISFIFO
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"
#include <algorithm>  // std::min, std::max
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <cstring>  // std::memset
#include <fstream>  // std::ifstream
#include <istream>  // std::istream
#include <limits>
#include <string>  // std::string
#include <vector>


constexpr auto blocksize = uint64_t{4096UL * 4096UL};

// The .udb binary format is read and written in host byte order with no
// byteswapping, so a database is portable and correctly parsed only on a
// little-endian host. Fail the build on a big-endian target rather than
// silently reading a .udb wrong (see the little-endian note in sff_convert.cc).
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
              "udb.cc assumes a little-endian host");
#endif


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto largeread(std::istream & input, void * buf, uint64_t nbyte, uint64_t offset, Progress & progress_bar) -> uint64_t
  {
    /* read the file in blocks and update progress */

    auto progress = offset;
    for (uint64_t i = 0; i < nbyte; i += blocksize)
      {
        auto const rem = std::min(blocksize, nbyte - i);
        input.read((static_cast<char *>(buf)) + i, static_cast<std::streamsize>(rem));
        if (static_cast<uint64_t>(input.gcount()) != rem)
          {
            fatal("Unable to read from UDB file or invalid UDB file");
          }

        progress += rem;
        progress_bar.update(progress);
      }
    return nbyte;
  }

}  // end of anonymous namespace


auto udb_detect_isudb(const char * filename) -> bool
{
  /*
    Detect whether the given filename seems to refer to an UDB file.
    It must be an uncompressed regular file, not a pipe.
  */

  constexpr static uint32_t udb_file_signature {0x55444246}; // 'FBDU UDBF'
  constexpr static uint64_t expected_n_bytes {sizeof(uint32_t)};

  /* Only a regular file can be probed here and then reopened from the
     start by the actual reader. A non-rewindable stream (a named pipe,
     or a character device such as FreeBSD's /dev/stdin and the /dev/fd/N
     entries created by shell process substitution) cannot be a UDB file,
     and reading its magic number would consume bytes the subsequent
     reader could not recover. Stat the open descriptor, not the path: on
     FreeBSD stat() of the path misreports such streams (/dev/stdin is a
     character device there, not a pipe), whereas fstat() on the opened
     descriptor reports the underlying pipe. open()+close() without a
     read() does not consume pipe data, so bailing out for anything that
     is not a regular file leaves the stream intact for the reader.
     open_input_file() also maps "-" to a duplicate of stdin (matching the
     reader), whereas stat() of the literal path "-" would fail. */

  auto const input = open_input_file(filename);
  if (not input)
    {
      fatal("Unable to open input file for reading (%s)", filename);
    }

  xstat_t fs;
  if (xfstat(fileno(input.get()), & fs) != 0)
    {
      fatal("Unable to get status for input file (%s)", filename);
    }

  if (not S_ISREG(fs.st_mode))
    {
      return false;
    }

  unsigned int magic = 0;
  auto const bytesread = std::fread(& magic, 1, static_cast<std::size_t>(expected_n_bytes), input.get());

  if ((static_cast<uint64_t>(bytesread) == expected_n_bytes) and (magic == udb_file_signature))
    {
      return true;
    }

  return false;
}


/* Validate-on-load helpers for untrusted UDB header fields.

   The values below are read verbatim from the file, so a crafted or
   corrupt UDB must be rejected with a clear error rather than allowed to
   drive an out-of-bounds allocation, index or write. There is no
   recoverable error channel (fatal() terminates the process), so a
   violation fatals. */

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto udb_checked_add(uint64_t const lhs, uint64_t const rhs) -> uint64_t
  {
    if (lhs > std::numeric_limits<uint64_t>::max() - rhs)
      {
        fatal("Invalid UDB file");
      }
    return lhs + rhs;
  }

}  // end of anonymous namespace


auto udb_read(const char * filename,
              bool const create_bitmaps,
              bool const parse_abundances,
              struct Dbindex & dbindex,
              struct Database & db,
              struct Parameters const & parameters) -> void
{
  /* read UDB as indexed database */

  auto seqcount = 0U;
  auto udb_wordlength = 0U;
  uint64_t nucleotides = 0;

  /* udb_read fills the reserved database buffers in place (it bypasses
     Database::add). These raw pointers are bound to the passed-in Database's
     vector storage right after udb_reserve() sizes it below; the buffers are
     not resized again during the load, so the pointers stay valid throughout. */
  char * datap = nullptr;
  seqinfo_t * seqindex = nullptr;

  xstat_t fs;
  if (xstat(filename, & fs) != 0)
    {
      fatal("Unable to get status for input file (%s)", filename);
    }

  auto const is_pipe = S_ISFIFO(fs.st_mode);
  if (is_pipe)
    {
      fatal("Cannot read UDB file from a pipe");
    }

  /* get file size */

  uint64_t const filesize = static_cast<uint64_t>(fs.st_size);

  /* open UDB file */

  std::ifstream in_stream(filename, std::ios::binary);
  if (not in_stream)
    {
      fatal("Unable to open UDB file for reading");
    }

  std::string const prompt = std::string("Reading UDB file ") + filename;


  /* header */

  std::array<unsigned int, 50> buffer {{}};
  uint64_t pos = 0;

  uint64_t longestheader = 0;
  auto shortest = std::numeric_limits<unsigned int>::max();
  auto longest = 0U;
  {
    Progress progress_bar(prompt.c_str(), filesize, parameters);
    pos += largeread(in_stream, buffer.data(), 4 * 50, pos, progress_bar);

    if ((buffer[0]  != 0x55444246) or
        (buffer[2] != 32) or
        (buffer[4] < 3) or
        (buffer[4] > 15) or
        (buffer[13] == 0) or
        (buffer[17] != 0x0000746e) or
        (buffer[49] != 0x55444266))
      {
        fatal("Invalid UDB file");
      }

    udb_wordlength = buffer[4];
    seqcount = buffer[13];
    dbindex.dbaccel = buffer[6];

    /* The per-sequence header-index and length tables each store 4 bytes
       per sequence, so a file cannot describe more than filesize/4
       sequences. Rejecting a larger seqcount also keeps it well clear of
       the seqcount + 1 wrap when the header index is sized below. */

    if (seqcount > filesize / 4)
      {
        fatal("Invalid UDB file");
      }

    /* The index is built at the UDB file's own word length. Publish it as the
       effective index width (read by the query-k-mer extractors) rather than
       mutating the opt_wordlength config global (E1); warn when it overrides the
       configured value. */
    if (udb_wordlength != static_cast<unsigned int>(parameters.opt_wordlength))
      {
        std::fprintf(stderr, "\nWARNING: Wordlength adjusted to %u as indicated in UDB file\n", udb_wordlength);
      }
    dbindex.wordlength = udb_wordlength;

    /* word match counts */

    dbindex.hashsize = 1U << (2 * udb_wordlength);
    dbindex.kmercount = static_cast<unsigned int *>(xmalloc(dbindex.hashsize * sizeof(unsigned int)));
    dbindex.kmerhash = static_cast<uint64_t *>(xmalloc(dbindex.hashsize * sizeof(uint64_t)));
    dbindex.kmerbitmap = std::vector<Bitmap>(dbindex.hashsize);

    pos += largeread(in_stream, dbindex.kmercount, 4 * dbindex.hashsize, pos, progress_bar);

    dbindex.indexsize = 0;
    for (uint64_t i = 0; i < dbindex.hashsize; i++)
      {
        dbindex.kmerhash[i] = dbindex.indexsize;
        dbindex.indexsize = udb_checked_add(dbindex.indexsize, dbindex.kmercount[i]);
      }

    /* The word-list section stores 4 bytes per index entry, so a file can
       hold at most filesize/4 entries; a larger total means the kmercount[]
       values do not match the on-disk section (padded/corrupt file). */

    if (dbindex.indexsize > filesize / 4)
      {
        fatal("Invalid UDB file");
      }

    /* signature */

    pos += largeread(in_stream, buffer.data(), 4, pos, progress_bar);

    if (buffer[0] != 0x55444233)
      {
        fatal("Invalid UDB file");
      }

    /* sequence numbers for word matches */

    dbindex.kmerindex = static_cast<unsigned int *>(xmalloc(dbindex.indexsize * 4));

    pos += largeread(in_stream, dbindex.kmerindex, 4 * dbindex.indexsize, pos, progress_bar);

    /* Every entry is a sequence number used both as a bit offset in the
       per-word bitmaps (Bitmap::set writes bitmap[value >> 3], no bounds
       check) and as an index into seqindex/dbindex_map during search. A
       value >= seqcount is therefore an out-of-bounds write or read, so
       reject it here rather than at use. */

    for (uint64_t i = 0; i < dbindex.indexsize; i++)
      {
        if (dbindex.kmerindex[i] >= seqcount)
          {
            fatal("Invalid UDB file");
          }
      }

    /* new header */

    pos += largeread(in_stream, buffer.data(), 4 * 8, pos, progress_bar);

    if ((buffer[0] != 0x55444234) or
        (buffer[1] != 0x005e0db3) or
        (buffer[2] != seqcount) or
        (buffer[7] != 0x005e0db4))
      {
        fatal("Invalid UDB file");
      }

    nucleotides = ((static_cast<uint64_t>(buffer[4])) << 32U) | buffer[3];
    auto const udb_headerchars = ((static_cast<uint64_t>(buffer[6])) << 32U) | buffer[5];

    /* allocate the two database buffers up front; udb_read fills them in place */

    uint64_t const datap_bytes =
      udb_checked_add(udb_checked_add(udb_headerchars, nucleotides), seqcount);
    db.udb_reserve(seqcount, datap_bytes);
    datap = db.data_.data();
    seqindex = db.seqindex_.data();

    /* header index */

    std::vector<unsigned int> header_index(seqcount + 1);

    pos += largeread(in_stream, header_index.data(), 4 * seqcount, pos, progress_bar);

    header_index[seqcount] = static_cast<unsigned int>(udb_headerchars);

    auto last = 0U;
    for (auto i = 0U; i < seqcount; i++)
      {
        unsigned int const current_index = header_index[i];
        if ((current_index < last) or (current_index >= udb_headerchars))
          {
            fatal("Invalid UDB file");
          }
        /* Header offsets must strictly increase: an equal (or smaller) next
           offset would make headerlen (next - current - 1) underflow. */
        if (header_index[i + 1] <= current_index)
          {
            fatal("Invalid UDB file");
          }
        seqindex[i].header_p = current_index;
        seqindex[i].headerlen = header_index[i + 1] - current_index - 1;
        if (static_cast<int64_t>(seqindex[i].headerlen) > std::numeric_limits<int>::max() - buffer_headroom)
          {
            fatal("UDB file contains a header too long for this version of vsearch");
          }
        seqindex[i].size = 1;
        last = current_index;
      }


    /* headers */

    pos += largeread(in_stream, datap, udb_headerchars, pos, progress_bar);

    for (auto i = 0U; i < seqcount; i++)
      {
        longestheader = std::max<uint64_t>(seqindex[i].headerlen, longestheader);
      }

    /* sequence lengths */

    std::vector<unsigned int> sequence_lengths(seqcount);

    pos += largeread(in_stream, sequence_lengths.data(), 4 * seqcount, pos, progress_bar);

    uint64_t sum = 0;

    for (auto i = 0U; i < seqcount; i++)
      {
        unsigned int const sequence_length = sequence_lengths[i];

        if (static_cast<int64_t>(sequence_length) > std::numeric_limits<int>::max() - buffer_headroom)
          {
            fatal("UDB file contains a sequence too long for this version of vsearch");
          }

        seqindex[i].seq_p = udb_headerchars + sum;
        seqindex[i].seqlen = sequence_length;
        seqindex[i].qual_p = 0;

        shortest = std::min(sequence_length, shortest);
        longest = std::max(sequence_length, longest);

        sum += sequence_length;

        if (sum > nucleotides)
          {
            fatal("Invalid UDB file");
          }
      }


    if (sum != nucleotides)
      {
        fatal("Invalid UDB file");
      }

    /* sequences */

    pos += largeread(in_stream, datap + udb_headerchars, nucleotides, pos, progress_bar);

    if (pos != filesize)
      {
        fatal("Incorrect UDB file size");
      }

    /* close UDB file */

    in_stream.close();
  }

  /* reorganize the sequences in memory and record the database statistics */

  db.udb_finalize(seqcount, nucleotides, longest, shortest, longestheader, parameters);

  /* Create bitmaps for the most frequent words */

  if (create_bitmaps)
    {
      auto const bitmap_mincount = seqcount / 8;
      {
        Progress progress("Creating bitmaps", dbindex.hashsize, parameters);
        for (auto i = 0U; i < dbindex.hashsize; i++)
          {
            if (dbindex.kmercount[i] >= bitmap_mincount)
              {
                dbindex.kmerbitmap[i] = Bitmap(seqcount + 127); // pad for xmm
                for (auto j = 0U; j < dbindex.kmercount[i]; j++)
                  {
                    dbindex.kmerbitmap[i].set(dbindex.kmerindex[dbindex.kmerhash[i]+j]);
                  }
              }
            progress.update(i + 1);
          }
      }
    }

  /* get abundances and longest header */

  if (parse_abundances)
    {
      {
        Progress progress("Parsing abundances", seqcount, parameters);
        for (auto i = 0U; i < seqcount; i++)
          {
            auto const size = header_get_size(datap + seqindex[i].header_p,
                                           static_cast<int>(seqindex[i].headerlen));
            if (size > 0)
              {
                seqindex[i].size = static_cast<uint64_t>(size);
              }
            else
              {
                seqindex[i].size = 1;
              }
            progress.update(i + 1);
          }
      }
    }

  /* the unique-kmer finder (dbindex.uhandle) is a Uniquer value member, ready to
     use as default-constructed; the UDB path does not build the index with it */

  /* make mapping from indexno to seqno */

  dbindex.map = static_cast<unsigned int *>(xmalloc(seqcount * sizeof(unsigned int)));
  dbindex.count = seqcount;

  for (auto i = 0U; i < seqcount; i++)
    {
      dbindex.map[i] = i;
    }

  /* done */

  /* some stats */

  if (not parameters.opt_quiet)
    {
      if (seqcount > 0)
        {
          std::fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, min %" PRIu64 ", max %" PRIu64 ", avg %.0f\n",
                  db.getnucleotidecount(),
                  db.getsequencecount(),
                  db.getshortestsequence(),
                  db.getlongestsequence(),
                  static_cast<double>(db.getnucleotidecount()) * 1.0 / static_cast<double>(db.getsequencecount()));
        }
      else
        {
          std::fprintf(stderr,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n",
                  db.getnucleotidecount(),
                  db.getsequencecount());
        }
    }

  if (parameters.opt_log != nullptr)
    {
      if (seqcount > 0)
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs, min %" PRIu64 ", max %" PRIu64 ", avg %.0f\n\n",
                  db.getnucleotidecount(),
                  db.getsequencecount(),
                  db.getshortestsequence(),
                  db.getlongestsequence(),
                  static_cast<double>(db.getnucleotidecount()) * 1.0 / static_cast<double>(db.getsequencecount()));
        }
      else
        {
          std::fprintf(parameters.fp_log,
                  "%" PRIu64 " nt in %" PRIu64 " seqs\n\n",
                  db.getnucleotidecount(),
                  db.getsequencecount());
        }
    }
}
