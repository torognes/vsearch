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
#include "core/udb.hpp"  // UdbUse, udb_read, udb_detect_isudb
#include "os/system.hpp"  // xstat_t, xstat, xfstat, S_ISREG, S_ISFIFO
#include "utils/decimal_digits.hpp"  // decimal::to_text
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/span.hpp"  // Span, make_span
#include "utils/warn.hpp"  // vsearch::warn
#include <algorithm>  // std::min, std::max
#include <array>
#include <cassert>  // assert
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <fstream>  // std::ifstream
#include <ios>
#include <istream>  // std::istream
#include <limits>
#include <string>  // std::string
#include <sys/stat.h>
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

  /* The UDB format's integer fields are 4 bytes by definition of the format,
     not because sizeof(unsigned int) happens to be 4 on this host. Taking the
     byte count from the span below makes that dependency implicit, so state
     it here: every cross-compilation target vsearch supports has a 32-bit
     int, so none of them would catch it. */
  static_assert(sizeof(unsigned int) == 4, "UDB stores 32-bit fields");

  /* Read buf.size() elements from the UDB file. A span rather than a
     (void *, nbyte) pair: the byte count and the destination used to be
     separate arguments, so a mismatch -- a count in bytes where elements were
     meant, a buffer sized seqcount + 1 read as seqcount, a scratch vector
     reused for a shorter field -- produced a file that still looked valid and
     desynchronised every later offset. Now the extent comes from the
     destination and the call site says how much of it is being filled. */
  template <typename Type>
  auto largeread(std::istream & input, Span<Type> const buf, uint64_t const offset,
                 Progress & progress_bar) -> uint64_t
  {
    /* read the file in blocks and update progress */

    /* the destination as raw bytes: std::istream reads chars, whatever the
       span's element type is. as_writable_bytes() is the one audited rebind
       (see utils/span.hpp), and the block below is a subspan of it, so the
       destination and its length can no longer drift apart within the loop
       either. The local is const -- it is the Span that does not move, not
       the bytes it writes through -- which is also what cppcheck asked for. */
    auto const bytes = buf.as_writable_bytes();
    auto const nbyte = static_cast<uint64_t>(bytes.size());
    auto progress = offset;
    for (uint64_t i = 0; i < nbyte; i += blocksize)
      {
        auto const rem = std::min(blocksize, nbyte - i);
        auto const block = bytes.subspan(static_cast<std::size_t>(i),
                                         static_cast<std::size_t>(rem));
        input.read(block.data(), static_cast<std::streamsize>(block.size()));
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
      fatal(std::string("Unable to open input file for reading (")
            + std::string(filename)
            + ")");
    }

  xstat_t fs;
  if (xfstat(fileno(input.get()), & fs) != 0)
    {
      fatal(std::string("Unable to get status for input file (")
            + std::string(filename)
            + ")");
    }

  if (not S_ISREG(fs.st_mode))
    {
      return false;
    }

  unsigned int magic = 0;
  auto const bytesread = std::fread(& magic, 1, static_cast<std::size_t>(expected_n_bytes), input.get());

  return (static_cast<uint64_t>(bytesread) == expected_n_bytes) and (magic == udb_file_signature);
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
              UdbUse const usage,
              struct Dbindex & dbindex,
              struct Database & db,
              struct Parameters const & parameters) -> void
{
  /* read UDB as indexed database */

  auto seqcount = 0U;
  auto udb_wordlength = 0U;
  uint64_t nucleotides = 0;

  xstat_t fs;
  if (xstat(filename, & fs) != 0)
    {
      fatal(std::string("Unable to get status for input file (")
            + std::string(filename)
            + ")");
    }

  auto const is_pipe = S_ISFIFO(fs.st_mode);
  if (is_pipe)
    {
      fatal("Cannot read UDB file from a pipe");
    }

  /* get file size */

  auto const filesize = static_cast<uint64_t>(fs.st_size);

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
    Progress progress_bar(prompt, filesize, parameters);
    pos += largeread(in_stream, make_span(buffer), pos, progress_bar);

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
        vsearch::warn("Wordlength adjusted to " + decimal::to_text(udb_wordlength)
                      + " as indicated in UDB file");
      }
    dbindex.wordlength = udb_wordlength;

    /* word match counts */

    dbindex.hashsize = 1U << (2 * udb_wordlength);
    dbindex.kmercount.resize(dbindex.hashsize);
    dbindex.bitmap_slots_reset(dbindex.hashsize);
    /* filled by push_back into reserved space below, not resized here: the loop
       writes every entry, so value-initialising them first is 8 bytes per slot
       of zeros nobody reads */
    dbindex.kmerhash.clear();
    dbindex.kmerhash.reserve(dbindex.hashsize);

    pos += largeread(in_stream, make_span(dbindex.kmercount).first(dbindex.hashsize), pos, progress_bar);

    dbindex.indexsize = 0;
    for (uint64_t i = 0; i < dbindex.hashsize; i++)
      {
        dbindex.kmerhash.push_back(dbindex.indexsize);
        dbindex.indexsize = udb_checked_add(dbindex.indexsize, dbindex.kmercount[i]);
      }
    /* one entry per slot: size() counts what the loop wrote, so this checks it
       covered every slot (this path stores no end marker, unlike prepare) */
    assert(dbindex.kmerhash.size() == dbindex.hashsize);

    /* The word-list section stores 4 bytes per index entry, so a file can
       hold at most filesize/4 entries; a larger total means the kmercount[]
       values do not match the on-disk section (padded/corrupt file). */

    if (dbindex.indexsize > filesize / 4)
      {
        fatal("Invalid UDB file");
      }

    /* signature */

    pos += largeread(in_stream, make_span(buffer).first(1), pos, progress_bar);

    if (buffer[0] != 0x55444233)
      {
        fatal("Invalid UDB file");
      }

    /* sequence numbers for word matches */

    dbindex.kmerindex.resize(dbindex.indexsize);

    pos += largeread(in_stream, make_span(dbindex.kmerindex).first(dbindex.indexsize), pos, progress_bar);

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

    pos += largeread(in_stream, make_span(buffer).first(8), pos, progress_bar);

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

    /* udb_read fills the reserved database buffers in place (it bypasses
       Database::add). These raw pointers are bound to the passed-in Database's
       vector storage right after udb_reserve() sizes it below; the buffers are
       not resized again during the load, so the pointers stay valid throughout. */
    char * datap = nullptr;
    seqinfo_t * seqindex = nullptr;

    uint64_t const datap_bytes =
      udb_checked_add(udb_checked_add(udb_headerchars, nucleotides), seqcount);
    db.udb_reserve(seqcount, datap_bytes);
    datap = db.data_.data();
    seqindex = db.seqindex_.data();

    /* The reserved sequence-data buffer, as one span the two reads below take
       their slices from. Both extents (udb_headerchars, nucleotides) come out
       of the file, and datap_bytes is their checked sum, so this is where a
       length the buffer cannot hold is caught -- at the read, not at the use. */
    auto const data_buffer = Span<char>{datap, static_cast<std::size_t>(datap_bytes)};

    /* header index */

    std::vector<unsigned int> header_index(seqcount + 1);

    /* .first(seqcount), not the whole vector: header_index holds
       seqcount + 1 entries and the last one is filled by hand below, from the
       header-section length rather than from the file. Reading the whole
       vector would take one element too many and shift every later offset. */
    pos += largeread(in_stream, make_span(header_index).first(seqcount), pos, progress_bar);

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

    pos += largeread(in_stream, data_buffer.first(static_cast<std::size_t>(udb_headerchars)), pos, progress_bar);

    for (auto i = 0U; i < seqcount; i++)
      {
        longestheader = std::max<uint64_t>(seqindex[i].headerlen, longestheader);
      }

    /* sequence lengths */

    std::vector<unsigned int> sequence_lengths(seqcount);

    pos += largeread(in_stream, make_span(sequence_lengths), pos, progress_bar);

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

    pos += largeread(in_stream,
                     data_buffer.subspan(static_cast<std::size_t>(udb_headerchars),
                                         static_cast<std::size_t>(nucleotides)),
                     pos, progress_bar);

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

  if (usage == UdbUse::search)
    {
      auto const bitmap_mincount = bitmap_min_matches(seqcount);
      dbindex.set_bitmap_width(seqcount);
      {
        Progress progress("Creating bitmaps", dbindex.hashsize, parameters);
        for (auto i = 0U; i < dbindex.hashsize; i++)
          {
            if (dbindex.kmercount[i] >= bitmap_mincount)
              {
                auto & bitmap = dbindex.bitmap_create(i);
                for (auto j = 0U; j < dbindex.kmercount[i]; j++)
                  {
                    bitmap.set(dbindex.kmerindex[dbindex.kmerhash[i]+j]);
                  }
              }
            progress.update(i + 1);
          }
      }

      /* Distinct k-mers per index element, for the low_kmer_targets list
         below. Counted from the stored k-mer lists, which are complete for
         every k-mer -- the bitmaps above are filled from them -- so this is
         one pass over the loaded index rather than a second pass over the
         sequences, whose masking at index build time this session does not
         know. Every entry was checked against seqcount when the index was
         read, so the counter below is always in range. Index element numbers
         and sequence numbers coincide for a UDB database (map[i] == i, set
         further down). */
      std::vector<unsigned int> kmers_per_target(seqcount, 0U);
      for (auto kmer = 0U; kmer < dbindex.hashsize; ++kmer)
        {
          auto const first = dbindex.kmerhash[kmer];
          for (auto j = 0U; j < dbindex.kmercount[kmer]; ++j)
            {
              ++kmers_per_target[dbindex.kmerindex[first + j]];
            }
        }

      assert(parameters.opt_minwordmatches >= 0);
      dbindex.minwordmatches = static_cast<unsigned int>(parameters.opt_minwordmatches);
      for (auto element = 0U; element < seqcount; ++element)
        {
          auto const kmers = kmers_per_target[element];
          if ((kmers != 0) and (kmers < dbindex.minwordmatches))
            {
              dbindex.low_kmer_targets.push_back(LowKmerTarget{element, kmers});
            }
        }
    }

  /* get abundances and longest header */

  if (usage == UdbUse::search)
    {
      {
        Progress progress("Parsing abundances", seqcount, parameters);
        for (auto i = 0U; i < seqcount; i++)
          {
            /* udb_finalize() above only moved the sequences and rewrote their
               seq_p, so the headers are where db's own accessor finds them */
            auto const size = header_get_size(db.header_view(i));
            db.set_abundance(i, (size > 0) ? size : 1);
            progress.update(i + 1);
          }
      }
    }

  /* the unique-kmer finder (dbindex.uhandle) is a Uniquer value member, ready to
     use as default-constructed; the UDB path does not build the index with it */

  /* make mapping from indexno to seqno */

  dbindex.map.resize(seqcount);
  dbindex.count = seqcount;

  for (auto i = 0U; i < seqcount; i++)
    {
      dbindex.map[i] = i;
    }

  /* done */

  /* some stats */

  /* print_database_size() (core/db.hpp) branches on getsequencecount(), where
     this used to branch on the local seqcount. The two are the same number:
     udb_finalize() above was handed seqcount, and it is what the database
     reports. */
  if (not parameters.opt_quiet)
    {
      print_database_size(stderr, db);
    }

  if (parameters.fp_log != nullptr)
    {
      print_database_size(parameters.fp_log, db);
      fprint(parameters.fp_log, '\n');
    }
}
