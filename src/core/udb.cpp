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
#include <algorithm>  // std::any_of, std::min, std::max
#include <array>
#include <cassert>  // assert
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t
#include <fstream>  // std::ifstream
#include <ios>
#include <istream>  // std::istream
#include <limits>
#include <numeric>  // std::accumulate, std::iota
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


auto udb_detect_isudb(char const * filename) -> bool
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
     read() does not consume the data of a stream that stays open
     elsewhere (stdin, or a /dev/fd/N entry duplicating an inherited
     descriptor), so bailing out for anything that is not a regular file
     leaves such a stream intact for the reader. A named FIFO is the
     exception, and is never opened here (see below).
     open_input_file() also maps "-" to a duplicate of stdin (matching the
     reader), whereas stat() of the literal path "-" would fail. */

  /* A named FIFO must not be opened here at all: this probe is a
     separate open() from the reader's, and a FIFO's buffered data is
     discarded when its last reader closes. If the writer has already
     written and exited when the probe closes, the reader's own open()
     then waits forever for a writer that will never come (a real named
     pipe, or <() on FreeBSD, where bash implements it with one). Only
     stat() of the path can tell without opening; a path that is not a
     FIFO, or that stat() cannot resolve (such as "-"), takes the
     descriptor route below. */
  xstat_t path_status;
  if ((xstat(filename, & path_status) == 0) and S_ISFIFO(path_status.st_mode))
    {
      return false;
    }

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


  /* A run of 4-byte fields in the file: how many there are, and the byte
     offset the first one starts at. Named because both are uint64_t, and a
     call site that swapped them would read a plausible number of fields from
     the wrong place -- the same class of mistake largeread()'s Span parameter
     removed for lengths. */
  struct FieldRun
  {
    uint64_t entries;
    uint64_t file_offset;
  };


  /* Read a section of 4-byte values through a scratch buffer, handing each
     block to `inspect` and keeping none of it.

     Every value in the k-mer sections is read for one of two reasons: because
     the session will use it, or because reading it is what proves the file is
     intact. A reporting session has only the second reason, and holding a
     section it will not use costs 4 bytes per entry -- 1.22 GB of word list
     for a 221 085-sequence reference, 268 MB of word counts at word length 13.

     The buffer is one largeread() block at most, so the progress bar advances
     at the same byte offsets it would with the whole section in memory, and no
     larger than the section, so a small database does not trade a 10 MB table
     for a 16 MB scratch buffer. */
  template <typename Inspect>
  auto udb_stream_section(std::istream & input, FieldRun const run,
                          Progress & progress_bar,
                          Inspect inspect) -> uint64_t
  {
    auto const block_entries = static_cast<std::size_t>(
      std::min<uint64_t>(blocksize / sizeof(unsigned int), run.entries));
    std::vector<unsigned int> block(block_entries);
    auto position = run.file_offset;

    for (uint64_t done = 0; done < run.entries; done += block_entries)
      {
        auto const wanted = static_cast<std::size_t>(
          std::min<uint64_t>(block_entries, run.entries - done));
        auto const chunk = make_span(block).first(wanted);
        position += largeread(input, chunk, position, progress_bar);
        inspect(chunk);
      }

    return position - run.file_offset;
  }


}  // end of anonymous namespace


auto udb_read(char const * filename,
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
       configured value.

       Only a search session has a configured value to override: --wordlength
       is not among the options --udb2fasta, --udbstats or --udbinfo accept
       (they reject it), so opt_wordlength is there the hardcoded default of 8
       and the warning told the user that a setting they could not have made
       had been adjusted -- on every UDB not built at word length 8. */
    if ((usage == UdbUse::search)
        and (udb_wordlength != static_cast<unsigned int>(parameters.opt_wordlength)))
      {
        vsearch::warn("Wordlength adjusted to " + decimal::to_text(udb_wordlength)
                      + " as indicated in UDB file");
      }
    dbindex.wordlength = udb_wordlength;

    /* word match counts */

    dbindex.hashsize = 1U << (2 * udb_wordlength);
    /* The k-mer -> bitmap lookup is one 4-byte slot per k-mer, allocated and
       zeroed here but filled only by the bitmap loop further down, which runs
       for UdbUse::search alone. A reporting session never asks a k-mer whether
       it has a bitmap, so it can leave the table empty rather than pay 4 bytes
       per slot for it -- 268 MB at word length 13. has_bitmap() asserts the
       k-mer is in range, so a session that starts reading it after all says so
       in a debug build rather than reading out of bounds. */
    if (usage == UdbUse::search)
      {
        dbindex.bitmap_slots_reset(dbindex.hashsize);
      }
    /* filled by push_back into reserved space below, not resized here: the loop
       writes every entry, so value-initialising them first is 8 bytes per slot
       of zeros nobody reads. Only search walks the whole offset table; the
       report needs a handful of its entries and recovers those from the counts,
       so it reserves nothing (537 MB at word length 13). */
    dbindex.kmerhash.clear();
    if (usage == UdbUse::search)
      {
        dbindex.kmerhash.reserve(dbindex.hashsize);
      }

    /* word counts */

    dbindex.indexsize = 0;

    if (usage == UdbUse::sequences)
      {
        /* The counts say nothing this session will report, but their sum is the
           length of the word list, which is what says where the sequence half
           of the file begins. So they are read for the total and dropped. */
        dbindex.kmercount.clear();
        dbindex.kmercount.shrink_to_fit();
        pos += udb_stream_section(in_stream, FieldRun{dbindex.hashsize, pos}, progress_bar,
                                  [&dbindex, seqcount](Span<unsigned int> const block) -> void
                                  {
                                    for (auto const count : block)
                                      {
                                        /* A word count is a number of
                                           database sequences; see the loop
                                           below for why a larger value cannot
                                           describe a valid word list. */
                                        if (count > seqcount)
                                          {
                                            fatal("Invalid UDB file");
                                          }
                                        dbindex.indexsize =
                                          udb_checked_add(dbindex.indexsize, count);
                                      }
                                  });
      }
    else
      {
        dbindex.kmercount.resize(dbindex.hashsize);
        pos += largeread(in_stream, make_span(dbindex.kmercount).first(dbindex.hashsize), pos, progress_bar);

        for (uint64_t i = 0; i < dbindex.hashsize; i++)
          {
            /* A word count is a number of database sequences: both index
               builders add a given sequence at most once per k-mer (its
               distinct k-mers come from Uniquer), so a count above seqcount
               cannot describe a valid word list, whatever wrote the file. The
               word list's own entries are already checked against seqcount
               further down; this is the same property one level up, and it is
               what lets --udbstats treat the counts as a distribution over
               0..seqcount rather than over an unbounded range. */
            if (dbindex.kmercount[i] > seqcount)
              {
                fatal("Invalid UDB file");
              }
            if (usage == UdbUse::search)
              {
                dbindex.kmerhash.push_back(dbindex.indexsize);
              }
            dbindex.indexsize = udb_checked_add(dbindex.indexsize, dbindex.kmercount[i]);
          }
        /* one entry per slot: size() counts what the loop wrote, so this checks
           it covered every slot (this path stores no end marker, unlike
           prepare) */
        assert((usage != UdbUse::search)
               or (dbindex.kmerhash.size() == dbindex.hashsize));
      }

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

    /* Every entry is a sequence number used both as a bit offset in the
       per-word bitmaps (Bitmap::set writes bitmap[value >> 3], no bounds
       check) and as an index into seqindex/dbindex_map during search. A
       value >= seqcount is therefore an out-of-bounds write or read, so
       reject it here rather than at use.

       Search needs the entries afterwards and keeps them. A reporting session
       does not: it validates them as they stream past and holds none, which is
       4 bytes per entry it no longer allocates -- 1.22 GB for a
       221 085-sequence reference at word length 8, where this section is 78 %
       of the file. The check is kept in both cases because it is a property of
       the file, not of what this session means to do with it: a truncated or
       rewritten word list is rejected by --udb2fasta and --udbstats exactly
       where it is rejected today. */

    if (usage == UdbUse::search)
      {
        dbindex.kmerindex.resize(dbindex.indexsize);

        pos += largeread(in_stream, make_span(dbindex.kmerindex).first(dbindex.indexsize), pos, progress_bar);

        /* every word-list entry is a sequence number, so none may reach
           seqcount; the whole list was just read, so this is the whole
           container */
        if (std::any_of(dbindex.kmerindex.cbegin(), dbindex.kmerindex.cend(),
                        [seqcount](unsigned int const seqno) -> bool
                        { return seqno >= seqcount; }))
          {
            fatal("Invalid UDB file");
          }
      }
    else
      {
        dbindex.kmerindex.clear();
        dbindex.kmerindex.shrink_to_fit();
        pos += udb_stream_section(in_stream, FieldRun{dbindex.indexsize, pos}, progress_bar,
                                  [seqcount](Span<unsigned int> const block) -> void
                                  {
                                    for (auto const entry : block)
                                      {
                                        if (entry >= seqcount)
                                          {
                                            fatal("Invalid UDB file");
                                          }
                                      }
                                  });
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

    longestheader = std::accumulate(seqindex, std::next(seqindex, seqcount), longestheader,
                                    [](uint64_t const longest_so_far, seqinfo_t const & record) -> uint64_t
                                    { return std::max<uint64_t>(longest_so_far, record.headerlen); });

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

  /* Search reads the abundances through db.getabundance() and gates their use
     on --sizein at each site (see msa.cpp, cluster.cpp, derep.cpp), so its
     parse is unconditional as before.

     --udb2fasta (UdbUse::sequences) parses them only when asked. A UDB stores
     headers verbatim, so the annotation is there to be read, but reading it is
     what --sizein means throughout vsearch -- without it an input's abundances
     are not read and every sequence counts as one. So --sizeout alone still
     writes size=1, matching every other command, and --sizein --sizeout now
     writes the stored value instead of overwriting it with 1. */

  auto const parse_abundances =
    (usage == UdbUse::search) or parameters.opt_sizein;

  if (parse_abundances)
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

  std::iota(dbindex.map.begin(), dbindex.map.end(), 0U);

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


auto udb_read_word_entries(char const * filename,
                           struct Dbindex const & dbindex,
                           unsigned int const seqcount,
                           uint64_t const first,
                           Span<unsigned int> const entries) -> void
{
  assert(dbindex.wordlength >= 3);
  assert(dbindex.wordlength <= 15);
  assert(dbindex.hashsize == (1U << (2 * dbindex.wordlength)));

  if (entries.empty())
    {
      return;
    }

  std::ifstream in_stream(filename, std::ios::binary);
  if (not in_stream)
    {
      fatal("Unable to open UDB file for reading");
    }

  /* Where entry number `first` sits: the 50-field header, then one count per
     k-mer slot, then the section signature, then the word list itself -- all of
     them 4-byte fields (see the static_assert in udb_read). Spelled out rather
     than as a literal so it stays tied to the layout above it. */
  auto const header_fields = uint64_t{50};
  auto const signature_fields = uint64_t{1};
  auto const slots = uint64_t{dbindex.hashsize};
  auto const field = uint64_t{sizeof(unsigned int)};
  auto const offset = field * (header_fields + slots + signature_fields + first);

  in_stream.seekg(static_cast<std::streamoff>(offset));
  if (not in_stream)
    {
      fatal("Unable to read from UDB file or invalid UDB file");
    }

  auto const bytes = entries.as_writable_bytes();
  in_stream.read(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  if (static_cast<std::size_t>(in_stream.gcount()) != bytes.size())
    {
      fatal("Unable to read from UDB file or invalid UDB file");
    }

  for (auto const entry : entries)
    {
      if (entry >= seqcount)
        {
          fatal("Invalid UDB file");
        }
    }
}
