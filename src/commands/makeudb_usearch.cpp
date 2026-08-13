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
#include "commands/makeudb_usearch.hpp"
#include "core/db.hpp"
#include "core/mask.hpp"
#include "core/dbindex.hpp"
#include "utils/fatal.hpp"
#include "utils/progress.hpp"
#include "utils/view.hpp"  // View, make_view
#include <algorithm>  // std::max, std::fill_n
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <fstream>  // std::ofstream
#include <iterator>  // std::next
#include <ios>
#include <ostream>  // std::ostream
#include <vector>


constexpr auto blocksize = uint64_t{4096UL * 4096UL};


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* The UDB format's integer fields are 4 bytes by definition of the format,
     not because sizeof(unsigned int) happens to be 4 on this host. The byte
     count now comes from the view below, which makes that dependency
     implicit, so state it -- the cross-compilation targets all have a 32-bit
     int and would not catch it. Same assertion as in core/udb.cpp; there is
     no shared UDB header to put it in. */
  static_assert(sizeof(unsigned int) == 4, "UDB stores 32-bit fields");

  /* Write buf.size() elements to the UDB file. A view rather than a
     (void const *, nbyte) pair, for the reason given at largeread() in
     core/udb.cpp: one scratch vector is reused here for six different fields
     of six different lengths, and the length used to travel beside the
     pointer instead of with it. */
  template <typename Type>
  auto largewrite(std::ostream & output, View<Type> const buf, uint64_t const offset,
                  Progress & progress_bar) -> uint64_t
  {
    /* call write multiple times and update progress */

    auto const nbyte = static_cast<uint64_t>(buf.size_bytes());
    auto const * const bytes = reinterpret_cast<char const *>(buf.data());
    auto progress = offset;
    for (uint64_t i = 0; i < nbyte; i += blocksize)
      {
        auto const rem = std::min(blocksize, nbyte - i);
        output.write(std::next(bytes, static_cast<std::ptrdiff_t>(i)),
                     static_cast<std::streamsize>(rem));
        if (not output)
          {
            fatal("Unable to write to UDB file");
          }

        progress += rem;
        progress_bar.update(progress);
      }
    return nbyte;
  }

}  // end of anonymous namespace


auto makeudb_usearch(struct Parameters const & parameters) -> void
{
  Database db;  /* the sequence database this run owns (RAII) */
  Dbindex dbindex;  /* the k-mer index this run owns (RAII) */

  if (parameters.opt_output == nullptr) {
    fatal("UDB output file must be specified with --output");
  }

  std::ofstream out_stream(parameters.opt_output, std::ios::binary | std::ios::trunc);
  if (not out_stream)
    {
      fatal(std::string("Unable to open output file for writing (")
            + std::string(parameters.opt_output)
            + ")");
    }

  db.read(parameters.opt_makeudb_usearch, 1, parameters);

  /* Every reader rejects a UDB whose header declares no sequences (the
     buffer[13] == 0 term of udb_read()'s "Invalid UDB file" test), so writing
     one produces a file vsearch cannot open again -- and the abort then blames
     the format of a file vsearch itself wrote, at whatever command later tries
     to use it rather than here. Reachable from an empty input, and more easily
     from an input whose every sequence was discarded by --minseqlength, which
     defaults to 32 for this command: a file of shorter sequences used to yield
     an unreadable database and exit 0, announced only by the informational
     "N sequences discarded" line. usearch aborts in the same situation, with
     "Empty database" (verified on 9.2.64, 10.0.240 and 11.0.667).

     Placed before the index is built, so an empty run also stops short of
     writing the 4^wordlength-entry k-mer table. --output is opened above, so
     an aborted run leaves it created but empty, as usearch does too. */
  if (db.getsequencecount() == 0)
    {
      fatal("Cannot write a UDB file for an empty database");
    }

  if (parameters.opt_dbmask == Masking::dust)
    {
      dust_all(db, parameters);
    }
  else if ((parameters.opt_dbmask == Masking::soft) and parameters.opt_hardmask)
    {
      hardmask_all(db);
    }

  dbindex.prepare(1, parameters.opt_dbmask, db, parameters);
  dbindex.add_all_sequences(parameters.opt_dbmask, db, parameters);

  auto const seqcount = static_cast<unsigned int>(db.getsequencecount());
  auto const ntcount = db.getnucleotidecount();

  uint64_t header_characters = 0;
  for (auto i = 0U; i < seqcount; i++)
    {
      header_characters += db.getheaderlen(i) + 1;
    }

  uint64_t const kmerhash_entries = uint64_t{1} << (2 * static_cast<uint64_t>(parameters.opt_wordlength));

  /* count word matches */
  uint64_t wordmatches = 0;
  for (auto i = 0U; i < kmerhash_entries; i++)
    {
      wordmatches += dbindex.kmercount[i];
    }

  uint64_t pos = 0;
  uint64_t const progress_all =
    (uint64_t{4} * 50) +
    (4 * kmerhash_entries) +
    (uint64_t{4} * 1) +
    (4 * wordmatches) +
    (uint64_t{4} * 8) +
    (uint64_t{4} * seqcount) +
    header_characters +
    (uint64_t{4} * seqcount) +
    ntcount;


  /* the UDB file header is 50 32-bit words */
  static constexpr auto header_words = std::size_t{50};
  uint64_t const buffersize = std::max(50U, seqcount);
  std::vector<unsigned int> buffer(buffersize);

  /* Header */
  buffer[0]  = 0x55444246; /* FBDU UDBF */
  buffer[2]  = 32; /* bits */
  buffer[4]  = static_cast<unsigned int>(parameters.opt_wordlength); /* default 8 */
  buffer[5]  = 1; /* dbstep */
  buffer[6]  = 100; /* dbaccelpct % */
  buffer[11] = 0; /* slots */
  buffer[13] = seqcount; /* number of sequences */
  buffer[17] = 0x0000746e; /* alphabet: "nt" */
  buffer[49] = 0x55444266; /* fBDU UDBf */
  {
    Progress progress_bar("Writing UDB file", progress_all, parameters);
    /* the 50-word file header, not the whole scratch vector: buffer is
       sized max(50, seqcount) and is reused below for five more fields */
    pos += largewrite(out_stream, make_view(buffer).first(header_words), 0, progress_bar);

    /* write 4^wordlength uint32_t's with word match counts */
    pos += largewrite(out_stream, make_view(dbindex.kmercount).first(kmerhash_entries), pos, progress_bar);

    /* 3BDU */
    buffer[0] = 0x55444233; /* 3BDU UDB3 */
    pos += largewrite(out_stream, make_view(buffer).first(1), pos, progress_bar);

    /* lists of sequence no's with matches for all words */
    for (auto i = 0U; i < kmerhash_entries; i++)
      {
        if (dbindex.has_bitmap(i))
          {
            auto const & bitmap = dbindex.bitmap_of(i);
            std::fill_n(buffer.begin(), dbindex.kmercount[i], 0U);
            auto elements = 0U;
            for (auto j = 0U; j < seqcount; j++)
              {
                if (bitmap.is_set(j))
                  {
                    buffer[elements++] = j;
                  }
              }
            pos += largewrite(out_stream, make_view(buffer).first(elements), pos, progress_bar);
          }
        else
          {
            if (dbindex.kmercount[i] > 0)
              {
                /* a window into kmerindex, not a pointer plus a count:
                   subspan checks the offset, which is a stored index value
                   rather than a constant */
                pos += largewrite(out_stream,
                                  make_view(dbindex.kmerindex)
                                    .subspan(dbindex.kmerhash[i], dbindex.kmercount[i]),
                                  pos,
                                  progress_bar);
              }
          }
      }

    /* New header */
    buffer[0] = 0x55444234; /* 4BDU UDB4 */
    /* 0x005e0db3 */
    buffer[1] = 0x005e0db3;
    /* number of sequences, uint32_t */
    buffer[2] = seqcount;
    /* total number of nucleotides, uint64_t */
    buffer[3] = static_cast<unsigned int>(ntcount & 0xffffffff);
    buffer[4] = static_cast<unsigned int>(ntcount >> 32U);
    /* total number of header characters, incl zero-terminator, uint64_t */
    buffer[5] = static_cast<unsigned int>(header_characters & 0xffffffff);
    buffer[6] = static_cast<unsigned int>(header_characters >> 32U);
    /* 0x005e0db4 */
    buffer[7] = 0x005e0db4;
    pos += largewrite(out_stream, make_view(buffer).first(8), pos, progress_bar);

    /* indices to headers (uint32_t) */
    auto sum = 0U;
    for (auto i = 0U; i < seqcount; i++)
      {
        buffer[i] = sum;
        sum += static_cast<unsigned int>(db.getheaderlen(i) + 1);
      }
    pos += largewrite(out_stream, make_view(buffer).first(seqcount), pos, progress_bar);

    /* headers (ascii, zero terminated, not padded) */
    /* the terminator is written here rather than read out of Database's
       buffer. That byte is real -- Database::add appends it -- but it is a
       property of the UDB format, not of how Database happens to store a
       header, and reading a len + 1'th byte through db.getheader(i) coupled
       this writer to that layout. It goes through largewrite so it keeps both
       things a bare out_stream.put() would drop: the stream check, and the
       pos accounting every later offset depends on. */
    static constexpr char terminator = '\0';
    for (auto i = 0U; i < seqcount; i++)
      {
        pos += largewrite(out_stream, db.header_view(i), pos, progress_bar);
        pos += largewrite(out_stream, View<char>{& terminator, 1}, pos, progress_bar);
      }

    /* sequence lengths (uint32_t) */
    for (auto i = 0U; i < seqcount; i++)
      {
        buffer[i] = static_cast<unsigned int>(db.getsequencelen(i));
      }
    pos += largewrite(out_stream, make_view(buffer).first(seqcount), pos, progress_bar);

    /* sequences (ascii, no term, no pad) */
    for (auto i = 0U; i < seqcount; i++)
      {
        pos += largewrite(out_stream, db.sequence_view(i), pos, progress_bar);
      }

    out_stream.close();
    if (not out_stream)
      {
        fatal("Unable to close UDB file");
      }
  }

  dbindex.clear();
  db.clear();
}
