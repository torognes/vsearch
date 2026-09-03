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
#include "commands/udbinfo.hpp"
#include "os/system.hpp"
#include "utils/fatal.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/span.hpp"  // Span, make_span
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::fprintf
#include <fstream>  // std::ifstream
#include <sys/stat.h>
#include <ios>


namespace {

  /* Twenty-four write calls, and they used to be spelled out twice. The
     header array is passed whole rather than the seven fields reported. */
  auto print_udb_header(std::FILE * output_stream,
                        std::array<unsigned int, 50> const & buffer) -> void
  {
    fprint(output_stream, "           Seqs  ");
    fprint_integer(output_stream, buffer[13]);
    fprint(output_stream, '\n');
    fprint(output_stream, "     SeqIx bits  ");
    fprint_integer(output_stream, buffer[2]);
    fprint(output_stream, '\n');
    fprint(output_stream, "          Alpha  nt (4)\n");
    fprint(output_stream, "     Word width  ");
    fprint_integer(output_stream, buffer[4]);
    fprint(output_stream, '\n');
    fprint(output_stream, "          Slots  ");
    fprint_integer(output_stream, buffer[11]);
    fprint(output_stream, '\n');
    fprint(output_stream, "      Dict size  ");
    fprint_integer(output_stream, (1U << (2 * buffer[4])));
    fprint(output_stream, " (");
    std::fprintf(output_stream, "%.1f", (1U << (2 * buffer[4])) * 1.0 / 1000.0);
    fprint(output_stream, "k)\n");
    fprint(output_stream, "         DBstep  ");
    fprint_integer(output_stream, buffer[5]);
    fprint(output_stream, '\n');
    fprint(output_stream, "        DBAccel  ");
    fprint_integer(output_stream, buffer[6]);
    fprint(output_stream, "%\n");
  }

}  // anonymous namespace


auto udbinfo(struct Parameters const & parameters) -> void
{
  /* Read UDB header and show basic info */

  char const * const filename = parameters.input_filename;

  /* Reject a pipe up front, matching udb_read(): --udbinfo needs the file
     size to validate the sequence count below and cannot obtain it from a
     stream, so a UDB must be a seekable file path rather than a pipe. */

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

  auto const filesize = static_cast<uint64_t>(fs.st_size);

  /* The UDB format's integer fields are 4 bytes by definition of the format,
     not because sizeof(unsigned int) happens to be 4 on this host. The read
     length below now comes from the buffer rather than from a literal, which
     makes that dependency implicit, so state it -- the cross-compilation
     targets all have a 32-bit int and would not catch it. Same assertion as
     in core/udb.cpp and commands/makeudb_usearch.cpp; there is no shared UDB
     header to put it in. */
  static_assert(sizeof(unsigned int) == 4, "UDB stores 32-bit fields");

  std::array<unsigned int, 50> buffer {{}};

  std::ifstream in_stream(filename, std::ios::binary);
  if (not in_stream)
    {
      fatal("Unable to open UDB file for reading");
    }

  /* the destination as raw bytes: std::istream reads chars, whatever the
     span's element type is. The length travels with the destination instead
     of beside it as the literal 4 * 50 it used to be, so the two cannot
     disagree; see largeread() in core/udb.cpp for the same reasoning. */
  auto const header_bytes = make_span(buffer).as_writable_bytes();

  in_stream.read(header_bytes.data(),
                 static_cast<std::streamsize>(header_bytes.size()));
  if (static_cast<std::size_t>(in_stream.gcount()) != header_bytes.size())
    {
      fatal("Unable to read from UDB file or invalid UDB file");
    }

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

  /* Reject an inflated sequence count, mirroring the guard in udb_read():
     the per-sequence header-index and length tables each store 4 bytes
     per sequence, so a file cannot describe more than filesize/4 of them.
     buffer[13] is the one file-derived field --udbinfo reports; without
     this check a corrupt UDB that every other reader rejects would print
     a garbage count and still exit 0. */

  if (buffer[13] > filesize / 4)
    {
      fatal("Invalid UDB file");
    }

  if (not parameters.opt_quiet)
    {
      print_udb_header(stderr, buffer);
    }

  if (parameters.fp_log != nullptr)
    {
      print_udb_header(parameters.fp_log, buffer);
    }
}
