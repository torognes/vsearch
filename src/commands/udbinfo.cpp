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
#include "commands/udbinfo.hpp"
#include "os/system.hpp"
#include "utils/fatal.hpp"
#include <array>
#include <cstdint>  // uint64_t
#include <cstdio>  // std::fprintf
#include <fstream>  // std::ifstream


auto udbinfo(struct Parameters const & parameters) -> void
{
  /* Read UDB header and show basic info */

  char const * const filename = parameters.opt_udbinfo;

  /* Reject a pipe up front, matching udb_read(): --udbinfo needs the file
     size to validate the sequence count below and cannot obtain it from a
     stream, so a UDB must be a seekable file path rather than a pipe. */

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

  uint64_t const filesize = static_cast<uint64_t>(fs.st_size);

  std::array<unsigned int, 50> buffer {{}};

  std::ifstream in_stream(filename, std::ios::binary);
  if (not in_stream)
    {
      fatal("Unable to open UDB file for reading");
    }

  in_stream.read(static_cast<char *>(static_cast<void *>(buffer.data())), 4 * 50);
  if (static_cast<uint64_t>(in_stream.gcount()) != 4 * 50)
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
      std::fprintf(stderr, "           Seqs  %u\n", buffer[13]);
      std::fprintf(stderr, "     SeqIx bits  %u\n", buffer[2]);
      std::fprintf(stderr, "          Alpha  nt (4)\n");
      std::fprintf(stderr, "     Word width  %u\n", buffer[4]);
      std::fprintf(stderr, "          Slots  %u\n", buffer[11]);
      std::fprintf(stderr, "      Dict size  %u (%.1fk)\n",
              (1U << (2 * buffer[4])),
              (1U << (2 * buffer[4])) * 1.0 / 1000.0);
      std::fprintf(stderr, "         DBstep  %u\n", buffer[5]);
      std::fprintf(stderr, "        DBAccel  %u%%\n", buffer[6]);
    }

  if (parameters.opt_log != nullptr)
    {
      std::fprintf(parameters.fp_log, "           Seqs  %u\n", buffer[13]);
      std::fprintf(parameters.fp_log, "     SeqIx bits  %u\n", buffer[2]);
      std::fprintf(parameters.fp_log, "          Alpha  nt (4)\n");
      std::fprintf(parameters.fp_log, "     Word width  %u\n", buffer[4]);
      std::fprintf(parameters.fp_log, "          Slots  %u\n", buffer[11]);
      std::fprintf(parameters.fp_log, "      Dict size  %u (%.1fk)\n",
              (1U << (2 * buffer[4])),
              (1U << (2 * buffer[4])) * 1.0 / 1000.0);
      std::fprintf(parameters.fp_log, "         DBstep  %u\n", buffer[5]);
      std::fprintf(parameters.fp_log, "        DBAccel  %u%%\n", buffer[6]);
    }
}
