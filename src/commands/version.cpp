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
#include "commands/version.hpp"
#include "os/dynlibs.hpp"
#include "utils/print_view.hpp"  // fprint
#include <cstdio>  // std::printf


auto show_publication() -> void
{
  fprint(stdout, "Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)\n"
                 "VSEARCH: a versatile open source tool for metagenomics\n"
                 "PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584\n"
                 "\n");
}


auto version(struct Parameters const & parameters) -> void
{
  if (parameters.opt_quiet) { return ; }

  show_publication();

  if (compression::gzip_supported)
    {
      fprint(stdout, "Compiled with support for gzip-compressed files,");
      if ((parameters.dyn_libs != nullptr) and parameters.dyn_libs->gzip_available())
        {
          fprint(stdout, " and the library is loaded.\n");

          char const * const gz_version = parameters.dyn_libs->gzip_version();
          unsigned long const flags = parameters.dyn_libs->gzip_compile_flags();

          fprint(stdout, "zlib version ");
          std::fputs(gz_version, stdout);
          fprint(stdout, ", compile flags ");
          /* the only hex conversion left in the tree beside fasta.cpp's
             "%02x": zlib reports its compile flags as a bit field, and a
             decimal rendering of them would not be readable. Deliberately
             kept, see TBD_20260804_c_style_elimination.md's Out of scope. */
          std::printf("%lx", flags);
          static constexpr auto check_10th_bit = 1024U; // 0x0400
          if ((flags & check_10th_bit) != 0U)
            {
              fprint(stdout, " (ZLIB_WINAPI)");
            }
          fprint(stdout, '\n');
        }
      else
        {
          fprint(stdout, " but the library was not found.\n");
        }
    }
  else
    {
      fprint(stdout, "Compiled without support for gzip-compressed files.\n");
    }

  if (compression::bzip2_supported)
    {
      fprint(stdout, "Compiled with support for bzip2-compressed files,");
      if ((parameters.dyn_libs != nullptr) and parameters.dyn_libs->bzip2_available())
        {
          fprint(stdout, " and the library is loaded.\n");
        }
      else
        {
          fprint(stdout, " but the library was not found.\n");
        }
    }
  else
    {
      fprint(stdout, "Compiled without support for bzip2-compressed files.\n");
    }
}
