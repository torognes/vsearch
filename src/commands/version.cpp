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
#include <cstdio>  // std::printf, std::fprintf


auto show_publication() -> void
{
  std::fprintf(stdout,
          "Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)\n"
          "VSEARCH: a versatile open source tool for metagenomics\n"
          "PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584\n"
          "\n");
}


auto version(struct Parameters const & parameters) -> void
{
  if (parameters.opt_quiet) { return ; }

  show_publication();

#ifdef HAVE_ZLIB_H
  std::printf("Compiled with support for gzip-compressed files,");
  if ((parameters.dyn_libs != nullptr) and parameters.dyn_libs->gzip_available())
    {
      std::printf(" and the library is loaded.\n");

      char const * const gz_version = parameters.dyn_libs->gzip_version();
      unsigned long const flags = parameters.dyn_libs->gzip_compile_flags();

      std::printf("zlib version %s, compile flags %lx", gz_version, flags);
      static constexpr auto check_10th_bit = 1024U; // 0x0400
      if ((flags & check_10th_bit) != 0U)
        {
          std::printf(" (ZLIB_WINAPI)");
        }
      std::printf("\n");
    }
  else
    {
      std::printf(" but the library was not found.\n");
    }
#else
  std::printf("Compiled without support for gzip-compressed files.\n");
#endif

#ifdef HAVE_BZLIB_H
  std::printf("Compiled with support for bzip2-compressed files,");
  if ((parameters.dyn_libs != nullptr) and parameters.dyn_libs->bzip2_available())
    {
      std::printf(" and the library is loaded.\n");
    }
  else
    {
      std::printf(" but the library was not found.\n");
    }
#else
  std::printf("Compiled without support for bzip2-compressed files.\n");
#endif
}
