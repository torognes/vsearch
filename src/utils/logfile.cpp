/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "logfile.hpp"
#include "vsearch.h"  // struct Parameters, fp_log
#include "util.h"  // open_optional_output, fclose_output
#include "arch.h"  // arch_get_memused
#include "timestamp.hpp"  // iso8601_local_timestamp
#include <chrono>  // std::chrono::steady_clock, std::chrono::duration
#include <cmath>  // std::floor
#include <cstdio>  // std::FILE, std::fprintf


LogFile::LogFile(struct Parameters & parameters)
{
  if (parameters.opt_log == nullptr) { return; }
  handle = open_optional_output(parameters.opt_log, "log");
  fp_log = handle;
  parameters.fp_log = handle;
  std::fprintf(handle, "%s\n", parameters.prog_header.c_str());
  std::fprintf(handle, "%s\n", parameters.command_line.c_str());

  start_time = std::chrono::steady_clock::now();
  std::fprintf(handle, "Started  %s\n", iso8601_local_timestamp().c_str());
}


LogFile::~LogFile()
{
  if (handle == nullptr) { return; }
  auto const finish_time = std::chrono::steady_clock::now();
  std::fprintf(handle, "\n");
  std::fprintf(handle, "Finished %s", iso8601_local_timestamp().c_str());

  double const time_diff =
    std::chrono::duration<double>(finish_time - start_time).count();
  std::fprintf(handle, "\n");
  std::fprintf(handle, "Elapsed time %02.0lf:%02.0lf\n",
          std::floor(time_diff / 60.0),
          std::floor(time_diff - (60.0 * std::floor(time_diff / 60.0))));
  double const maxmem = static_cast<double>(arch_get_memused()) / 1048576.0;
  if (maxmem < 1024.0)
    {
      std::fprintf(handle, "Max memory %.1lfMB\n", maxmem);
    }
  else
    {
      std::fprintf(handle, "Max memory %.1lfGB\n", maxmem / 1024.0);
    }
  fclose_output(handle);
}
