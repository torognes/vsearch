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

#include "open_file.hpp"  // OutputFileHandle
#include <chrono>  // std::chrono::steady_clock
#include <cstdio>  // std::FILE

struct Parameters;


/* Accessor for the optional --log file handle, used by the process-wide
   error/warning reporters (fatal(), fastx.cc's warn()) that run with no
   Parameters in scope. The handle is owned and published by the LogFile RAII
   object below; handle() returns nullptr when no --log file is open. Code that
   already holds a Parameters should read parameters.fp_log directly instead. */
namespace log_file
{
  auto handle() noexcept -> std::FILE *;
  auto set_handle(std::FILE * new_handle) noexcept -> void;
}


/* RAII owner of the optional --log file. When --log is given, the constructor
   opens the file, publishes it through parameters.fp_log (for code that holds a
   Parameters) and log_file::set_handle() (for the Parameters-less reporters)
   so the rest of the program logs to it, and writes the program header, the
   command line and the "Started" timestamp; the destructor writes the
   "Finished"/elapsed-time/max-memory footer and closes the file. Co-locating
   open and close keeps the two halves of the log lifecycle together and emits
   the footer on every exit path out of the enclosing scope. Without --log the
   object owns no file and both halves are no-ops. */
class LogFile
{
public:
  explicit LogFile(struct Parameters & parameters);
  ~LogFile();

  LogFile(LogFile const &) = delete;
  LogFile(LogFile &&) = delete;
  auto operator=(LogFile const &) -> LogFile & = delete;
  auto operator=(LogFile &&) -> LogFile & = delete;

private:
  OutputFileHandle handle;
  std::chrono::steady_clock::time_point start_time;
};
