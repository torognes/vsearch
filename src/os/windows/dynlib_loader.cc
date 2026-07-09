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

#include "../../utils/dynlib_loader.hpp"
#include <windows.h>  // LoadLibraryExA, GetProcAddress, FreeLibrary, HMODULE


/*
  A bare DLL name passed to LoadLibrary is resolved through the Windows DLL
  search order, which can include the current working directory — so running
  vsearch from a directory an attacker can write to lets a planted
  zlib1.dll / libbz2.dll be loaded (DLL planting). LoadLibraryEx with
  LOAD_LIBRARY_SEARCH_DEFAULT_DIRS restricts the search to the application
  directory, any AddDllDirectory dirs, and System32 — never the current
  directory. The flag is defined on Windows 8+ (and Vista/7 with KB2533623);
  define it here in case the build SDK headers predate it.
*/
#ifndef LOAD_LIBRARY_SEARCH_DEFAULT_DIRS
#define LOAD_LIBRARY_SEARCH_DEFAULT_DIRS 0x00001000
#endif


namespace dynlib {

  auto open(const char * library_name) -> handle
  {
    return reinterpret_cast<handle>(
      LoadLibraryExA(library_name, nullptr, LOAD_LIBRARY_SEARCH_DEFAULT_DIRS));
  }

  auto symbol(handle library, const char * symbol_name) -> symbol_ptr
  {
    return reinterpret_cast<symbol_ptr>(
      GetProcAddress(reinterpret_cast<HMODULE>(library), symbol_name));
  }

  auto close(handle library) -> void
  {
    FreeLibrary(reinterpret_cast<HMODULE>(library));
  }

}  // namespace dynlib
