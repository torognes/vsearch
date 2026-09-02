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

#include "fatal.hpp"
#include "open_file.hpp"
#include <unistd.h>  // dup, STDIN_FILENO, STDOUT_FILENO
#include <cassert>
#include <cerrno>  // errno
#include <cstddef>  // std::size_t
#include <cstdio>  // std::fopen, fdopen
#include <cstring>  // std::strcmp
#include <map>  // std::map (the registry of open output streams)
#include <memory>  // std::weak_ptr
#include <string>  // std::string (registry keys, open-failure message)
#include <utility>  // std::move


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // C++17 refactoring:
  // constexpr std::string_view a_dash = "-";
  // simpler string comparisons: if (filename == a_dash) {

  // type-safe string mode wrapper
  struct ModeString {
    explicit constexpr ModeString(char const * str) noexcept
      : mode(str) {
    }
    char const * mode;
  };


  auto is_dash(char const * filename) -> bool {
    assert(filename != nullptr);
    return std::strcmp(filename, "-") == 0;
  }


  auto check_file_descriptor(int const file_descriptor) -> void {
    assert(file_descriptor >= -1);
    if (file_descriptor != -1) {
      return;
    }
    if (errno == EBADF) {
      fatal("original fd is not an open file descriptor.");
    }
    if (errno == EMFILE) {
      fatal("too many open file descriptors.");
    }
    fatal("cannot duplicate input or output stream.");
  }


  // Open a stdio input stream and return a raw FILE * for the caller to wrap in
  // the matching RAII handle. A filename of '-' is served by a duplicate of
  // stdin. Outputs no longer come through here: they go through the registry
  // below, which has to look a name up before opening it, so it cannot hand the
  // job to a function that opens first.
  auto open_input_stream(char const * filename,
                         ModeString const & mode) -> std::FILE * {
    assert(filename != nullptr);
    assert(mode.mode != nullptr);
    if (is_dash(filename)) {
      auto const file_descriptor = dup(STDIN_FILENO);
      check_file_descriptor(file_descriptor);
      return fdopen(file_descriptor, mode.mode);
    }
    return std::fopen(filename, mode.mode);
  }


  // Both output openers report an open failure the same way; only the
  // null-filename policy differs between them.
  auto fatal_output_open_failed(char const * filename,
                                OutputOption const option) -> void {
    assert(filename != nullptr);
    std::string const message = std::string("unable to open output file for ")
      + option.name + " (" + filename + ")";
    fatal(message);
  }

  /* The output streams currently open, keyed by the target the user named, so
     that two options naming the same target write through one std::FILE -- one
     buffer, one file offset, one close -- instead of clobbering each other (two
     identical paths) or splicing their buffer flushes together mid-record (two
     streams on stdout).

     The values are weak, so a stream's lifetime is exactly the union of the
     handles that alias it: nothing outlives the command that opened it, and no
     static destructor can run the checked close -- and hence fatal() -- after
     main() has returned. Same "encapsulated owner" shape as log_handle in
     utils/logfile.cpp, but behind an accessor so that the table is built on
     first use rather than during static initialisation, and so that nothing
     else in this file can name it directly.

     Not synchronised, and it does not need to be: every output is opened from a
     command's prologue, on the thread that runs the command, before that
     command starts any workers. A worker thread that opened an output would
     need a lock here. */
  using OutputRegistry = std::map<std::string, std::weak_ptr<std::FILE>>;

  auto open_outputs() -> OutputRegistry & {
    static OutputRegistry registry;
    return registry;
  }


  /* The registry key for an output name.

     Only transforms that cannot make two *different* files look like the same
     one are allowed here, because the two kinds of mistake are not symmetric: a
     missed alias merely leaves the pre-existing clobbering behaviour in place,
     whereas a false alias would merge two streams into one file and never write
     the other at all. So: strip a leading "./", collapse runs of '/' that are
     not at the start of the name, and stop there.

     In particular ".." is never resolved -- removing it textually is wrong
     whenever the component before it is a symlink -- and a leading run of
     slashes is kept verbatim, because POSIX leaves a name beginning with
     exactly two slashes implementation-defined. Symlinks, hard links, a
     relative name spelled against an absolute one, and two mount points onto
     one file are therefore not deduplicated. This is best effort, not a
     guarantee.

     A real resolver (realpath() on the dirname plus the basename on POSIX,
     GetFullPathName on Windows) would close the relative/absolute gap, at the
     cost of a new src/os split; deferred until someone asks for it. */
  auto output_key(char const * filename) -> std::string {
    assert(filename != nullptr);
    std::string name{filename};
    while (name.compare(0, 2, "./") == 0) {
      name.erase(0, 2);
    }
    std::string key;
    key.reserve(name.size());
    auto leading_slashes = true;
    for (auto const character : name) {
      if (character != '/') {
        leading_slashes = false;
      }
      auto const duplicate_separator = (character == '/')
        and (not leading_slashes)
        and (not key.empty())
        and (key.back() == '/');
      if (duplicate_separator) {
        continue;
      }
      key.push_back(character);
    }
    return key;
  }


  /* Expired entries are swept lazily rather than erased by each handle's
     deleter. That keeps handle destruction free of any dependency on this
     table -- a deleter reaching back into global state would have to outlive
     it -- and keeps OutputFileHandle's deleter the plain checked close.

     A bound is needed because a command can open output files in a loop:
     --clusters opens one per cluster (see core/cluster.cpp), so a table that
     only ever grew would carry one expired entry per cluster for the rest of
     the run. Sweeping once the table passes a small bound is amortised O(1)
     per open -- each sweep costs O(bound) and reclaims everything but the
     handful of streams still open, and no command opens more than about
     fifteen outputs at a time. */
  constexpr auto max_stale_outputs = std::size_t{64};

  auto sweep_expired_outputs(OutputRegistry & registry) -> void {
    for (auto entry = registry.begin(); entry != registry.end(); ) {
      if (entry->second.expired()) {
        entry = registry.erase(entry);
      }
      else {
        ++entry;
      }
    }
  }


  /* Wrap a freshly opened output stream and record it under its key. A null
     stream (the open failed) is returned as an empty handle and not recorded,
     so a later attempt on the same name can retry rather than inherit the
     failure. The key is either absent from the table or present but expired --
     open_output_file() has already asked for a live one -- so assigning over
     whatever insert() finds is right either way. */
  auto register_output(std::string key, std::FILE * const stream) -> OutputFileHandle {
    if (stream == nullptr) {
      return OutputFileHandle{};
    }
    auto & registry = open_outputs();
    if (registry.size() > max_stale_outputs) {
      sweep_expired_outputs(registry);
    }
    OutputFileHandle handle{stream, CheckedCloseOutputHandle{}};
    auto const entry = registry.insert(
        OutputRegistry::value_type{std::move(key), std::weak_ptr<std::FILE>{}}).first;
    entry->second = handle;
    return handle;
  }


  /* The stream already open under this key, or an empty handle. */
  auto registered_output(std::string const & key) -> OutputFileHandle {
    auto const & registry = open_outputs();
    auto const entry = registry.find(key);
    if (entry == registry.end()) {
      return OutputFileHandle{};
    }
    return entry->second.lock();
  }


}  // end of anonymous namespace


// read_file, file to read, open_input_file, open_istream
auto open_input_file(char const * filename) -> FileHandle {
  if (filename == nullptr) {
    return FileHandle{nullptr};
  }
  auto const mode = ModeString{"rb"};  // r: reading, b: non-UNIX environments
  /* open the input stream given by filename, but if name is '-' then
     use a duplicate of stdin (fd = STDIN_FILENO = 0) */
  return FileHandle{open_input_stream(filename, mode)};
}


auto CheckedCloseOutputHandle::operator()(std::FILE * file_handle) const noexcept -> void {
  if (file_handle == nullptr) {
    return;
  }
  /* A write error (full disk, quota, broken pipe) is often deferred by stdio
     until the buffer is flushed, so check fflush and the error flag before
     closing; fclose also flushes and can report the same error. Fail loudly
     rather than leave a silently truncated output file. */
  if ((std::fflush(file_handle) != 0) or (std::ferror(file_handle) != 0)) {
    fatal("Unable to write to output file (disk full, quota exceeded, or broken pipe?)");
  }
  if (std::fclose(file_handle) != 0) {
    fatal("Unable to close output file (disk full or quota exceeded?)");
  }
}


// write_file, file to write, open_output_file, open_ostream
auto open_output_file(char const * filename) -> OutputFileHandle {
  if (filename == nullptr) {
    return OutputFileHandle{};
  }
  /* The registry lookup has to come before the open, not after: mode "wb"
     truncates, so opening a second stream on a file another handle is already
     writing would destroy what has been written so far. That ordering is also
     why the key can only be derived from the name -- identifying the file by
     (device, inode) would need it open, which is already too late. */
  auto const key = output_key(filename);
  if (auto const already_open = registered_output(key)) {
    return already_open;
  }
  auto const mode = ModeString{"wb"};  // w: writing, b: binary (no \n->\r\n on non-UNIX), matches input "rb"
  /* open the output stream given by filename, but if name is '-' then
     use a duplicate of stdout (fd = STDOUT_FILENO = 1) */
  if (is_dash(filename)) {
    auto const file_descriptor = dup(STDOUT_FILENO);
    check_file_descriptor(file_descriptor);
    return register_output(key, fdopen(file_descriptor, mode.mode));
  }
  return register_output(key, std::fopen(filename, mode.mode));
}


auto open_mandatory_output_file(char const * filename,
                                OutputOption const option) -> OutputFileHandle {
  if (filename == nullptr) {
    fatal(std::string("output file must be specified with ")
          + std::string(option.name));
  }
  auto output_handle = open_output_file(filename);
  if (not output_handle) {
    fatal_output_open_failed(filename, option);
  }
  return output_handle;
}


auto open_optional_output_file(char const * filename,
                               OutputOption const option) -> OutputFileHandle {
  if (filename == nullptr) {
    return OutputFileHandle{};
  }
  auto output_handle = open_output_file(filename);
  if (not output_handle) {
    fatal_output_open_failed(filename, option);
  }
  return output_handle;
}
