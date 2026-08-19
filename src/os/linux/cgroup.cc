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
#include "cgroup.hpp"
#include "utils/open_file.hpp"  // FileHandle
#include <algorithm>  // std::max
#include <array>
#include <cassert>  // assert
#include <cerrno>  // errno, ERANGE
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fopen, std::fgets
#include <cstdlib>  // std::strtoull
#include <numeric>  // std::accumulate
#include <string>
#include <vector>


namespace vsearch
{
  namespace cgroup
  {

    namespace
    {

      /* One line of 'input', however long, in 'line', which is cleared first;
         the newline is removed. False at end of input. */
      auto read_line(std::FILE * const input, std::string & line) -> bool
      {
        line.clear();
        static constexpr std::size_t chunk_size = 256;
        std::array<char, chunk_size> buffer {{}};
        while (std::fgets(buffer.data(), static_cast<int>(buffer.size()), input) != nullptr)
          {
            line.append(buffer.data());
            if ((not line.empty()) and (line.back() == '\n'))
              {
                line.pop_back();
                return true;
              }
          }
        return not line.empty();
      }


      /* True when 'controllers', the comma-separated list in the middle field
         of a cgroup v1 /proc line ("cpu,cpuacct"), contains 'controller' as a
         whole entry. A plain substring search would let "cpu" be answered by
         "cpuacct" or "cpuset", which name different hierarchies. */
      auto lists_controller(std::string const & controllers,
                            std::string const & controller) -> bool
      {
        std::string::size_type start = 0;
        while (start <= controllers.size())
          {
            auto const comma = controllers.find(',', start);
            auto const end = (comma == std::string::npos) ? controllers.size() : comma;
            if (controllers.compare(start, end - start, controller) == 0) { return true; }
            if (comma == std::string::npos) { break; }
            start = comma + 1;
          }
        return false;
      }


      /* This process's path within the hierarchy of 'controller', as the
         /proc file gives it, and whether that hierarchy is the unified one.

         The unified hierarchy (cgroup v2) writes a single line, "0::/path",
         for every controller at once. cgroup v1 writes one line per
         hierarchy, "id:controllers:/path", and the relevant one is the line
         whose controller list names ours. */
      auto own_path(Locations const & locations,
                    std::string const & controller,
                    bool & unified) -> std::string
      {
        unified = false;
        FileHandle const input {std::fopen(locations.proc_file.c_str(), "r")};
        if (not input) { return {}; }

        std::string line;
        std::string controller_path;
        while (read_line(input.get(), line))
          {
            static constexpr auto unified_prefix = "0::";
            static constexpr std::string::size_type unified_prefix_length = 3;
            if (line.compare(0, unified_prefix_length, unified_prefix) == 0)
              {
                unified = true;
                return line.substr(unified_prefix_length);
              }
            auto const first = line.find(':');
            if (first == std::string::npos) { continue; }
            auto const second = line.find(':', first + 1);
            if (second == std::string::npos) { continue; }
            if (lists_controller(line.substr(first + 1, second - first - 1), controller))
              {
                controller_path = line.substr(second + 1);
              }
          }
        return controller_path;
      }


      /* The CPU quota that one cgroup directory states, in cores, or no_limit
         when it states none. The two spellings and the rounding are described
         at smallest_core_quota() in the header. */
      auto directory_core_quota(std::string const & directory,
                                bool const unified) -> uint64_t
      {
        std::string line;
        auto quota = no_limit;
        auto period = no_limit;

        if (unified)
          {
            if (not read_first_line(directory + "/cpu.max", line)) { return no_limit; }
            auto const space = line.find(' ');
            if (space == std::string::npos) { return no_limit; }
            quota = parse_limit(line.substr(0, space));
            period = parse_limit(line.substr(space + 1));
          }
        else
          {
            if (read_first_line(directory + "/cpu.cfs_quota_us", line))
              {
                quota = parse_limit(line);  // "-1" when there is no quota
              }
            if (read_first_line(directory + "/cpu.cfs_period_us", line))
              {
                period = parse_limit(line);
              }
          }

        if (quota == no_limit) { return no_limit; }
        if (period == no_limit) { return no_limit; }

        /* parse_limit() reports a zero as no_limit, so the guard above is also
           what makes the division safe. */
        assert(period != 0);
        auto const cores = (quota + period - 1) / period;  // rounded up
        static constexpr uint64_t at_least_one_core {1};
        return std::max(cores, at_least_one_core);
      }

    }  // namespace


    auto default_locations() -> Locations
    {
      return Locations {"/sys/fs/cgroup", "/proc/self/cgroup"};
    }


    auto own_hierarchy(Locations const & locations,
                       std::string const & controller) -> Hierarchy
    {
      bool unified {false};
      auto path = own_path(locations, controller, unified);

      /* Under cgroup v1 each controller is mounted in its own subdirectory of
         the mount point, named after it; the unified hierarchy is one tree
         holding them all. A v1 host that mounts a controller under its full
         co-mounted name only ("cpu,cpuacct") rather than providing the usual
         per-controller symlink beside it is not found, and then reports no
         limit -- the pre-cgroup behaviour, which is the failure this whole
         module degrades to everywhere. */
      auto const base = unified
        ? locations.mount_point
        : locations.mount_point + "/" + controller;

      Hierarchy hierarchy {std::vector<std::string>(), unified};

      /* An unreadable /proc file, or a cgroup v1 host not mounting this
         controller, leaves nothing to walk. A genuine path always begins
         with a '/', even for a process in the root cgroup. */
      if (path.empty()) { return hierarchy; }

      /* Leaf first, then each ancestor in turn, ending with the hierarchy
         root itself (what remains once the path is empty). */
      while (true)
        {
          hierarchy.directories.push_back(base + path);
          auto const slash = path.rfind('/');
          if (slash == std::string::npos) { break; }
          path.erase(slash);
        }
      return hierarchy;
    }


    auto read_first_line(std::string const & path, std::string & line) -> bool
    {
      line.clear();
      FileHandle const input {std::fopen(path.c_str(), "r")};
      if (not input) { return false; }
      return read_line(input.get(), line);
    }


    auto parse_limit(std::string const & text) -> uint64_t
    {
      /* std::strtoull accepts a leading sign and wraps the result, which would
         turn cgroup v1's "no quota" marker -1 into 2^64-1. Rejected here
         rather than left to the magnitude test below, so that the reason is
         visible where it applies. */
      if (text.empty()) { return no_limit; }
      if ((text.front() == '-') or (text.front() == '+')) { return no_limit; }

      static constexpr int base_ten = 10;
      char * end_ptr = nullptr;
      errno = 0;
      auto const value = std::strtoull(text.c_str(), &end_ptr, base_ten);
      if (errno == ERANGE) { return no_limit; }
      if (end_ptr == text.c_str()) { return no_limit; }  // "max", or not a number
      if (*end_ptr != '\0') { return no_limit; }  // trailing characters

      /* cgroup v1's stand-in for "unlimited"; see the header. */
      static constexpr uint64_t implausible_limit {uint64_t{1} << 62U};
      if (value >= implausible_limit) { return no_limit; }

      return value;
    }


    auto tighten(uint64_t const limit, uint64_t const candidate) -> uint64_t
    {
      if (candidate == no_limit) { return limit; }
      if (limit == no_limit) { return candidate; }
      return (candidate < limit) ? candidate : limit;
    }


    /* Both walks are the same fold: start from "no opinion" and tighten() in
       what each directory says, leaf first. */

    auto smallest_limit(Hierarchy const & hierarchy,
                        std::string const & filename) -> uint64_t
    {
      return std::accumulate(hierarchy.directories.cbegin(),
                             hierarchy.directories.cend(),
                             no_limit,
                             [&filename](uint64_t const limit,
                                         std::string const & directory) -> uint64_t {
                               std::string line;
                               if (not read_first_line(directory + "/" + filename, line))
                                 {
                                   return limit;
                                 }
                               return tighten(limit, parse_limit(line));
                             });
    }


    auto smallest_core_quota(Hierarchy const & hierarchy) -> uint64_t
    {
      return std::accumulate(hierarchy.directories.cbegin(),
                             hierarchy.directories.cend(),
                             no_limit,
                             [&hierarchy](uint64_t const quota,
                                          std::string const & directory) -> uint64_t {
                               return tighten(quota,
                                              directory_core_quota(directory,
                                                                   hierarchy.unified));
                             });
    }

  }  // namespace cgroup
}  // namespace vsearch
