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
#pragma once

#include <cstdint>  // uint64_t
#include <string>
#include <vector>


/* Reading the control-group (cgroup) limits that apply to this process.
   Linux-only, and used by both os/linux/system_memory.cc (the memory limit)
   and os/linux/system_cpu.cc (the CPU quota), which is why the hierarchy walk
   lives here instead of being written twice.

   The problem these limits solve for us: a container shares the host's
   kernel, so sysconf(_SC_PHYS_PAGES) and sysconf(_SC_NPROCESSORS_ONLN) report
   the *machine's* memory and cores and know nothing of the cgroup that
   actually caps the run. Inside a 4 GB container on a 512 GB host they report
   512 GB (issue #584: a shell session confined to 50 MB and one core was told
   it had 500+ GB and 128 cores, and then crashed). Slurm with cgroups,
   Docker, podman, Kubernetes and Apptainer are where large vsearch runs live.

   A virtual machine needs none of this and is already correct: it has its own
   kernel, and sysconf() returns what the hypervisor gave the guest. The
   problem is cgroup isolation, not virtualisation.

   Nothing here calls fatal(). An absent, unreadable or unparseable file is
   not an error, it is "no limit here": a kernel without cgroups, a restricted
   /proc, and a future format change all degrade to the pre-cgroup behaviour
   of reporting the machine's figures rather than failing the run. */
namespace vsearch
{
  namespace cgroup
  {

    /* "No limit is set here", the value parse_limit() returns for every
       reading that does not yield a usable number. Zero cannot be a genuine
       limit -- a cgroup allowed no memory and no CPU time at all could not
       run a process -- so it needs no separate flag. */
    constexpr uint64_t no_limit {0};

    /* Where the kernel publishes the two things this module reads. Bundled
       rather than passed as two adjacent string parameters that a caller
       could hand over in either order, and parameters at all rather than
       constants so that the parsing and the walk can be exercised against a
       synthetic hierarchy. That matters here: neither vsearch-tests nor an
       unprivileged CI job can create a real cgroup, and the case a real host
       almost always presents is the one that finds no limit at all. */
    struct Locations
    {
      std::string mount_point;  // where the cgroup filesystem is mounted
      std::string proc_file;    // this process's cgroup membership

      /* No default member initializers: under C++11 they would make this a
         non-aggregate and break the brace-initialization at its call sites.
         C++14 refactoring: add them, aggregates may have them. */
    };

    /* The real ones, "/sys/fs/cgroup" and "/proc/self/cgroup". */
    auto default_locations() -> Locations;

    /* Where this process sits in the cgroup hierarchy. */
    struct Hierarchy
    {
      /* Every directory whose limits apply to this process: the cgroup it
         sits in first, then each ancestor in turn, up to the hierarchy root.
         Empty when this process is in no cgroup, or when the /proc file
         cannot be read.

         The ancestors are the point, not a refinement: Kubernetes sets the
         limit on the *pod's* cgroup, an ancestor of the container's, whose
         own memory.max reads "max". Reading the leaf alone finds nothing and
         concludes there is no limit, which is the bug this module exists to
         fix. The tightest limit anywhere on the path is the one that binds,
         so a caller reads the same file in each directory and keeps the
         smallest answer, folding the readings together with tighten(). */
      std::vector<std::string> directories;

      /* True under the unified hierarchy (cgroup v2), false under cgroup v1.
         The two spell every limit differently -- memory.max against
         memory.limit_in_bytes, cpu.max against cpu.cfs_quota_us -- so the
         caller needs this to know which file to ask for. */
      bool unified;

      /* No default member initializers, as in Locations above. */
    };

    /* The hierarchy of 'controller' ("memory", "cpu") that this process
       belongs to. The controller is what cgroup v1 needs and cgroup v2
       ignores: v1 mounts every controller separately, in its own
       subdirectory of the mount point and with its own line in the /proc
       file, so a process can sit at different paths in each; v2 has a single
       tree that carries them all. */
    auto own_hierarchy(Locations const & locations,
                       std::string const & controller) -> Hierarchy;

    /* The first line of the file at 'path', without its newline, in 'line',
       which is cleared first. False when the file cannot be opened or is
       empty, in which case 'line' is left empty.

       The line is accumulated until the newline rather than read into a
       fixed-size buffer: cgroup paths under Kubernetes run to a couple of
       hundred characters, and a truncated path would silently name a
       directory that does not exist. */
    auto read_first_line(std::string const & path, std::string & line) -> bool;

    /* The unsigned decimal integer that is the whole of 'text', or no_limit
       when 'text' is not one: cgroup v2's literal "max", cgroup v1's "-1",
       an empty field, trailing characters, or a value out of range.

       A value at or above 2^62 is also no_limit. cgroup v1 has no word for
       "unlimited" and writes 9223372036854771712 instead (INT64_MAX rounded
       down to a page multiple), which read literally is an 8 EiB budget --
       worse than no cap at all, since it would be believed. No real limit
       comes near that, in bytes or in microseconds. */
    auto parse_limit(std::string const & text) -> uint64_t;

    /* The smaller of 'limit' and 'candidate', with no_limit meaning "no
       opinion" rather than zero, so that a walk over a Hierarchy can fold
       every directory's reading into one running answer -- most directories
       on a real host have no limit set. */
    auto tighten(uint64_t limit, uint64_t candidate) -> uint64_t;

    /* The tightest reading of 'filename' -- a file holding one number, such
       as "memory.max" -- over every directory of 'hierarchy', or no_limit
       when none of them sets it. This is the whole of what a caller wanting a
       simple limit has to do:

         auto const own = own_hierarchy(default_locations(), "memory");
         auto const limit = smallest_limit(own, own.unified ? "memory.max"
                                                           : "memory.limit_in_bytes");

       The caller picks the filename because only it knows which limit it
       wants, and cgroup v1 and v2 never spell one the same way. */
    auto smallest_limit(Hierarchy const & hierarchy,
                        std::string const & filename) -> uint64_t;

    /* The tightest CPU quota over every directory of 'hierarchy', as a number
       of cores, or no_limit when none of them sets one. This is what
       "docker run --cpus", a Kubernetes CPU limit and a Slurm cgroup with a
       CPU quota write; confinement to a *set* of CPUs (taskset, cpuset,
       Slurm's --cpu-bind) is not here, it is in the affinity mask, which
       os/linux/system_cpu.cc reads instead.

       A quota is a ratio and not a single number, which is why it needs its
       own function rather than a smallest_limit() call: cgroup v2 states both
       halves in cpu.max as "<quota> <period>", cgroup v1 splits them over
       cpu.cfs_quota_us and cpu.cfs_period_us, and both are microseconds of
       CPU time per period. A quota of 150000 against a period of 100000 is
       1.5 cores, which is what "docker run --cpus=1.5" writes.

       The ratio is rounded *up*, and never below one core: 1.5 cores of
       throughput can be spent by two threads and only 1.0 of it by one, so
       rounding down would leave a third of the allocation unused. Cores are
       also what is compared across the hierarchy, rather than quotas, since
       two cgroups on one path may state their quotas over different
       periods. */
    auto smallest_core_quota(Hierarchy const & hierarchy) -> uint64_t;

  }  // namespace cgroup
}  // namespace vsearch
