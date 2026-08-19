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
#include "os/system.hpp"
#include <algorithm>  // std::min, std::max

/* sched_getaffinity() and CPU_COUNT are Linux interfaces, and this file also
   serves as the generic-Unix backend (the build's default when the host is
   neither Windows, macOS nor FreeBSD), so the affinity mask is read only where
   it exists. Everything else here is portable: the cgroup files are simply
   absent on a host without them, which cgroup.hpp reports as "no limit". */
#ifdef __linux__
#include <sched.h>  // sched_getaffinity, cpu_set_t, CPU_ZERO, CPU_COUNT
#endif


auto system_get_available_cores() -> long
{
  auto const online = system_get_cores();

  /* sysconf() could not say, and there is then nothing to reduce. Passed
     through rather than corrected, so that this reports exactly what
     system_get_cores() has always reported on such a host. */
  if (online < 1) { return online; }

  auto available = online;

#if defined(__linux__) && defined(CPU_COUNT)
  /* The CPUs this process may actually be scheduled on. This covers taskset,
     a cpuset cgroup, Slurm's --cpu-bind and "docker run --cpuset-cpus" in one
     reading, and it is authoritative: the kernel will not schedule us
     elsewhere however many threads we start.

     A failure is not an error. In particular sched_getaffinity() reports
     EINVAL on a machine with more CPUs than a cpu_set_t can hold (1024), and
     such a host is not one whose core count we need to reduce. */
  cpu_set_t mask;
  CPU_ZERO(&mask);
  if (sched_getaffinity(0, sizeof(mask), &mask) == 0)
    {
      auto const allowed = CPU_COUNT(&mask);
      if (allowed > 0)
        {
          available = std::min(available, static_cast<long>(allowed));
        }
    }
#endif

  /* The tightest CPU quota on this process's cgroup path, in cores: what
     "docker run --cpus", a Kubernetes CPU limit and a Slurm cgroup with a CPU
     quota set. Confinement to a set of CPUs is not here, it is in the affinity
     mask read above. See os/linux/cgroup.hpp. */
  auto const quota = vsearch::cgroup::smallest_core_quota(
      vsearch::cgroup::own_hierarchy(vsearch::cgroup::default_locations(), "cpu"));
  if (quota != vsearch::cgroup::no_limit)
    {
      /* The quota is a core count derived from microseconds and cannot exceed
         what a long holds on any host that could also report 'online'. */
      available = std::min(available, static_cast<long>(quota));
    }

  static constexpr long at_least_one_core {1};
  return std::max(available, at_least_one_core);
}
