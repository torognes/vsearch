/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include <cstdint> // uint64_t
#include <cstdio>  // std::FILE, std::size_t


constexpr auto md5_digest_length = 16;
constexpr auto sha1_digest_length = 20;

constexpr auto len_hex_dig_md5 = (2 * md5_digest_length) + 1;
constexpr int len_hex_dig_sha1 = (2 * sha1_digest_length) + 1;

auto fatal(const char * msg) -> void;
auto fatal(const char * format, const char * message) -> void;
auto xstrdup(const char * src) -> char *;
auto xstrchrnul(char * str, int target) -> char *;
auto xsprintf(char * * ret, const char * format, ...) -> int;
auto hash_cityhash64(char * sequence, uint64_t length) -> uint64_t;
auto hash_cityhash128(char * sequence, uint64_t length) -> uint128;
auto show_rusage() -> void;

auto progress_init(const char * prompt, uint64_t size) -> void;
auto progress_update(uint64_t progress) -> void;
auto progress_done() -> void;

auto random_init() -> void;
auto random_int(int64_t upper_limit) -> int64_t;
auto random_ulong(uint64_t upper_limit) -> uint64_t;

auto string_normalize(char * normalized, char * raw_seq, unsigned int len) -> void;

auto reverse_complement(char * rc_seq, char * seq, int64_t len) -> void;

auto fprint_hex(std::FILE * output_handle, unsigned char * data, int len) -> void;

auto get_hex_seq_digest_sha1(char * hex, char * seq, int seqlen) -> void;
auto get_hex_seq_digest_md5(char * hex, char * seq, int seqlen) -> void;

auto fprint_seq_digest_sha1(std::FILE * output_handle, char * seq, int seqlen) -> void;
auto fprint_seq_digest_md5(std::FILE * output_handle, char * seq, int seqlen) -> void;

auto fopen_input(const char * filename) -> std::FILE *;
auto fopen_output(const char * filename) -> std::FILE *;

inline auto xpthread_attr_init(pthread_attr_t * attr) -> void
{
  if (pthread_attr_init(attr) != 0)
    {
      fatal("Unable to init thread attributes");
    }
}

inline auto xpthread_attr_destroy(pthread_attr_t * attr) -> void
{
  if (pthread_attr_destroy(attr) != 0)
    {
      fatal("Unable to destroy thread attributes");
    }
}

inline auto xpthread_attr_setdetachstate(pthread_attr_t * attr, int detachstate) -> void
{
  if (pthread_attr_setdetachstate(attr, detachstate) != 0)
    {
      fatal("Unable to set thread attributes detach state");
    }
}

inline auto xpthread_create(pthread_t * thread, const pthread_attr_t * attr,
                            void *(*start_routine)(void *), void * arg) -> void
{
  if (pthread_create(thread, attr, start_routine, arg) != 0)
    {
      fatal("Unable to create thread");
    }
}

inline auto xpthread_join(pthread_t thread, void ** value_ptr) -> void
{
  if (pthread_join(thread, value_ptr) != 0)
    {
      fatal("Unable to join thread");
    }
}

inline auto xpthread_mutex_init(pthread_mutex_t * mutex,
                                const pthread_mutexattr_t * attr) -> void
{
  if (pthread_mutex_init(mutex, attr) != 0)
    {
      fatal("Unable to init mutex");
    }
}

inline auto xpthread_mutex_destroy(pthread_mutex_t * mutex) -> void
{
  if (pthread_mutex_destroy(mutex) != 0)
    {
      fatal("Unable to destroy mutex");
    }
}

inline auto xpthread_mutex_lock(pthread_mutex_t * mutex) -> void
{
  if (pthread_mutex_lock(mutex) != 0)
    {
      fatal("Unable to lock mutex");
    }
}

inline auto xpthread_mutex_unlock(pthread_mutex_t * mutex) -> void
{
  if (pthread_mutex_unlock(mutex) != 0)
    {
      fatal("Unable to unlock mutex");
    }
}

inline auto xpthread_cond_init(pthread_cond_t * cond,
                               const pthread_condattr_t * attr) -> void
{
  if (pthread_cond_init(cond, attr) != 0)
    {
      fatal("Unable to init condition variable");
    }
}

inline auto xpthread_cond_destroy(pthread_cond_t * cond) -> void
{
  if (pthread_cond_destroy(cond) != 0)
    {
      fatal("Unable to destroy condition variable");
    }
}

inline auto xpthread_cond_wait(pthread_cond_t * cond, pthread_mutex_t * mutex) -> void
{
  if (pthread_cond_wait(cond, mutex) != 0)
    {
      fatal("Unable to wait on condition variable");
    }
}

inline auto xpthread_cond_signal(pthread_cond_t * cond) -> void
{
  if (pthread_cond_signal(cond) != 0)
    {
      fatal("Unable to signal condition variable");
    }
}

inline auto xpthread_cond_broadcast(pthread_cond_t * cond) -> void
{
  if (pthread_cond_broadcast(cond) != 0)
    {
      fatal("Unable to broadcast condition variable");
    }
}
