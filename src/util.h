/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef exp10
#define exp10(x) (pow(10.0,(x)))
#endif

#define SHA_DIGEST_LENGTH SHA1_DIGEST_SIZE

constexpr int MD5_DIGEST_LENGTH {16};
#define LEN_DIG_SHA1 SHA_DIGEST_LENGTH

constexpr int LEN_HEX_DIG_MD5 {2 * MD5_DIGEST_LENGTH + 1};
#define LEN_HEX_DIG_SHA1 (2*LEN_DIG_SHA1+1)

void fatal(const char * msg);
void fatal(const char * format, const char * message);
char * xstrdup(const char *s);
char * xstrchrnul(char *s, int c);
int xsprintf(char * * ret, const char * format, ...);
uint64_t hash_cityhash64(char * s, uint64_t n);
uint128 hash_cityhash128(char * s, uint64_t n);
int64_t getusec();
void show_rusage();

void progress_init(const char * prompt, uint64_t size);
void progress_update(uint64_t progress);
void progress_done();

void random_init();
int64_t random_int(int64_t n);
uint64_t random_ulong(uint64_t n);

void string_normalize(char * normalized, char * s, unsigned int len);

void reverse_complement(char * rc, char * seq, int64_t len);

void fprint_hex(FILE * fp, unsigned char * data, int len);

void get_hex_seq_digest_sha1(char * hex, char * seq, int seqlen);
void get_hex_seq_digest_md5(char * hex, char * seq, int seqlen);

void fprint_seq_digest_sha1(FILE * fp, char * seq, int seqlen);
void fprint_seq_digest_md5(FILE * fp, char * seq, int seqlen);

FILE * fopen_input(const char * filename);
FILE * fopen_output(const char * filename);

void inline xpthread_attr_init(pthread_attr_t *attr)
{
  if (pthread_attr_init(attr))
    {
      fatal("Unable to init thread attributes");
    }
}

void inline xpthread_attr_destroy(pthread_attr_t *attr)
{
  if (pthread_attr_destroy(attr))
    {
      fatal("Unable to destroy thread attributes");
    }
}

void inline xpthread_attr_setdetachstate(pthread_attr_t *attr, int detachstate)
{
  if (pthread_attr_setdetachstate(attr, detachstate))
    {
      fatal("Unable to set thread attributes detach state");
    }
}

void inline xpthread_create(pthread_t *thread, const pthread_attr_t *attr,
                            void *(*start_routine)(void *), void *arg)
{
  if (pthread_create(thread, attr, start_routine, arg))
    {
      fatal("Unable to create thread");
    }
}

void inline xpthread_join(pthread_t thread, void **value_ptr)
{
  if (pthread_join(thread, value_ptr))
    {
      fatal("Unable to join thread");
    }
}

void inline xpthread_mutex_init(pthread_mutex_t *mutex,
                                const pthread_mutexattr_t *attr)
{
  if (pthread_mutex_init(mutex, attr))
    {
      fatal("Unable to init mutex");
    }
}

void inline xpthread_mutex_destroy(pthread_mutex_t *mutex)
{
  if (pthread_mutex_destroy(mutex))
    {
      fatal("Unable to destroy mutex");
    }
}

void inline xpthread_mutex_lock(pthread_mutex_t *mutex)
{
  if (pthread_mutex_lock(mutex))
    {
      fatal("Unable to lock mutex");
    }
}

void inline xpthread_mutex_unlock(pthread_mutex_t *mutex)
{
  if (pthread_mutex_unlock(mutex))
    {
      fatal("Unable to unlock mutex");
    }
}

void inline xpthread_cond_init(pthread_cond_t *cond,
                               const pthread_condattr_t *attr)
{
  if (pthread_cond_init(cond, attr))
    {
      fatal("Unable to init condition variable");
    }
}

void inline xpthread_cond_destroy(pthread_cond_t *cond)
{
  if (pthread_cond_destroy(cond))
    {
      fatal("Unable to destroy condition variable");
    }
}

void inline xpthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex)
{
  if (pthread_cond_wait(cond, mutex))
    {
      fatal("Unable to wait on condition variable");
    }
}

void inline xpthread_cond_signal(pthread_cond_t *cond)
{
  if (pthread_cond_signal(cond))
    {
      fatal("Unable to signal condition variable");
    }
}

void inline xpthread_cond_broadcast(pthread_cond_t *cond)
{
  if (pthread_cond_broadcast(cond))
    {
      fatal("Unable to broadcast condition variable");
    }
}
