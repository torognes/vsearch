/*

  VSEARCH5D: a modified version of VSEARCH

  Copyright (C) 2016, Akifumi S. Tanabe

  Contact: Akifumi S. Tanabe
  https://github.com/astanabe/vsearch5d

  Original version of VSEARCH
  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri

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

#define MD5_DIGEST_LENGTH 16
#define SHA_DIGEST_LENGTH SHA1_DIGEST_SIZE

#define LEN_DIG_MD5 MD5_DIGEST_LENGTH
#define LEN_DIG_SHA1 SHA_DIGEST_LENGTH

#define LEN_HEX_DIG_MD5 (2*LEN_DIG_MD5+1)
#define LEN_HEX_DIG_SHA1 (2*LEN_DIG_SHA1+1)

long gcd(long a, long b);
void fatal(const char * msg);
void fatal(const char * format, const char * message);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
char * xstrdup(const char *s);
char * xstrchrnul(char *s, int c);
unsigned long hash_cityhash64(char * s, unsigned long n);
long getusec(void);
void show_rusage();

void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done();

void random_init();
long random_int(long n);
unsigned long random_ulong(unsigned long n);

void string_normalize(char * normalized, char * s, unsigned int len);

void reverse_complement(char * rc, char * seq, long len);

void fprint_hex(FILE * fp, unsigned char * data, int len);

void get_hex_seq_digest_sha1(char * hex, char * seq, int seqlen);
void get_hex_seq_digest_md5(char * hex, char * seq, int seqlen);

void fprint_seq_digest_sha1(FILE * fp, char * seq, int seqlen);
void fprint_seq_digest_md5(FILE * fp, char * seq, int seqlen);
