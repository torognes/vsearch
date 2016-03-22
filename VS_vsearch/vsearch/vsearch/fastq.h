/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

struct fastq_buffer_s
{
  char * data;
  unsigned long length;
  unsigned long alloc;
  unsigned long position;
};

struct fastq_s
{
  FILE * fp;

#ifdef HAVE_ZLIB_H
  gzFile fp_gz;
#endif

#ifdef HAVE_BZLIB_H
  BZFILE * fp_bz;
#endif

  struct fastq_buffer_s file_buffer;

  struct fastq_buffer_s header_buffer;
  struct fastq_buffer_s sequence_buffer;
  struct fastq_buffer_s quality_buffer;

  unsigned long file_size;
  unsigned long file_position;

  unsigned long lineno;
  unsigned long lineno_start;
  long seqno;

  unsigned long stripped_all;
  unsigned long stripped[256];

  int format;

};

typedef struct fastq_s * fastq_handle;

fastq_handle fastq_open(const char * filename);
void fastq_close(fastq_handle h);
bool fastq_next(fastq_handle h,
                bool truncateatspace,
                char * char_mapping);
unsigned long fastq_get_position(fastq_handle h);
unsigned long fastq_get_size(fastq_handle h);
unsigned long fastq_get_lineno(fastq_handle h);
unsigned long fastq_get_seqno(fastq_handle h);
char * fastq_get_header(fastq_handle h);
char * fastq_get_sequence(fastq_handle h);
char * fastq_get_quality(fastq_handle h);
long fastq_get_abundance(fastq_handle h);
unsigned long fastq_get_header_length(fastq_handle h);
unsigned long fastq_get_sequence_length(fastq_handle h);
unsigned long fastq_get_quality_length(fastq_handle h);

void fastq_print(FILE * fp, char * header, char * sequence, char * quality);

void fastq_print_with_ee(FILE * fp, char * header, char * sequence,
                         char * quality, double ee);

void fastq_print_relabel(FILE * fp,
                         char * seq,
                         int len,
                         char * header,
                         int header_len,
                         char * quality,
                         int abundance,
                         int ordinal);


void fastq_print_db(FILE * fp, unsigned long seqno);
