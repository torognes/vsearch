/*
    Copyright (C) 2014-2015 Torbjorn Rognes and Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

struct fastq_buffer_s
{
  char * data;
  size_t length;
  size_t alloc;
  size_t position;
};

struct fastq_s
{
  FILE * fp;

#ifdef HAVE_ZLIB
  gzFile fp_gz;
#endif

#ifdef HAVE_BZLIB
  BZFILE * fp_bz;
#endif

  int format;

  struct fastq_buffer_s file_buffer;

  struct fastq_buffer_s header_buffer;
  struct fastq_buffer_s sequence_buffer;
  struct fastq_buffer_s quality_buffer;

  size_t file_size;
  size_t file_position;

  size_t lineno;
  long seqno;

  size_t stripped_all;
  size_t stripped[256];
};

typedef struct fastq_s * fastq_handle;

fastq_handle fastq_open(const char * filename);
void fastq_close(fastq_handle h);
bool fastq_next(fastq_handle h,
                bool truncateatspace,
                unsigned int * char_action,
                char * char_mapping);
size_t fastq_get_position(fastq_handle h);
size_t fastq_get_size(fastq_handle h);
size_t fastq_get_lineno(fastq_handle h);
size_t fastq_get_seqno(fastq_handle h);
char * fastq_get_header(fastq_handle h);
char * fastq_get_sequence(fastq_handle h);
char * fastq_get_quality(fastq_handle h);
size_t fastq_get_header_length(fastq_handle h);
size_t fastq_get_sequence_length(fastq_handle h);
size_t fastq_get_quality_length(fastq_handle h);
