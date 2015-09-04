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

struct fasta_buffer_s
{
  char * data;
  unsigned long length;
  unsigned long alloc;
  unsigned long position;
};

struct fasta_s
{
  FILE * fp;

#ifdef HAVE_ZLIB
  gzFile fp_gz;
#endif

#ifdef HAVE_BZLIB
  BZFILE * fp_bz;
#endif

  regex_t size_regexp;

  struct fasta_buffer_s file_buffer;
  struct fasta_buffer_s header_buffer;
  struct fasta_buffer_s sequence_buffer;

  unsigned long file_size;
  unsigned long file_position;

  unsigned long lineno;
  long seqno;

  unsigned long stripped_all;
  unsigned long stripped[256];

  int format;
};

typedef struct fasta_s * fasta_handle;

fasta_handle fasta_open(const char * filename);
void fasta_close(fasta_handle h);
bool fasta_next(fasta_handle h,
                bool truncateatspace,
                char * char_mapping);
unsigned long fasta_get_position(fasta_handle h);
unsigned long fasta_get_size(fasta_handle h);
unsigned long fasta_get_lineno(fasta_handle h);
unsigned long fasta_get_seqno(fasta_handle h);
char * fasta_get_header(fasta_handle h);
char * fasta_get_sequence(fasta_handle h);
unsigned long fasta_get_header_length(fasta_handle h);
unsigned long fasta_get_sequence_length(fasta_handle h);
long fasta_get_abundance(fasta_handle h);
