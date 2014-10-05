/*
    Copyright (C) 2014 Torbjorn Rognes & Tomas Flouri

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

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

long gcd(long a, long b);
void fatal(const char * msg);
void fatal(const char * format, const char * message);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
char * xstrchrnul(char *s, int c);
unsigned long hash_cityhash64(char * s, unsigned long n);
long getusec(void);
void show_rusage();
void fprint_fasta_hdr_only(FILE * fp, char * hdr);
void fprint_fasta_seq_only(FILE * fp, char * seq, unsigned long len, int width);
void db_fprint_fasta(FILE * fp, unsigned long seqno);
void db_fprint_fasta_with_size(FILE * fp, unsigned long seqno, unsigned long size);
void reverse_complement(char * rc, char * seq, long len);

void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done();

int detect_compress_format (const char * filename);

#ifdef HAVE_BZLIB
char * bz_fgets (char * s, int size, BZFILE * stream, long linealloc,
                 int * bz_error_ptr, char * buf_internal, long * buf_internal_len);
#endif
