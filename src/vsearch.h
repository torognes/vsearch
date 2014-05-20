/*
    Copyright (C) 2014 Torbjorn Rognes

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

#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <x86intrin.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>


/* constants */

#define PROG_NAME "vsearch"
#define PROG_VERSION "v0.0.4"

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#define KMERLENGTH 4
#define KMERVECTORBITS (1<<(2*KMERLENGTH))
#define KMERVECTORBYTES (KMERVECTORBITS/8)

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

struct seqinfo_s
{
  /* 16 byte alignment of kmervector required for SSE4 version of kmerdiff */
  unsigned char kmervector[KMERVECTORBYTES];

  char * header;
  char * seq;
  unsigned long headerlen;
  unsigned long headeridlen;
  unsigned long seqlen;
  unsigned long dummy;

  /* 32 + 6*8 = 80 bytes */
};

typedef struct seqinfo_s seqinfo_t;

extern seqinfo_t * seqindex;

struct queryinfo
{
  unsigned long qno;
  long len;
  char * seq;
};

typedef struct queryinfo queryinfo_t;


/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))


/* common data */

extern char * queryname;
extern char * matrixname;

extern char * databasefilename;
extern char * alnoutfilename;

extern int wordlength;
extern int maxrejects;
extern int maxaccepts;
extern int match_score;
extern int mismatch_score;
extern double identity;
extern int mismatch_cost;
extern int gapopen_cost;
extern int gapextend_cost;
extern int rowlen;

extern FILE * alnoutfile;

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;

extern char sym_nt[];

extern unsigned long longestdbsequence;

extern queryinfo_t query;


/* functions in util.cc */

long gcd(long a, long b);
void fatal(const char * msg);
void fatal(const char * format, const char * message);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
char * xstrchrnul(char *s, int c);
unsigned long hash_fnv_1a_64(unsigned char * s, unsigned long n);
long getusec(void);
void show_rusage();


/* functions in kmer.cc */

void printkmers(unsigned char * kmervector);

void findkmers(unsigned char * seq, unsigned long seqlen,
                unsigned char * kmervector);

unsigned long comparekmervectors(unsigned char * a, unsigned char * b);

unsigned long kmer_diff(unsigned long a, unsigned long b);

void kmer_diff_parallel(unsigned long seed,
                         unsigned long listlen,
                         unsigned long * amplist,
                         unsigned long * difflist);


/* functions in db.cc */

void db_read(const char * filename);

unsigned long db_getsequencecount();
unsigned long db_getnucleotidecount();

unsigned long db_getlongestheader();
unsigned long db_getlongestsequence();
unsigned long db_getshortestsequence();

seqinfo_t * db_getseqinfo(unsigned long seqno);

char * db_getsequence(unsigned long seqno);
unsigned long db_getsequencelen(unsigned long seqno);

void db_getsequenceandlength(unsigned long seqno,
                             char ** address,
                             long * length);

char * db_getheader(unsigned long seqno);
unsigned long db_getheaderlen(unsigned long seqno);

unsigned long db_getabundance(unsigned long seqno);

void db_showsequence(unsigned long seqno);
void db_showall();
void db_free();

void db_putseq(long seqno);

inline unsigned char * db_getkmervector(unsigned long seqno)
{
  return seqindex[seqno].kmervector;
}


/* functions in query.cc */

void query_open(const char * filename);

int query_getnext(char ** header, long * header_length,
                  char ** seq, long * seq_length, long * query_no);

void query_close();


/* functions in popcount.cc */

unsigned long popcount(unsigned long x);
unsigned long popcount_128(__m128i x);


/* functions in nw.cc */

void nw_init();

void nw_exit();

void nw_align(char * dseq,
              char * dend,
              char * qseq,
              char * qend,
              long * score_matrix,
              unsigned long gapopen,
              unsigned long gapextend,
              unsigned long * nwscore,
              unsigned long * nwdiff,
              unsigned long * nwgaps,
              unsigned long * nwindels,
              unsigned long * nwalignmentlength,
              char ** nwalignment,
              unsigned long queryno,
              unsigned long dbseqno);


/* functions in kmercount.cc */

struct kmercountelem
{
  unsigned int kmer;
  unsigned int count;
};

extern struct kmercountelem * kmercounthash;

void count_kmers_init();
void count_kmers_exit();
unsigned int count_kmers_gethashsize();
void count_kmers(unsigned int k, char * seq, unsigned int seqlen);
unsigned int count_kmers_unique();

unsigned int count_kmers_getcount(unsigned int wordlength, unsigned kmer);


/* functions in dbindex.cc */

extern unsigned int * kmerhash;
extern unsigned int * kmerindex;

extern unsigned int kmerhashsize;
extern unsigned int kmerindexsize;

void fprint_kmer(FILE * f, unsigned int k, unsigned long kmer);
void dbindex_build();
void dbindex_free();
int dbindex_getkmermatchcount(int kmer);
int dbindex_getkmermatch(int kmer, int matchno);


/* functions in search.cc */

void search();


/* functions in showalign.cc */

void showalign(FILE * f,
               char * seq1,
               long seq1len,
               const char * seq1name,
               char * seq2,
               long seq2len,
               const char * seq2name,
               char * cigar,
               int numwidth,
               int namewidth,
               int alignwidth);
