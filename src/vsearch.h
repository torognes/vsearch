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

#ifdef HAVE_BZLIB_H
#include <bzlib.h>
#endif


/* constants */

#define PROG_NAME "vsearch"
#define PROG_VERSION "v0.0.5"

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#ifdef HAVE_BZLIB_H
#define BZ_VERBOSE_0 0
#define BZ_VERBOSE_1 1
#define BZ_VERBOSE_2 2
#define BZ_VERBOSE_3 3
#define BZ_VERBOSE_4 4

#define BZ_MORE_MEM 0  /* faster decompression using more memory */
#define BZ_LESS_MEM 1  /* slower decompression but requires less memory */
#endif

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

struct seqinfo_s
{
  char * header;
  char * seq;
  unsigned long headerlen;
  unsigned long headeridlen;
  unsigned long seqlen;
  unsigned long dummy;
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


struct hit
{
  long target;
  char strand;
  long count; /* number of word matches */

  long matches;
  long mismatches;

  /* info for global alignment, entire sequences */

  long nwscore;
  long nwdiff; /* indels and mismatches in global alignment */
  long nwgaps; /* gaps in global alignment */
  long nwindels; /* indels in global alignment */
  long nwalignmentlength; /* length of global alignment */
  double nwid; /* percent identity of global alignment */
  char * nwalignment; /* alignment string (cigar) of global alignment */

  /* info for semi-global alignment, excluding gaps at ends */

  long internal_alignmentlength;
  long internal_gaps;
  long internal_indels;
  double internal_id;

  long trim_q_left;
  long trim_q_right;
  long trim_t_left;
  long trim_t_right;
};



/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))


/* common data */

extern char * queryname;
extern char * matrixname;

extern char * databasefilename;
extern char * alnoutfilename;
extern char * useroutfilename;
extern char * blast6outfilename;
extern char * ucfilename;

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
extern int opt_self;

extern FILE * alnoutfile;
extern FILE * useroutfile;
extern FILE * blast6outfile;
extern FILE * ucfile;

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


/* functions in query.cc */

void query_open(const char * filename);

int query_getnext(char ** header, long * header_length,
                  char ** seq, long * seq_length, long * query_no);

void query_close();

long query_getfilesize();

long query_getfilepos();


/* functions in bzquery.cc */

void query_bz_open(const char * filename);

int query_bz_getnext(char ** header, long * header_length,
                  char ** seq, long * seq_length, long * query_no);

void query_bz_close();


/* functions in nw.cc */

void nw_init();

void nw_exit();

void nw_align(char * dseq,
              char * dend,
              char * qseq,
              char * qend,
              long * score_matrix,
	      unsigned long gapopen_q_left,
	      unsigned long gapopen_q_internal,
	      unsigned long gapopen_q_right,
	      unsigned long gapopen_t_left,
	      unsigned long gapopen_t_internal,
	      unsigned long gapopen_t_right,
	      unsigned long gapextend_q_left,
	      unsigned long gapextend_q_internal,
	      unsigned long gapextend_q_right,
	      unsigned long gapextend_t_left,
	      unsigned long gapextend_t_internal,
	      unsigned long gapextend_t_right,
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

char * align_getrow(char * seq, char * cigar, int alignlen, int origin);

void align_show(FILE * f,
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

/* functions in userfields.cc */

extern int * userfields_requested;
extern int userfields_requested_count;

int parse_userfields_arg(char * arg);

/* functions in results.cc */

void results_show_blast6out(struct hit * hits, int accepts,
			   char * query_head,
			   char * qsequence, long qseqlen);

void results_show_uc(struct hit * hits, int accepts,
		    char * query_head,
		    char * qsequence, long qseqlen,
		    int allhits);

void results_show_userout(struct hit * hits, int accepts,
			 char * query_head,
			 char * qsequence, long qseqlen);

void results_show_alnout(struct hit * hits, int accepts,
                        char * query_head,
                        char * qsequence, long qseqlen);
