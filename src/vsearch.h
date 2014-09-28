/*
    Copyright (C) 2014 Torbjorn Rognes and Tomas Flouri

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

#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
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
#include <ctype.h>
#include <regex.h>
#include <fcntl.h>
#include <unistd.h>
#include <city.h>

#ifndef __APPLE__
#include <sys/sysinfo.h>
#endif

#ifdef HAVE_BZLIB
#include <bzlib.h>
#endif

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

/* constants */

#define PROG_NAME "vsearch"
#define PROG_VERSION "v0.0.14"

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#ifdef HAVE_BZLIB
#define BZ_VERBOSE_0 0
#define BZ_VERBOSE_1 1
#define BZ_VERBOSE_2 2
#define BZ_VERBOSE_3 3
#define BZ_VERBOSE_4 4

#define BZ_MORE_MEM 0  /* faster decompression using more memory */
#define BZ_LESS_MEM 1  /* slower decompression but requires less memory */
#endif

#define FORMAT_PLAIN 1
#define FORMAT_BZIP  2
#define FORMAT_GZIP  3

#define MASK_ERROR -1
#define MASK_NONE 0
#define MASK_DUST 1
#define MASK_SOFT 2

/* structures and data types */

struct seqinfo_s
{
  char * header;
  char * seq;
  unsigned long headerlen;
  unsigned long headeridlen;
  unsigned long seqlen;
  unsigned long size;
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
  int strand;
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
  long trim_aln_left;
  long trim_aln_right;
};


/* macros */

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif


/* common data */

extern char * queryname;

extern char * opt_usearch_global;
extern char * opt_derep_fulllength;
extern char * opt_sortbysize;
extern char * opt_sortbylength;
extern char * opt_shuffle;
extern char * opt_mask;

extern long opt_hardmask;
extern long opt_qmask;
extern long opt_dbmask;

extern char * alnoutfilename;
extern char * useroutfilename;
extern char * blast6outfilename;
extern char * ucfilename;
extern char * opt_fastapairs;
extern char * opt_matched;
extern char * opt_notmatched;
extern char * opt_dbmatched;
extern char * opt_dbnotmatched;

extern char * databasefilename;
extern char * opt_output;
extern char * opt_relabel;

extern long opt_threads;
extern long opt_fulldp;

extern long wordlength;
extern long maxrejects;
extern long maxaccepts;
extern long match_score;
extern long mismatch_score;
extern long rowlen;
extern long opt_seed;

extern double opt_id;
extern double opt_query_cov;
extern double opt_target_cov;
extern long opt_idprefix;
extern long opt_idsuffix;
extern double opt_minqt;
extern double opt_maxqt;
extern double opt_minsl;
extern double opt_maxsl;
extern long opt_leftjust;
extern long opt_rightjust;
extern long opt_self;
extern long opt_selfid;
extern double opt_maxid;
extern double opt_minsizeratio;
extern double opt_maxsizeratio;
extern long opt_maxdiffs;
extern long opt_maxsubs;
extern long opt_maxgaps;
extern long opt_mincols;
extern long opt_maxqsize;
extern long opt_mintsize;
extern double opt_mid;

extern double opt_weak_id;
extern long opt_strand;
extern long opt_uc_allhits;
extern long opt_output_no_hits;
extern long opt_notrunclabels;
extern long opt_minsize;
extern long opt_maxsize;
extern long opt_sizein;
extern long opt_sizeout;
extern long opt_minseqlength;
extern long opt_maxseqlength;
extern long opt_minuniquesize;
extern long opt_topn;
extern long opt_maxhits;
extern long opt_top_hits_only;
extern long opt_fasta_width;

extern int opt_gap_open_query_left;
extern int opt_gap_open_target_left;
extern int opt_gap_open_query_interior;
extern int opt_gap_open_target_interior;
extern int opt_gap_open_query_right;
extern int opt_gap_open_target_right;
extern int opt_gap_extension_query_left;
extern int opt_gap_extension_target_left;
extern int opt_gap_extension_query_interior;
extern int opt_gap_extension_target_interior;
extern int opt_gap_extension_query_right;
extern int opt_gap_extension_target_right;

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

extern unsigned long longestdbsequence;

extern queryinfo_t query;

extern regex_t db_regexp;


/* search16.cc */

typedef signed short CELL;
typedef unsigned short WORD;
typedef unsigned char BYTE;
struct s16info_s;

struct s16info_s *
search16_init(CELL score_match,
              CELL score_mismatch,
              CELL penalty_gap_open_query_left,
              CELL penalty_gap_open_target_left,
              CELL penalty_gap_open_query_interior,
              CELL penalty_gap_open_target_interior,
              CELL penalty_gap_open_query_right,
              CELL penalty_gap_open_target_right,
              CELL penalty_gap_extension_query_left,
              CELL penalty_gap_extension_target_left,
              CELL penalty_gap_extension_query_interior,
              CELL penalty_gap_extension_target_interior,
              CELL penalty_gap_extension_query_right,
              CELL penalty_gap_extension_target_right);

void
search16_exit(s16info_s * s);

void
search16_qprep(s16info_s * s, char * qseq, int qlen);

void
search16(s16info_s * s,
         unsigned int sequences,
         unsigned int * seqnos,
         CELL * pscores,
         unsigned short * paligned,
         unsigned short * pmatches,
         unsigned short * pmismatches,
         unsigned short * pgaps,
         char * * pcigar);

  
/* maps.cc */

extern char sym_nt_2bit[5];
extern char sym_nt_4bit[17];

extern unsigned int chrstatus[256];
extern unsigned int chrmap_2bit[256];
extern unsigned int chrmap_4bit[256];
extern char chrmap_complement[256];


/* functions in arch.cc */

unsigned long arch_get_memused();
unsigned long arch_get_memtotal();
void arch_srandom_init();


/* functions in util.cc */

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


/* functions in db.cc */

inline char * db_getheader(unsigned long seqno)
{
  return seqindex[seqno].header;
}

inline char * db_getsequence(unsigned long seqno)
{
  return seqindex[seqno].seq;
}

inline unsigned long db_getabundance(unsigned long seqno)
{
  return seqindex[seqno].size;
}

inline unsigned long db_getsequencelen(unsigned long seqno)
{
  return seqindex[seqno].seqlen;
}

inline unsigned long db_getheaderlen(unsigned long seqno)
{
  return seqindex[seqno].headerlen;
}

void db_read(const char * filename, int upcase);
void db_free();

unsigned long db_getsequencecount();
unsigned long db_getnucleotidecount();
unsigned long db_getlongestheader();
unsigned long db_getlongestsequence();
unsigned long db_getshortestsequence();


/* functions in query.cc */

void query_open(const char * filename);

int query_getnext(char ** head, long * head_len,
                  char ** seq, long * seq_len, long * qno,
                  long * qsize, int upcase);

void query_close();

long query_getfilesize();

long query_getfilepos();


/* functions in nw.cc */

struct nwinfo_s;

struct nwinfo_s * nw_init();

void nw_exit(struct nwinfo_s * nw);

void nw_align(char * dseq,
              char * dend,
              char * qseq,
              char * qend,
              long * score_matrix,
              long gapopen_q_left,
              long gapopen_q_internal,
              long gapopen_q_right,
              long gapopen_t_left,
              long gapopen_t_internal,
              long gapopen_t_right,
              long gapextend_q_left,
              long gapextend_q_internal,
              long gapextend_q_right,
              long gapextend_t_left,
              long gapextend_t_internal,
              long gapextend_t_right,
              long * nwscore,
              long * nwdiff,
              long * nwgaps,
              long * nwindels,
              long * nwalignmentlength,
              char ** nwalignment,
              long queryno,
              long dbseqno,
              struct nwinfo_s * nw);


/* functions in unique.cc */

struct bucket_s;
struct uhandle_s;

struct uhandle_s * unique_init();

void unique_exit(struct uhandle_s * u);

void unique_count(struct uhandle_s * uh, 
                  int k,
                  int seqlen,
                  char * seq,
                  unsigned int * listlen,
                  unsigned int * * list);


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

void search(char * cmdline, char * progheader);


/* functions in showalign.cc */

char * align_getrow(char * seq, char * cigar, int alignlen, int origin);

void align_fprint_uncompressed_alignment(FILE * f, char * cigar);

void align_show(FILE * f,
                char * seq1,
                long seq1len,
                long seq1off,
                const char * seq1name,
                char * seq2,
                long seq2len,
                long seq2off,
                const char * seq2name,
                char * cigar,
                long cigarlen,
                int numwidth,
                int namewidth,
                int alignwidth,
                int strand);


/* functions in userfields.cc */

extern int * userfields_requested;
extern int userfields_requested_count;

int parse_userfields_arg(char * arg);

/* functions in results.cc */

void results_show_alnout(FILE * fp,
                         struct hit * hits,
                         int hitcount,
                         char * query_head,
                         char * qsequence,
                         long qseqlen,
                         char * rc);

void results_show_blast6out_one(FILE * fp,
                                struct hit * hp,
                                char * query_head,
                                char * qsequence,
                                long qseqlen,
                                char * rc);

void results_show_uc_one(FILE * fp,
                         struct hit * hp,
                         char * query_head,
                         char * qsequence,
                         long qseqlen,
                         char * rc);

void results_show_userout_one(FILE * fp,
                              struct hit * hp,
                              char * query_head,
                              char * qsequence,
                              long qseqlen,
                              char * rc);

void results_show_fastapairs_one(FILE * fp,
                                 struct hit * hp,
                                 char * query_head,
                                 char * qsequence,
                                 long qseqlen,
                                 char * rc);

/* functions in sortbysize.cc */

void sortbysize();


/* functions in sortbylength.cc */

void sortbylength();


/* functions in derep.cc */

void derep_fulllength();


/* functions in shuffle.cc */

void shuffle();


/* functions in mask.cc */

void mask();
void dust(char * m, int len);
void hardmask(char * m, int len);
void dust_all();
void hardmask_all();


/* minheap.cc */

typedef struct topscore
{
  unsigned int count;
  unsigned int seqno;
  unsigned int length;
} elem_t;


typedef struct minheap_s
{
  int alloc;
  int count;
  elem_t * array;
} minheap_t;

inline int minheap_isempty(minheap_t * m)
{
  return !m->count;
}

inline void minheap_empty(minheap_t * m)
{
  m->count = 0;
}

elem_t minheap_poplast(minheap_t * m);
void minheap_sort(minheap_t * m);
minheap_t * minheap_init(int size);
void minheap_exit(minheap_t * m);
void minheap_add(minheap_t * m, elem_t * n);
elem_t minheap_pop(minheap_t * m);
