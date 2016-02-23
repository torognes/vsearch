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

#include "vsearch.h"

/* options */

bool opt_fastq_allowmergestagger;
bool opt_fastq_nostagger;
bool opt_fastq_eeout;
bool opt_clusterout_id;
bool opt_clusterout_sort;
bool opt_eeout;
bool opt_quiet;
bool opt_relabel_keep;
bool opt_relabel_md5;
bool opt_relabel_sha1;
bool opt_samheader;
bool opt_sizeorder;
bool opt_xsize;
char * opt_eetabbedout;
char * opt_fastaout_notmerged_fwd;
char * opt_fastaout_notmerged_rev;
char * opt_fastq_mergepairs;
char * opt_fastqout_notmerged_fwd;
char * opt_fastqout_notmerged_rev;
char * opt_allpairs_global;
char * opt_alnout;
char * opt_blast6out;
char * opt_borderline;
char * opt_centroids;
char * opt_chimeras;
char * opt_cluster_fast;
char * opt_cluster_size;
char * opt_cluster_smallmem;
char * opt_clusters;
char * opt_consout;
char * opt_db;
char * opt_dbmatched;
char * opt_dbnotmatched;
char * opt_derep_fulllength;
char * opt_derep_prefix;
char * opt_fastaout;
char * opt_fastaout_discarded;
char * opt_fastapairs;
char * opt_fastq_chars;
char * opt_fastq_convert;
char * opt_fastq_filter;
char * opt_fastq_stats;
char * opt_fastqout;
char * opt_fastqout_discarded;
char * opt_fastx_mask;
char * opt_fastx_revcomp;
char * opt_fastx_subsample;
char * opt_label_suffix;
char * opt_log;
char * opt_maskfasta;
char * opt_matched;
char * opt_msaout;
char * opt_nonchimeras;
char * opt_notmatched;
char * opt_output;
char * opt_pattern;
char * opt_profile;
char * opt_relabel;
char * opt_samout;
char * opt_search_exact;
char * opt_shuffle;
char * opt_sortbylength;
char * opt_sortbysize;
char * opt_uc;
char * opt_uchime_denovo;
char * opt_uchime_ref;
char * opt_uchimealns;
char * opt_uchimeout;
char * opt_usearch_global;
char * opt_userout;
char * opt_reverse;
double opt_abskew;
double opt_dn;
double opt_fastq_maxee;
double opt_fastq_maxee_rate;
double opt_id;
double opt_max_unmasked_pct;
double opt_maxid;
double opt_maxqt;
double opt_maxsizeratio;
double opt_maxsl;
double opt_mid;
double opt_min_unmasked_pct;
double opt_mindiv;
double opt_minh;
double opt_minqt;
double opt_minsizeratio;
double opt_minsl;
double opt_query_cov;
double opt_sample_pct;
double opt_target_cov;
double opt_weak_id;
double opt_xn;
int opt_acceptall;
int opt_alignwidth;
int opt_cons_truncate;
int opt_gap_extension_query_interior;
int opt_gap_extension_query_left;
int opt_gap_extension_query_right;
int opt_gap_extension_target_interior;
int opt_gap_extension_target_left;
int opt_gap_extension_target_right;
int opt_gap_open_query_interior;
int opt_gap_open_query_left;
int opt_gap_open_query_right;
int opt_gap_open_target_interior;
int opt_gap_open_target_left;
int opt_gap_open_target_right;
int opt_help;
int opt_mindiffs;
int opt_slots;
int opt_uchimeout5;
int opt_usersort;
int opt_version;
long opt_fastq_maxdiffs;
long opt_fastq_maxmergelen;
long opt_fastq_minmergelen;
long opt_fastq_minovlen;
long opt_dbmask;
long opt_fasta_width;
long opt_fastq_ascii;
long opt_fastq_asciiout;
long opt_fastq_maxns;
long opt_fastq_minlen;
long opt_fastq_qmax;
long opt_fastq_qmaxout;
long opt_fastq_qmin;
long opt_fastq_qminout;
long opt_fastq_stripleft;
long opt_fastq_tail;
long opt_fastq_trunclen;
long opt_fastq_truncqual;
long opt_fulldp;
long opt_hardmask;
long opt_iddef;
long opt_idprefix;
long opt_idsuffix;
long opt_leftjust;
long opt_match;
long opt_maxaccepts;
long opt_maxdiffs;
long opt_maxgaps;
long opt_maxhits;
long opt_maxqsize;
long opt_maxrejects;
long opt_maxseqlength;
long opt_maxsize;
long opt_maxsubs;
long opt_maxuniquesize;
long opt_mincols;
long opt_minseqlength;
long opt_minsize;
long opt_mintsize;
long opt_minuniquesize;
long opt_minwordmatches;
long opt_mismatch;
long opt_notrunclabels;
long opt_output_no_hits;
long opt_qmask;
long opt_randseed;
long opt_rightjust;
long opt_rowlen;
long opt_sample_size;
long opt_self;
long opt_selfid;
long opt_sizein;
long opt_sizeout;
long opt_strand;
long opt_threads;
long opt_top_hits_only;
long opt_topn;
long opt_uc_allhits;
long opt_wordlength;
long opt_idoffset;

/* Other variables */

/* cpu features available */

long mmx_present = 0;
long sse_present = 0;
long sse2_present = 0;
long sse3_present = 0;
long ssse3_present = 0;
long sse41_present = 0;
long sse42_present = 0;
long popcnt_present = 0;
long avx_present = 0;
long avx2_present = 0;

static char * progname;
static char progheader[80];
static char * cmdline;
static time_t time_start;
static time_t time_finish;

FILE * fp_log = 0;

abundance_t * global_abundance;

#define cpuid(f1, f2, a, b, c, d)                                \
  __asm__ __volatile__ ("cpuid"                                  \
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d) \
                        : "a" (f1), "c" (f2));

void cpu_features_detect()
{
  unsigned int a, b, c, d;

  cpuid(0, 0, a, b, c, d);
  unsigned int maxlevel = a & 0xff;

  if (maxlevel >= 1)
  {
    cpuid(1, 0, a, b, c, d);
    mmx_present    = (d >> 23) & 1;
    sse_present    = (d >> 25) & 1;
    sse2_present   = (d >> 26) & 1;
    sse3_present   = (c >>  0) & 1;
    ssse3_present  = (c >>  9) & 1;
    sse41_present  = (c >> 19) & 1;
    sse42_present  = (c >> 20) & 1;
    popcnt_present = (c >> 23) & 1;
    avx_present    = (c >> 28) & 1;
    
    if (maxlevel >= 7)
    {
      cpuid(7, 0, a, b, c, d);
      avx2_present = (b >>  5) & 1;
    }
  }
}

void cpu_features_show()
{
  fprintf(stderr, "CPU features:");
  if (mmx_present)
    fprintf(stderr, " mmx");
  if (sse_present)
    fprintf(stderr, " sse");
  if (sse2_present)
    fprintf(stderr, " sse2");
  if (sse3_present)
    fprintf(stderr, " sse3");
  if (ssse3_present)
    fprintf(stderr, " ssse3");
  if (sse41_present)
    fprintf(stderr, " sse4.1");
  if (sse42_present)
    fprintf(stderr, " sse4.2");
  if (popcnt_present)
    fprintf(stderr, " popcnt");
  if (avx_present)
    fprintf(stderr, " avx");
  if (avx2_present)
    fprintf(stderr, " avx2");
  fprintf(stderr, "\n");
}


void args_get_gap_penalty_string(char * arg, int is_open)
{
  /* See http://www.drive5.com/usearch/manual/aln_params.html
     
     --gapopen *E/10I/1E/2L/3RQ/4RT/1IQ
     --gapext *E/10I/1E/2L/3RQ/4RT/1IQ
     
     integer or *
     followed by I, E, L, R, Q or T characters
     separated by /
     * means infinitely high (disallow)        
     E=end
     I=interior
     L=left
     R=right
     Q=query
     T=target

     E cannot be combined with L or R
     
     We do not support floating point values. Therefore,
     all default score and penalties are multiplied by 2.

  */
  
  char *p = arg;

  while (*p)
    {
      int skip = 0;
      int pen = 0;

      if (sscanf(p, "%d%n", &pen, &skip) == 1)
        {
          p += skip;
        }
      else if (*p == '*')
        {
          pen = 1000;
          p++;
        }
      else
        fatal("Invalid gap penalty argument (%s)", p);

      char * q = p;

      int set_E = 0;
      int set_I = 0;
      int set_L = 0;
      int set_R = 0;
      int set_Q = 0;
      int set_T = 0;

      while((*p) && (*p != '/'))
        {
          switch(*p)
            {
            case 'E':
              set_E = 1;
              break;
            case 'I':
              set_I = 1;
              break;
            case 'L':
              set_L = 1;
              break;
            case 'R':
              set_R = 1;
              break;
            case 'Q':
              set_Q = 1;
              break;
            case 'T':
              set_T = 1;
              break;
            default:
              fatal("Invalid char '%.1s' in gap penalty string", p);
              break;
            }
          p++;
        }
      
      if (*p == '/')
        p++;

      if (set_E && (set_L || set_R))
        fatal("Invalid gap penalty string (E and L or R) '%s'", q);
      
      if (set_E)
        {
          set_L = 1;
          set_R = 1;
        }

      /* if neither L, I, R nor E is specified, it applies to all */

      if ((!set_L) && (!set_I) && (!set_R))
        {
          set_L = 1;
          set_I = 1;
          set_R = 1;
        }

      /* if neither Q nor T is specified, it applies to both */

      if ((!set_Q) && (!set_T))
        {
          set_Q = 1;
          set_T = 1;
        }

      if (is_open)
        {
          if (set_Q)
            {
              if (set_L)
                opt_gap_open_query_left = pen;
              if (set_I)
                opt_gap_open_query_interior = pen;
              if (set_R)
                opt_gap_open_query_right = pen;
            }
          if (set_T)
            {
              if (set_L)
                opt_gap_open_target_left = pen;
              if (set_I)
                opt_gap_open_target_interior = pen;
              if (set_R)
                opt_gap_open_target_right = pen;
            }     
        }
      else
        {
          if (set_Q)
            {
              if (set_L)
                opt_gap_extension_query_left = pen;
              if (set_I)
                opt_gap_extension_query_interior = pen;
              if (set_R)
                opt_gap_extension_query_right = pen;
            }
          if (set_T)
            {
              if (set_L)
                opt_gap_extension_target_left = pen;
              if (set_I)
                opt_gap_extension_target_interior = pen;
              if (set_R)
                opt_gap_extension_target_right = pen;
            }     
        }
    }
}


long args_getlong(char * arg)
{
  int len = 0;
  long temp = 0;
  int ret = sscanf(arg, "%ld%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}

double args_getdouble(char * arg)
{
  int len = 0;
  double temp = 0;
  int ret = sscanf(arg, "%lf%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}

void args_init(int argc, char **argv)
{
  /* Set defaults */

  progname = argv[0];

  opt_abskew = 2.0;
  opt_acceptall = 0;
  opt_alignwidth = 80;
  opt_allpairs_global = 0;
  opt_alnout = 0;
  opt_blast6out = 0;
  opt_borderline = 0;
  opt_centroids = 0;
  opt_chimeras = 0;
  opt_cluster_fast = 0;
  opt_cluster_size = 0;
  opt_cluster_smallmem = 0;
  opt_clusterout_id = 0;
  opt_clusterout_sort = 0;
  opt_clusters = 0;
  opt_cons_truncate = 0;
  opt_consout = 0;
  opt_db = 0;
  opt_dbmask = MASK_DUST;
  opt_dbmatched = 0;
  opt_dbnotmatched = 0;
  opt_derep_fulllength = 0;
  opt_derep_prefix = 0;
  opt_dn = 1.4;
  opt_eeout = 0;
  opt_eetabbedout = 0;
  opt_fastaout_notmerged_fwd = 0;
  opt_fastaout_notmerged_rev = 0;
  opt_fasta_width = 80;
  opt_fastaout = 0;
  opt_fastaout_discarded = 0;
  opt_fastapairs = 0;
  opt_fastq_allowmergestagger = 0;
  opt_fastq_ascii = 33;
  opt_fastq_asciiout = 33;
  opt_fastq_chars = 0;
  opt_fastq_convert = 0;
  opt_fastq_eeout = 0;
  opt_fastq_filter = 0;
  opt_fastq_maxdiffs = 1000000;
  opt_fastq_maxee = DBL_MAX;
  opt_fastq_maxee_rate = DBL_MAX;
  opt_fastq_maxmergelen  = 1000000;
  opt_fastq_maxns = LONG_MAX;
  opt_fastq_mergepairs = 0;
  opt_fastq_minlen = 1;
  opt_fastq_minmergelen = 0;
  opt_fastq_minovlen = 16;
  opt_fastq_nostagger = 1;
  opt_fastqout_notmerged_fwd = 0;
  opt_fastqout_notmerged_rev = 0;
  opt_fastq_qmax = 41;
  opt_fastq_qmaxout = 41;
  opt_fastq_qmin = 0;
  opt_fastq_qminout = 0;
  opt_fastq_stats = 0;
  opt_fastq_stripleft = 0;
  opt_fastq_tail = 4;
  opt_fastq_trunclen = 0;
  opt_fastq_truncqual = LONG_MIN;
  opt_fastqout = 0;
  opt_fastqout_discarded = 0;
  opt_fastx_mask = 0;
  opt_fastx_revcomp = 0;
  opt_fastx_subsample = 0;
  opt_fulldp = 0;
  opt_gap_extension_query_interior=2;
  opt_gap_extension_query_left=1;
  opt_gap_extension_query_right=1;
  opt_gap_extension_target_interior=2;
  opt_gap_extension_target_left=1;
  opt_gap_extension_target_right=1;
  opt_gap_open_query_interior=20;
  opt_gap_open_query_left=2;
  opt_gap_open_query_right=2;
  opt_gap_open_target_interior=20;
  opt_gap_open_target_left=2;
  opt_gap_open_target_right=2;
  opt_hardmask = 0;
  opt_help = 0;
  opt_id = -1.0;
  opt_iddef = 2;
  opt_idoffset = 0;
  opt_idprefix = 0;
  opt_idsuffix = 0;
  opt_label_suffix = 0;
  opt_leftjust = 0;
  opt_log = 0;
  opt_maskfasta = 0;
  opt_match = 2;
  opt_matched = 0;
  opt_max_unmasked_pct = 100.0;
  opt_maxaccepts = 1;
  opt_maxdiffs = INT_MAX;
  opt_maxgaps = INT_MAX;
  opt_maxhits = LONG_MAX;
  opt_maxid = 1.0;
  opt_maxqsize = INT_MAX;
  opt_maxqt = DBL_MAX;
  opt_maxrejects = -1;
  opt_maxseqlength = 50000;
  opt_maxsize = LONG_MAX;
  opt_maxsizeratio = DBL_MAX;
  opt_maxsl = DBL_MAX;
  opt_maxsubs = INT_MAX;
  opt_maxuniquesize = LONG_MAX;
  opt_mid = 0.0;
  opt_min_unmasked_pct = 0.0;
  opt_mincols = 0;
  opt_mindiffs = 3;
  opt_mindiv = 0.8;
  opt_minh = 0.28;
  opt_minqt = 0.0;
  opt_minseqlength = 0;
  opt_minsize = 0;
  opt_minsizeratio = 0.0;
  opt_minsl = 0.0;
  opt_mintsize = 0;
  opt_minuniquesize = 0;
  opt_minwordmatches = 0;
  opt_mismatch = -4;
  opt_msaout = 0;
  opt_nonchimeras = 0;
  opt_notmatched = 0;
  opt_notrunclabels = 0;
  opt_output = 0;
  opt_output_no_hits = 0;
  opt_pattern = 0;
  opt_profile = 0;
  opt_qmask = MASK_DUST;
  opt_query_cov = 0.0;
  opt_quiet = false;
  opt_randseed = 0;
  opt_relabel = 0;
  opt_relabel_keep = 0;
  opt_relabel_md5 = 0;
  opt_relabel_sha1 = 0;
  opt_reverse = 0;
  opt_rightjust = 0;
  opt_rowlen = 64;
  opt_samheader = 0;
  opt_samout = 0;
  opt_sample_pct = 0;
  opt_sample_size = 0;
  opt_search_exact = 0;
  opt_self = 0;
  opt_selfid = 0;
  opt_shuffle = 0;
  opt_sizein = 0;
  opt_sizeorder = 0;
  opt_sizeout = 0;
  opt_slots = 0;
  opt_sortbylength = 0;
  opt_sortbysize = 0;
  opt_strand = 1;
  opt_target_cov = 0.0;
  opt_threads = 0;
  opt_top_hits_only = 0;
  opt_topn = LONG_MAX;
  opt_uc = 0;
  opt_uc_allhits = 0;
  opt_uchime_denovo = 0;
  opt_uchime_ref = 0;
  opt_uchimealns = 0;
  opt_uchimeout = 0;
  opt_uchimeout5 = 0;
  opt_usearch_global = 0;
  opt_userout = 0;
  opt_usersort = 0;
  opt_version = 0;
  opt_weak_id = 10.0;
  opt_wordlength = 8;
  opt_xn = 8.0;
  opt_xsize = 0;

  opterr = 1;

  static struct option long_options[] =

  {
    {"help",                  no_argument,       0, 0 },
    {"version",               no_argument,       0, 0 },
    {"alnout",                required_argument, 0, 0 },
    {"usearch_global",        required_argument, 0, 0 },
    {"db",                    required_argument, 0, 0 },
    {"id",                    required_argument, 0, 0 },
    {"maxaccepts",            required_argument, 0, 0 },
    {"maxrejects",            required_argument, 0, 0 },
    {"wordlength",            required_argument, 0, 0 },
    {"match",                 required_argument, 0, 0 },
    {"mismatch",              required_argument, 0, 0 },
    {"fulldp",                no_argument,       0, 0 },
    {"strand",                required_argument, 0, 0 },
    {"threads",               required_argument, 0, 0 },
    {"gapopen",               required_argument, 0, 0 },
    {"gapext",                required_argument, 0, 0 },
    {"rowlen",                required_argument, 0, 0 },
    {"userfields",            required_argument, 0, 0 },
    {"userout",               required_argument, 0, 0 },
    {"self",                  no_argument,       0, 0 },
    {"blast6out",             required_argument, 0, 0 },
    {"uc",                    required_argument, 0, 0 },
    {"weak_id",               required_argument, 0, 0 },
    {"uc_allhits",            no_argument,       0, 0 },
    {"notrunclabels",         no_argument,       0, 0 },
    {"sortbysize",            required_argument, 0, 0 },
    {"output",                required_argument, 0, 0 },
    {"minsize",               required_argument, 0, 0 },
    {"maxsize",               required_argument, 0, 0 },
    {"relabel",               required_argument, 0, 0 },
    {"sizeout",               no_argument,       0, 0 },
    {"derep_fulllength",      required_argument, 0, 0 },
    {"minseqlength",          required_argument, 0, 0 },
    {"minuniquesize",         required_argument, 0, 0 },
    {"topn",                  required_argument, 0, 0 },
    {"maxseqlength",          required_argument, 0, 0 },
    {"sizein",                no_argument,       0, 0 },
    {"sortbylength",          required_argument, 0, 0 },
    {"matched",               required_argument, 0, 0 },
    {"notmatched",            required_argument, 0, 0 },
    {"dbmatched",             required_argument, 0, 0 },
    {"dbnotmatched",          required_argument, 0, 0 },
    {"fastapairs",            required_argument, 0, 0 },
    {"output_no_hits",        no_argument,       0, 0 },
    {"maxhits",               required_argument, 0, 0 },
    {"top_hits_only",         no_argument,       0, 0 },
    {"fasta_width",           required_argument, 0, 0 },
    {"query_cov",             required_argument, 0, 0 },
    {"target_cov",            required_argument, 0, 0 },
    {"idprefix",              required_argument, 0, 0 },
    {"idsuffix",              required_argument, 0, 0 },
    {"minqt",                 required_argument, 0, 0 },
    {"maxqt",                 required_argument, 0, 0 },
    {"minsl",                 required_argument, 0, 0 },
    {"maxsl",                 required_argument, 0, 0 },
    {"leftjust",              no_argument,       0, 0 },
    {"rightjust",             no_argument,       0, 0 },
    {"selfid",                no_argument,       0, 0 },
    {"maxid",                 required_argument, 0, 0 },
    {"minsizeratio",          required_argument, 0, 0 },
    {"maxsizeratio",          required_argument, 0, 0 },
    {"maxdiffs",              required_argument, 0, 0 },
    {"maxsubs",               required_argument, 0, 0 },
    {"maxgaps",               required_argument, 0, 0 },
    {"mincols",               required_argument, 0, 0 },
    {"maxqsize",              required_argument, 0, 0 },
    {"mintsize",              required_argument, 0, 0 },
    {"mid",                   required_argument, 0, 0 },
    {"shuffle",               required_argument, 0, 0 },
    {"randseed",              required_argument, 0, 0 },
    {"maskfasta",             required_argument, 0, 0 },
    {"hardmask",              no_argument,       0, 0 },
    {"qmask",                 required_argument, 0, 0 },
    {"dbmask",                required_argument, 0, 0 },
    {"cluster_smallmem",      required_argument, 0, 0 },
    {"cluster_fast",          required_argument, 0, 0 },
    {"centroids",             required_argument, 0, 0 },
    {"clusters",              required_argument, 0, 0 },
    {"consout",               required_argument, 0, 0 },
    {"cons_truncate",         no_argument,       0, 0 },
    {"msaout",                required_argument, 0, 0 },
    {"usersort",              no_argument,       0, 0 },
    {"xn",                    required_argument, 0, 0 },
    {"iddef",                 required_argument, 0, 0 },
    {"slots",                 required_argument, 0, 0 },
    {"pattern",               required_argument, 0, 0 },
    {"maxuniquesize",         required_argument, 0, 0 },
    {"abskew",                required_argument, 0, 0 },
    {"chimeras",              required_argument, 0, 0 },
    {"dn",                    required_argument, 0, 0 },
    {"mindiffs",              required_argument, 0, 0 },
    {"mindiv",                required_argument, 0, 0 },
    {"minh",                  required_argument, 0, 0 },
    {"nonchimeras",           required_argument, 0, 0 },
    {"uchime_denovo",         required_argument, 0, 0 },
    {"uchime_ref",            required_argument, 0, 0 },
    {"uchimealns",            required_argument, 0, 0 },
    {"uchimeout",             required_argument, 0, 0 },
    {"uchimeout5",            no_argument,       0, 0 },
    {"alignwidth",            required_argument, 0, 0 },
    {"allpairs_global",       required_argument, 0, 0 },
    {"acceptall",             no_argument,       0, 0 },
    {"cluster_size",          required_argument, 0, 0 },
    {"samout",                required_argument, 0, 0 },
    {"log",                   required_argument, 0, 0 },
    {"quiet",                 no_argument,       0, 0 },
    {"fastx_subsample",       required_argument, 0, 0 },
    {"sample_pct",            required_argument, 0, 0 },
    {"fastq_chars",           required_argument, 0, 0 },
    {"profile",               required_argument, 0, 0 },
    {"sample_size",           required_argument, 0, 0 },
    {"fastaout",              required_argument, 0, 0 },
    {"xsize",                 no_argument,       0, 0 },
    {"clusterout_id",         no_argument,       0, 0 },
    {"clusterout_sort",       no_argument,       0, 0 },
    {"borderline",            required_argument, 0, 0 },
    {"relabel_sha1",          no_argument,       0, 0 },
    {"relabel_md5",           no_argument,       0, 0 },
    {"derep_prefix",          required_argument, 0, 0 },
    {"fastq_filter",          required_argument, 0, 0 },
    {"fastqout",              required_argument, 0, 0 },
    {"fastaout_discarded",    required_argument, 0, 0 },
    {"fastqout_discarded",    required_argument, 0, 0 },
    {"fastq_truncqual",       required_argument, 0, 0 },
    {"fastq_maxee",           required_argument, 0, 0 },
    {"fastq_trunclen",        required_argument, 0, 0 },
    {"fastq_minlen",          required_argument, 0, 0 },
    {"fastq_stripleft",       required_argument, 0, 0 },
    {"fastq_maxee_rate",      required_argument, 0, 0 },
    {"fastq_maxns",           required_argument, 0, 0 },
    {"eeout",                 no_argument,       0, 0 },
    {"fastq_ascii",           required_argument, 0, 0 },
    {"fastq_qmin",            required_argument, 0, 0 },
    {"fastq_qmax",            required_argument, 0, 0 },
    {"fastq_qmaxout",         required_argument, 0, 0 },
    {"fastq_stats",           required_argument, 0, 0 },
    {"fastq_tail",            required_argument, 0, 0 },
    {"fastx_revcomp",         required_argument, 0, 0 },
    {"label_suffix",          required_argument, 0, 0 },
    {"h",                     no_argument,       0, 0 },
    {"samheader",             no_argument,       0, 0 },
    {"sizeorder",             no_argument,       0, 0 },
    {"minwordmatches",        required_argument, 0, 0 },
    {"v",                     no_argument,       0, 0 },
    {"relabel_keep",          no_argument,       0, 0 },
    {"search_exact",          required_argument, 0, 0 },
    {"fastx_mask",            required_argument, 0, 0 },
    {"min_unmasked_pct",      required_argument, 0, 0 },
    {"max_unmasked_pct",      required_argument, 0, 0 },
    {"fastq_convert",         required_argument, 0, 0 },
    {"fastq_asciiout",        required_argument, 0, 0 },
    {"fastq_qminout",         required_argument, 0, 0 },
    {"fastq_mergepairs",      required_argument, 0, 0 },
    {"fastq_eeout",           no_argument,       0, 0 },
    {"fastqout_notmerged_fwd",required_argument, 0, 0 },
    {"fastqout_notmerged_rev",required_argument, 0, 0 },
    {"fastq_minovlen",        required_argument, 0, 0 },
    {"fastq_minmergelen",     required_argument, 0, 0 },
    {"fastq_maxmergelen",     required_argument, 0, 0 },
    {"fastq_nostagger",       no_argument,       0, 0 },
    {"fastq_allowmergestagger", no_argument,     0, 0 },
    {"fastq_maxdiffs",        required_argument, 0, 0 },
    {"fastaout_notmerged_fwd",required_argument, 0, 0 },
    {"fastaout_notmerged_rev",required_argument, 0, 0 },
    {"reverse",               required_argument, 0, 0 },
    {"eetabbedout",           required_argument, 0, 0 },
    {"idoffset",              required_argument, 0, 0 },
    { 0, 0, 0, 0 }
  };
  
  int option_count = sizeof(long_options) / sizeof(struct option);
  bool options_selected[option_count];

  memset(options_selected, 0, sizeof(options_selected));

  int option_index = 0;
  int c;
  
  while ((c = getopt_long_only(argc, argv, "", long_options, 
                               &option_index)) == 0)
    {
      if (option_index < option_count)
        options_selected[option_index] = 1;

      switch(option_index)
        {
        case 0:
          opt_help = 1;
          break;
              
        case 1:
          opt_version = 1;
          break;

        case 2:
          opt_alnout = optarg;
          break;
          
        case 3:
          opt_usearch_global = optarg;
          break;

        case 4:
          opt_db = optarg;
          break;

        case 5:
          opt_id = args_getdouble(optarg);
          break;

        case 6:
          opt_maxaccepts = args_getlong(optarg);
          break;

        case 7:
          opt_maxrejects = args_getlong(optarg);
          break;

        case 8:
          opt_wordlength = args_getlong(optarg);
          break;

        case 9:
          opt_match = args_getlong(optarg);
          break;

        case 10:
          opt_mismatch = args_getlong(optarg);
          break;

        case 11:
          opt_fulldp = 1;
          break;

        case 12:
          if (strcasecmp(optarg, "plus") == 0)
            opt_strand = 1;
          else if (strcasecmp(optarg, "both") == 0)
            opt_strand = 2;
          else
            opt_strand = 0;
          break;

        case 13:
          opt_threads = (long) args_getdouble(optarg);
          break;

        case 14:
          args_get_gap_penalty_string(optarg, 1);
          break;

        case 15:
          args_get_gap_penalty_string(optarg, 0);
          break;

        case 16:
          opt_rowlen = args_getlong(optarg);
          break;

        case 17:
          if (!parse_userfields_arg(optarg))
            fatal("Unrecognized userfield argument");
          break;

        case 18:
          opt_userout = optarg;
          break;
      
        case 19:
          opt_self = 1;
          break;
      
        case 20:
          opt_blast6out = optarg;
          break;
      
        case 21:
          opt_uc = optarg;
          break;
      
        case 22:
          opt_weak_id = args_getdouble(optarg);
          break;

        case 23:
          opt_uc_allhits = 1;
          break;

        case 24:
          opt_notrunclabels = 1;
          break;

        case 25:
          opt_sortbysize = optarg;
          break;

        case 26:
          opt_output = optarg;
          break;

        case 27:
          opt_minsize = args_getlong(optarg);
          break;

        case 28:
          opt_maxsize = args_getlong(optarg);
          break;

        case 29:
          opt_relabel = optarg;
          break;

        case 30:
          opt_sizeout = 1;
          break;

        case 31:
          opt_derep_fulllength = optarg;
          break;

        case 32:
          opt_minseqlength = args_getlong(optarg);
          break;

        case 33:
          opt_minuniquesize = args_getlong(optarg);
          break;

        case 34:
          opt_topn = args_getlong(optarg);
          break;

        case 35:
          opt_maxseqlength = args_getlong(optarg);
          break;

        case 36:
          opt_sizein = 1;
          break;

        case 37:
          opt_sortbylength = optarg;
          break;

        case 38:
          opt_matched = optarg;
          break;

        case 39:
          opt_notmatched = optarg;
          break;

        case 40:
          opt_dbmatched = optarg;
          break;

        case 41:
          opt_dbnotmatched = optarg;
          break;

        case 42:
          opt_fastapairs = optarg;
          break;

        case 43:
          opt_output_no_hits = 1;
          break;

        case 44:
          opt_maxhits = args_getlong(optarg);
          break;

        case 45:
          opt_top_hits_only = 1;
          break;

        case 46:
          opt_fasta_width = args_getlong(optarg);
          break;

        case 47:
          opt_query_cov = args_getdouble(optarg);
          break;

        case 48:
          opt_target_cov = args_getdouble(optarg);
          break;

        case 49:
          opt_idprefix = args_getlong(optarg);
          break;

        case 50:
          opt_idsuffix = args_getlong(optarg);
          break;

        case 51:
          opt_minqt = args_getdouble(optarg);
          break;

        case 52:
          opt_maxqt = args_getdouble(optarg);
          break;

        case 53:
          opt_minsl = args_getdouble(optarg);
          break;

        case 54:
          opt_maxsl = args_getdouble(optarg);
          break;

        case 55:
          opt_leftjust = 1;
          break;

        case 56:
          opt_rightjust = 1;
          break;

        case 57:
          opt_selfid = 1;
          break;

        case 58:
          opt_maxid = args_getdouble(optarg);
          break;

        case 59:
          opt_minsizeratio = args_getdouble(optarg);
          break;

        case 60:
          opt_maxsizeratio = args_getdouble(optarg);
          break;

        case 61:
          opt_maxdiffs = args_getlong(optarg);
          break;

        case 62:
          opt_maxsubs = args_getlong(optarg);
          break;

        case 63:
          opt_maxgaps = args_getlong(optarg);
          break;

        case 64:
          opt_mincols = args_getlong(optarg);
          break;

        case 65:
          opt_maxqsize = args_getlong(optarg);
          break;

        case 66:
          opt_mintsize = args_getlong(optarg);
          break;

        case 67:
          opt_mid = args_getdouble(optarg);
          break;

        case 68:
          opt_shuffle = optarg;
          break;

        case 69:
          opt_randseed = args_getlong(optarg);
          break;

        case 70:
          opt_maskfasta = optarg;
          break;

        case 71:
          opt_hardmask = 1;
          break;

        case 72:
          if (strcasecmp(optarg, "none") == 0)
            opt_qmask = MASK_NONE;
          else if (strcasecmp(optarg, "dust") == 0)
            opt_qmask = MASK_DUST;
          else if (strcasecmp(optarg, "soft") == 0)
            opt_qmask = MASK_SOFT;
          else
            opt_qmask = MASK_ERROR;
          break;

        case 73:
          if (strcasecmp(optarg, "none") == 0)
            opt_dbmask = MASK_NONE;
          else if (strcasecmp(optarg, "dust") == 0)
            opt_dbmask = MASK_DUST;
          else if (strcasecmp(optarg, "soft") == 0)
            opt_dbmask = MASK_SOFT;
          else
            opt_dbmask = MASK_ERROR;
          break;

        case 74:
          opt_cluster_smallmem = optarg;
          break;

        case 75:
          opt_cluster_fast = optarg;
          break;

        case 76:
          opt_centroids = optarg;
          break;

        case 77:
          opt_clusters = optarg;
          break;

        case 78:
          opt_consout = optarg;
          break;

        case 79:
          fprintf(stderr, "WARNING: Option --cons_truncate is ignored\n");
          opt_cons_truncate = 1;
          break;

        case 80:
          opt_msaout = optarg;
          break;

        case 81:
          opt_usersort = 1;
          break;

        case 82:
          opt_xn = args_getdouble(optarg);
          break;

        case 83:
          opt_iddef = args_getlong(optarg);
          break;

        case 84:
          fprintf(stderr, "WARNING: Option --slots is ignored\n");
          opt_slots = args_getlong(optarg);
          break;

        case 85:
          fprintf(stderr, "WARNING: Option --pattern is ignored\n");
          opt_pattern = optarg;
          break;

        case 86:
          opt_maxuniquesize = args_getlong(optarg);
          break;

        case 87:
          opt_abskew = args_getdouble(optarg);
          break;
          
        case 88:
          opt_chimeras = optarg;
          break;
          
        case 89:
          opt_dn = args_getdouble(optarg);
          break;
          
        case 90:
          opt_mindiffs = args_getlong(optarg);
          break;
          
        case 91:
          opt_mindiv = args_getdouble(optarg);
          break;
          
        case 92:
          opt_minh = args_getdouble(optarg);
          break;
          
        case 93:
          opt_nonchimeras = optarg;
          break;
          
        case 94:
          opt_uchime_denovo = optarg;
          break;
          
        case 95:
          opt_uchime_ref = optarg;
          break;
          
        case 96:
          opt_uchimealns = optarg;
          break;
          
        case 97:
          opt_uchimeout = optarg;
          break;
          
        case 98:
          opt_uchimeout5 = 1;
          break;
          
        case 99:
          opt_alignwidth = args_getlong(optarg);
          break;
          
        case 100:
          opt_allpairs_global = optarg;
          break;
          
        case 101:
          opt_acceptall = 1;
          break;
          
        case 102:
          opt_cluster_size = optarg;
          break;

        case 103:
          opt_samout = optarg;
          break;

        case 104:
          opt_log = optarg;
          break;

        case 105:
          opt_quiet = true;
          break;

        case 106:
          opt_fastx_subsample = optarg;
          break;

        case 107:
          opt_sample_pct = args_getdouble(optarg);
          break;

        case 108:
          opt_fastq_chars = optarg;
          break;

        case 109:
          opt_profile = optarg;
          break;

        case 110:
          opt_sample_size = args_getlong(optarg);
          break;
          
        case 111:
          opt_fastaout = optarg;
          break;

        case 112:
          opt_xsize = 1;
          break;

        case 113:
          opt_clusterout_id = 1;
          break;

        case 114:
          opt_clusterout_sort = 1;
          break;

        case 115:
          opt_borderline = optarg;
          break;

        case 116:
          opt_relabel_sha1 = 1;
          break;

        case 117:
          opt_relabel_md5 = 1;
          break;

        case 118:
          opt_derep_prefix = optarg;
          break;

        case 119:
          opt_fastq_filter = optarg;
          break;

        case 120:
          opt_fastqout = optarg;
          break;

        case 121:
          opt_fastaout_discarded = optarg;
          break;

        case 122:
          opt_fastqout_discarded = optarg;
          break;

        case 123:
          opt_fastq_truncqual = args_getlong(optarg);
          break;

        case 124:
          opt_fastq_maxee = args_getdouble(optarg);
          break;

        case 125:
          opt_fastq_trunclen = args_getlong(optarg);
          break;

        case 126:
          opt_fastq_minlen = args_getlong(optarg);
          break;

        case 127:
          opt_fastq_stripleft = args_getlong(optarg);
          break;

        case 128:
          opt_fastq_maxee_rate = args_getdouble(optarg);
          break;

        case 129:
          opt_fastq_maxns = args_getlong(optarg);
          break;

        case 130:
          opt_eeout = 1;
          break;

        case 131:
          opt_fastq_ascii = args_getlong(optarg);
          break;

        case 132:
          opt_fastq_qmin = args_getlong(optarg);
          break;

        case 133:
          opt_fastq_qmax = args_getlong(optarg);
          break;

        case 134:
          opt_fastq_qmaxout = args_getlong(optarg);
          break;

        case 135:
          opt_fastq_stats = optarg;
          break;

        case 136:
          opt_fastq_tail = args_getlong(optarg);
          break;

        case 137:
          opt_fastx_revcomp = optarg;
          break;

        case 138:
          opt_label_suffix = optarg;
          break;

        case 139:
          opt_help = 1;
          break;

        case 140:
          opt_samheader = 1;
          break;

        case 141:
          opt_sizeorder = 1;
          break;

        case 142:
          opt_minwordmatches = args_getlong(optarg);
          break;

        case 143:
          opt_version = 1;
          break;

        case 144:
          opt_relabel_keep = 1;
          break;

        case 145:
          opt_search_exact = optarg;
          break;

        case 146:
          opt_fastx_mask = optarg;
          break;

        case 147:
          opt_min_unmasked_pct = args_getdouble(optarg);
          break;

        case 148:
          opt_max_unmasked_pct = args_getdouble(optarg);
          break;

        case 149:
          opt_fastq_convert = optarg;
          break;

        case 150:
          opt_fastq_asciiout = args_getlong(optarg);
          break;

        case 151:
          opt_fastq_qminout = args_getlong(optarg);
          break;

        case 152:
          opt_fastq_mergepairs = optarg;
          break;

        case 153:
          opt_fastq_eeout = 1;
          break;

        case 154:
          opt_fastqout_notmerged_fwd = optarg;
          break;

        case 155:
          opt_fastqout_notmerged_rev = optarg;
          break;

        case 156:
          opt_fastq_minovlen = args_getlong(optarg);
          break;

        case 157:
          opt_fastq_minmergelen = args_getlong(optarg);
          break;

        case 158:
          opt_fastq_maxmergelen = args_getlong(optarg);
          break;

        case 159:
          opt_fastq_nostagger = optarg;
          break;

        case 160:
          opt_fastq_allowmergestagger = 1;
          break;

        case 161:
          opt_fastq_maxdiffs = args_getlong(optarg);
          break;

        case 162:
          opt_fastaout_notmerged_fwd = optarg;
          break;

        case 163:
          opt_fastaout_notmerged_rev = optarg;
          break;

        case 164:
          opt_reverse = optarg;
          break;

        case 165:
          opt_eetabbedout = optarg;
          break;

        case 166:
          opt_idoffset = args_getlong(optarg);
          break;

        default:
          fatal("Internal error in option parsing");
        }
    }

  /* Terminate if ambiguous or illegal options have been detected */
  if (c != -1)
    exit(EXIT_FAILURE);

  /* Terminate after reporting any extra non-option arguments */
  if (optind < argc)
    fatal("Unrecognized string on command line (%s)", argv[optind]);

  int commands = 0;

  if (opt_fastq_chars)
    commands++;
  if (opt_fastq_filter)
    commands++;
  if (opt_fastq_stats)
    commands++;
  if (opt_usearch_global)
    commands++;
  if (opt_sortbysize)
    commands++;
  if (opt_sortbylength)
    commands++;
  if (opt_derep_fulllength)
    commands++;
  if (opt_derep_prefix)
    commands++;
  if (opt_help)
    commands++;
  if (opt_version)
    commands++;
  if (opt_shuffle)
    commands++;
  if (opt_fastx_subsample)
    commands++;
  if (opt_maskfasta)
    commands++;
  if (opt_cluster_smallmem)
    commands++;
  if (opt_cluster_fast)
    commands++;
  if (opt_cluster_size)
    commands++;
  if (opt_uchime_denovo)
    commands++;
  if (opt_uchime_ref)
    commands++;
  if (opt_allpairs_global)
    commands++;
  if (opt_fastx_revcomp)
    commands++;
  if (opt_search_exact)
    commands++;
  if (opt_fastx_mask)
    commands++;
  if (opt_fastq_convert)
    commands++;
  if (opt_fastq_mergepairs)
    commands++;
  
  if (commands > 1)
    fatal("More than one command specified");

  if (opt_weak_id > opt_id)
    opt_weak_id = opt_id;

  if (opt_maxrejects == -1)
    {
      if (opt_cluster_fast)
        opt_maxrejects = 8;
      else
        opt_maxrejects = 32;
    }

  if (opt_minseqlength < 0)
    fatal("The argument to --minseqlength must be positive");

  if (opt_maxaccepts < 0)
    fatal("The argument to --maxaccepts must not be negative");

  if (opt_maxrejects < 0)
    fatal("The argument to --maxrejects must not be negative");

  if ((opt_threads < 0) || (opt_threads > 1024))
    fatal("The argument to --threads must be in the range 0 (default) to 1024");

  if ((opt_wordlength < 7) || (opt_wordlength > 15))
    fatal("The argument to --wordlength must be in the range 7 to 15");

  if ((opt_iddef < 0) || (opt_iddef > 4))
    fatal("The argument to --iddef must in the range 0 to 4");

  if ((opt_idoffset < 0) || (opt_idoffset > 16))
    fatal("The argument to --idoffset must in the range 0 to 16");

  if (opt_match <= 0)
    fatal("The argument to --match must be positive");

  if (opt_mismatch >= 0)
    fatal("The argument to --mismatch must be negative");

  if (opt_alignwidth < 0)
    fatal("The argument to --alignwidth must not be negative");

  if (opt_rowlen < 0)
    fatal("The argument to --rowlen must not be negative");

  if (opt_strand < 1)
    fatal("The argument to --strand must be plus or both");
  
  if (opt_qmask == MASK_ERROR)
    fatal("The argument to --qmask must be none, dust or soft");
  
  if (opt_dbmask == MASK_ERROR)
    fatal("The argument to --dbmask must be none, dust or soft");

  if ((opt_sample_pct < 0.0) || (opt_sample_pct > 100.0))
    fatal("The argument to --sample_pct must be in the range 0.0 to 100.0");
  
  if (opt_sample_size < 0)
    fatal("The argument to --sample_size must not be negative");

  if (opt_relabel_sha1 && opt_relabel_md5)
    fatal("Specify either --relabel_sha1 or --relabel_md5, not both");

  if (opt_fastq_tail < 1)
    fatal("The argument to --fastq_tail must be positive");

  if (opt_minwordmatches < 0)
    fatal("The argument to --minwordmatches must not be negative");

  if ((opt_min_unmasked_pct < 0.0) && (opt_min_unmasked_pct > 100.0))
    fatal("The argument to --min_unmasked_pct must be between 0.0 and 100.0");

  if ((opt_max_unmasked_pct < 0.0) && (opt_max_unmasked_pct > 100.0))
    fatal("The argument to --max_unmasked_pct must be between 0.0 and 100.0");

  if (opt_min_unmasked_pct > opt_max_unmasked_pct)
    fatal("The argument to --min_unmasked_pct cannot be larger than to --max_unmasked_pct");

  if (opt_fastq_qmin > opt_fastq_qmax)
    fatal("The argument to --fastq_qmin cannot be larger than to --fastq_qmax");
  
  if (opt_fastq_qminout > opt_fastq_qmaxout)
    fatal("The argument to --fastq_qminout cannot be larger than to --fastq_qmaxout");
  
  /* TODO: check valid range of gap penalties */

  /* adapt/adjust parameters */

#if 1

  /* 
     Adjust gap open penalty according to convention.

     The specified gap open penalties include the penalty for
     a single nucleotide gap:

     gap penalty = gap open penalty + (gap length - 1) * gap extension penalty

     The rest of the code assumes the first nucleotide gap penalty is not
     included in the gap opening penalty.
  */

  opt_gap_open_query_left -= opt_gap_extension_query_left;
  opt_gap_open_target_left -= opt_gap_extension_target_left;
  opt_gap_open_query_interior -= opt_gap_extension_query_interior;
  opt_gap_open_target_interior -= opt_gap_extension_target_interior;
  opt_gap_open_query_right -= opt_gap_extension_query_right;
  opt_gap_open_target_right -= opt_gap_extension_target_right;

#endif

  if (opt_minwordmatches == 0)
    opt_minwordmatches = minwordmatches_defaults[opt_wordlength];

  if (opt_threads == 0)
    opt_threads = sysconf(_SC_NPROCESSORS_ONLN);

  /* set default opt_minseqlength depending on command */

  if (opt_minseqlength == 0)
    {
      if (opt_cluster_smallmem || opt_cluster_fast || opt_cluster_size ||
          opt_usearch_global || opt_derep_fulllength || opt_derep_prefix )
        opt_minseqlength = 32;
      else
        opt_minseqlength = 1;
    }

  if (opt_idoffset >= opt_minseqlength)
    fatal("The argument to --idoffset must be smaller than to --minseqlength");
}


void cmd_help()
{
  /*       0         1         2         3         4         5         6         7          */
  /*       01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  if (! opt_quiet)
    {
      fprintf(stdout, 
              "Usage: %s [OPTIONS]\n", progname);
      
      fprintf(stdout, 
              "\n"
              "General options\n"
              "  --fasta_width INT           width of FASTA seq lines, 0 for no wrap (80)\n"
              "  --help | --h                display help information\n"
              "  --log FILENAME              write messages, timing and memory info to file\n"
              "  --maxseqlength INT          maximum sequence length (50000)\n"
              "  --minseqlength INT          min seq length (clust/derep/search: 32, other:1)\n"
              "  --notrunclabels             do not truncate labels at first space\n"
              "  --quiet                     output just warnings and fatal errors to stderr\n"
              "  --threads INT               number of threads to use, zero for all cores (0)\n"
              "  --version                   display version information\n"
              "\n"
              "Chimera detection\n"
              "  --uchime_denovo FILENAME    detect chimeras de novo\n"
              "  --uchime_ref FILENAME       detect chimeras using a reference database\n"
              "Options\n"
              "  --abskew REAL               min abundance ratio of parent vs chimera (2.0)\n"
              "  --alignwidth INT            width of alignment in uchimealn output (80)\n"
              "  --borderline FILENAME       output borderline chimeric sequences to file\n"
              "  --chimeras FILENAME         output chimeric sequences to file\n"
              "  --db FILENAME               reference database for --uchime_ref\n"
              "  --dn REAL                   'no' vote pseudo-count (1.4)\n"
              "  --mindiffs INT              minimum number of differences in segment (3)\n"
              "  --mindiv REAL               minimum divergence from closest parent (0.8)\n"
              "  --minh REAL                 minimum score (0.28)\n"
              "  --nonchimeras FILENAME      output non-chimeric sequences to file\n"
              "  --relabel STRING            relabel nonchimeras with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --self                      exclude identical labels for --uchime_ref\n"
              "  --selfid                    exclude identical sequences for --uchime_ref\n"
              "  --sizeout                   include abundance information when relabelling\n"
              "  --uchimealns FILENAME       output chimera alignments to file\n"
              "  --uchimeout FILENAME        output to chimera info to tab-separated file\n"
              "  --uchimeout5                make output compatible with uchime version 5\n"
              "  --xn REAL                   'no' vote weight (8.0)\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "Clustering\n"
              "  --cluster_fast FILENAME     cluster sequences after sorting by length\n"
              "  --cluster_size FILENAME     cluster sequences after sorting by abundance\n"
              "  --cluster_smallmem FILENAME cluster already sorted sequences (see -usersort)\n"
              "Options (most searching options also apply)\n"
              "  --centroids FILENAME        output centroid sequences to FASTA file\n"
              "  --clusterout_id             add cluster id info to consout and profile files\n"
              "  --clusterout_sort           order msaout, consout, profile by decr abundance\n"
              "  --clusters STRING           output each cluster to a separate FASTA file\n"
              "  --consout FILENAME          output cluster consensus sequences to FASTA file\n"
              "  --cons_truncate             do not ignore terminal gaps in MSA for consensus\n"
              "  --id REAL                   reject if identity lower\n"
              "  --iddef INT                 id definition, 0-4=CD-HIT,all,int,MBL,BLAST (2)\n"
              "  --idoffset INT              id offset (0)\n"
              "  --msaout FILENAME           output multiple seq. alignments to FASTA file\n"
              "  --profile FILENAME          output sequence profile of each cluster to file\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              "  --relabel STRING            relabel centroids with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --sizeorder                 sort accepted centroids by abundance (AGC)\n"
              "  --sizeout                   write cluster abundances to centroid file\n"
              "  --strand plus|both          cluster using plus or both strands (plus)\n"
              "  --uc FILENAME               specify filename for UCLUST-like output\n"
              "  --usersort                  indicate sequences not pre-sorted by length\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "Dereplication\n"
              "  --derep_fulllength FILENAME dereplicate sequences in the given FASTA file\n"
              "  --derep_prefix FILENAME     dereplicate sequences in file based on prefixes\n"
              "Options\n"
              "  --maxuniquesize INT         maximum abundance for output from dereplication\n"
              "  --minuniquesize INT         minimum abundance for output from dereplication\n"
              "  --output FILENAME           output FASTA file\n"
              "  --relabel STRING            relabel with this prefix string after derep.\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --sizeout                   write abundance annotation to output\n"
              "  --strand plus|both          dereplicate plus or both strands (plus)\n"
              "  --topn INT                  output just the n most abundant sequences\n"
              "  --uc FILENAME               filename for UCLUST-like output\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "FASTQ filtering\n"
              "  --fastq_filter FILENAME     filter FASTQ file, output to FASTQ or FASTA file\n"
              "Options\n"
              "  --eeout                     include expected errors in FASTQ filter output\n"
              "  --fastaout FILENAME         FASTA output filename for passed sequences\n"
              "  --fastaout_discarded FNAME  FASTA filename for discarded sequences\n"
              "  --fastqout FILENAME         FASTQ output filename for passed sequences\n"
              "  --fastqout_discarded FNAME  FASTQ filename for discarded sequences\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_maxee REAL          maximum expected error value for FASTQ filter\n"
              "  --fastq_maxee_rate REAL     maximum expected error rate for FASTQ filter\n"
              "  --fastq_maxns INT           maximum number of N's for FASTQ filter\n"
              "  --fastq_minlen INT          minimum length for FASTQ filter\n"
              "  --fastq_stripleft INT       bases on the left to delete for FASTQ filter\n"
              "  --fastq_trunclen INT        read length for FASTQ filter truncation\n"
              "  --fastq_truncqual INT       base quality value for FASTQ filter truncation\n"
              "  --relabel STRING            relabel filtered sequences with given prefix\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel filtered sequences with md5 digest\n"
              "  --relabel_sha1              relabel filtered sequences with sha1 digest\n"
              "  --sizeout                   include abundance information when relabelling\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "FASTQ format conversion\n"
              "  --fastq_convert FILENAME    convert between FASTQ file formats\n"
              "Options\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_asciiout INT        FASTQ output quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmaxout INT         maximum base quality value for FASTQ output (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastq_qminout INT         minimum base quality value for FASTQ output (0)\n"
              "\n"
              "FASTQ format detection and quality analysis\n"
              "  --fastq_chars FILENAME      analyse FASTQ file for version and quality range\n"
              "Options\n"
              "  --fastq_tail INT            min length of tails to count for fastq_chars (4)\n"
              "\n"
              "FASTQ paired-end reads merging\n"
              "  --fastq_mergepairs FILENAME merge paired-end reads into one sequence\n"
              "Options:\n"
              "  --eetabbedout FILENAME      output error statistics to specified file\n"
              "  --fastaout FILENAME         FASTA output filename for merged sequences\n"
              "  --fastaout_notmerged_fwd FN FASTA filename for non-merged forward sequences\n"
              "  --fastaout_notmerged_rev FN FASTA filename for non-merged reverse sequences\n"
              "  --fastq_allowmergestagger   Allow merging of staggered reads\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_eeout               include expected errors in FASTQ output\n"
              "  --fastq_maxdiffs            maximum number of different bases in overlap\n"
              "  --fastq_maxee REAL          maximum expected error value for merged sequence\n"
              "  --fastq_maxmergelen         maximum length of entire merged sequence\n"
              "  --fastq_maxns INT           maximum number of N's\n"
              "  --fastq_minlen INT          minimum input read length after truncation (1)\n"
              "  --fastq_minmergelen         minimum length of entire merged sequence\n"
              "  --fastq_minovlen            minimum length of overlap between reads\n"
              "  --fastq_nostagger           disallow merging of staggered reads (default)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmaxout INT         maximum base quality value for FASTQ output (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastq_qminout INT         minimum base quality value for FASTQ output (0)\n"
              "  --fastq_truncqual INT       base quality value for truncation\n"
              "  --fastqout FILENAME         FASTQ output filename for merged sequences\n"
              "  --fastqout_notmerged_fwd  F FASTQ filename for non-merged forward sequences\n"
              "  --fastqout_notmerged_rev  F FASTQ filename for non-merged reverse sequences\n"
              "  --label_suffix              suffix to append to label of merged sequences\n"
              "  --reverse FILENAME          specify FASTQ file with reverse reads\n"
              "\n"
              "FASTQ quality statistics\n"
              "  --fastq_stats FILENAME      report FASTQ file statistics\n"
              "Options\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "\n"
              "Masking (new)\n"
              "  --fastx_mask FILENAME       mask sequences in the given FASTA or FASTQ file\n"
              "Options\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastaout FILENAME         output to specified FASTA file\n"
              "  --fastqout FILENAME         output to specified FASTQ file\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --max_unmasked_pct          max unmasked %% of sequences to keep (100.0)\n"
              "  --min_unmasked_pct          min unmasked %% of sequences to keep (0.0)\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              "\n"
              "Masking (old)\n"
              "  --maskfasta FILENAME        mask sequences in the given FASTA file\n"
              "Options\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --output FILENAME           output to specified FASTA file\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              "\n"
              "Pairwise alignment\n"
              "  --allpairs_global FILENAME  perform global alignment of all sequence pairs\n"
              "Options (most searching options also apply)\n"
              "  --alnout FILENAME           filename for human-readable alignment output\n"
              "  --acceptall                 output all pairwise alignments\n"
              "\n"
              "Reverse complementation\n"
              "  --fastx_revcomp FILENAME    Reverse-complement seqs in FASTA or FASTQ file\n"
              "Options\n"
              "  --fastaout FILENAME         FASTA output filename\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastqout FILENAME         FASTQ output filename\n"
              "  --label_suffix STRING       Label to append to identifier in the output\n"
              "\n"
              "Searching\n"
              "  --search_exact FILENAME     filename of queries for exact match search\n"
              "  --usearch_global FILENAME   filename of queries for global alignment search\n"
              "Options\n"
              "  --alnout FILENAME           filename for human-readable alignment output\n"
              "  --blast6out FILENAME        filename for blast-like tab-separated output\n"
              "  --db FILENAME               filename for FASTA formatted database for search\n"
              "  --dbmask none|dust|soft     mask db with dust, soft or no method (dust)\n"
              "  --dbmatched FILENAME        FASTA file for matching database sequences\n"
              "  --dbnotmatched FILENAME     FASTA file for non-matching database sequences\n"
              "  --fastapairs FILENAME       FASTA file with pairs of query and target\n"
              "  --fulldp                    full dynamic programming alignment (always on)\n"
              "  --gapext STRING             penalties for gap extension (2I/1E)\n"
              "  --gapopen STRING            penalties for gap opening (20I/2E)\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --id REAL                   reject if identity lower\n"
              "  --iddef INT                 id definition, 0-4=CD-HIT,all,int,MBL,BLAST (2)\n"
              "  --idoffset INT              id offset (0)\n"
              "  --idprefix INT              reject if first n nucleotides do not match\n"
              "  --idsuffix INT              reject if last n nucleotides do not match\n"
              "  --leftjust                  reject if terminal gaps at alignment left end\n"
              "  --match INT                 score for match (2)\n"
              "  --matched FILENAME          FASTA file for matching query sequences\n"
              "  --maxaccepts INT            number of hits to accept and show per strand (1)\n"
              "  --maxdiffs INT              reject if more substitutions or indels\n"
              "  --maxgaps INT               reject if more indels\n"
              "  --maxhits INT               maximum number of hits to show (unlimited)\n"
              "  --maxid REAL                reject if identity higher\n"
              "  --maxqsize INT              reject if query abundance larger\n"
              "  --maxqt REAL                reject if query/target length ratio higher\n"
              "  --maxrejects INT            number of non-matching hits to consider (32)\n"
              "  --maxsizeratio REAL         reject if query/target abundance ratio higher\n"
              "  --maxsl REAL                reject if shorter/longer length ratio higher\n"
              "  --maxsubs INT               reject if more substitutions\n"
              "  --mid REAL                  reject if percent identity lower, ignoring gaps\n"
              "  --mincols INT               reject if alignment length shorter\n"
              "  --minqt REAL                reject if query/target length ratio lower\n"
              "  --minsizeratio REAL         reject if query/target abundance ratio lower\n"
              "  --minsl REAL                reject if shorter/longer length ratio lower\n"
              "  --mintsize INT              reject if target abundance lower\n"
              "  --minwordmatches INT        minimum number of word matches required (10)\n"
              "  --mismatch INT              score for mismatch (-4)\n"
              "  --notmatched FILENAME       FASTA file for non-matching query sequences\n"
              "  --output_no_hits            output non-matching queries to output files\n"
              "  --pattern STRING            option is ignored\n"
              "  --qmask none|dust|soft      mask query with dust, soft or no method (dust)\n"
              "  --query_cov REAL            reject if fraction of query seq. aligned lower\n"
              "  --rightjust                 reject if terminal gaps at alignment right end\n"
              "  --rowlen INT                width of alignment lines in alnout output (64)\n"
              "  --samheader                 include a header in the SAM output file\n"
              "  --samout FILENAME           filename for SAM format output\n"
              "  --self                      reject if labels identical\n"
              "  --selfid                    reject if sequences identical\n"
              "  --sizeout                   write abundance annotation to dbmatched file\n"
              "  --slots INT                 option is ignored\n"
              "  --strand plus|both          search plus or both strands (plus)\n"
              "  --target_cov REAL           reject if fraction of target seq. aligned lower\n"
              "  --top_hits_only             output only hits with identity equal to the best\n"
              "  --uc FILENAME               filename for UCLUST-like output\n"
              "  --uc_allhits                show all, not just top hit with uc output\n"
              "  --userfields STRING         fields to output in userout file\n"
              "  --userout FILENAME          filename for user-defined tab-separated output\n"
              "  --weak_id REAL              include aligned hits with >= id; continue search\n"
              "  --wordlength INT            length of words for database index 3-15 (8)\n"
              "\n"
              "Shuffling and sorting\n"
              "  --shuffle FILENAME          shuffle order of sequences in FASTA file randomly\n"
              "  --sortbylength FILENAME     sort sequences by length in given FASTA file\n"
              "  --sortbysize FILENAME       abundance sort sequences in given FASTA file\n"
              "Options\n"
              "  --maxsize INT               maximum abundance for sortbysize\n"
              "  --minsize INT               minimum abundance for sortbysize\n"
              "  --output FILENAME           output to specified FASTA file\n"
              "  --randseed INT              seed for PRNG, zero to use random data source (0)\n"
              "  --relabel STRING            relabel sequences with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizeout                   include abundance information when relabelling\n"
              "  --topn INT                  output just first n sequences\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "Subsampling\n"
              "  --fastx_subsample FILENAME  subsample sequences from given FASTA/FASTQ file\n"
              "Options\n"
              "  --fastaout FILENAME         output FASTA file for subsamples\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastqout FILENAME         output FASTQ file for subsamples\n"
              "  --randseed INT              seed for PRNG, zero to use random data source (0)\n"
              "  --relabel STRING            relabel sequences with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sample_pct REAL           sampling percentage between 0.0 and 100.0\n"
              "  --sample_size INT           sampling size\n"
              "  --sizein                    consider abundance info from input, do not ignore\n"
              "  --sizeout                   update abundance information in output\n"
              "  --xsize                     strip abundance information in output\n"
          );
    }
}

void cmd_allpairs_global()
{
  /* check options */

  if ((!opt_alnout) && (!opt_userout) &&
      (!opt_uc) && (!opt_blast6out) &&
      (!opt_matched) && (!opt_notmatched) &&
      (!opt_samout))
    fatal("No output files specified");
  
  if (! (opt_acceptall || ((opt_id >= 0.0) && (opt_id <= 1.0)))) 
    fatal("Specify either --acceptall or --id with an identity from 0.0 to 1.0");

  allpairs_global(cmdline, progheader);
}

void cmd_usearch_global()
{
  /* check options */

  if ((!opt_alnout) && (!opt_userout) &&
      (!opt_uc) && (!opt_blast6out) &&
      (!opt_matched) && (!opt_notmatched) &&
      (!opt_dbmatched) && (!opt_dbnotmatched) &&
      (!opt_samout))
    fatal("No output files specified");
  
  if (!opt_db)
    fatal("Database filename not specified with --db");
  
  if ((opt_id < 0.0) || (opt_id > 1.0))
    fatal("Identity between 0.0 and 1.0 must be specified with --id");

  usearch_global(cmdline, progheader);
}

void cmd_search_exact()
{
  /* check options */

  if ((!opt_alnout) && (!opt_userout) &&
      (!opt_uc) && (!opt_blast6out) &&
      (!opt_matched) && (!opt_notmatched) &&
      (!opt_dbmatched) && (!opt_dbnotmatched) &&
      (!opt_samout))
    fatal("No output files specified");

  if (!opt_db)
    fatal("Database filename not specified with --db");

  search_exact(cmdline, progheader);
}

void cmd_sortbysize()
{
  if (!opt_output)
    fatal("FASTA output file for sortbysize must be specified with --output");

  sortbysize();
}

void cmd_sortbylength()
{
  if (!opt_output)
    fatal("FASTA output file for sortbylength must be specified with --output");

  sortbylength();
}

void cmd_derep()
{
  if ((!opt_output) && (!opt_uc))
    fatal("Output file for derepl_fulllength must be specified with --output or --uc");
  
  if (opt_derep_fulllength)
    derep_fulllength();
  else
    {
      if (opt_strand > 1)
        fatal("Option '--strand both' not supported with --derep_prefix");
      else
        derep_prefix();
    }
}

void cmd_shuffle()
{
  if (!opt_output)
    fatal("Output file for shuffling must be specified with --output");
  
  shuffle();
}

void cmd_subsample()
{
  if ((!opt_fastaout) && (!opt_fastqout))
    fatal("Specifiy output files for subsampling with --fastaout and/or --fastqout");

  if ((opt_sample_pct > 0) == (opt_sample_size > 0))
    fatal("Specify either --sample_pct or --sample_size, not both");

  subsample();
}

void cmd_maskfasta()
{
  if (!opt_output)
    fatal("Output file for masking must be specified with --output");
  
  maskfasta();
}

void cmd_fastx_mask()
{
  if ((!opt_fastaout) && (!opt_fastqout))
    fatal("Specifiy output files for masking with --fastaout and/or --fastqout");

  fastx_mask();
}

void cmd_none()
{
  if (! opt_quiet)
    fprintf(stderr,
            "For help, please enter: %s --help\n"
            "\n"
            "For further details, please see the manual by entering: man vsearch\n"
            "\n"
            "Example commands:\n"
            "\n"
            "vsearch --allpairs_global FILENAME --id 0.5 --alnout FILENAME\n"
            "vsearch --cluster_fast FILENAME --id 0.97 --centroids FILENAME\n"
            "vsearch --cluster_size FILENAME --id 0.97 --centroids FILENAME\n"
            "vsearch --cluster_smallmem FILENAME --usersort --id 0.97 --centroids FILENAME\n"
            "vsearch --derep_fulllength FILENAME --output FILENAME\n"
            "vsearch --derep_prefix FILENAME --output FILENAME\n"
            "vsearch --fastq_chars FILENAME\n"
            "vsearch --fastq_convert FILENAME --fastqout FILENAME --fastq_ascii 64\n"
            "vsearch --fastq_filter FILENAME --fastqout FILENAME --fastq_truncqual 20\n"
            "vsearch --fastq_mergepairs FILENAME --reverse FILENAME --fastqout FILENAME\n"
            "vsearch --fastq_stats FILENAME --log FILENAME\n"
            "vsearch --fastx_mask FILENAME --fastaout FILENAME\n"
            "vsearch --fastx_revcomp FILENAME --fastqout FILENAME\n"
            "vsearch --fastx_subsample FILENAME --fastaout FILENAME --sample_pct 1\n"
            "vsearch --search_exact FILENAME --db FILENAME --alnout FILENAME\n"
            "vsearch --shuffle FILENAME --output FILENAME\n"
            "vsearch --sortbylength FILENAME --output FILENAME\n"
            "vsearch --sortbysize FILENAME --output FILENAME\n"
            "vsearch --uchime_denovo FILENAME --nonchimeras FILENAME\n"
            "vsearch --uchime_ref FILENAME --db FILENAME --nonchimeras FILENAME\n"
            "vsearch --usearch_global FILENAME --db FILENAME --id 0.97 --alnout FILENAME\n"
            "\n",
            progname);
}

void cmd_fastx_revcomp()
{
  if ((!opt_fastaout) && (!opt_fastqout))
    fatal("No output files specified");
  
  fastx_revcomp();
}

void cmd_fastq_convert()
{
  if (! opt_fastqout)
    fatal("No output file specified with --fastqout");
  
  fastq_convert();
}

void cmd_cluster()
{
  if ((!opt_alnout) && (!opt_userout) &&
      (!opt_uc) && (!opt_blast6out) &&
      (!opt_matched) && (!opt_notmatched) &&
      (!opt_centroids) && (!opt_clusters) &&
      (!opt_consout) && (!opt_msaout) &&
      (!opt_samout) && (!opt_profile))
    fatal("No output files specified");
  
  if ((opt_id < 0.0) || (opt_id > 1.0))
    fatal("Identity between 0.0 and 1.0 must be specified with --id");

  if (opt_cluster_fast)
    cluster_fast(cmdline, progheader);
  else if (opt_cluster_smallmem)
    cluster_smallmem(cmdline, progheader);
  else if (opt_cluster_size)
    cluster_size(cmdline, progheader);
}

void cmd_uchime()
{
  if ((!opt_chimeras)  && (!opt_nonchimeras) &&
      (!opt_uchimeout) && (!opt_uchimealns))
    fatal("No output files specified");

  if (opt_uchime_ref && ! opt_db)
    fatal("Database filename not specified with --db");

  if (opt_xn <= 1.0)
    fatal("Argument to --xn must be > 1");

  if (opt_dn <= 0.0)
    fatal("Argument to --dn must be > 0");
  
  if (opt_mindiffs <= 0)
    fatal("Argument to --mindiffs must be > 0");

  if (opt_mindiv <= 0.0)
    fatal("Argument to --mindiv must be > 0");

  if (opt_minh <= 0.0)
    fatal("Argument to --minh must be > 0");

  if (opt_abskew <= 1.0)
    fatal("Argument to --abskew must be > 1");

  chimera();
}

void cmd_fastq_filter()
{
  if ((!opt_fastqout) && (!opt_fastaout) &&
      (!opt_fastqout_discarded) && (!opt_fastaout_discarded))
    fatal("No output files specified");
  fastq_filter();
}

void cmd_fastq_mergepairs()
{
  if (!opt_reverse)
    fatal("No reverse reads file specified with --reverse");
  if ((!opt_fastqout) &&
      (!opt_fastaout) &&
      (!opt_fastqout_notmerged_fwd) &&
      (!opt_fastqout_notmerged_rev) &&
      (!opt_fastaout_notmerged_fwd) &&
      (!opt_fastaout_notmerged_rev) &&
      (!opt_eetabbedout))
    fatal("No output files specified");
  fastq_mergepairs();
}

void fillheader()
{
  snprintf(progheader, 80, 
           "%s v%s_%s, %.1fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           sysconf(_SC_NPROCESSORS_ONLN));
}

void getentirecommandline(int argc, char** argv)
{
  int len = 0;
  for (int i=0; i<argc; i++)
    len += strlen(argv[i]);

  cmdline = (char*) xmalloc(len+argc+1);
  cmdline[0] = 0;

  for (int i=0; i<argc; i++)
    {
      strcat(cmdline, argv[i]);
      strcat(cmdline, " ");
    }
}

void show_header()
{
  if (! opt_quiet)
    {
      fprintf(stdout, "%s\n", progheader);
      fprintf(stdout, "https://github.com/torognes/vsearch\n");
      fprintf(stdout, "\n");
    }
}

int main(int argc, char** argv)
{
  fillheader();
  getentirecommandline(argc, argv);

  cpu_features_detect();

  args_init(argc, argv);

  dynlibs_open();

  if (opt_log)
    {
      fp_log = fopen(opt_log, "w");
      if (!fp_log)
        fatal("Unable to open log file for writing");
      fprintf(fp_log, "%s\n", progheader);
      fprintf(fp_log, "%s\n", cmdline);

      char time_string[26];
      time_start = time(0);
      struct tm tm_start;
      localtime_r(& time_start, & tm_start);
      strftime(time_string, 26, "%c", & tm_start);
      fprintf(fp_log, "Started  %s", time_string);
    }

  show_header();

  if (!sse2_present)
    fatal("Sorry, this program requires a cpu with SSE2.");

  global_abundance = abundance_init();

  if (opt_help)
    cmd_help();
  else if (opt_allpairs_global)
    cmd_allpairs_global();
  else if (opt_usearch_global)
    cmd_usearch_global();
  else if (opt_sortbysize)
    cmd_sortbysize();
  else if (opt_sortbylength)
    cmd_sortbylength();
  else if (opt_derep_fulllength || opt_derep_prefix)
    cmd_derep();
  else if (opt_shuffle)
    cmd_shuffle();
  else if (opt_fastx_subsample)
    cmd_subsample();
  else if (opt_maskfasta)
    cmd_maskfasta();
  else if (opt_cluster_smallmem || opt_cluster_fast || opt_cluster_size)
    cmd_cluster();
  else if (opt_uchime_denovo || opt_uchime_ref)
    cmd_uchime();
  else if (opt_fastq_chars)
    fastq_chars();
  else if (opt_fastq_stats)
    fastq_stats();
  else if (opt_fastq_filter)
    cmd_fastq_filter();
  else if (opt_fastx_revcomp)
    cmd_fastx_revcomp();
  else if (opt_search_exact)
    cmd_search_exact();
  else if (opt_fastx_mask)
    cmd_fastx_mask();
  else if (opt_fastq_convert)
    cmd_fastq_convert();
  else if (opt_fastq_mergepairs)
    cmd_fastq_mergepairs();
  else if (opt_version)
    {
    }
  else
    cmd_none();
  
  if (opt_log)
    {

      time_finish = time(0);
      struct tm tm_finish;
      localtime_r(& time_finish, & tm_finish);
      char time_string[26];
      strftime(time_string, 26, "%c", & tm_finish);
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Finished %s", time_string);

      time_t time_diff = time_finish - time_start;
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Elapsed time %02lu:%02lu\n", 
              time_diff / 60, time_diff % 60);
      double maxmem = arch_get_memused() / 1048576.0;
      if (maxmem < 1024.0)
        fprintf(fp_log, "Max memory %.1lfMB\n", maxmem);
      else
        fprintf(fp_log, "Max memory %.1lfGB\n", maxmem/1024.0);
      fclose(fp_log);
    }

  free(cmdline);

  abundance_exit(global_abundance);

  dynlibs_close();
}
