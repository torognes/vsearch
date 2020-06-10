/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2020, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

bool opt_bzip2_decompress;
bool opt_clusterout_id;
bool opt_clusterout_sort;
bool opt_eeout;
bool opt_fasta_score;
bool opt_fastq_allowmergestagger;
bool opt_fastq_eeout;
bool opt_fastq_nostagger;
bool opt_gzip_decompress;
bool opt_label_substr_match;
bool opt_no_progress;
bool opt_quiet;
bool opt_relabel_keep;
bool opt_relabel_md5;
bool opt_relabel_self;
bool opt_relabel_sha1;
bool opt_samheader;
bool opt_sff_clip;
bool opt_sizeorder;
bool opt_xee;
bool opt_xsize;
char * opt_allpairs_global;
char * opt_alnout;
char * opt_biomout;
char * opt_blast6out;
char * opt_borderline;
char * opt_centroids;
char * opt_chimeras;
char * opt_cluster_fast;
char * opt_cluster_size;
char * opt_cluster_smallmem;
char * opt_cluster_unoise;
char * opt_clusters;
char * opt_consout;
char * opt_cut;
char * opt_cut_pattern;
char * opt_db;
char * opt_dbmatched;
char * opt_dbnotmatched;
char * opt_derep_fulllength;
char * opt_derep_prefix;
char * opt_eetabbedout;
char * opt_fastaout;
char * opt_fastaout_discarded;
char * opt_fastaout_discarded_rev;
char * opt_fastaout_notmerged_fwd;
char * opt_fastaout_notmerged_rev;
char * opt_fastaout_rev;
char * opt_fastapairs;
char * opt_fastq_chars;
char * opt_fastq_convert;
char * opt_fastq_eestats;
char * opt_fastq_eestats2;
char * opt_fastq_filter;
char * opt_fastq_join;
char * opt_fastq_mergepairs;
char * opt_fastq_stats;
char * opt_fastqout;
char * opt_fastqout_discarded;
char * opt_fastqout_discarded_rev;
char * opt_fastqout_notmerged_fwd;
char * opt_fastqout_notmerged_rev;
char * opt_fastqout_rev;
char * opt_fastx_filter;
char * opt_fastx_getseq;
char * opt_fastx_getseqs;
char * opt_fastx_getsubseq;
char * opt_fastx_mask;
char * opt_fastx_revcomp;
char * opt_fastx_subsample;
char * opt_join_padgap;
char * opt_join_padgapq;
char * opt_label;
char * opt_labels;
char * opt_label_suffix;
char * opt_label_word;
char * opt_label_words;
char * opt_label_field;
char * opt_log;
char * opt_makeudb_usearch;
char * opt_maskfasta;
char * opt_matched;
char * opt_mothur_shared_out;
char * opt_msaout;
char * opt_nonchimeras;
char * opt_notmatched;
char * opt_notmatchedfq;
char * opt_otutabout;
char * opt_output;
char * opt_pattern;
char * opt_profile;
char * opt_relabel;
char * opt_rereplicate;
char * opt_reverse;
char * opt_samout;
char * opt_search_exact;
char * opt_sff_convert;
char * opt_shuffle;
char * opt_sintax;
char * opt_sortbylength;
char * opt_sortbysize;
char * opt_tabbedout;
char * opt_udb2fasta;
char * opt_udbinfo;
char * opt_udbstats;
char * opt_uc;
char * opt_uchime_denovo;
char * opt_uchime2_denovo;
char * opt_uchime3_denovo;
char * opt_uchime_ref;
char * opt_uchimealns;
char * opt_uchimeout;
char * opt_usearch_global;
char * opt_userout;
double * opt_ee_cutoffs_values;
double opt_abskew;
double opt_dn;
double opt_fastq_maxdiffpct;
double opt_fastq_maxee;
double opt_fastq_maxee_rate;
double opt_fastq_truncee;
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
double opt_sintax_cutoff;
double opt_target_cov;
double opt_unoise_alpha;
double opt_weak_id;
double opt_xn;
int opt_acceptall;
int opt_alignwidth;
int opt_cons_truncate;
int opt_ee_cutoffs_count;
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
int opt_length_cutoffs_shortest;
int opt_length_cutoffs_longest;
int opt_length_cutoffs_increment;
int opt_mindiffs;
int opt_slots;
int opt_uchimeout5;
int opt_usersort;
int opt_version;
int64_t opt_dbmask;
int64_t opt_fasta_width;
int64_t opt_fastq_ascii;
int64_t opt_fastq_asciiout;
int64_t opt_fastq_maxdiffs;
int64_t opt_fastq_maxlen;
int64_t opt_fastq_maxmergelen;
int64_t opt_fastq_maxns;
int64_t opt_fastq_minlen;
int64_t opt_fastq_minmergelen;
int64_t opt_fastq_minovlen;
int64_t opt_fastq_qmax;
int64_t opt_fastq_qmaxout;
int64_t opt_fastq_qmin;
int64_t opt_fastq_qminout;
int64_t opt_fastq_stripleft;
int64_t opt_fastq_stripright;
int64_t opt_fastq_tail;
int64_t opt_fastq_trunclen;
int64_t opt_fastq_trunclen_keep;
int64_t opt_fastq_truncqual;
int64_t opt_fulldp;
int64_t opt_hardmask;
int64_t opt_iddef;
int64_t opt_idprefix;
int64_t opt_idsuffix;
int64_t opt_leftjust;
int64_t opt_match;
int64_t opt_maxaccepts;
int64_t opt_maxdiffs;
int64_t opt_maxgaps;
int64_t opt_maxhits;
int64_t opt_maxqsize;
int64_t opt_maxrejects;
int64_t opt_maxseqlength;
int64_t opt_maxsize;
int64_t opt_maxsubs;
int64_t opt_maxuniquesize;
int64_t opt_mincols;
int64_t opt_minseqlength;
int64_t opt_minsize;
int64_t opt_mintsize;
int64_t opt_minuniquesize;
int64_t opt_minwordmatches;
int64_t opt_mismatch;
int64_t opt_notrunclabels;
int64_t opt_output_no_hits;
int64_t opt_qmask;
int64_t opt_randseed;
int64_t opt_rightjust;
int64_t opt_rowlen;
int64_t opt_sample_size;
int64_t opt_self;
int64_t opt_selfid;
int64_t opt_sizein;
int64_t opt_sizeout;
int64_t opt_strand;
int64_t opt_subseq_start;
int64_t opt_subseq_end;
int64_t opt_threads;
int64_t opt_top_hits_only;
int64_t opt_topn;
int64_t opt_uc_allhits;
int64_t opt_wordlength;

/* Other variables */

/* cpu features available */

int64_t altivec_present = 0;
int64_t neon_present = 0;
int64_t mmx_present = 0;
int64_t sse_present = 0;
int64_t sse2_present = 0;
int64_t sse3_present = 0;
int64_t ssse3_present = 0;
int64_t sse41_present = 0;
int64_t sse42_present = 0;
int64_t popcnt_present = 0;
int64_t avx_present = 0;
int64_t avx2_present = 0;

static char * progname;
static char progheader[80];
static char * cmdline;
static time_t time_start;
static time_t time_finish;

FILE * fp_log = 0;

char * STDIN_NAME = (char*) "/dev/stdin";
char * STDOUT_NAME = (char*) "/dev/stdout";

#ifdef __x86_64__
#define cpuid(f1, f2, a, b, c, d)                                \
  __asm__ __volatile__ ("cpuid"                                  \
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d) \
                        : "a" (f1), "c" (f2));
#endif

void cpu_features_detect()
{
#ifdef __aarch64__
#ifdef __ARM_NEON
  /* may check /proc/cpuinfo for asimd or neon */
  neon_present = 1;
#else
#error ARM Neon not present
#endif
#elif __PPC__
  altivec_present = 1;
#elif __x86_64__
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
#else
#error Unknown architecture
#endif
}

void cpu_features_show()
{
  fprintf(stderr, "CPU features:");
  if (neon_present)
    fprintf(stderr, " neon");
  if (altivec_present)
    fprintf(stderr, " altivec");
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

void args_get_ee_cutoffs(char * arg)
{
  /* get comma-separated list of floating point numbers */
  /* save in ee_cutoffs_count and ee_cutoffs_values */

  int commas = 0;
  for (size_t i=0; i<strlen(arg); i++)
    if (arg[i] == ',')
      commas++;

  opt_ee_cutoffs_count = 0;
  opt_ee_cutoffs_values = (double*) xrealloc(opt_ee_cutoffs_values, (commas+1) * sizeof(double));

  char * s = arg;
  while(1)
    {
      double val = 0;
      int skip = 0;

      if ((sscanf(s, "%lf%n", &val, &skip) != 1) || (val <= 0.0))
        fatal("Invalid arguments to ee_cutoffs");

      opt_ee_cutoffs_values[opt_ee_cutoffs_count++] = val;

      s += skip;

      if (*s == ',')
        s++;
      else if (*s == 0)
        break;
      else
        fatal("Invalid arguments to ee_cutoffs");
    }
}

void args_get_length_cutoffs(char * arg)
{
  /* get comma-separated list of 3 integers: */
  /* smallest, largest and increment. */
  /* second value may be * indicating no limit */
  /* save in length_cutoffs_{smallest,largest,increment} */

  int skip = 0;
  if (sscanf(arg, "%d,%d,%d%n", &opt_length_cutoffs_shortest, &opt_length_cutoffs_longest, &opt_length_cutoffs_increment, & skip) == 3)
    {
      if ((size_t)skip < strlen(arg))
        fatal("Invalid arguments to length_cutoffs");
    }
  else if (sscanf(arg, "%d,*,%d%n", &opt_length_cutoffs_shortest, &opt_length_cutoffs_increment, &skip) == 2)
    {
      if ((size_t)skip < strlen(arg))
        fatal("Invalid arguments to length_cutoffs");
      opt_length_cutoffs_longest = INT_MAX;
    }
  else
    fatal("Invalid arguments to length_cutoffs");

  if ((opt_length_cutoffs_shortest < 1) ||
      (opt_length_cutoffs_shortest > opt_length_cutoffs_longest) ||
      (opt_length_cutoffs_increment < 1))
    fatal("Invalid arguments to length_cutoffs");

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


int64_t args_getlong(char * arg)
{
  int len = 0;
  int64_t temp = 0;
  int ret = sscanf(arg, "%" PRId64 "%n", &temp, &len);
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

  opt_abskew = -1.0;
  opt_acceptall = 0;
  opt_alignwidth = 80;
  opt_allpairs_global = 0;
  opt_alnout = 0;
  opt_blast6out = 0;
  opt_biomout = 0;
  opt_borderline = 0;
  opt_bzip2_decompress = 0;
  opt_centroids = 0;
  opt_chimeras = 0;
  opt_cluster_fast = 0;
  opt_cluster_size = 0;
  opt_cluster_smallmem = 0;
  opt_cluster_unoise = 0;
  opt_clusterout_id = 0;
  opt_clusterout_sort = 0;
  opt_clusters = 0;
  opt_cons_truncate = 0;
  opt_consout = 0;
  opt_cut = 0;
  opt_cut_pattern = 0;
  opt_db = 0;
  opt_dbmask = MASK_DUST;
  opt_dbmatched = 0;
  opt_dbnotmatched = 0;
  opt_derep_fulllength = 0;
  opt_derep_prefix = 0;
  opt_dn = 1.4;
  opt_ee_cutoffs_count = 3;
  opt_ee_cutoffs_values = (double*) xmalloc(opt_ee_cutoffs_count * sizeof(double));
  opt_ee_cutoffs_values[0] = 0.5;
  opt_ee_cutoffs_values[1] = 1.0;
  opt_ee_cutoffs_values[2] = 2.0;
  opt_eeout = 0;
  opt_eetabbedout = 0;
  opt_fastaout_notmerged_fwd = 0;
  opt_fastaout_notmerged_rev = 0;
  opt_fasta_score = 0;
  opt_fasta_width = 80;
  opt_fastaout = 0;
  opt_fastaout_discarded = 0;
  opt_fastaout_discarded_rev = 0;
  opt_fastaout_rev = 0;
  opt_fastapairs = 0;
  opt_fastq_allowmergestagger = 0;
  opt_fastq_ascii = 33;
  opt_fastq_asciiout = 33;
  opt_fastq_chars = 0;
  opt_fastq_convert = 0;
  opt_fastq_eeout = 0;
  opt_fastq_eestats = 0;
  opt_fastq_eestats2 = 0;
  opt_fastq_filter = 0;
  opt_fastq_join = 0;
  opt_fastq_maxdiffpct = 100.0;
  opt_fastq_maxdiffs = 10;
  opt_fastq_maxee = DBL_MAX;
  opt_fastq_maxee_rate = DBL_MAX;
  opt_fastq_maxlen = LONG_MAX;
  opt_fastq_maxmergelen  = 1000000;
  opt_fastq_maxns = LONG_MAX;
  opt_fastq_mergepairs = 0;
  opt_fastq_minlen = 1;
  opt_fastq_minmergelen = 0;
  opt_fastq_minovlen = 10;
  opt_fastq_nostagger = 1;
  opt_fastqout_notmerged_fwd = 0;
  opt_fastqout_notmerged_rev = 0;
  opt_fastq_qmax = 41;
  opt_fastq_qmaxout = 41;
  opt_fastq_qmin = 0;
  opt_fastq_qminout = 0;
  opt_fastq_stats = 0;
  opt_fastq_stripleft = 0;
  opt_fastq_stripright = 0;
  opt_fastq_tail = 4;
  opt_fastq_truncee = DBL_MAX;
  opt_fastq_trunclen = -1;
  opt_fastq_trunclen_keep = -1;
  opt_fastq_truncqual = LONG_MIN;
  opt_fastqout = 0;
  opt_fastqout_discarded = 0;
  opt_fastqout_discarded_rev = 0;
  opt_fastqout_rev = 0;
  opt_fastx_filter = 0;
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
  opt_fastx_getseq = 0;
  opt_fastx_getseqs = 0;
  opt_fastx_getsubseq = 0;
  opt_gzip_decompress = 0;
  opt_hardmask = 0;
  opt_help = 0;
  opt_id = -1.0;
  opt_iddef = 2;
  opt_idprefix = 0;
  opt_idsuffix = 0;
  opt_join_padgap = 0;
  opt_join_padgapq = 0;
  opt_label = 0;
  opt_label_substr_match = false;
  opt_label_suffix = 0;
  opt_labels = 0;
  opt_label_field = 0;
  opt_label_word = 0;
  opt_label_words = 0;
  opt_leftjust = 0;
  opt_length_cutoffs_increment = 50;
  opt_length_cutoffs_longest = INT_MAX;
  opt_length_cutoffs_shortest = 50;
  opt_log = 0;
  opt_makeudb_usearch = 0;
  opt_maskfasta = 0;
  opt_match = 2;
  opt_matched = 0;
  opt_max_unmasked_pct = 100.0;
  opt_maxaccepts = 1;
  opt_maxdiffs = INT_MAX;
  opt_maxgaps = INT_MAX;
  opt_maxhits = 0;
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
  opt_minseqlength = -1;
  opt_minsize = 0;
  opt_minsizeratio = 0.0;
  opt_minsl = 0.0;
  opt_mintsize = 0;
  opt_minuniquesize = 1;
  opt_minwordmatches = -1;
  opt_mismatch = -4;
  opt_mothur_shared_out = 0;
  opt_msaout = 0;
  opt_no_progress = 0;
  opt_nonchimeras = 0;
  opt_notmatched = 0;
  opt_notmatched = 0;
  opt_notrunclabels = 0;
  opt_otutabout = 0;
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
  opt_relabel_self = 0;
  opt_relabel_sha1 = 0;
  opt_rereplicate = 0;
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
  opt_sff_convert = 0;
  opt_sff_clip = 0;
  opt_shuffle = 0;
  opt_sintax = 0;
  opt_sintax_cutoff = 0.0;
  opt_sizein = 0;
  opt_sizeorder = 0;
  opt_sizeout = 0;
  opt_slots = 0;
  opt_sortbylength = 0;
  opt_sortbysize = 0;
  opt_strand = 1;
  opt_subseq_start = 1;
  opt_subseq_end = LONG_MAX;
  opt_tabbedout = 0;
  opt_target_cov = 0.0;
  opt_threads = 0;
  opt_top_hits_only = 0;
  opt_topn = LONG_MAX;
  opt_udb2fasta = 0;
  opt_udbinfo = 0;
  opt_udbstats = 0;
  opt_uc = 0;
  opt_uc_allhits = 0;
  opt_uchime_denovo = 0;
  opt_uchime2_denovo = 0;
  opt_uchime3_denovo = 0;
  opt_uchime_ref = 0;
  opt_uchimealns = 0;
  opt_uchimeout = 0;
  opt_uchimeout5 = 0;
  opt_unoise_alpha = 2.0;
  opt_usearch_global = 0;
  opt_userout = 0;
  opt_usersort = 0;
  opt_version = 0;
  opt_weak_id = 10.0;
  opt_wordlength = 8;
  opt_xn = 8.0;
  opt_xsize = 0;
  opt_xee = 0;

  opterr = 1;

  enum
  {
    option_abskew,
    option_acceptall,
    option_alignwidth,
    option_allpairs_global,
    option_alnout,
    option_band,
    option_biomout,
    option_blast6out,
    option_borderline,
    option_bzip2_decompress,
    option_centroids,
    option_chimeras,
    option_cluster_fast,
    option_cluster_size,
    option_cluster_smallmem,
    option_cluster_unoise,
    option_clusterout_id,
    option_clusterout_sort,
    option_clusters,
    option_cons_truncate,
    option_consout,
    option_cut,
    option_cut_pattern,
    option_db,
    option_dbmask,
    option_dbmatched,
    option_dbnotmatched,
    option_derep_fulllength,
    option_derep_prefix,
    option_dn,
    option_ee_cutoffs,
    option_eeout,
    option_eetabbedout,
    option_fasta_score,
    option_fasta_width,
    option_fastaout,
    option_fastaout_discarded,
    option_fastaout_discarded_rev,
    option_fastaout_notmerged_fwd,
    option_fastaout_notmerged_rev,
    option_fastaout_rev,
    option_fastapairs,
    option_fastq_allowmergestagger,
    option_fastq_ascii,
    option_fastq_asciiout,
    option_fastq_chars,
    option_fastq_convert,
    option_fastq_eeout,
    option_fastq_eestats,
    option_fastq_eestats2,
    option_fastq_filter,
    option_fastq_join,
    option_fastq_maxdiffpct,
    option_fastq_maxdiffs,
    option_fastq_maxee,
    option_fastq_maxee_rate,
    option_fastq_maxlen,
    option_fastq_maxmergelen,
    option_fastq_maxns,
    option_fastq_mergepairs,
    option_fastq_minlen,
    option_fastq_minmergelen,
    option_fastq_minovlen,
    option_fastq_nostagger,
    option_fastq_qmax,
    option_fastq_qmaxout,
    option_fastq_qmin,
    option_fastq_qminout,
    option_fastq_stats,
    option_fastq_stripleft,
    option_fastq_stripright,
    option_fastq_tail,
    option_fastq_truncee,
    option_fastq_trunclen,
    option_fastq_trunclen_keep,
    option_fastq_truncqual,
    option_fastqout,
    option_fastqout_discarded,
    option_fastqout_discarded_rev,
    option_fastqout_notmerged_fwd,
    option_fastqout_notmerged_rev,
    option_fastqout_rev,
    option_fastx_filter,
    option_fastx_getseq,
    option_fastx_getseqs,
    option_fastx_getsubseq,
    option_fastx_mask,
    option_fastx_revcomp,
    option_fastx_subsample,
    option_fulldp,
    option_gapext,
    option_gapopen,
    option_gzip_decompress,
    option_h,
    option_hardmask,
    option_help,
    option_hspw,
    option_id,
    option_iddef,
    option_idprefix,
    option_idsuffix,
    option_join_padgap,
    option_join_padgapq,
    option_label,
    option_label_field,
    option_label_substr_match,
    option_label_suffix,
    option_label_word,
    option_label_words,
    option_labels,
    option_leftjust,
    option_length_cutoffs,
    option_log,
    option_makeudb_usearch,
    option_maskfasta,
    option_match,
    option_matched,
    option_max_unmasked_pct,
    option_maxaccepts,
    option_maxdiffs,
    option_maxgaps,
    option_maxhits,
    option_maxid,
    option_maxqsize,
    option_maxqt,
    option_maxrejects,
    option_maxseqlength,
    option_maxsize,
    option_maxsizeratio,
    option_maxsl,
    option_maxsubs,
    option_maxuniquesize,
    option_mid,
    option_min_unmasked_pct,
    option_mincols,
    option_mindiffs,
    option_mindiv,
    option_minh,
    option_minhsp,
    option_minqt,
    option_minseqlength,
    option_minsize,
    option_minsizeratio,
    option_minsl,
    option_mintsize,
    option_minuniquesize,
    option_minwordmatches,
    option_mismatch,
    option_mothur_shared_out,
    option_msaout,
    option_no_progress,
    option_nonchimeras,
    option_notmatched,
    option_notmatchedfq,
    option_notrunclabels,
    option_otutabout,
    option_output,
    option_output_no_hits,
    option_pattern,
    option_profile,
    option_qmask,
    option_query_cov,
    option_quiet,
    option_randseed,
    option_relabel,
    option_relabel_keep,
    option_relabel_md5,
    option_relabel_self,
    option_relabel_sha1,
    option_rereplicate,
    option_reverse,
    option_rightjust,
    option_rowlen,
    option_samheader,
    option_samout,
    option_sample_pct,
    option_sample_size,
    option_search_exact,
    option_self,
    option_selfid,
    option_sff_clip,
    option_sff_convert,
    option_shuffle,
    option_sintax,
    option_sintax_cutoff,
    option_sizein,
    option_sizeorder,
    option_sizeout,
    option_slots,
    option_sortbylength,
    option_sortbysize,
    option_strand,
    option_subseq_end,
    option_subseq_start,
    option_tabbedout,
    option_target_cov,
    option_threads,
    option_top_hits_only,
    option_topn,
    option_uc,
    option_uc_allhits,
    option_uchime2_denovo,
    option_uchime3_denovo,
    option_uchime_denovo,
    option_uchime_ref,
    option_uchimealns,
    option_uchimeout,
    option_uchimeout5,
    option_udb2fasta,
    option_udbinfo,
    option_udbstats,
    option_unoise_alpha,
    option_usearch_global,
    option_userfields,
    option_userout,
    option_usersort,
    option_v,
    option_version,
    option_weak_id,
    option_wordlength,
    option_xdrop_nw,
    option_xee,
    option_xn,
    option_xsize
  };

  static struct option long_options[] =
  {
    {"abskew",                required_argument, 0, 0 },
    {"acceptall",             no_argument,       0, 0 },
    {"alignwidth",            required_argument, 0, 0 },
    {"allpairs_global",       required_argument, 0, 0 },
    {"alnout",                required_argument, 0, 0 },
    {"band",                  required_argument, 0, 0 },
    {"biomout",               required_argument, 0, 0 },
    {"blast6out",             required_argument, 0, 0 },
    {"borderline",            required_argument, 0, 0 },
    {"bzip2_decompress",      no_argument,       0, 0 },
    {"centroids",             required_argument, 0, 0 },
    {"chimeras",              required_argument, 0, 0 },
    {"cluster_fast",          required_argument, 0, 0 },
    {"cluster_size",          required_argument, 0, 0 },
    {"cluster_smallmem",      required_argument, 0, 0 },
    {"cluster_unoise",        required_argument, 0, 0 },
    {"clusterout_id",         no_argument,       0, 0 },
    {"clusterout_sort",       no_argument,       0, 0 },
    {"clusters",              required_argument, 0, 0 },
    {"cons_truncate",         no_argument,       0, 0 },
    {"consout",               required_argument, 0, 0 },
    {"cut",                   required_argument, 0, 0 },
    {"cut_pattern",           required_argument, 0, 0 },
    {"db",                    required_argument, 0, 0 },
    {"dbmask",                required_argument, 0, 0 },
    {"dbmatched",             required_argument, 0, 0 },
    {"dbnotmatched",          required_argument, 0, 0 },
    {"derep_fulllength",      required_argument, 0, 0 },
    {"derep_prefix",          required_argument, 0, 0 },
    {"dn",                    required_argument, 0, 0 },
    {"ee_cutoffs",            required_argument, 0, 0 },
    {"eeout",                 no_argument,       0, 0 },
    {"eetabbedout",           required_argument, 0, 0 },
    {"fasta_score",           no_argument,       0, 0 },
    {"fasta_width",           required_argument, 0, 0 },
    {"fastaout",              required_argument, 0, 0 },
    {"fastaout_discarded",    required_argument, 0, 0 },
    {"fastaout_discarded_rev",required_argument, 0, 0 },
    {"fastaout_notmerged_fwd",required_argument, 0, 0 },
    {"fastaout_notmerged_rev",required_argument, 0, 0 },
    {"fastaout_rev",          required_argument, 0, 0 },
    {"fastapairs",            required_argument, 0, 0 },
    {"fastq_allowmergestagger", no_argument,     0, 0 },
    {"fastq_ascii",           required_argument, 0, 0 },
    {"fastq_asciiout",        required_argument, 0, 0 },
    {"fastq_chars",           required_argument, 0, 0 },
    {"fastq_convert",         required_argument, 0, 0 },
    {"fastq_eeout",           no_argument,       0, 0 },
    {"fastq_eestats",         required_argument, 0, 0 },
    {"fastq_eestats2",        required_argument, 0, 0 },
    {"fastq_filter",          required_argument, 0, 0 },
    {"fastq_join",            required_argument, 0, 0 },
    {"fastq_maxdiffpct",      required_argument, 0, 0 },
    {"fastq_maxdiffs",        required_argument, 0, 0 },
    {"fastq_maxee",           required_argument, 0, 0 },
    {"fastq_maxee_rate",      required_argument, 0, 0 },
    {"fastq_maxlen",          required_argument, 0, 0 },
    {"fastq_maxmergelen",     required_argument, 0, 0 },
    {"fastq_maxns",           required_argument, 0, 0 },
    {"fastq_mergepairs",      required_argument, 0, 0 },
    {"fastq_minlen",          required_argument, 0, 0 },
    {"fastq_minmergelen",     required_argument, 0, 0 },
    {"fastq_minovlen",        required_argument, 0, 0 },
    {"fastq_nostagger",       no_argument,       0, 0 },
    {"fastq_qmax",            required_argument, 0, 0 },
    {"fastq_qmaxout",         required_argument, 0, 0 },
    {"fastq_qmin",            required_argument, 0, 0 },
    {"fastq_qminout",         required_argument, 0, 0 },
    {"fastq_stats",           required_argument, 0, 0 },
    {"fastq_stripleft",       required_argument, 0, 0 },
    {"fastq_stripright",      required_argument, 0, 0 },
    {"fastq_tail",            required_argument, 0, 0 },
    {"fastq_truncee",         required_argument, 0, 0 },
    {"fastq_trunclen",        required_argument, 0, 0 },
    {"fastq_trunclen_keep",   required_argument, 0, 0 },
    {"fastq_truncqual",       required_argument, 0, 0 },
    {"fastqout",              required_argument, 0, 0 },
    {"fastqout_discarded",    required_argument, 0, 0 },
    {"fastqout_discarded_rev",required_argument, 0, 0 },
    {"fastqout_notmerged_fwd",required_argument, 0, 0 },
    {"fastqout_notmerged_rev",required_argument, 0, 0 },
    {"fastqout_rev",          required_argument, 0, 0 },
    {"fastx_filter",          required_argument, 0, 0 },
    {"fastx_getseq",          required_argument, 0, 0 },
    {"fastx_getseqs",         required_argument, 0, 0 },
    {"fastx_getsubseq",       required_argument, 0, 0 },
    {"fastx_mask",            required_argument, 0, 0 },
    {"fastx_revcomp",         required_argument, 0, 0 },
    {"fastx_subsample",       required_argument, 0, 0 },
    {"fulldp",                no_argument,       0, 0 },
    {"gapext",                required_argument, 0, 0 },
    {"gapopen",               required_argument, 0, 0 },
    {"gzip_decompress",       no_argument,       0, 0 },
    {"h",                     no_argument,       0, 0 },
    {"hardmask",              no_argument,       0, 0 },
    {"help",                  no_argument,       0, 0 },
    {"hspw",                  required_argument, 0, 0 },
    {"id",                    required_argument, 0, 0 },
    {"iddef",                 required_argument, 0, 0 },
    {"idprefix",              required_argument, 0, 0 },
    {"idsuffix",              required_argument, 0, 0 },
    {"join_padgap",           required_argument, 0, 0 },
    {"join_padgapq",          required_argument, 0, 0 },
    {"label",                 required_argument, 0, 0 },
    {"label_field",           required_argument, 0, 0 },
    {"label_substr_match",    no_argument,       0, 0 },
    {"label_suffix",          required_argument, 0, 0 },
    {"label_word",            required_argument, 0, 0 },
    {"label_words",           required_argument, 0, 0 },
    {"labels",                required_argument, 0, 0 },
    {"leftjust",              no_argument,       0, 0 },
    {"length_cutoffs",        required_argument, 0, 0 },
    {"log",                   required_argument, 0, 0 },
    {"makeudb_usearch",       required_argument, 0, 0 },
    {"maskfasta",             required_argument, 0, 0 },
    {"match",                 required_argument, 0, 0 },
    {"matched",               required_argument, 0, 0 },
    {"max_unmasked_pct",      required_argument, 0, 0 },
    {"maxaccepts",            required_argument, 0, 0 },
    {"maxdiffs",              required_argument, 0, 0 },
    {"maxgaps",               required_argument, 0, 0 },
    {"maxhits",               required_argument, 0, 0 },
    {"maxid",                 required_argument, 0, 0 },
    {"maxqsize",              required_argument, 0, 0 },
    {"maxqt",                 required_argument, 0, 0 },
    {"maxrejects",            required_argument, 0, 0 },
    {"maxseqlength",          required_argument, 0, 0 },
    {"maxsize",               required_argument, 0, 0 },
    {"maxsizeratio",          required_argument, 0, 0 },
    {"maxsl",                 required_argument, 0, 0 },
    {"maxsubs",               required_argument, 0, 0 },
    {"maxuniquesize",         required_argument, 0, 0 },
    {"mid",                   required_argument, 0, 0 },
    {"min_unmasked_pct",      required_argument, 0, 0 },
    {"mincols",               required_argument, 0, 0 },
    {"mindiffs",              required_argument, 0, 0 },
    {"mindiv",                required_argument, 0, 0 },
    {"minh",                  required_argument, 0, 0 },
    {"minhsp",                required_argument, 0, 0 },
    {"minqt",                 required_argument, 0, 0 },
    {"minseqlength",          required_argument, 0, 0 },
    {"minsize",               required_argument, 0, 0 },
    {"minsizeratio",          required_argument, 0, 0 },
    {"minsl",                 required_argument, 0, 0 },
    {"mintsize",              required_argument, 0, 0 },
    {"minuniquesize",         required_argument, 0, 0 },
    {"minwordmatches",        required_argument, 0, 0 },
    {"mismatch",              required_argument, 0, 0 },
    {"mothur_shared_out",     required_argument, 0, 0 },
    {"msaout",                required_argument, 0, 0 },
    {"no_progress",           no_argument,       0, 0 },
    {"nonchimeras",           required_argument, 0, 0 },
    {"notmatched",            required_argument, 0, 0 },
    {"notmatchedfq",          required_argument, 0, 0 },
    {"notrunclabels",         no_argument,       0, 0 },
    {"otutabout",             required_argument, 0, 0 },
    {"output",                required_argument, 0, 0 },
    {"output_no_hits",        no_argument,       0, 0 },
    {"pattern",               required_argument, 0, 0 },
    {"profile",               required_argument, 0, 0 },
    {"qmask",                 required_argument, 0, 0 },
    {"query_cov",             required_argument, 0, 0 },
    {"quiet",                 no_argument,       0, 0 },
    {"randseed",              required_argument, 0, 0 },
    {"relabel",               required_argument, 0, 0 },
    {"relabel_keep",          no_argument,       0, 0 },
    {"relabel_md5",           no_argument,       0, 0 },
    {"relabel_self",          no_argument,       0, 0 },
    {"relabel_sha1",          no_argument,       0, 0 },
    {"rereplicate",           required_argument, 0, 0 },
    {"reverse",               required_argument, 0, 0 },
    {"rightjust",             no_argument,       0, 0 },
    {"rowlen",                required_argument, 0, 0 },
    {"samheader",             no_argument,       0, 0 },
    {"samout",                required_argument, 0, 0 },
    {"sample_pct",            required_argument, 0, 0 },
    {"sample_size",           required_argument, 0, 0 },
    {"search_exact",          required_argument, 0, 0 },
    {"self",                  no_argument,       0, 0 },
    {"selfid",                no_argument,       0, 0 },
    {"sff_clip",              no_argument,       0, 0 },
    {"sff_convert",           required_argument, 0, 0 },
    {"shuffle",               required_argument, 0, 0 },
    {"sintax",                required_argument, 0, 0 },
    {"sintax_cutoff",         required_argument, 0, 0 },
    {"sizein",                no_argument,       0, 0 },
    {"sizeorder",             no_argument,       0, 0 },
    {"sizeout",               no_argument,       0, 0 },
    {"slots",                 required_argument, 0, 0 },
    {"sortbylength",          required_argument, 0, 0 },
    {"sortbysize",            required_argument, 0, 0 },
    {"strand",                required_argument, 0, 0 },
    {"subseq_end",            required_argument, 0, 0 },
    {"subseq_start",          required_argument, 0, 0 },
    {"tabbedout",             required_argument, 0, 0 },
    {"target_cov",            required_argument, 0, 0 },
    {"threads",               required_argument, 0, 0 },
    {"top_hits_only",         no_argument,       0, 0 },
    {"topn",                  required_argument, 0, 0 },
    {"uc",                    required_argument, 0, 0 },
    {"uc_allhits",            no_argument,       0, 0 },
    {"uchime2_denovo",        required_argument, 0, 0 },
    {"uchime3_denovo",        required_argument, 0, 0 },
    {"uchime_denovo",         required_argument, 0, 0 },
    {"uchime_ref",            required_argument, 0, 0 },
    {"uchimealns",            required_argument, 0, 0 },
    {"uchimeout",             required_argument, 0, 0 },
    {"uchimeout5",            no_argument,       0, 0 },
    {"udb2fasta",             required_argument, 0, 0 },
    {"udbinfo",               required_argument, 0, 0 },
    {"udbstats",              required_argument, 0, 0 },
    {"unoise_alpha",          required_argument, 0, 0 },
    {"usearch_global",        required_argument, 0, 0 },
    {"userfields",            required_argument, 0, 0 },
    {"userout",               required_argument, 0, 0 },
    {"usersort",              no_argument,       0, 0 },
    {"v",                     no_argument,       0, 0 },
    {"version",               no_argument,       0, 0 },
    {"weak_id",               required_argument, 0, 0 },
    {"wordlength",            required_argument, 0, 0 },
    {"xdrop_nw",              required_argument, 0, 0 },
    {"xee",                   no_argument,       0, 0 },
    {"xn",                    required_argument, 0, 0 },
    {"xsize",                 no_argument,       0, 0 },
    { 0,                      0,                 0, 0 }
  };

  const int options_count = (sizeof(long_options) / sizeof(struct option)) - 1;

  bool options_selected[options_count];
  memset(options_selected, 0, sizeof(options_selected));

  int options_index = 0;
  int c;

  while ((c = getopt_long_only(argc, argv, "", long_options,
                               &options_index)) == 0)
    {
      if (options_index < options_count)
        options_selected[options_index] = 1;

      switch(options_index)
        {
        case option_help:
          opt_help = 1;
          break;

        case option_version:
          opt_version = 1;
          break;

        case option_alnout:
          opt_alnout = optarg;
          break;

        case option_usearch_global:
          opt_usearch_global = optarg;
          break;

        case option_db:
          opt_db = optarg;
          break;

        case option_id:
          opt_id = args_getdouble(optarg);
          break;

        case option_maxaccepts:
          opt_maxaccepts = args_getlong(optarg);
          break;

        case option_maxrejects:
          opt_maxrejects = args_getlong(optarg);
          break;

        case option_wordlength:
          opt_wordlength = args_getlong(optarg);
          break;

        case option_match:
          opt_match = args_getlong(optarg);
          break;

        case option_mismatch:
          opt_mismatch = args_getlong(optarg);
          break;

        case option_fulldp:
          opt_fulldp = 1;
          fprintf(stderr, "WARNING: Option --fulldp is ignored\n");
          break;

        case option_strand:
          if (strcasecmp(optarg, "plus") == 0)
            opt_strand = 1;
          else if (strcasecmp(optarg, "both") == 0)
            opt_strand = 2;
          else
            fatal("The argument to --strand must be plus or both");
          break;

        case option_threads:
          opt_threads = (int64_t) args_getdouble(optarg);
          break;

        case option_gapopen:
          args_get_gap_penalty_string(optarg, 1);
          break;

        case option_gapext:
          args_get_gap_penalty_string(optarg, 0);
          break;

        case option_rowlen:
          opt_rowlen = args_getlong(optarg);
          break;

        case option_userfields:
          if (!parse_userfields_arg(optarg))
            fatal("Unrecognized userfield argument");
          break;

        case option_userout:
          opt_userout = optarg;
          break;

        case option_self:
          opt_self = 1;
          break;

        case option_blast6out:
          opt_blast6out = optarg;
          break;

        case option_uc:
          opt_uc = optarg;
          break;

        case option_weak_id:
          opt_weak_id = args_getdouble(optarg);
          break;

        case option_uc_allhits:
          opt_uc_allhits = 1;
          break;

        case option_notrunclabels:
          opt_notrunclabels = 1;
          break;

        case option_sortbysize:
          opt_sortbysize = optarg;
          break;

        case option_output:
          opt_output = optarg;
          break;

        case option_minsize:
          opt_minsize = args_getlong(optarg);
          if (opt_minsize <= 0)
            fatal("The argument to --minsize must be at least 1");
          break;

        case option_maxsize:
          opt_maxsize = args_getlong(optarg);
          break;

        case option_relabel:
          opt_relabel = optarg;
          break;

        case option_sizeout:
          opt_sizeout = 1;
          break;

        case option_derep_fulllength:
          opt_derep_fulllength = optarg;
          break;

        case option_minseqlength:
          opt_minseqlength = args_getlong(optarg);
          if (opt_minseqlength < 0)
            fatal("The argument to --minseqlength must not be negative");
          break;

        case option_minuniquesize:
          opt_minuniquesize = args_getlong(optarg);
          break;

        case option_topn:
          opt_topn = args_getlong(optarg);
          break;

        case option_maxseqlength:
          opt_maxseqlength = args_getlong(optarg);
          break;

        case option_sizein:
          opt_sizein = 1;
          break;

        case option_sortbylength:
          opt_sortbylength = optarg;
          break;

        case option_matched:
          opt_matched = optarg;
          break;

        case option_notmatched:
          opt_notmatched = optarg;
          break;

        case option_dbmatched:
          opt_dbmatched = optarg;
          break;

        case option_dbnotmatched:
          opt_dbnotmatched = optarg;
          break;

        case option_fastapairs:
          opt_fastapairs = optarg;
          break;

        case option_output_no_hits:
          opt_output_no_hits = 1;
          break;

        case option_maxhits:
          opt_maxhits = args_getlong(optarg);
          break;

        case option_top_hits_only:
          opt_top_hits_only = 1;
          break;

        case option_fasta_width:
          opt_fasta_width = args_getlong(optarg);
          break;

        case option_query_cov:
          opt_query_cov = args_getdouble(optarg);
          break;

        case option_target_cov:
          opt_target_cov = args_getdouble(optarg);
          break;

        case option_idprefix:
          opt_idprefix = args_getlong(optarg);
          break;

        case option_idsuffix:
          opt_idsuffix = args_getlong(optarg);
          break;

        case option_minqt:
          opt_minqt = args_getdouble(optarg);
          break;

        case option_maxqt:
          opt_maxqt = args_getdouble(optarg);
          break;

        case option_minsl:
          opt_minsl = args_getdouble(optarg);
          break;

        case option_maxsl:
          opt_maxsl = args_getdouble(optarg);
          break;

        case option_leftjust:
          opt_leftjust = 1;
          break;

        case option_rightjust:
          opt_rightjust = 1;
          break;

        case option_selfid:
          opt_selfid = 1;
          break;

        case option_maxid:
          opt_maxid = args_getdouble(optarg);
          break;

        case option_minsizeratio:
          opt_minsizeratio = args_getdouble(optarg);
          break;

        case option_maxsizeratio:
          opt_maxsizeratio = args_getdouble(optarg);
          break;

        case option_maxdiffs:
          opt_maxdiffs = args_getlong(optarg);
          break;

        case option_maxsubs:
          opt_maxsubs = args_getlong(optarg);
          break;

        case option_maxgaps:
          opt_maxgaps = args_getlong(optarg);
          break;

        case option_mincols:
          opt_mincols = args_getlong(optarg);
          break;

        case option_maxqsize:
          opt_maxqsize = args_getlong(optarg);
          break;

        case option_mintsize:
          opt_mintsize = args_getlong(optarg);
          break;

        case option_mid:
          opt_mid = args_getdouble(optarg);
          break;

        case option_shuffle:
          opt_shuffle = optarg;
          break;

        case option_randseed:
          opt_randseed = args_getlong(optarg);
          break;

        case option_maskfasta:
          opt_maskfasta = optarg;
          break;

        case option_hardmask:
          opt_hardmask = 1;
          break;

        case option_qmask:
          if (strcasecmp(optarg, "none") == 0)
            opt_qmask = MASK_NONE;
          else if (strcasecmp(optarg, "dust") == 0)
            opt_qmask = MASK_DUST;
          else if (strcasecmp(optarg, "soft") == 0)
            opt_qmask = MASK_SOFT;
          else
            opt_qmask = MASK_ERROR;
          break;

        case option_dbmask:
          if (strcasecmp(optarg, "none") == 0)
            opt_dbmask = MASK_NONE;
          else if (strcasecmp(optarg, "dust") == 0)
            opt_dbmask = MASK_DUST;
          else if (strcasecmp(optarg, "soft") == 0)
            opt_dbmask = MASK_SOFT;
          else
            opt_dbmask = MASK_ERROR;
          break;

        case option_cluster_smallmem:
          opt_cluster_smallmem = optarg;
          break;

        case option_cluster_fast:
          opt_cluster_fast = optarg;
          break;

        case option_centroids:
          opt_centroids = optarg;
          break;

        case option_clusters:
          opt_clusters = optarg;
          break;

        case option_consout:
          opt_consout = optarg;
          break;

        case option_cons_truncate:
          fprintf(stderr, "WARNING: Option --cons_truncate is ignored\n");
          opt_cons_truncate = 1;
          break;

        case option_msaout:
          opt_msaout = optarg;
          break;

        case option_usersort:
          opt_usersort = 1;
          break;

        case option_xn:
          opt_xn = args_getdouble(optarg);
          break;

        case option_iddef:
          opt_iddef = args_getlong(optarg);
          break;

        case option_slots:
          fprintf(stderr, "WARNING: Option --slots is ignored\n");
          opt_slots = args_getlong(optarg);
          break;

        case option_pattern:
          fprintf(stderr, "WARNING: Option --pattern is ignored\n");
          opt_pattern = optarg;
          break;

        case option_maxuniquesize:
          opt_maxuniquesize = args_getlong(optarg);
          break;

        case option_abskew:
          opt_abskew = args_getdouble(optarg);
          break;

        case option_chimeras:
          opt_chimeras = optarg;
          break;

        case option_dn:
          opt_dn = args_getdouble(optarg);
          break;

        case option_mindiffs:
          opt_mindiffs = args_getlong(optarg);
          break;

        case option_mindiv:
          opt_mindiv = args_getdouble(optarg);
          break;

        case option_minh:
          opt_minh = args_getdouble(optarg);
          break;

        case option_nonchimeras:
          opt_nonchimeras = optarg;
          break;

        case option_uchime_denovo:
          opt_uchime_denovo = optarg;
          break;

        case option_uchime_ref:
          opt_uchime_ref = optarg;
          break;

        case option_uchimealns:
          opt_uchimealns = optarg;
          break;

        case option_uchimeout:
          opt_uchimeout = optarg;
          break;

        case option_uchimeout5:
          opt_uchimeout5 = 1;
          break;

        case option_alignwidth:
          opt_alignwidth = args_getlong(optarg);
          break;

        case option_allpairs_global:
          opt_allpairs_global = optarg;
          break;

        case option_acceptall:
          opt_acceptall = 1;
          break;

        case option_cluster_size:
          opt_cluster_size = optarg;
          break;

        case option_samout:
          opt_samout = optarg;
          break;

        case option_log:
          opt_log = optarg;
          break;

        case option_quiet:
          opt_quiet = true;
          break;

        case option_fastx_subsample:
          opt_fastx_subsample = optarg;
          break;

        case option_sample_pct:
          opt_sample_pct = args_getdouble(optarg);
          break;

        case option_fastq_chars:
          opt_fastq_chars = optarg;
          break;

        case option_profile:
          opt_profile = optarg;
          break;

        case option_sample_size:
          opt_sample_size = args_getlong(optarg);
          break;

        case option_fastaout:
          opt_fastaout = optarg;
          break;

        case option_xsize:
          opt_xsize = 1;
          break;

        case option_clusterout_id:
          opt_clusterout_id = 1;
          break;

        case option_clusterout_sort:
          opt_clusterout_sort = 1;
          break;

        case option_borderline:
          opt_borderline = optarg;
          break;

        case option_relabel_sha1:
          opt_relabel_sha1 = 1;
          break;

        case option_relabel_md5:
          opt_relabel_md5 = 1;
          break;

        case option_derep_prefix:
          opt_derep_prefix = optarg;
          break;

        case option_fastq_filter:
          opt_fastq_filter = optarg;
          break;

        case option_fastqout:
          opt_fastqout = optarg;
          break;

        case option_fastaout_discarded:
          opt_fastaout_discarded = optarg;
          break;

        case option_fastqout_discarded:
          opt_fastqout_discarded = optarg;
          break;

        case option_fastq_truncqual:
          opt_fastq_truncqual = args_getlong(optarg);
          break;

        case option_fastq_maxee:
          opt_fastq_maxee = args_getdouble(optarg);
          break;

        case option_fastq_trunclen:
          opt_fastq_trunclen = args_getlong(optarg);
          break;

        case option_fastq_minlen:
          opt_fastq_minlen = args_getlong(optarg);
          break;

        case option_fastq_stripleft:
          opt_fastq_stripleft = args_getlong(optarg);
          break;

        case option_fastq_maxee_rate:
          opt_fastq_maxee_rate = args_getdouble(optarg);
          break;

        case option_fastq_maxns:
          opt_fastq_maxns = args_getlong(optarg);
          break;

        case option_eeout:
          opt_eeout = 1;
          break;

        case option_fastq_ascii:
          opt_fastq_ascii = args_getlong(optarg);
          break;

        case option_fastq_qmin:
          opt_fastq_qmin = args_getlong(optarg);
          break;

        case option_fastq_qmax:
          opt_fastq_qmax = args_getlong(optarg);
          break;

        case option_fastq_qmaxout:
          opt_fastq_qmaxout = args_getlong(optarg);
          break;

        case option_fastq_stats:
          opt_fastq_stats = optarg;
          break;

        case option_fastq_tail:
          opt_fastq_tail = args_getlong(optarg);
          break;

        case option_fastx_revcomp:
          opt_fastx_revcomp = optarg;
          break;

        case option_label_suffix:
          opt_label_suffix = optarg;
          break;

        case option_h:
          opt_help = 1;
          break;

        case option_samheader:
          opt_samheader = 1;
          break;

        case option_sizeorder:
          opt_sizeorder = 1;
          break;

        case option_minwordmatches:
          opt_minwordmatches = args_getlong(optarg);
          if (opt_minwordmatches < 0)
            fatal("The argument to --minwordmatches must not be negative");
          break;

        case option_v:
          opt_version = 1;
          break;

        case option_relabel_keep:
          opt_relabel_keep = 1;
          break;

        case option_search_exact:
          opt_search_exact = optarg;
          break;

        case option_fastx_mask:
          opt_fastx_mask = optarg;
          break;

        case option_min_unmasked_pct:
          opt_min_unmasked_pct = args_getdouble(optarg);
          break;

        case option_max_unmasked_pct:
          opt_max_unmasked_pct = args_getdouble(optarg);
          break;

        case option_fastq_convert:
          opt_fastq_convert = optarg;
          break;

        case option_fastq_asciiout:
          opt_fastq_asciiout = args_getlong(optarg);
          break;

        case option_fastq_qminout:
          opt_fastq_qminout = args_getlong(optarg);
          break;

        case option_fastq_mergepairs:
          opt_fastq_mergepairs = optarg;
          break;

        case option_fastq_eeout:
          opt_fastq_eeout = 1;
          break;

        case option_fastqout_notmerged_fwd:
          opt_fastqout_notmerged_fwd = optarg;
          break;

        case option_fastqout_notmerged_rev:
          opt_fastqout_notmerged_rev = optarg;
          break;

        case option_fastq_minovlen:
          opt_fastq_minovlen = args_getlong(optarg);
          break;

        case option_fastq_minmergelen:
          opt_fastq_minmergelen = args_getlong(optarg);
          break;

        case option_fastq_maxmergelen:
          opt_fastq_maxmergelen = args_getlong(optarg);
          break;

        case option_fastq_nostagger:
          opt_fastq_nostagger = optarg;
          break;

        case option_fastq_allowmergestagger:
          opt_fastq_allowmergestagger = 1;
          break;

        case option_fastq_maxdiffs:
          opt_fastq_maxdiffs = args_getlong(optarg);
          break;

        case option_fastaout_notmerged_fwd:
          opt_fastaout_notmerged_fwd = optarg;
          break;

        case option_fastaout_notmerged_rev:
          opt_fastaout_notmerged_rev = optarg;
          break;

        case option_reverse:
          opt_reverse = optarg;
          break;

        case option_eetabbedout:
          opt_eetabbedout = optarg;
          break;

        case option_fasta_score:
          opt_fasta_score = 1;
          break;

        case option_fastq_eestats:
          opt_fastq_eestats = optarg;
          break;

        case option_rereplicate:
          opt_rereplicate = optarg;
          break;

        case option_xdrop_nw:
          /* xdrop_nw ignored */
          fprintf(stderr, "WARNING: Option --xdrop_nw is ignored\n");
          break;

        case option_minhsp:
          /* minhsp ignored */
          fprintf(stderr, "WARNING: Option --minhsp is ignored\n");
          break;

        case option_band:
          /* band ignored */
          fprintf(stderr, "WARNING: Option --band is ignored\n");
          break;

        case option_hspw:
          /* hspw ignored */
          fprintf(stderr, "WARNING: Option --hspw is ignored\n");
          break;

        case option_gzip_decompress:
          opt_gzip_decompress = 1;
          break;

        case option_bzip2_decompress:
          opt_bzip2_decompress = 1;
          break;

        case option_fastq_maxlen:
          opt_fastq_maxlen = args_getlong(optarg);
          break;

        case option_fastq_truncee:
          opt_fastq_truncee = args_getdouble(optarg);
          break;

        case option_fastx_filter:
          opt_fastx_filter = optarg;
          break;

        case option_otutabout:
          opt_otutabout = optarg;
          break;

        case option_mothur_shared_out:
          opt_mothur_shared_out = optarg;
          break;

        case option_biomout:
          opt_biomout = optarg;
          break;

        case option_fastq_trunclen_keep:
          opt_fastq_trunclen_keep = args_getlong(optarg);
          break;

        case option_fastq_stripright:
          opt_fastq_stripright = args_getlong(optarg);
          break;

        case option_no_progress:
          opt_no_progress = 1;
          break;

        case option_fastq_eestats2:
          opt_fastq_eestats2 = optarg;
          break;

        case option_ee_cutoffs:
          args_get_ee_cutoffs(optarg);
          break;

        case option_length_cutoffs:
          args_get_length_cutoffs(optarg);
          break;

        case option_makeudb_usearch:
          opt_makeudb_usearch = optarg;
          break;

        case option_udb2fasta:
          opt_udb2fasta = optarg;
          break;

        case option_udbinfo:
          opt_udbinfo = optarg;
          break;

        case option_udbstats:
          opt_udbstats = optarg;
          break;

        case option_cluster_unoise:
          opt_cluster_unoise = optarg;
          break;

        case option_unoise_alpha:
          opt_unoise_alpha = args_getdouble(optarg);
          break;

        case option_uchime2_denovo:
          opt_uchime2_denovo = optarg;
          break;

        case option_uchime3_denovo:
          opt_uchime3_denovo = optarg;
          break;

        case option_sintax:
          opt_sintax = optarg;
          break;

        case option_sintax_cutoff:
          opt_sintax_cutoff = args_getdouble(optarg);
          break;

        case option_tabbedout:
          opt_tabbedout = optarg;
          break;

        case option_fastq_maxdiffpct:
          opt_fastq_maxdiffpct = args_getdouble(optarg);
          break;

        case option_fastq_join:
          opt_fastq_join = optarg;
          break;

        case option_join_padgap:
          opt_join_padgap = optarg;
          break;

        case option_join_padgapq:
          opt_join_padgapq = optarg;
          break;

        case option_sff_convert:
          opt_sff_convert = optarg;
          break;

        case option_sff_clip:
          opt_sff_clip = 1;
          break;

        case option_fastaout_rev:
          opt_fastaout_rev = optarg;
          break;

        case option_fastaout_discarded_rev:
          opt_fastaout_discarded_rev = optarg;
          break;

        case option_fastqout_rev:
          opt_fastqout_rev = optarg;
          break;

        case option_fastqout_discarded_rev:
          opt_fastqout_discarded_rev = optarg;
          break;

        case option_xee:
          opt_xee = 1;
          break;

        case option_fastx_getseq:
          opt_fastx_getseq = optarg;
          break;

        case option_fastx_getseqs:
          opt_fastx_getseqs = optarg;
          break;

        case option_fastx_getsubseq:
          opt_fastx_getsubseq = optarg;
          break;

        case option_label_substr_match:
          opt_label_substr_match = true;
          break;

        case option_label:
          opt_label = optarg;
          break;

        case option_subseq_start:
          opt_subseq_start = args_getlong(optarg);
          break;

        case option_subseq_end:
          opt_subseq_end = args_getlong(optarg);
          break;

        case option_notmatchedfq:
          opt_notmatchedfq = optarg;
          break;

        case option_label_field:
          opt_label_field = optarg;
          break;

        case option_label_word:
          opt_label_word = optarg;
          break;

        case option_label_words:
          opt_label_words = optarg;
          break;

        case option_labels:
          opt_labels = optarg;
          break;

        case option_cut:
          opt_cut = optarg;
          break;

        case option_cut_pattern:
          opt_cut_pattern = optarg;
          break;

        case option_relabel_self:
          opt_relabel_self = 1;
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

  /* Below is a list of all command names, in alphabetical order. */

  int command_options[] =
    {
      option_allpairs_global,
      option_cluster_fast,
      option_cluster_size,
      option_cluster_smallmem,
      option_cluster_unoise,
      option_cut,
      option_derep_fulllength,
      option_derep_prefix,
      option_fastq_chars,
      option_fastq_convert,
      option_fastq_eestats,
      option_fastq_eestats2,
      option_fastq_filter,
      option_fastq_join,
      option_fastq_mergepairs,
      option_fastq_stats,
      option_fastx_filter,
      option_fastx_getseq,
      option_fastx_getseqs,
      option_fastx_getsubseq,
      option_fastx_mask,
      option_fastx_revcomp,
      option_fastx_subsample,
      option_h,
      option_help,
      option_makeudb_usearch,
      option_maskfasta,
      option_rereplicate,
      option_search_exact,
      option_sff_convert,
      option_shuffle,
      option_sintax,
      option_sortbylength,
      option_sortbysize,
      option_uchime2_denovo,
      option_uchime3_denovo,
      option_uchime_denovo,
      option_uchime_ref,
      option_udb2fasta,
      option_udbinfo,
      option_udbstats,
      option_usearch_global,
      option_v,
      option_version
    };

  const int commands_count = sizeof(command_options) / sizeof(int);

  /*
    Below is a list of all the options that are valid for each command.
    The first line is the command and the lines below are the valid options.
  */

  const int valid_options[][92] =
    {
      {
        option_allpairs_global,
        option_acceptall,
        option_alnout,
        option_band,
        option_blast6out,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastapairs,
        option_fulldp,
        option_gapext,
        option_gapopen,
        option_gzip_decompress,
        option_hardmask,
        option_hspw,
        option_id,
        option_iddef,
        option_idprefix,
        option_idsuffix,
        option_leftjust,
        option_log,
        option_match,
        option_matched,
        option_maxaccepts,
        option_maxdiffs,
        option_maxgaps,
        option_maxhits,
        option_maxid,
        option_maxqsize,
        option_maxqt,
        option_maxrejects,
        option_maxseqlength,
        option_maxsizeratio,
        option_maxsl,
        option_maxsubs,
        option_mid,
        option_mincols,
        option_minhsp,
        option_minqt,
        option_minseqlength,
        option_minsizeratio,
        option_minsl,
        option_mintsize,
        option_minwordmatches,
        option_mismatch,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_output_no_hits,
        option_pattern,
        option_qmask,
        option_query_cov,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rightjust,
        option_rowlen,
        option_samheader,
        option_samout,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeout,
        option_slots,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_uc,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xsize,
        -1 },

      { option_cluster_fast,
        option_alnout,
        option_band,
        option_biomout,
        option_blast6out,
        option_bzip2_decompress,
        option_centroids,
        option_clusterout_id,
        option_clusterout_sort,
        option_clusters,
        option_cons_truncate,
        option_consout,
        option_fasta_width,
        option_fastapairs,
        option_fulldp,
        option_gapext,
        option_gapopen,
        option_gzip_decompress,
        option_hardmask,
        option_hspw,
        option_id,
        option_iddef,
        option_idprefix,
        option_idsuffix,
        option_leftjust,
        option_log,
        option_match,
        option_matched,
        option_maxaccepts,
        option_maxdiffs,
        option_maxgaps,
        option_maxhits,
        option_maxid,
        option_maxqsize,
        option_maxqt,
        option_maxrejects,
        option_maxseqlength,
        option_maxsizeratio,
        option_maxsl,
        option_maxsubs,
        option_mid,
        option_mincols,
        option_minhsp,
        option_minqt,
        option_minseqlength,
        option_minsizeratio,
        option_minsl,
        option_mintsize,
        option_minwordmatches,
        option_mismatch,
        option_mothur_shared_out,
        option_msaout,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_profile,
        option_qmask,
        option_query_cov,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rightjust,
        option_rowlen,
        option_samheader,
        option_samout,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeorder,
        option_sizeout,
        option_slots,
        option_strand,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_uc,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xsize,
        -1 },

      { option_cluster_size,
        option_alnout,
        option_band,
        option_biomout,
        option_blast6out,
        option_bzip2_decompress,
        option_centroids,
        option_clusterout_id,
        option_clusterout_sort,
        option_clusters,
        option_cons_truncate,
        option_consout,
        option_fasta_width,
        option_fastapairs,
        option_fulldp,
        option_gapext,
        option_gapopen,
        option_gzip_decompress,
        option_hardmask,
        option_hspw,
        option_id,
        option_iddef,
        option_idprefix,
        option_idsuffix,
        option_leftjust,
        option_log,
        option_match,
        option_matched,
        option_maxaccepts,
        option_maxdiffs,
        option_maxgaps,
        option_maxhits,
        option_maxid,
        option_maxqsize,
        option_maxqt,
        option_maxrejects,
        option_maxseqlength,
        option_maxsizeratio,
        option_maxsl,
        option_maxsubs,
        option_mid,
        option_mincols,
        option_minhsp,
        option_minqt,
        option_minseqlength,
        option_minsizeratio,
        option_minsl,
        option_mintsize,
        option_minwordmatches,
        option_mismatch,
        option_mothur_shared_out,
        option_msaout,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_profile,
        option_qmask,
        option_query_cov,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rightjust,
        option_rowlen,
        option_samheader,
        option_samout,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeorder,
        option_sizeout,
        option_slots,
        option_strand,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_uc,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xsize,
        -1 },

      { option_cluster_smallmem,
        option_alnout,
        option_band,
        option_biomout,
        option_blast6out,
        option_bzip2_decompress,
        option_centroids,
        option_clusterout_id,
        option_clusterout_sort,
        option_clusters,
        option_cons_truncate,
        option_consout,
        option_fasta_width,
        option_fastapairs,
        option_fulldp,
        option_gapext,
        option_gapopen,
        option_gzip_decompress,
        option_hardmask,
        option_hspw,
        option_id,
        option_iddef,
        option_idprefix,
        option_idsuffix,
        option_leftjust,
        option_log,
        option_match,
        option_matched,
        option_maxaccepts,
        option_maxdiffs,
        option_maxgaps,
        option_maxhits,
        option_maxid,
        option_maxqsize,
        option_maxqt,
        option_maxrejects,
        option_maxseqlength,
        option_maxsizeratio,
        option_maxsl,
        option_maxsubs,
        option_mid,
        option_mincols,
        option_minhsp,
        option_minqt,
        option_minseqlength,
        option_minsizeratio,
        option_minsl,
        option_mintsize,
        option_minwordmatches,
        option_mismatch,
        option_mothur_shared_out,
        option_msaout,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_profile,
        option_qmask,
        option_query_cov,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rightjust,
        option_rowlen,
        option_samheader,
        option_samout,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeorder,
        option_sizeout,
        option_slots,
        option_strand,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_uc,
        option_userfields,
        option_userout,
        option_usersort,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xsize,
        -1 },

      { option_cluster_unoise,
        option_alnout,
        option_band,
        option_biomout,
        option_blast6out,
        option_bzip2_decompress,
        option_centroids,
        option_clusterout_id,
        option_clusterout_sort,
        option_clusters,
        option_cons_truncate,
        option_consout,
        option_fasta_width,
        option_fastapairs,
        option_fulldp,
        option_gapext,
        option_gapopen,
        option_gzip_decompress,
        option_hardmask,
        option_hspw,
        option_id,
        option_iddef,
        option_idprefix,
        option_idsuffix,
        option_leftjust,
        option_log,
        option_match,
        option_matched,
        option_maxaccepts,
        option_maxdiffs,
        option_maxgaps,
        option_maxhits,
        option_maxid,
        option_maxqsize,
        option_maxqt,
        option_maxrejects,
        option_maxseqlength,
        option_maxsizeratio,
        option_maxsl,
        option_maxsubs,
        option_mid,
        option_mincols,
        option_minhsp,
        option_minqt,
        option_minseqlength,
        option_minsizeratio,
        option_minsize,
        option_minsl,
        option_mintsize,
        option_minwordmatches,
        option_mismatch,
        option_mothur_shared_out,
        option_msaout,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_profile,
        option_qmask,
        option_query_cov,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rightjust,
        option_rowlen,
        option_samheader,
        option_samout,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeorder,
        option_sizeout,
        option_slots,
        option_strand,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_uc,
        option_unoise_alpha,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xsize,
        -1 },

      { option_cut,
        option_bzip2_decompress,
        option_cut_pattern,
        option_fasta_width,
        option_fastaout,
        option_fastaout_discarded,
        option_fastaout_discarded_rev,
        option_fastaout_rev,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_xee,
        option_xsize,
        -1 },

      { option_derep_fulllength,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_log,
        option_maxseqlength,
        option_maxuniquesize,
        option_minseqlength,
        option_minuniquesize,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_topn,
        option_uc,
        option_xee,
        option_xsize,
        -1 },

      { option_derep_prefix,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_log,
        option_maxseqlength,
        option_maxuniquesize,
        option_minseqlength,
        option_minuniquesize,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_topn,
        option_uc,
        option_xee,
        option_xsize,
        -1 },

      { option_fastq_chars,
        option_bzip2_decompress,
        option_fastq_tail,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_quiet,
        option_threads,
        -1 },

      { option_fastq_convert,
        option_bzip2_decompress,
        option_fastq_ascii,
        option_fastq_asciiout,
        option_fastq_qmax,
        option_fastq_qmaxout,
        option_fastq_qmin,
        option_fastq_qminout,
        option_fastqout,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastq_eestats,
        option_bzip2_decompress,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_output,
        option_quiet,
        option_threads,
        -1 },

      { option_fastq_eestats2,
        option_bzip2_decompress,
        option_ee_cutoffs,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_length_cutoffs,
        option_log,
        option_no_progress,
        option_output,
        option_quiet,
        option_threads,
        -1 },

      { option_fastq_filter,
        option_bzip2_decompress,
        option_eeout,
        option_fasta_width,
        option_fastaout,
        option_fastaout_discarded,
        option_fastaout_discarded_rev,
        option_fastaout_rev,
        option_fastq_ascii,
        option_fastq_eeout,
        option_fastq_maxee,
        option_fastq_maxee_rate,
        option_fastq_maxlen,
        option_fastq_maxns,
        option_fastq_minlen,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastq_stripleft,
        option_fastq_stripright,
        option_fastq_truncee,
        option_fastq_trunclen,
        option_fastq_trunclen_keep,
        option_fastq_truncqual,
        option_fastqout,
        option_fastqout_discarded,
        option_fastqout_discarded_rev,
        option_fastqout_rev,
        option_gzip_decompress,
        option_log,
        option_maxsize,
        option_minsize,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_reverse,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastq_join,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_gzip_decompress,
        option_join_padgap,
        option_join_padgapq,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_reverse,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastq_mergepairs,
        option_bzip2_decompress,
        option_eeout,
        option_eetabbedout,
        option_fasta_width,
        option_fastaout,
        option_fastaout_notmerged_fwd,
        option_fastaout_notmerged_rev,
        option_fastq_allowmergestagger,
        option_fastq_ascii,
        option_fastq_eeout,
        option_fastq_maxdiffpct,
        option_fastq_maxdiffs,
        option_fastq_maxee,
        option_fastq_maxlen,
        option_fastq_maxmergelen,
        option_fastq_maxns,
        option_fastq_minlen,
        option_fastq_minmergelen,
        option_fastq_minovlen,
        option_fastq_nostagger,
        option_fastq_qmax,
        option_fastq_qmaxout,
        option_fastq_qmin,
        option_fastq_truncqual,
        option_fastqout,
        option_fastqout_notmerged_fwd,
        option_fastqout_notmerged_rev,
        option_gzip_decompress,
        option_label_suffix,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_reverse,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastq_stats,
        option_bzip2_decompress,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_output,
        option_quiet,
        option_threads,
        -1 },

      { option_fastx_filter,
        option_bzip2_decompress,
        option_eeout,
        option_fasta_width,
        option_fastaout,
        option_fastaout_discarded,
        option_fastaout_discarded_rev,
        option_fastaout_rev,
        option_fastq_ascii,
        option_fastq_eeout,
        option_fastq_maxee,
        option_fastq_maxee_rate,
        option_fastq_maxlen,
        option_fastq_maxns,
        option_fastq_minlen,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastq_stripleft,
        option_fastq_stripright,
        option_fastq_truncee,
        option_fastq_trunclen,
        option_fastq_trunclen_keep,
        option_fastq_truncqual,
        option_fastqout,
        option_fastqout_discarded,
        option_fastqout_discarded_rev,
        option_fastqout_rev,
        option_gzip_decompress,
        option_log,
        option_maxsize,
        option_minsize,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_reverse,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastx_getseq,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_gzip_decompress,
        option_label,
        option_label_substr_match,
        option_log,
        option_no_progress,
        option_notmatched,
        option_notmatchedfq,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastx_getseqs,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_gzip_decompress,
        option_label,
        option_label_field,
        option_label_substr_match,
        option_label_word,
        option_label_words,
        option_labels,
        option_log,
        option_no_progress,
        option_notmatched,
        option_notmatchedfq,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastx_getsubseq,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_gzip_decompress,
        option_label,
        option_label_substr_match,
        option_log,
        option_no_progress,
        option_notmatched,
        option_notmatchedfq,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_subseq_end,
        option_subseq_start,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastx_mask,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_gzip_decompress,
        option_hardmask,
        option_log,
        option_max_unmasked_pct,
        option_min_unmasked_pct,
        option_no_progress,
        option_notrunclabels,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastx_revcomp,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_gzip_decompress,
        option_label_suffix,
        option_log,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_fastx_subsample,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastaout_discarded,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_fastqout_discarded,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_randseed,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample_pct,
        option_sample_size,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_h,
        option_log,
        option_quiet,
        option_threads,
        -1 },

      { option_help,
        option_log,
        option_quiet,
        option_threads,
        -1 },

      { option_makeudb_usearch,
        option_bzip2_decompress,
        option_dbmask,
        option_gzip_decompress,
        option_hardmask,
        option_log,
        option_minseqlength,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_quiet,
        option_threads,
        option_wordlength,
        -1 },

      { option_maskfasta,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_hardmask,
        option_log,
        option_max_unmasked_pct,
        option_maxseqlength,
        option_min_unmasked_pct,
        option_minseqlength,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_rereplicate,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_search_exact,
        option_alnout,
        option_biomout,
        option_blast6out,
        option_bzip2_decompress,
        option_db,
        option_dbmask,
        option_dbmatched,
        option_dbnotmatched,
        option_fasta_width,
        option_fastapairs,
        option_gzip_decompress,
        option_hardmask,
        option_log,
        option_match,
        option_matched,
        option_maxhits,
        option_maxqsize,
        option_maxqt,
        option_maxseqlength,
        option_maxsizeratio,
        option_maxsl,
        option_mincols,
        option_minqt,
        option_minseqlength,
        option_minsizeratio,
        option_minsl,
        option_mintsize,
        option_mismatch,
        option_mothur_shared_out,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rowlen,
        option_samheader,
        option_samout,
        option_self,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_top_hits_only,
        option_uc,
        option_uc_allhits,
        option_userfields,
        option_userout,
        option_xee,
        option_xsize,
        -1 },

      { option_sff_convert,
        option_fastq_asciiout,
        option_fastq_qmaxout,
        option_fastq_qminout,
        option_fastqout,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sff_clip,
        option_sizeout,
        option_threads,
        -1 },

      { option_shuffle,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_log,
        option_maxseqlength,
        option_minseqlength,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_quiet,
        option_randseed,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_topn,
        option_xee,
        option_xsize,
        -1 },

      { option_sintax,
        option_bzip2_decompress,
        option_db,
        option_dbmask,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_sintax_cutoff,
        option_strand,
        option_tabbedout,
        option_threads,
        option_wordlength,
        -1 },

      { option_sortbylength,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_log,
        option_maxseqlength,
        option_minseqlength,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_topn,
        option_xee,
        option_xsize,
        -1 },

      { option_sortbysize,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_log,
        option_maxseqlength,
        option_maxsize,
        option_minseqlength,
        option_minsize,
        option_no_progress,
        option_notrunclabels,
        option_output,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_topn,
        option_xee,
        option_xsize,
        -1 },

      { option_uchime2_denovo,
        option_abskew,
        option_alignwidth,
        option_borderline,
        option_chimeras,
        option_dn,
        option_fasta_score,
        option_fasta_width,
        option_gapext,
        option_gapopen,
        option_hardmask,
        option_log,
        option_match,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_mismatch,
        option_no_progress,
        option_nonchimeras,
        option_notrunclabels,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_uchimealns,
        option_uchimeout,
        option_uchimeout5,
        option_xee,
        option_xn,
        option_xsize,
        -1 },

      { option_uchime3_denovo,
        option_abskew,
        option_alignwidth,
        option_borderline,
        option_chimeras,
        option_dn,
        option_fasta_score,
        option_fasta_width,
        option_gapext,
        option_gapopen,
        option_hardmask,
        option_log,
        option_match,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_mismatch,
        option_no_progress,
        option_nonchimeras,
        option_notrunclabels,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_uchimealns,
        option_uchimeout,
        option_uchimeout5,
        option_xee,
        option_xn,
        option_xsize,
        -1 },

      { option_uchime_denovo,
        option_abskew,
        option_alignwidth,
        option_borderline,
        option_chimeras,
        option_dn,
        option_fasta_score,
        option_fasta_width,
        option_gapext,
        option_gapopen,
        option_hardmask,
        option_log,
        option_match,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_mismatch,
        option_no_progress,
        option_nonchimeras,
        option_notrunclabels,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_uchimealns,
        option_uchimeout,
        option_uchimeout5,
        option_xee,
        option_xn,
        option_xsize,
        -1 },

      { option_uchime_ref,
        option_abskew,
        option_alignwidth,
        option_borderline,
        option_chimeras,
        option_db,
        option_dbmask,
        option_dn,
        option_fasta_score,
        option_fasta_width,
        option_gapext,
        option_gapopen,
        option_hardmask,
        option_log,
        option_match,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_mismatch,
        option_no_progress,
        option_nonchimeras,
        option_notrunclabels,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_uchimealns,
        option_uchimeout,
        option_uchimeout5,
        option_xee,
        option_xn,
        option_xsize,
        -1 },

      { option_udb2fasta,
        option_fasta_width,
        option_log,
        option_no_progress,
        option_output,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xsize,
        -1 },

      { option_udbinfo,
        option_log,
        option_quiet,
        option_threads,
        -1 },

      { option_udbstats,
        option_log,
        option_no_progress,
        option_quiet,
        option_threads,
        -1 },

      { option_usearch_global,
        option_alnout,
        option_band,
        option_biomout,
        option_blast6out,
        option_bzip2_decompress,
        option_db,
        option_dbmask,
        option_dbmatched,
        option_dbnotmatched,
        option_fasta_width,
        option_fastapairs,
        option_fulldp,
        option_gapext,
        option_gapopen,
        option_gzip_decompress,
        option_hardmask,
        option_hspw,
        option_id,
        option_iddef,
        option_idprefix,
        option_idsuffix,
        option_leftjust,
        option_log,
        option_match,
        option_matched,
        option_maxaccepts,
        option_maxdiffs,
        option_maxgaps,
        option_maxhits,
        option_maxid,
        option_maxqsize,
        option_maxqt,
        option_maxrejects,
        option_maxseqlength,
        option_maxsizeratio,
        option_maxsl,
        option_maxsubs,
        option_mid,
        option_mincols,
        option_minhsp,
        option_minqt,
        option_minseqlength,
        option_minsizeratio,
        option_minsl,
        option_mintsize,
        option_minwordmatches,
        option_mismatch,
        option_mothur_shared_out,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_qmask,
        option_query_cov,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rightjust,
        option_rowlen,
        option_samheader,
        option_samout,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeout,
        option_slots,
        option_strand,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_uc,
        option_uc_allhits,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xsize,
        -1 },

      { option_v,
        option_log,
        option_quiet,
        option_threads,
        -1 },

      { option_version,
        option_log,
        option_quiet,
        option_threads,
        -1 }
    };

  /* check that only one commmand is specified */
  int commands = 0;
  int k = -1;
  for (int i = 0; i < commands_count; i++)
    if (options_selected[command_options[i]])
      {
        commands++;
        k = i;
      }
  if (commands > 1)
    fatal("More than one command specified");

  /* check that only valid options are specified */
  int invalid_options = 0;

  if (commands == 0)
    {
      /* check if any options are specified */
      bool any_options = false;
      for (int i = 0; i < options_count; i++)
        if (options_selected[i])
          any_options = true;
      if (any_options)
        fprintf(stderr, "WARNING: Options given, but no valid command specified.\n");
    }
  else
    {
      for (int i = 0; i < options_count; i++)
        {
          if (options_selected[i])
            {
              int j = 0;
              bool ok = false;
              while (valid_options[k][j] >= 0)
                {
                  if (valid_options[k][j] == i)
                    {
                      ok = true;
                      break;
                    }
                  j++;
                }
              if (! ok)
                {
                  invalid_options++;

                  if (invalid_options == 1)
                    {
                      fprintf(stderr,
                              "Fatal error: Invalid options to command %s\n",
                              long_options[command_options[k]].name);
                      fprintf(stderr,
                              "Invalid option(s):");
                    }
                  fprintf(stderr, " --%s",
                          long_options[i].name);
                }
            }
        }

      if (invalid_options > 0)
        {
          fprintf(stderr, "\nThe valid options for the %s command are:",
                  long_options[command_options[k]].name);
          int count = 0;
          for(int j = 1; valid_options[k][j] >= 0; j++)
            {
              fprintf(stderr, " --%s", long_options[valid_options[k][j]].name);
              count++;
            }
          if (! count)
            fprintf(stderr, " (none)");
          fprintf(stderr, "\n");
          exit(EXIT_FAILURE);
        }
    }

  /* multi-threaded commands */

  if ((opt_threads < 0) || (opt_threads > 1024))
    fatal("The argument to --threads must be in the range 0 (default) to 1024");

  if (opt_allpairs_global || opt_cluster_fast || opt_cluster_size ||
      opt_cluster_smallmem || opt_cluster_unoise || opt_fastq_mergepairs ||
      opt_fastx_mask || opt_maskfasta || opt_search_exact || opt_sintax ||
      opt_uchime_ref || opt_usearch_global)
    {
      if (opt_threads == 0)
        opt_threads = arch_get_cores();
    }
  else
    {
      if (opt_threads > 1)
        fprintf(stderr, "WARNING: The %s command does not support multithreading.\nOnly 1 thread used.\n", long_options[command_options[k]].name);
      opt_threads = 1;
    }

  if (opt_cluster_unoise)
    opt_weak_id = 0.90;
  else
    if (opt_weak_id > opt_id)
      opt_weak_id = opt_id;

  if (opt_maxrejects == -1)
    {
      if (opt_cluster_fast)
        opt_maxrejects = 8;
      else
        opt_maxrejects = 32;
    }

  if (opt_maxaccepts < 0)
    fatal("The argument to --maxaccepts must not be negative");

  if (opt_maxrejects < 0)
    fatal("The argument to --maxrejects must not be negative");

  if ((opt_wordlength < 3) || (opt_wordlength > 15))
    fatal("The argument to --wordlength must be in the range 3 to 15");

  if ((opt_iddef < 0) || (opt_iddef > 4))
    fatal("The argument to --iddef must in the range 0 to 4");

#if 0

  if (opt_match <= 0)
    fatal("The argument to --match must be positive");

  if (opt_mismatch >= 0)
    fatal("The argument to --mismatch must be negative");

#endif


  if (opt_alignwidth < 0)
    fatal("The argument to --alignwidth must not be negative");

  if (opt_rowlen < 0)
    fatal("The argument to --rowlen must not be negative");

  if (opt_qmask == MASK_ERROR)
    fatal("The argument to --qmask must be none, dust or soft");

  if (opt_dbmask == MASK_ERROR)
    fatal("The argument to --dbmask must be none, dust or soft");

  if ((opt_sample_pct < 0.0) || (opt_sample_pct > 100.0))
    fatal("The argument to --sample_pct must be in the range 0.0 to 100.0");

  if (opt_sample_size < 0)
    fatal("The argument to --sample_size must not be negative");

  if (((opt_relabel ? 1 : 0) +
       opt_relabel_md5 + opt_relabel_self + opt_relabel_sha1) > 1)
    fatal("Specify only one of --relabel, --relabel_self, --relabel_sha1, or --relabel_md5");

  if (opt_fastq_tail < 1)
    fatal("The argument to --fastq_tail must be positive");

  if ((opt_min_unmasked_pct < 0.0) && (opt_min_unmasked_pct > 100.0))
    fatal("The argument to --min_unmasked_pct must be between 0.0 and 100.0");

  if ((opt_max_unmasked_pct < 0.0) && (opt_max_unmasked_pct > 100.0))
    fatal("The argument to --max_unmasked_pct must be between 0.0 and 100.0");

  if (opt_min_unmasked_pct > opt_max_unmasked_pct)
    fatal("The argument to --min_unmasked_pct cannot be larger than to --max_unmasked_pct");

  if ((opt_fastq_ascii != 33) && (opt_fastq_ascii != 64))
    fatal("The argument to --fastq_ascii must be 33 or 64");

  if (opt_fastq_qmin > opt_fastq_qmax)
    fatal("The argument to --fastq_qmin cannot be larger than to --fastq_qmax");

  if (opt_fastq_ascii + opt_fastq_qmin < 33)
    fatal("Sum of arguments to --fastq_ascii and --fastq_qmin must be no less than 33");

  if (opt_fastq_ascii + opt_fastq_qmax > 126)
    fatal("Sum of arguments to --fastq_ascii and --fastq_qmax must be no more than 126");

  if (opt_fastq_qminout > opt_fastq_qmaxout)
    fatal("The argument to --fastq_qminout cannot be larger than to --fastq_qmaxout");

  if ((opt_fastq_asciiout != 33) && (opt_fastq_asciiout != 64))
    fatal("The argument to --fastq_asciiout must be 33 or 64");

  if (opt_fastq_asciiout + opt_fastq_qminout < 33)
    fatal("Sum of arguments to --fastq_asciiout and --fastq_qminout must be no less than 33");

  if (opt_fastq_asciiout + opt_fastq_qmaxout > 126)
    fatal("Sum of arguments to --fastq_asciiout and --fastq_qmaxout must be no more than 126");

  if (opt_gzip_decompress && opt_bzip2_decompress)
    fatal("Specify either --gzip_decompress or --bzip2_decompress, not both");

  if ((opt_sintax_cutoff < 0.0) || (opt_sintax_cutoff > 1.0))
    fatal("The argument to sintax_cutoff must be in the range 0.0 to 1.0");

  if (opt_minuniquesize < 1)
    fatal("The argument to minuniquesize must be at least 1");

  if (opt_maxuniquesize < 1)
    fatal("The argument to maxuniquesize must be at least 1");

  if (opt_maxsize < 1)
    fatal("The argument to maxsize must be at least 1");

  if (opt_maxhits < 0)
    fatal("The argument to maxhits cannot be negative");


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

  /* set defaults parameters, if not specified */

  if (opt_maxhits == 0)
    opt_maxhits = LONG_MAX;

  if (opt_minwordmatches < 0)
    opt_minwordmatches = minwordmatches_defaults[opt_wordlength];

  /* set default opt_minsize depending on command */
  if (opt_minsize == 0)
    {
      if (opt_cluster_unoise)
        opt_minsize = 8;
      else
        opt_minsize = 1;
    }

  /* set default opt_abskew depending on command */
  if (opt_abskew < 0.0)
    {
      if (opt_uchime3_denovo)
        opt_abskew = 16.0;
      else
        opt_abskew = 2.0;
    }

  /* set default opt_minseqlength depending on command */

  if (opt_minseqlength < 0)
    {
      if (opt_cluster_smallmem || opt_cluster_fast || opt_cluster_size ||
          opt_usearch_global || opt_derep_fulllength || opt_derep_prefix ||
          opt_makeudb_usearch || opt_cluster_unoise || opt_sintax)
        opt_minseqlength = 32;
      else
        opt_minseqlength = 1;
    }

  if (opt_sintax)
    opt_notrunclabels = 1;
}

void show_publication()
{
  fprintf(stdout,
          "Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)\n"
          "VSEARCH: a versatile open source tool for metagenomics\n"
          "PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584\n"
          "\n");
}

void cmd_version()
{
  if (! opt_quiet)
    {
      show_publication();

#ifdef HAVE_ZLIB_H
      printf("Compiled with support for gzip-compressed files,");
      if (gz_lib)
        printf(" and the library is loaded.\n");
      else
        printf(" but the library was not found.\n");
#else
      printf("Compiled without support for gzip-compressed files.\n");
#endif

#ifdef HAVE_BZLIB_H
      printf("Compiled with support for bzip2-compressed files,");
      if (bz2_lib)
        printf(" and the library is loaded.\n");
      else
        printf(" but the library was not found.\n");
#else
      printf("Compiled without support for bzip2-compressed files.\n");
#endif
    }
}

void cmd_help()
{
  /*       0         1         2         3         4         5         6         7          */
  /*       01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  if (! opt_quiet)
    {
      show_publication();

      fprintf(stdout,
              "Usage: %s [OPTIONS]\n", progname);

      fprintf(stdout,
              "\n"
              "General options\n"
              "  --bzip2_decompress          decompress input with bzip2 (required if pipe)\n"
              "  --fasta_width INT           width of FASTA seq lines, 0 for no wrap (80)\n"
              "  --gzip_decompress           decompress input with gzip (required if pipe)\n"
              "  --help | -h                 display help information\n"
              "  --log FILENAME              write messages, timing and memory info to file\n"
              "  --maxseqlength INT          maximum sequence length (50000)\n"
              "  --minseqlength INT          min seq length (clust/derep/search: 32, other:1)\n"
              "  --no_progress               do not show progress indicator\n"
              "  --notrunclabels             do not truncate labels at first space\n"
              "  --quiet                     output just warnings and fatal errors to stderr\n"
              "  --threads INT               number of threads to use, zero for all cores (0)\n"
              "  --version | -v              display version information\n"
              "\n"
              "Chimera detection\n"
              "  --uchime_denovo FILENAME    detect chimeras de novo\n"
              "  --uchime2_denovo FILENAME   detect chimeras de novo in denoised amplicons\n"
              "  --uchime3_denovo FILENAME   detect chimeras de novo in denoised amplicons\n"
              "  --uchime_ref FILENAME       detect chimeras using a reference database\n"
              " Data\n"
              "  --db FILENAME               reference database for --uchime_ref\n"
              " Parameters\n"
              "  --abskew REAL               minimum abundance ratio (2.0, 16.0 for uchime3)\n"
              "  --dn REAL                   'no' vote pseudo-count (1.4)\n"
              "  --mindiffs INT              minimum number of differences in segment (3) *\n"
              "  --mindiv REAL               minimum divergence from closest parent (0.8) *\n"
              "  --minh REAL                 minimum score (0.28) * ignored in uchime2/3\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --self                      exclude identical labels for --uchime_ref\n"
              "  --selfid                    exclude identical sequences for --uchime_ref\n"
              "  --xn REAL                   'no' vote weight (8.0)\n"
              " Output\n"
              "  --alignwidth INT            width of alignment in uchimealn output (80)\n"
              "  --borderline FILENAME       output borderline chimeric sequences to file\n"
              "  --chimeras FILENAME         output chimeric sequences to file\n"
              "  --fasta_score               include chimera score in fasta output\n"
              "  --nonchimeras FILENAME      output non-chimeric sequences to file\n"
              "  --relabel STRING            relabel nonchimeras with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_self              relabel with the sequence itself as label\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizeout                   include abundance information when relabelling\n"
              "  --uchimealns FILENAME       output chimera alignments to file\n"
              "  --uchimeout FILENAME        output to chimera info to tab-separated file\n"
              "  --uchimeout5                make output compatible with uchime version 5\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "Clustering\n"
              "  --cluster_fast FILENAME     cluster sequences after sorting by length\n"
              "  --cluster_size FILENAME     cluster sequences after sorting by abundance\n"
              "  --cluster_smallmem FILENAME cluster already sorted sequences (see -usersort)\n"
              "  --cluster_unoise FILENAME   denoise Illumina amplicon reads\n"
              " Parameters (most searching options also apply)\n"
              "  --cons_truncate             do not ignore terminal gaps in MSA for consensus\n"
              "  --id REAL                   reject if identity lower, accepted values: 0-1.0\n"
              "  --iddef INT                 id definition, 0-4=CD-HIT,all,int,MBL,BLAST (2)\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --strand plus|both          cluster using plus or both strands (plus)\n"
              "  --usersort                  indicate sequences not pre-sorted by length\n"
              "  --minsize INT               minimum abundance (unoise only) (8)\n"
              "  --unoise_alpha REAL         alpha parameter (unoise only) (2.0)\n"
              " Output\n"
              "  --biomout FILENAME          filename for OTU table output in biom 1.0 format\n"
              "  --centroids FILENAME        output centroid sequences to FASTA file\n"
              "  --clusterout_id             add cluster id info to consout and profile files\n"
              "  --clusterout_sort           order msaout, consout, profile by decr abundance\n"
              "  --clusters STRING           output each cluster to a separate FASTA file\n"
              "  --consout FILENAME          output cluster consensus sequences to FASTA file\n"
              "  --mothur_shared_out FN      filename for OTU table output in mothur format\n"
              "  --msaout FILENAME           output multiple seq. alignments to FASTA file\n"
              "  --otutabout FILENAME        filename for OTU table output in classic format\n"
              "  --profile FILENAME          output sequence profile of each cluster to file\n"
              "  --relabel STRING            relabel centroids with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_self              relabel with the sequence itself as label\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizeorder                 sort accepted centroids by abundance, AGC\n"
              "  --sizeout                   write cluster abundances to centroid file\n"
              "  --uc FILENAME               specify filename for UCLUST-like output\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "Convert SFF to FASTQ\n"
              "  --sff_convert FILENAME      convert given SFF file to FASTQ format\n"
              " Parameters\n"
              "  --sff_clip                  clip ends of sequences as indicated in file (no)\n"
              "  --fastq_asciiout INT        FASTQ output quality score ASCII base char (33)\n"
              "  --fastq_qmaxout INT         maximum base quality value for FASTQ output (41)\n"
              "  --fastq_qminout INT         minimum base quality value for FASTQ output (0)\n"
              " Output\n"
              "  --fastqout FILENAME         output converted sequences to given FASTQ file\n"
              "\n"
              "Dereplication and rereplication\n"
              "  --derep_fulllength FILENAME dereplicate sequences in the given FASTA file\n"
              "  --derep_prefix FILENAME     dereplicate sequences in file based on prefixes\n"
              "  --rereplicate FILENAME      rereplicate sequences in the given FASTA file\n"
              " Parameters\n"
              "  --maxuniquesize INT         maximum abundance for output from dereplication\n"
              "  --minuniquesize INT         minimum abundance for output from dereplication\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --strand plus|both          dereplicate plus or both strands (plus)\n"
              " Output\n"
              "  --output FILENAME           output FASTA file\n"
              "  --relabel STRING            relabel with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_self              relabel with the sequence itself as label\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizeout                   write abundance annotation to output\n"
              "  --topn INT                  output only n most abundant sequences after derep\n"
              "  --uc FILENAME               filename for UCLUST-like dereplication output\n"
              "  --xsize                     strip abundance information in derep output\n"
              "\n"
              "FASTQ format conversion\n"
              "  --fastq_convert FILENAME    convert between FASTQ file formats\n"
              " Parameters\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_asciiout INT        FASTQ output quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmaxout INT         maximum base quality value for FASTQ output (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastq_qminout INT         minimum base quality value for FASTQ output (0)\n"
              " Output\n"
              "  --fastqout FILENAME         FASTQ output filename for converted sequences\n"
              "\n"
              "FASTQ format detection and quality analysis\n"
              "  --fastq_chars FILENAME      analyse FASTQ file for version and quality range\n"
              " Parameters\n"
              "  --fastq_tail INT            min length of tails to count for fastq_chars (4)\n"
              "\n"
              "FASTQ quality statistics\n"
              "  --fastq_stats FILENAME      report statistics on FASTQ file\n"
              "  --fastq_eestats FILENAME    quality score and expected error statistics\n"
              "  --fastq_eestats2 FILENAME   expected error and length cutoff statistics\n"
              " Parameters\n"
              "  --ee_cutoffs REAL,...       fastq_eestats2 expected error cutoffs (0.5,1,2)\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --length_cutoffs INT,INT,INT fastq_eestats2 length (min,max,incr) (50,*,50)\n"
              " Output\n"
              "  --log FILENAME              output file for fastq_stats statistics\n"
              "  --output FILENAME           output file for fastq_eestats(2) statistics\n"
              "\n"
              "Masking (new)\n"
              "  --fastx_mask FILENAME       mask sequences in the given FASTA or FASTQ file\n"
              " Parameters\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --max_unmasked_pct          max unmasked %% of sequences to keep (100.0)\n"
              "  --min_unmasked_pct          min unmasked %% of sequences to keep (0.0)\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              " Output\n"
              "  --fastaout FILENAME         output to specified FASTA file\n"
              "  --fastqout FILENAME         output to specified FASTQ file\n"
              "\n"
              "Masking (old)\n"
              "  --maskfasta FILENAME        mask sequences in the given FASTA file\n"
              " Parameters\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              " Output\n"
              "  --output FILENAME           output to specified FASTA file\n"
              "\n"
              "Paired-end reads joining\n"
              "  --fastq_join FILENAME       join paired-end reads into one sequence with gap\n"
              " Data\n"
              "  --reverse FILENAME          specify FASTQ file with reverse reads\n"
              "  --join_padgap STRING        sequence string used for padding (NNNNNNNN)\n"
              "  --join_padgapq STRING       quality string used for padding (IIIIIIII)\n"
              " Output\n"
              "  --fastaout FILENAME         FASTA output filename for joined sequences\n"
              "  --fastqout FILENAME         FASTQ output filename for joined sequences\n"
              "\n"
              "Paired-end reads merging\n"
              "  --fastq_mergepairs FILENAME merge paired-end reads into one sequence\n"
              " Data\n"
              "  --reverse FILENAME          specify FASTQ file with reverse reads\n"
              " Parameters\n"
              "  --fastq_allowmergestagger   allow merging of staggered reads\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_maxdiffpct REAL     maximum percentage diff. bases in overlap (100.0)\n"
              "  --fastq_maxdiffs INT        maximum number of different bases in overlap (10)\n"
              "  --fastq_maxee REAL          maximum expected error value for merged sequence\n"
              "  --fastq_maxmergelen         maximum length of entire merged sequence\n"
              "  --fastq_maxns INT           maximum number of N's\n"
              "  --fastq_minlen INT          minimum input read length after truncation (1)\n"
              "  --fastq_minmergelen         minimum length of entire merged sequence\n"
              "  --fastq_minovlen            minimum length of overlap between reads (10)\n"
              "  --fastq_nostagger           disallow merging of staggered reads (default)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmaxout INT         maximum base quality value for FASTQ output (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastq_qminout INT         minimum base quality value for FASTQ output (0)\n"
              "  --fastq_truncqual INT       base quality value for truncation\n"
              " Output\n"
              "  --eetabbedout FILENAME      output error statistics to specified file\n"
              "  --fastaout FILENAME         FASTA output filename for merged sequences\n"
              "  --fastaout_notmerged_fwd FN FASTA filename for non-merged forward sequences\n"
              "  --fastaout_notmerged_rev FN FASTA filename for non-merged reverse sequences\n"
              "  --fastq_eeout               include expected errors (ee) in FASTQ output\n"
              "  --fastqout FILENAME         FASTQ output filename for merged sequences\n"
              "  --fastqout_notmerged_fwd FN FASTQ filename for non-merged forward sequences\n"
              "  --fastqout_notmerged_rev FN FASTQ filename for non-merged reverse sequences\n"
              "  --label_suffix              suffix to append to label of merged sequences\n"
              "  --xee                       remove expected errors (ee) info from output\n"
              "\n"
              "Pairwise alignment\n"
              "  --allpairs_global FILENAME  perform global alignment of all sequence pairs\n"
              " Output (most searching options also apply)\n"
              "  --alnout FILENAME           filename for human-readable alignment output\n"
              "  --acceptall                 output all pairwise alignments\n"
              "\n"
              "Restriction site cutting\n"
              "  --cut FILENAME              filename of FASTA formatted input sequences\n"
              " Parameters\n"
              "  --cut_pattern STRING        pattern to match with ^ and _ at cut sites\n"
              " Output\n"
              "  --fastaout FILENAME         FASTA filename for fragments on forward strand\n"
              "  --fastaout_rev FILENAME     FASTA filename for fragments on reverse strand\n"
              "  --fastaout_discarded FN     FASTA filename for non-matching sequences\n"
              "  --fastaout_discarded_rev FN FASTA filename for non-matching, reverse compl.\n"
              "\n"
              "Reverse complementation\n"
              "  --fastx_revcomp FILENAME    reverse-complement seqs in FASTA or FASTQ file\n"
              " Parameters\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              " Output\n"
              "  --fastaout FILENAME         FASTA output filename\n"
              "  --fastqout FILENAME         FASTQ output filename\n"
              "  --label_suffix STRING       label to append to identifier in the output\n"
              "\n"
              "Searching\n"
              "  --search_exact FILENAME     filename of queries for exact match search\n"
              "  --usearch_global FILENAME   filename of queries for global alignment search\n"
              " Data\n"
              "  --db FILENAME               name of UDB or FASTA database for search\n"
              " Parameters\n"
              "  --dbmask none|dust|soft     mask db with dust, soft or no method (dust)\n"
              "  --fulldp                    full dynamic programming alignment (always on)\n"
              "  --gapext STRING             penalties for gap extension (2I/1E)\n"
              "  --gapopen STRING            penalties for gap opening (20I/2E)\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --id REAL                   reject if identity lower\n"
              "  --iddef INT                 id definition, 0-4=CD-HIT,all,int,MBL,BLAST (2)\n"
              "  --idprefix INT              reject if first n nucleotides do not match\n"
              "  --idsuffix INT              reject if last n nucleotides do not match\n"
              "  --leftjust                  reject if terminal gaps at alignment left end\n"
              "  --match INT                 score for match (2)\n"
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
              "  --minwordmatches INT        minimum number of word matches required (12)\n"
              "  --mismatch INT              score for mismatch (-4)\n"
              "  --pattern STRING            option is ignored\n"
              "  --qmask none|dust|soft      mask query with dust, soft or no method (dust)\n"
              "  --query_cov REAL            reject if fraction of query seq. aligned lower\n"
              "  --rightjust                 reject if terminal gaps at alignment right end\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --self                      reject if labels identical\n"
              "  --selfid                    reject if sequences identical\n"
              "  --slots INT                 option is ignored\n"
              "  --strand plus|both          search plus or both strands (plus)\n"
              "  --target_cov REAL           reject if fraction of target seq. aligned lower\n"
              "  --weak_id REAL              include aligned hits with >= id; continue search\n"
              "  --wordlength INT            length of words for database index 3-15 (8)\n"
              " Output\n"
              "  --alnout FILENAME           filename for human-readable alignment output\n"
              "  --biomout FILENAME          filename for OTU table output in biom 1.0 format\n"
              "  --blast6out FILENAME        filename for blast-like tab-separated output\n"
              "  --dbmatched FILENAME        FASTA file for matching database sequences\n"
              "  --dbnotmatched FILENAME     FASTA file for non-matching database sequences\n"
              "  --fastapairs FILENAME       FASTA file with pairs of query and target\n"
              "  --matched FILENAME          FASTA file for matching query sequences\n"
              "  --mothur_shared_out FN      filename for OTU table output in mothur format\n"
              "  --notmatched FILENAME       FASTA file for non-matching query sequences\n"
              "  --otutabout FILENAME        filename for OTU table output in classic format\n"
              "  --output_no_hits            output non-matching queries to output files\n"
              "  --rowlen INT                width of alignment lines in alnout output (64)\n"
              "  --samheader                 include a header in the SAM output file\n"
              "  --samout FILENAME           filename for SAM format output\n"
              "  --sizeout                   write abundance annotation to dbmatched file\n"
              "  --top_hits_only             output only hits with identity equal to the best\n"
              "  --uc FILENAME               filename for UCLUST-like output\n"
              "  --uc_allhits                show all, not just top hit with uc output\n"
              "  --userfields STRING         fields to output in userout file\n"
              "  --userout FILENAME          filename for user-defined tab-separated output\n"
              "\n"
              "Shuffling and sorting\n"
              "  --shuffle FILENAME          shuffle order of sequences in FASTA file randomly\n"
              "  --sortbylength FILENAME     sort sequences by length in given FASTA file\n"
              "  --sortbysize FILENAME       abundance sort sequences in given FASTA file\n"
              " Parameters\n"
              "  --maxsize INT               maximum abundance for sortbysize\n"
              "  --minsize INT               minimum abundance for sortbysize\n"
              "  --randseed INT              seed for PRNG, zero to use random data source (0)\n"
              "  --sizein                    propagate abundance annotation from input\n"
              " Output\n"
              "  --output FILENAME           output to specified FASTA file\n"
              "  --relabel STRING            relabel sequences with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_self              relabel with the sequence itself as label\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizeout                   include abundance information when relabelling\n"
              "  --topn INT                  output just first n sequences\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "Subsampling\n"
              "  --fastx_subsample FILENAME  subsample sequences from given FASTA/FASTQ file\n"
              " Parameters\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --randseed INT              seed for PRNG, zero to use random data source (0)\n"
              "  --sample_pct REAL           sampling percentage between 0.0 and 100.0\n"
              "  --sample_size INT           sampling size\n"
              "  --sizein                    consider abundance info from input, do not ignore\n"
              " Output\n"
              "  --fastaout FILENAME         output subsampled sequences to FASTA file\n"
              "  --fastaout_discarded FILE   output non-subsampled sequences to FASTA file\n"
              "  --fastqout FILENAME         output subsampled sequences to FASTQ file\n"
              "  --fastqout_discarded        output non-subsampled sequences to FASTQ file\n"
              "  --relabel STRING            relabel sequences with this prefix string\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
              "  --relabel_self              relabel with the sequence itself as label\n"
              "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
              "  --sizeout                   update abundance information in output\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "Taxonomic classification\n"
              "  --sintax FILENAME           classify sequences in given FASTA/FASTQ file\n"
              " Parameters\n"
              "  --db FILENAME               taxonomic reference db in given FASTA or UDB file\n"
              "  --sintax_cutoff REAL        confidence value cutoff level (0.0)\n"
              " Output\n"
              "  --tabbedout FILENAME        write results to given tab-delimited file\n"
              "\n"
              "Trimming and filtering\n"
              "  --fastx_filter FILENAME     trim and filter sequences in FASTA/FASTQ file\n"
              "  --fastq_filter FILENAME     trim and filter sequences in FASTQ file\n"
              "  --reverse FILENAME          FASTQ file with other end of paired-end reads\n"
              " Parameters\n"
              "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
              "  --fastq_maxee REAL          discard if expected error value is higher\n"
              "  --fastq_maxee_rate REAL     discard if expected error rate is higher\n"
              "  --fastq_maxlen INT          discard if length of sequence is longer\n"
              "  --fastq_maxns INT           discard if number of N's is higher\n"
              "  --fastq_minlen INT          discard if length of sequence is shorter\n"
              "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
              "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
              "  --fastq_stripleft INT       delete given number of bases from the 5' end\n"
              "  --fastq_stripright INT      delete given number of bases from the 3' end\n"
              "  --fastq_truncee REAL        truncate to given maximum expected error\n"
              "  --fastq_trunclen INT        truncate to given length (discard if shorter)\n"
              "  --fastq_trunclen_keep INT   truncate to given length (keep if shorter)\n"
              "  --fastq_truncqual INT       truncate to given minimum base quality\n"
              "  --maxsize INT               discard if abundance of sequence is above\n"
              "  --minsize INT               discard if abundance of sequence is below\n"
              " Output\n"
              "  --eeout                     include expected errors in output\n"
              "  --fastaout FN               FASTA filename for passed sequences\n"
              "  --fastaout_discarded FN     FASTA filename for discarded sequences\n"
              "  --fastaout_discarded_rev FN FASTA filename for discarded reverse sequences\n"
              "  --fastaout_rev FN           FASTA filename for passed reverse sequences\n"
              "  --fastqout FN               FASTQ filename for passed sequences\n"
              "  --fastqout_discarded FN     FASTQ filename for discarded sequences\n"
              "  --fastqout_discarded_rev FN FASTQ filename for discarded reverse sequences\n"
              "  --fastqout_rev FN           FASTQ filename for passed reverse sequences\n"
              "  --relabel STRING            relabel filtered sequences with given prefix\n"
              "  --relabel_keep              keep the old label after the new when relabelling\n"
              "  --relabel_md5               relabel filtered sequences with md5 digest\n"
              "  --relabel_self              relabel with the sequence itself as label\n"
              "  --relabel_sha1              relabel filtered sequences with sha1 digest\n"
              "  --sizeout                   include abundance information when relabelling\n"
              "  --xee                       remove expected errors (ee) info from output\n"
              "  --xsize                     strip abundance information in output\n"
              "\n"
              "UDB files\n"
              "  --makeudb_usearch FILENAME  make UDB file from given FASTA file\n"
              "  --udb2fasta FILENAME        output FASTA file from given UDB file\n"
              "  --udbinfo FILENAME          show information about UDB file\n"
              "  --udbstats FILENAME         report statistics about indexed words in UDB file\n"
              " Parameters\n"
              "  --dbmask none|dust|soft     mask db with dust, soft or no method (dust)\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --wordlength INT            length of words for database index 3-15 (8)\n"
              " Output\n"
              "  --output FILENAME           UDB or FASTA output file\n"
          );
    }
}

void cmd_allpairs_global()
{
  /* check options */

  if ((!opt_alnout) && (!opt_userout) &&
      (!opt_uc) && (!opt_blast6out) &&
      (!opt_matched) && (!opt_notmatched) &&
      (!opt_samout) && (!opt_fastapairs))
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
      (!opt_samout) && (!opt_otutabout) &&
      (!opt_biomout) && (!opt_mothur_shared_out) &&
      (!opt_fastapairs))
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
      (!opt_samout) && (!opt_otutabout) &&
      (!opt_biomout) && (!opt_mothur_shared_out) &&
      (!opt_fastapairs))
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

void cmd_rereplicate()
{
  if (!opt_output)
    fatal("FASTA output file for rereplicate must be specified with --output");

  rereplicate();
}

void cmd_derep()
{
  if ((!opt_output) && (!opt_uc))
    fatal("Output file for dereplication must be specified with --output or --uc");

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

void cmd_fastq_eestats()
{
  if (!opt_output)
    fatal("Output file for fastq_eestats must be specified with --output");

  fastq_eestats();
}

void cmd_fastq_eestats2()
{
  if (!opt_output)
    fatal("Output file for fastq_eestats2 must be specified with --output");

  fastq_eestats2();
}

void cmd_subsample()
{
  if ((!opt_fastaout) && (!opt_fastqout))
    fatal("Specify output files for subsampling with --fastaout and/or --fastqout");

  if ((opt_sample_pct > 0) == (opt_sample_size > 0))
    fatal("Specify either --sample_pct or --sample_size");

  subsample();
}

void cmd_maskfasta()
{
  if (!opt_output)
    fatal("Output file for masking must be specified with --output");

  maskfasta();
}

void cmd_makeudb_usearch()
{
  if (!opt_output)
    fatal("UDB output file must be specified with --output");
  udb_make();
}

void cmd_udb2fasta()
{
  if (!opt_output)
    fatal("FASTA output file must be specified with --output");
  udb_fasta();
}

void cmd_fastx_mask()
{
  if ((!opt_fastaout) && (!opt_fastqout))
    fatal("Specify output files for masking with --fastaout and/or --fastqout");

  fastx_mask();
}

void cmd_none()
{
  if (! opt_quiet)
    fprintf(stderr,
            "For help, please enter: %s --help | less\n"
            "For further details, please consult the manual by entering: man vsearch\n"
            "\n"
            "Selected commands:\n"
            "\n"
            "vsearch --allpairs_global FILENAME --id 0.5 --alnout FILENAME\n"
            "vsearch --cluster_size FILENAME --id 0.97 --centroids FILENAME\n"
            "vsearch --cut FILENAME --cut_pattern G^AATT_C --fastaout FILENAME\n"
            "vsearch --derep_fulllength FILENAME --output FILENAME\n"
            "vsearch --fastq_chars FILENAME\n"
            "vsearch --fastq_convert FILENAME --fastqout FILENAME --fastq_ascii 64\n"
            "vsearch --fastq_eestats FILENAME --output FILENAME\n"
            "vsearch --fastq_eestats2 FILENAME --output FILENAME\n"
            "vsearch --fastq_mergepairs FILENAME --reverse FILENAME --fastqout FILENAME\n"
            "vsearch --fastq_stats FILENAME --log FILENAME\n"
            "vsearch --fastx_filter FILENAME --fastaout FILENAME --fastq_trunclen 100\n"
            "vsearch --fastx_getseq FILENAME --label LABEL --fastaout FILENAME\n"
            "vsearch --fastx_mask FILENAME --fastaout FILENAME\n"
            "vsearch --fastx_revcomp FILENAME --fastqout FILENAME\n"
            "vsearch --fastx_subsample FILENAME --fastaout FILENAME --sample_pct 1\n"
            "vsearch --makeudb_usearch FILENAME --output FILENAME\n"
            "vsearch --search_exact FILENAME --db FILENAME --alnout FILENAME\n"
            "vsearch --sff_convert FILENAME --output FILENAME --sff_clip\n"
            "vsearch --shuffle FILENAME --output FILENAME\n"
            "vsearch --sintax FILENAME --db FILENAME --tabbedout FILENAME\n"
            "vsearch --sortbylength FILENAME --output FILENAME\n"
            "vsearch --sortbysize FILENAME --output FILENAME\n"
            "vsearch --uchime_denovo FILENAME --nonchimeras FILENAME\n"
            "vsearch --uchime_ref FILENAME --db FILENAME --nonchimeras FILENAME\n"
            "vsearch --usearch_global FILENAME --db FILENAME --id 0.97 --alnout FILENAME\n"
            "\n"
            "Other commands: cluster_fast, cluster_smallmem, cluster_unoise, derep_prefix,\n"
            "                fastq_filter, fastq_join, fastx_getseqs, fastx_getsubseqs,\n"
            "                maskfasta, rereplicate, uchime2_denovo, uchime3_denovo,\n"
            "                udb2fasta, udbinfo, udbstats, version\n"
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
      (!opt_samout) && (!opt_profile) &&
      (!opt_otutabout) && (!opt_biomout) &&
      (!opt_mothur_shared_out))
    fatal("No output files specified");

  if (!opt_cluster_unoise)
    if ((opt_id < 0.0) || (opt_id > 1.0))
      fatal("Identity between 0.0 and 1.0 must be specified with --id");

  if (opt_cluster_fast)
    cluster_fast(cmdline, progheader);
  else if (opt_cluster_smallmem)
    cluster_smallmem(cmdline, progheader);
  else if (opt_cluster_size)
    cluster_size(cmdline, progheader);
  else if (opt_cluster_unoise)
    cluster_unoise(cmdline, progheader);
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

  if ((!opt_uchime2_denovo) && (!opt_uchime3_denovo))
  {
    if (opt_mindiffs <= 0)
      fatal("Argument to --mindiffs must be > 0");

    if (opt_mindiv <= 0.0)
      fatal("Argument to --mindiv must be > 0");

    if (opt_minh <= 0.0)
      fatal("Argument to --minh must be > 0");
  }

#if 0
  if (opt_abskew <= 1.0)
    fatal("Argument to --abskew must be > 1");
#endif

  chimera();
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
           arch_get_cores());
}

void getentirecommandline(int argc, char** argv)
{
  int len = 0;
  for (int i=0; i<argc; i++)
    len += strlen(argv[i]);

  cmdline = (char*) xmalloc(len+argc);
  cmdline[0] = 0;

  for (int i=0; i<argc; i++)
    {
      if (i>0)
        strcat(cmdline, " ");
      strcat(cmdline, argv[i]);
    }
}

void show_header()
{
  if (! opt_quiet)
    {
      fprintf(stderr, "%s\n", progheader);
      fprintf(stderr, "https://github.com/torognes/vsearch\n");
      fprintf(stderr, "\n");
    }
}

int main(int argc, char** argv)
{
  fillheader();

  getentirecommandline(argc, argv);

  cpu_features_detect();

  args_init(argc, argv);

  if (opt_log)
    {
      fp_log = fopen_output(opt_log);
      if (!fp_log)
        fatal("Unable to open log file for writing");
      fprintf(fp_log, "%s\n", progheader);
      fprintf(fp_log, "%s\n", cmdline);

      char time_string[26];
      time_start = time(0);
      struct tm tm_start;
      localtime_r(& time_start, & tm_start);
      strftime(time_string, 26, "%c", & tm_start);
      fprintf(fp_log, "Started  %s\n", time_string);
    }

  random_init();

  show_header();

  dynlibs_open();

#ifdef __x86_64__
  if (!sse2_present)
    fatal("Sorry, this program requires a cpu with SSE2.");
#endif

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
  else if (opt_cluster_smallmem || opt_cluster_fast || opt_cluster_size || opt_cluster_unoise)
    cmd_cluster();
  else if (opt_uchime_denovo || opt_uchime_ref || opt_uchime2_denovo || opt_uchime3_denovo)
    cmd_uchime();
  else if (opt_fastq_chars)
    fastq_chars();
  else if (opt_fastq_stats)
    fastq_stats();
  else if (opt_fastq_filter)
    fastq_filter();
  else if (opt_fastx_filter)
    fastx_filter();
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
  else if (opt_fastq_eestats)
    cmd_fastq_eestats();
  else if (opt_fastq_eestats2)
    cmd_fastq_eestats2();
  else if (opt_fastq_join)
    fastq_join();
  else if (opt_rereplicate)
    cmd_rereplicate();
  else if (opt_version)
    cmd_version();
  else if (opt_makeudb_usearch)
    cmd_makeudb_usearch();
  else if (opt_udb2fasta)
    cmd_udb2fasta();
  else if (opt_udbinfo)
    udb_info();
  else if (opt_udbstats)
    udb_stats();
  else if (opt_sintax)
    sintax();
  else if (opt_sff_convert)
    sff_convert();
  else if (opt_fastx_getseq)
    fastx_getseq();
  else if (opt_fastx_getseqs)
    fastx_getseqs();
  else if (opt_fastx_getsubseq)
    fastx_getsubseq();
  else if (opt_cut)
    cut();
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

      double time_diff = difftime(time_finish, time_start);
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Elapsed time %02.0lf:%02.0lf\n",
              floor(time_diff / 60.0),
              floor(time_diff - 60.0 * floor(time_diff / 60.0)));
      double maxmem = arch_get_memused() / 1048576.0;
      if (maxmem < 1024.0)
        fprintf(fp_log, "Max memory %.1lfMB\n", maxmem);
      else
        fprintf(fp_log, "Max memory %.1lfGB\n", maxmem/1024.0);
      fclose(fp_log);
    }

  if (opt_ee_cutoffs_values)
    xfree(opt_ee_cutoffs_values);
  opt_ee_cutoffs_values = 0;

  xfree(cmdline);
  dynlibs_close();
}
