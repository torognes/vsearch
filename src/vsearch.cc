/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "allpairs.h"
#include "chimera.h"
#include "cluster.h"
#include "cut.h"
#include "derep.h"
#include "derep_prefix.h"
#include "derep_smallmem.h"
#include "dynlibs.h"
#include "eestats.h"
#include "fasta2fastq.h"
#include "fastq_chars.h"
#include "fastq_join.h"
#include "fastq_stats.h"
#include "fastqops.h"
#include "filter.h"
#include "getseq.h"
#include "mask.h"
#include "mergepairs.h"
#include "orient.h"
#include "rereplicate.h"
#include "search.h"
#include "search_exact.h"
#include "sff_convert.h"
#include "shuffle.h"
#include "sintax.h"
#include "sortbylength.h"
#include "sortbysize.h"
#include "subsample.h"
#include "udb.h"
#include "userfields.h"
#include <algorithm>  // std::count
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::floor
#include <ctime>  // std::strftime, std::localtime, std::time, std::time_t, std::tm, std::difftime
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t, std::sscanf, std::fclose, std::snprintf, std::printf, std::strcat
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::strlen, std::memset
#include <getopt.h>  // getopt_long_only, optarg, optind, opterr, struct
                     // option (no_argument, required_argument)
#include <limits>
#include <string.h>  // strcasecmp
#include <string>
#include <vector>


constexpr int64_t n_threads_max = 1024;
constexpr auto max_line_length = std::size_t{80};

/* options */

bool opt_bzip2_decompress = false;
bool opt_clusterout_id = false;
bool opt_clusterout_sort = false;
bool opt_eeout;
bool opt_fasta_score;
bool opt_fastq_allowmergestagger;
bool opt_fastq_eeout;
bool opt_fastq_nostagger;
bool opt_gzip_decompress;
bool opt_label_substr_match;
bool opt_lengthout;
bool opt_n_mismatch;
bool opt_no_progress;
bool opt_quiet;
bool opt_relabel_keep;
bool opt_relabel_md5;
bool opt_relabel_self;
bool opt_relabel_sha1;
bool opt_samheader;
bool opt_sintax_random;
bool opt_sizein;
bool opt_sizeorder;
bool opt_sizeout;
bool opt_xee;
bool opt_xlength;
bool opt_xsize;
char * opt_alnout;
char * opt_biomout;
char * opt_blast6out;
char * opt_borderline;
char * opt_centroids;
char * opt_chimeras;
char * opt_chimeras_alnout;
char * opt_chimeras_denovo;
char * opt_cluster_fast;
char * opt_cluster_size;
char * opt_cluster_smallmem;
char * opt_cluster_unoise;
char * opt_clusters;
char * opt_consout;
char * opt_db;
char * opt_dbmatched;
char * opt_dbnotmatched;
char * opt_eetabbedout;
char * opt_fastaout;
char * opt_fastaout_discarded;
char * opt_fastaout_discarded_rev;
char * opt_fastaout_notmerged_fwd;
char * opt_fastaout_notmerged_rev;
char * opt_fastaout_rev;
char * opt_fastapairs;
char * opt_fastq_convert;
char * opt_fastqout;
char * opt_fastqout_discarded;
char * opt_fastqout_discarded_rev;
char * opt_fastqout_notmerged_fwd;
char * opt_fastqout_notmerged_rev;
char * opt_fastqout_rev;
char * opt_label;
char * opt_label_field;
char * opt_label_suffix;
char * opt_label_word;
char * opt_label_words;
char * opt_labels;
char * opt_lcaout;
char * opt_log;
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
char * opt_qsegout;
char * opt_relabel;
char * opt_reverse;
char * opt_samout;
char * opt_sample;
char * opt_tabbedout;
char * opt_tsegout;
char * opt_uc;
char * opt_uchime2_denovo;
char * opt_uchime3_denovo;
char * opt_uchime_denovo;
char * opt_uchime_ref;
char * opt_uchimealns;
char * opt_uchimeout;
char * opt_userout;
double * opt_ee_cutoffs_values;
double opt_abskew;
double opt_chimeras_diff_pct;
double opt_dn;
double opt_fastq_maxdiffpct;
double opt_fastq_maxee;
double opt_fastq_maxee_rate;
double opt_fastq_truncee;
double opt_fastq_truncee_rate;
double opt_id;
double opt_lca_cutoff;
double opt_maxid;
double opt_maxqt;
double opt_maxsizeratio;
double opt_maxsl;
double opt_mid;
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
int opt_chimeras_length_min;
int opt_chimeras_parents_max;
int opt_chimeras_parts;
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
int opt_length_cutoffs_increment;
int opt_length_cutoffs_longest;
int opt_length_cutoffs_shortest;
int opt_mindiffs;
int opt_slots;
int opt_uchimeout5;
int opt_usersort;
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
int64_t opt_fastq_minqual;
int64_t opt_fastq_qmax;
int64_t opt_fastq_qmaxout;
int64_t opt_fastq_qmin;
int64_t opt_fastq_qminout;
int64_t opt_fastq_stripleft;
int64_t opt_fastq_stripright;
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
int64_t opt_strand;
int64_t opt_subseq_end;
int64_t opt_subseq_start;
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

static std::array<char, max_line_length> prog_header {{}};
static char * cmdline;
static time_t time_start;
static time_t time_finish;

std::FILE * fp_log = nullptr;


#ifdef __x86_64__
#define cpuid(f1, f2, a, b, c, d)                                       \
  __asm__ __volatile__ ("cpuid"                                         \
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d)        \
                        : "a" (f1), "c" (f2));
#endif


auto cpu_features_detect() -> void
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
  unsigned int a = 0;
  unsigned int b = 0;
  unsigned int c = 0;
  unsigned int d = 0;

  cpuid(0, 0, a, b, c, d);
  unsigned int const maxlevel = a & 0xff;

  if (maxlevel >= 1)
    {
      cpuid(1, 0, a, b, c, d);
      mmx_present    = (d >> 23U) & 1U;
      sse_present    = (d >> 25U) & 1U;
      sse2_present   = (d >> 26U) & 1U;
      sse3_present   = (c >>  0U) & 1U;
      ssse3_present  = (c >>  9U) & 1U;
      sse41_present  = (c >> 19U) & 1U;
      sse42_present  = (c >> 20U) & 1U;
      popcnt_present = (c >> 23U) & 1U;
      avx_present    = (c >> 28U) & 1U;

      if (maxlevel >= 7)
        {
          cpuid(7, 0, a, b, c, d);
          avx2_present = (b >>  5U) & 1U;
        }
    }
#else
    // simde
#endif
}


auto cpu_features_show() -> void
{
  fprintf(stderr, "CPU features:");
  if (neon_present != 0)
    {
      fprintf(stderr, " neon");
    }
  if (altivec_present != 0)
    {
      fprintf(stderr, " altivec");
    }
  if (mmx_present != 0)
    {
      fprintf(stderr, " mmx");
    }
  if (sse_present != 0)
    {
      fprintf(stderr, " sse");
    }
  if (sse2_present != 0)
    {
      fprintf(stderr, " sse2");
    }
  if (sse3_present != 0)
    {
      fprintf(stderr, " sse3");
    }
  if (ssse3_present != 0)
    {
      fprintf(stderr, " ssse3");
    }
  if (sse41_present != 0)
    {
      fprintf(stderr, " sse4.1");
    }
  if (sse42_present != 0)
    {
      fprintf(stderr, " sse4.2");
    }
  if (popcnt_present != 0)
    {
      fprintf(stderr, " popcnt");
    }
  if (avx_present != 0)
    {
      fprintf(stderr, " avx");
    }
  if (avx2_present != 0)
    {
      fprintf(stderr, " avx2");
    }
  fprintf(stderr, "\n");
}


auto args_get_ee_cutoffs(char * arg) -> void
{
  /* get comma-separated list of floating point numbers */
  /* save in ee_cutoffs_count and ee_cutoffs_values */

  auto const commas = std::count(arg, arg + std::strlen(arg), ',');

  opt_ee_cutoffs_count = 0;
  opt_ee_cutoffs_values = (double *) xrealloc(opt_ee_cutoffs_values, (commas + 1) * sizeof(double));

  char * cursor = arg;
  while (true)
    {
      double val = 0;
      int skip = 0;

      if ((std::sscanf(cursor, "%lf%n", &val, &skip) != 1) or (val <= 0.0))
        {
          fatal("Invalid arguments to ee_cutoffs");
        }

      opt_ee_cutoffs_values[opt_ee_cutoffs_count] = val;
      ++opt_ee_cutoffs_count;

      cursor += skip;

      if (*cursor == ',')
        {
          ++cursor;
        }
      else if (*cursor == '\0')
        {
          break;
        }
      else
        {
          fatal("Invalid arguments to ee_cutoffs");
        }
    }
}


auto args_get_length_cutoffs(char * arg) -> void
{
  /* get comma-separated list of 3 integers: */
  /* smallest, largest and increment. */
  /* second value may be * indicating no limit */
  /* save in length_cutoffs_{smallest,largest,increment} */

  // refactoring: std::stoi(), faster than sscanf()
  static constexpr auto n_of_expected_assignments= 3;
  int skip = 0;  // receives the number of characters read so far ('%n')
  if (std::sscanf(arg, "%d,%d,%d%n", &opt_length_cutoffs_shortest, &opt_length_cutoffs_longest, &opt_length_cutoffs_increment, & skip) == n_of_expected_assignments)
    {
      if ((size_t) skip < std::strlen(arg))
        {
          fatal("Invalid arguments to length_cutoffs");
        }
    }
  else if (std::sscanf(arg, "%d,*,%d%n", &opt_length_cutoffs_shortest, &opt_length_cutoffs_increment, &skip) == 2)
    {
      if ((size_t) skip < std::strlen(arg))
        {
          fatal("Invalid arguments to length_cutoffs");
        }
      opt_length_cutoffs_longest = std::numeric_limits<int>::max();
    }
  else
    {
      fatal("Invalid arguments to length_cutoffs");
    }

  if ((opt_length_cutoffs_shortest < 1) or
      (opt_length_cutoffs_shortest > opt_length_cutoffs_longest) or
      (opt_length_cutoffs_increment < 1))
    {
      fatal("Invalid arguments to length_cutoffs");
    }
}


auto args_get_gap_penalty_string(char * arg, bool const is_open) -> void
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

  char * cursor = arg;

  while (*cursor != '\0')
    {
      int skip = 0;
      int pen = 0;

      if (std::sscanf(cursor, "%d%n", &pen, &skip) == 1)
        {
          cursor += skip;
        }
      else if (*cursor == '*')
        {
          pen = 1000;
          ++cursor;
        }
      else
        {
          fatal("Invalid gap penalty argument (%s)", cursor);
        }

      char * q = cursor;

      int set_E = 0;
      int set_I = 0;
      int set_L = 0;
      int set_R = 0;
      int set_Q = 0;
      int set_T = 0;

      while ((*cursor != '\0') and (*cursor != '/'))
        {
          switch(*cursor)
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
              fatal("Invalid char '%.1s' in gap penalty string", cursor);
              break;
            }
          ++cursor;
        }

      if (*cursor == '/')
        {
          ++cursor;
        }

      if ((set_E != 0) and ((set_L != 0) or (set_R != 0)))
        {
          fatal("Invalid gap penalty string (E and L or R) '%s'", q);
        }

      if (set_E != 0)
        {
          set_L = 1;
          set_R = 1;
        }

      /* if neither L, I, R nor E is specified, it applies to all */

      if ((set_L == 0) and (set_I == 0) and (set_R == 0))
        {
          set_L = 1;
          set_I = 1;
          set_R = 1;
        }

      /* if neither Q nor T is specified, it applies to both */

      if ((set_Q == 0) and (set_T == 0))
        {
          set_Q = 1;
          set_T = 1;
        }

      if (is_open)
        {
          if (set_Q != 0)
            {
              if (set_L != 0)
                {
                  opt_gap_open_query_left = pen;
                }
              if (set_I != 0)
                {
                  opt_gap_open_query_interior = pen;
                }
              if (set_R != 0)
                {
                  opt_gap_open_query_right = pen;
                }
            }
          if (set_T != 0)
            {
              if (set_L != 0)
                {
                  opt_gap_open_target_left = pen;
                }
              if (set_I != 0)
                {
                  opt_gap_open_target_interior = pen;
                }
              if (set_R != 0)
                {
                  opt_gap_open_target_right = pen;
                }
            }
        }
      else
        {
          if (set_Q != 0)
            {
              if (set_L != 0)
                {
                  opt_gap_extension_query_left = pen;
                }
              if (set_I != 0)
                {
                  opt_gap_extension_query_interior = pen;
                }
              if (set_R != 0)
                {
                  opt_gap_extension_query_right = pen;
                }
            }
          if (set_T != 0)
            {
              if (set_L != 0)
                {
                  opt_gap_extension_target_left = pen;
                }
              if (set_I != 0)
                {
                  opt_gap_extension_target_interior = pen;
                }
              if (set_R != 0)
                {
                  opt_gap_extension_target_right = pen;
                }
            }
        }
    }
}


auto args_getlong(char * arg) -> int64_t
{
  int len = 0;
  int64_t temp = 0;
  auto const ret = std::sscanf(arg, "%" PRId64 "%n", &temp, &len);
  if ((ret == 0) or (((unsigned int) (len)) < std::strlen(arg)))
    {
      fatal("Illegal option argument");
    }
  return temp;
}


auto args_getdouble(char * arg) -> double
{
  int len = 0;
  double temp = 0;
  auto const ret = std::sscanf(arg, "%lf%n", &temp, &len);

  if ((ret == 0) or (((unsigned int) (len)) < std::strlen(arg)))
    {
      fatal("Illegal option argument");
    }
  return temp;
}


auto args_init(int argc, char ** argv, struct Parameters & parameters) -> void
{
  /* Set defaults */
  static constexpr auto dbl_max = std::numeric_limits<double>::max();
  static constexpr auto int_max = std::numeric_limits<int>::max();
  static constexpr auto long_min = std::numeric_limits<long>::min();
  static constexpr auto number_of_commands = std::size_t{50};
  static constexpr auto max_number_of_options = std::size_t{99};

  parameters.progname = argv[0];

  opt_abskew = 0.0;
  opt_acceptall = 0;
  opt_alignwidth = 80;
  opt_alnout = nullptr;
  opt_biomout = nullptr;
  opt_blast6out = nullptr;
  opt_borderline = nullptr;
  opt_centroids = nullptr;
  opt_chimeras = nullptr;
  opt_chimeras_denovo = nullptr;
  opt_chimeras_diff_pct = 0.0;
  opt_chimeras_length_min = 10;
  opt_chimeras_parents_max = 3;
  opt_chimeras_parts = 0;
  opt_cluster_fast = nullptr;
  opt_cluster_size = nullptr;
  opt_cluster_smallmem = nullptr;
  opt_cluster_unoise = nullptr;
  opt_clusters = nullptr;
  opt_cons_truncate = 0;
  opt_consout = nullptr;
  opt_db = nullptr;
  opt_dbmask = MASK_DUST;
  opt_dbmatched = nullptr;
  opt_dbnotmatched = nullptr;
  opt_dn = 1.4;
  opt_ee_cutoffs_count = 3;
  opt_ee_cutoffs_values = (double *) xmalloc(opt_ee_cutoffs_count * sizeof(double));
  opt_ee_cutoffs_values[0] = 0.5;
  opt_ee_cutoffs_values[1] = 1.0;
  opt_ee_cutoffs_values[2] = 2.0;
  opt_eeout = false;
  opt_eetabbedout = nullptr;
  opt_fasta_score = false;
  opt_fasta_width = 80;
  opt_fastaout = nullptr;
  opt_fastaout_discarded = nullptr;
  opt_fastaout_discarded_rev = nullptr;
  opt_fastaout_notmerged_fwd = nullptr;
  opt_fastaout_notmerged_rev = nullptr;
  opt_fastaout_rev = nullptr;
  opt_fastapairs = nullptr;
  opt_fastq_allowmergestagger = false;
  opt_fastq_ascii = 33;
  opt_fastq_asciiout = 33;
  opt_fastq_convert = nullptr;
  opt_fastq_eeout = false;
  opt_fastq_maxdiffpct = 100.0;
  opt_fastq_maxdiffs = 10;
  opt_fastq_maxee = dbl_max;
  opt_fastq_maxee_rate = dbl_max;
  opt_fastq_maxlen = int64_max;
  opt_fastq_maxmergelen  = 1000000;
  opt_fastq_maxns = int64_max;
  opt_fastq_minlen = 1;
  opt_fastq_minmergelen = 0;
  opt_fastq_minovlen = 10;
  opt_fastq_minqual = 0;
  opt_fastq_nostagger = true;
  opt_fastq_qmax = 41;
  opt_fastq_qmaxout = 41;
  opt_fastq_qmin = 0;
  opt_fastq_qminout = 0;
  opt_fastq_stripleft = 0;
  opt_fastq_stripright = 0;
  opt_fastq_truncee = dbl_max;
  opt_fastq_truncee_rate = dbl_max;
  opt_fastq_trunclen = -1;
  opt_fastq_trunclen_keep = -1;
  opt_fastq_truncqual = long_min;
  opt_fastqout = nullptr;
  opt_fastqout_discarded = nullptr;
  opt_fastqout_discarded_rev = nullptr;
  opt_fastqout_notmerged_fwd = nullptr;
  opt_fastqout_notmerged_rev = nullptr;
  opt_fastqout_rev = nullptr;
  opt_fulldp = 0;
  opt_gap_extension_query_interior = 2;
  opt_gap_extension_query_left = 1;
  opt_gap_extension_query_right = 1;
  opt_gap_extension_target_interior = 2;
  opt_gap_extension_target_left = 1;
  opt_gap_extension_target_right = 1;
  opt_gap_open_query_interior = 20;
  opt_gap_open_query_left = 2;
  opt_gap_open_query_right = 2;
  opt_gap_open_target_interior = 20;
  opt_gap_open_target_left = 2;
  opt_gap_open_target_right = 2;
  opt_gzip_decompress = false;
  opt_hardmask = 0;
  opt_id = -1.0;
  opt_iddef = 2;
  opt_idprefix = 0;
  opt_idsuffix = 0;
  opt_label = nullptr;
  opt_label_field = nullptr;
  opt_label_substr_match = false;
  opt_label_suffix = nullptr;
  opt_label_word = nullptr;
  opt_label_words = nullptr;
  opt_labels = nullptr;
  opt_lca_cutoff = 1.0;
  opt_lcaout = nullptr;
  opt_leftjust = 0;
  opt_length_cutoffs_increment = 50;
  opt_length_cutoffs_longest = int_max;
  opt_length_cutoffs_shortest = 50;
  opt_lengthout = false;
  opt_log = nullptr;
  opt_match = 2;
  opt_matched = nullptr;
  opt_maxaccepts = 1;
  opt_maxdiffs = int_max;
  opt_maxgaps = int_max;
  opt_maxhits = 0;
  opt_maxid = 1.0;
  opt_maxqsize = int_max;
  opt_maxqt = dbl_max;
  opt_maxrejects = -1;
  opt_maxseqlength = default_maxseqlength;
  opt_maxsize = int64_max;
  opt_maxsizeratio = dbl_max;
  opt_maxsl = dbl_max;
  opt_maxsubs = int_max;
  opt_maxuniquesize = int64_max;
  opt_mid = 0.0;
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
  opt_mothur_shared_out = nullptr;
  opt_msaout = nullptr;
  opt_n_mismatch = false;
  opt_no_progress = false;
  opt_nonchimeras = nullptr;
  opt_notmatched = nullptr;
  opt_notmatched = nullptr;
  opt_notrunclabels = 0;
  opt_otutabout = nullptr;
  opt_output = nullptr;
  opt_output_no_hits = 0;
  opt_pattern = nullptr;
  opt_profile = nullptr;
  opt_qmask = MASK_DUST;
  opt_qsegout = nullptr;
  opt_query_cov = 0.0;
  opt_quiet = false;
  opt_randseed = 0;
  opt_relabel = nullptr;
  opt_relabel_keep = false;
  opt_relabel_md5 = false;
  opt_relabel_self = false;
  opt_relabel_sha1 = false;
  opt_reverse = nullptr;
  opt_rightjust = 0;
  opt_rowlen = 64;
  opt_samheader = false;
  opt_samout = nullptr;
  opt_sample = nullptr;
  opt_sample_pct = 0;
  opt_sample_size = 0;
  opt_self = 0;
  opt_selfid = 0;
  opt_sintax_cutoff = 0.0;
  opt_sintax_random = false;
  opt_sizein = false;
  opt_sizeorder = false;
  opt_sizeout = false;
  opt_slots = 0;
  opt_strand = 1;
  opt_subseq_end = int64_max;
  opt_subseq_start = 1;
  opt_tabbedout = nullptr;
  opt_target_cov = 0.0;
  opt_threads = 0;
  opt_top_hits_only = 0;
  opt_topn = int64_max;
  opt_tsegout = nullptr;
  opt_uc = nullptr;
  opt_uc_allhits = 0;
  opt_uchime2_denovo = nullptr;
  opt_uchime3_denovo = nullptr;
  opt_uchime_denovo = nullptr;
  opt_uchime_ref = nullptr;
  opt_uchimealns = nullptr;
  opt_uchimeout = nullptr;
  opt_uchimeout5 = 0;
  opt_unoise_alpha = 2.0;
  opt_userout = nullptr;
  opt_usersort = 0;
  opt_weak_id = 10.0;
  opt_wordlength = 0;
  opt_xee = false;
  opt_xlength = false;
  opt_xn = 8.0;
  opt_xsize = false;

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
      option_chimeras_denovo,
      option_chimeras_diff_pct,
      option_chimeras_length_min,
      option_chimeras_parents_max,
      option_chimeras_parts,
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
      option_derep_id,
      option_derep_prefix,
      option_derep_smallmem,
      option_dn,
      option_ee_cutoffs,
      option_eeout,
      option_eetabbedout,
      option_fasta2fastq,
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
      option_fastq_minqual,
      option_fastq_nostagger,
      option_fastq_qmax,
      option_fastq_qmaxout,
      option_fastq_qmin,
      option_fastq_qminout,
      option_fastq_qout_max,
      option_fastq_stats,
      option_fastq_stripleft,
      option_fastq_stripright,
      option_fastq_tail,
      option_fastq_truncee,
      option_fastq_truncee_rate,
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
      option_fastx_uniques,
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
      option_lca_cutoff,
      option_lcaout,
      option_leftjust,
      option_length_cutoffs,
      option_lengthout,
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
      option_n_mismatch,
      option_no_progress,
      option_nonchimeras,
      option_notmatched,
      option_notmatchedfq,
      option_notrunclabels,
      option_orient,
      option_otutabout,
      option_output,
      option_output_no_hits,
      option_pattern,
      option_profile,
      option_qmask,
      option_qsegout,
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
      option_sample,
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
      option_sintax_random,
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
      option_tsegout,
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
      option_xlength,
      option_xn,
      option_xsize
    };

  static constexpr std::array<struct option, 247> long_options =
    {{
      {"abskew",                required_argument, nullptr, 0 },
      {"acceptall",             no_argument,       nullptr, 0 },
      {"alignwidth",            required_argument, nullptr, 0 },
      {"allpairs_global",       required_argument, nullptr, 0 },
      {"alnout",                required_argument, nullptr, 0 },
      {"band",                  required_argument, nullptr, 0 },
      {"biomout",               required_argument, nullptr, 0 },
      {"blast6out",             required_argument, nullptr, 0 },
      {"borderline",            required_argument, nullptr, 0 },
      {"bzip2_decompress",      no_argument,       nullptr, 0 },
      {"centroids",             required_argument, nullptr, 0 },
      {"chimeras",              required_argument, nullptr, 0 },
      {"chimeras_denovo",       required_argument, nullptr, 0 },
      {"chimeras_diff_pct",     required_argument, nullptr, 0 },
      {"chimeras_length_min",   required_argument, nullptr, 0 },
      {"chimeras_parents_max",  required_argument, nullptr, 0 },
      {"chimeras_parts",        required_argument, nullptr, 0 },
      {"cluster_fast",          required_argument, nullptr, 0 },
      {"cluster_size",          required_argument, nullptr, 0 },
      {"cluster_smallmem",      required_argument, nullptr, 0 },
      {"cluster_unoise",        required_argument, nullptr, 0 },
      {"clusterout_id",         no_argument,       nullptr, 0 },
      {"clusterout_sort",       no_argument,       nullptr, 0 },
      {"clusters",              required_argument, nullptr, 0 },
      {"cons_truncate",         no_argument,       nullptr, 0 },
      {"consout",               required_argument, nullptr, 0 },
      {"cut",                   required_argument, nullptr, 0 },
      {"cut_pattern",           required_argument, nullptr, 0 },
      {"db",                    required_argument, nullptr, 0 },
      {"dbmask",                required_argument, nullptr, 0 },
      {"dbmatched",             required_argument, nullptr, 0 },
      {"dbnotmatched",          required_argument, nullptr, 0 },
      {"derep_fulllength",      required_argument, nullptr, 0 },
      {"derep_id",              required_argument, nullptr, 0 },
      {"derep_prefix",          required_argument, nullptr, 0 },
      {"derep_smallmem",        required_argument, nullptr, 0 },
      {"dn",                    required_argument, nullptr, 0 },
      {"ee_cutoffs",            required_argument, nullptr, 0 },
      {"eeout",                 no_argument,       nullptr, 0 },
      {"eetabbedout",           required_argument, nullptr, 0 },
      {"fasta2fastq",           required_argument, nullptr, 0 },
      {"fasta_score",           no_argument,       nullptr, 0 },
      {"fasta_width",           required_argument, nullptr, 0 },
      {"fastaout",              required_argument, nullptr, 0 },
      {"fastaout_discarded",    required_argument, nullptr, 0 },
      {"fastaout_discarded_rev",required_argument, nullptr, 0 },
      {"fastaout_notmerged_fwd",required_argument, nullptr, 0 },
      {"fastaout_notmerged_rev",required_argument, nullptr, 0 },
      {"fastaout_rev",          required_argument, nullptr, 0 },
      {"fastapairs",            required_argument, nullptr, 0 },
      {"fastq_allowmergestagger", no_argument,     nullptr, 0 },
      {"fastq_ascii",           required_argument, nullptr, 0 },
      {"fastq_asciiout",        required_argument, nullptr, 0 },
      {"fastq_chars",           required_argument, nullptr, 0 },
      {"fastq_convert",         required_argument, nullptr, 0 },
      {"fastq_eeout",           no_argument,       nullptr, 0 },
      {"fastq_eestats",         required_argument, nullptr, 0 },
      {"fastq_eestats2",        required_argument, nullptr, 0 },
      {"fastq_filter",          required_argument, nullptr, 0 },
      {"fastq_join",            required_argument, nullptr, 0 },
      {"fastq_maxdiffpct",      required_argument, nullptr, 0 },
      {"fastq_maxdiffs",        required_argument, nullptr, 0 },
      {"fastq_maxee",           required_argument, nullptr, 0 },
      {"fastq_maxee_rate",      required_argument, nullptr, 0 },
      {"fastq_maxlen",          required_argument, nullptr, 0 },
      {"fastq_maxmergelen",     required_argument, nullptr, 0 },
      {"fastq_maxns",           required_argument, nullptr, 0 },
      {"fastq_mergepairs",      required_argument, nullptr, 0 },
      {"fastq_minlen",          required_argument, nullptr, 0 },
      {"fastq_minmergelen",     required_argument, nullptr, 0 },
      {"fastq_minovlen",        required_argument, nullptr, 0 },
      {"fastq_minqual",         required_argument, nullptr, 0 },
      {"fastq_nostagger",       no_argument,       nullptr, 0 },
      {"fastq_qmax",            required_argument, nullptr, 0 },
      {"fastq_qmaxout",         required_argument, nullptr, 0 },
      {"fastq_qmin",            required_argument, nullptr, 0 },
      {"fastq_qminout",         required_argument, nullptr, 0 },
      {"fastq_qout_max",        no_argument,       nullptr, 0 },
      {"fastq_stats",           required_argument, nullptr, 0 },
      {"fastq_stripleft",       required_argument, nullptr, 0 },
      {"fastq_stripright",      required_argument, nullptr, 0 },
      {"fastq_tail",            required_argument, nullptr, 0 },
      {"fastq_truncee",         required_argument, nullptr, 0 },
      {"fastq_truncee_rate",    required_argument, nullptr, 0 },
      {"fastq_trunclen",        required_argument, nullptr, 0 },
      {"fastq_trunclen_keep",   required_argument, nullptr, 0 },
      {"fastq_truncqual",       required_argument, nullptr, 0 },
      {"fastqout",              required_argument, nullptr, 0 },
      {"fastqout_discarded",    required_argument, nullptr, 0 },
      {"fastqout_discarded_rev",required_argument, nullptr, 0 },
      {"fastqout_notmerged_fwd",required_argument, nullptr, 0 },
      {"fastqout_notmerged_rev",required_argument, nullptr, 0 },
      {"fastqout_rev",          required_argument, nullptr, 0 },
      {"fastx_filter",          required_argument, nullptr, 0 },
      {"fastx_getseq",          required_argument, nullptr, 0 },
      {"fastx_getseqs",         required_argument, nullptr, 0 },
      {"fastx_getsubseq",       required_argument, nullptr, 0 },
      {"fastx_mask",            required_argument, nullptr, 0 },
      {"fastx_revcomp",         required_argument, nullptr, 0 },
      {"fastx_subsample",       required_argument, nullptr, 0 },
      {"fastx_uniques",         required_argument, nullptr, 0 },
      {"fulldp",                no_argument,       nullptr, 0 },
      {"gapext",                required_argument, nullptr, 0 },
      {"gapopen",               required_argument, nullptr, 0 },
      {"gzip_decompress",       no_argument,       nullptr, 0 },
      {"h",                     no_argument,       nullptr, 0 },
      {"hardmask",              no_argument,       nullptr, 0 },
      {"help",                  no_argument,       nullptr, 0 },
      {"hspw",                  required_argument, nullptr, 0 },
      {"id",                    required_argument, nullptr, 0 },
      {"iddef",                 required_argument, nullptr, 0 },
      {"idprefix",              required_argument, nullptr, 0 },
      {"idsuffix",              required_argument, nullptr, 0 },
      {"join_padgap",           required_argument, nullptr, 0 },
      {"join_padgapq",          required_argument, nullptr, 0 },
      {"label",                 required_argument, nullptr, 0 },
      {"label_field",           required_argument, nullptr, 0 },
      {"label_substr_match",    no_argument,       nullptr, 0 },
      {"label_suffix",          required_argument, nullptr, 0 },
      {"label_word",            required_argument, nullptr, 0 },
      {"label_words",           required_argument, nullptr, 0 },
      {"labels",                required_argument, nullptr, 0 },
      {"lca_cutoff",            required_argument, nullptr, 0 },
      {"lcaout",                required_argument, nullptr, 0 },
      {"leftjust",              no_argument,       nullptr, 0 },
      {"length_cutoffs",        required_argument, nullptr, 0 },
      {"lengthout",             no_argument,       nullptr, 0 },
      {"log",                   required_argument, nullptr, 0 },
      {"makeudb_usearch",       required_argument, nullptr, 0 },
      {"maskfasta",             required_argument, nullptr, 0 },
      {"match",                 required_argument, nullptr, 0 },
      {"matched",               required_argument, nullptr, 0 },
      {"max_unmasked_pct",      required_argument, nullptr, 0 },
      {"maxaccepts",            required_argument, nullptr, 0 },
      {"maxdiffs",              required_argument, nullptr, 0 },
      {"maxgaps",               required_argument, nullptr, 0 },
      {"maxhits",               required_argument, nullptr, 0 },
      {"maxid",                 required_argument, nullptr, 0 },
      {"maxqsize",              required_argument, nullptr, 0 },
      {"maxqt",                 required_argument, nullptr, 0 },
      {"maxrejects",            required_argument, nullptr, 0 },
      {"maxseqlength",          required_argument, nullptr, 0 },
      {"maxsize",               required_argument, nullptr, 0 },
      {"maxsizeratio",          required_argument, nullptr, 0 },
      {"maxsl",                 required_argument, nullptr, 0 },
      {"maxsubs",               required_argument, nullptr, 0 },
      {"maxuniquesize",         required_argument, nullptr, 0 },
      {"mid",                   required_argument, nullptr, 0 },
      {"min_unmasked_pct",      required_argument, nullptr, 0 },
      {"mincols",               required_argument, nullptr, 0 },
      {"mindiffs",              required_argument, nullptr, 0 },
      {"mindiv",                required_argument, nullptr, 0 },
      {"minh",                  required_argument, nullptr, 0 },
      {"minhsp",                required_argument, nullptr, 0 },
      {"minqt",                 required_argument, nullptr, 0 },
      {"minseqlength",          required_argument, nullptr, 0 },
      {"minsize",               required_argument, nullptr, 0 },
      {"minsizeratio",          required_argument, nullptr, 0 },
      {"minsl",                 required_argument, nullptr, 0 },
      {"mintsize",              required_argument, nullptr, 0 },
      {"minuniquesize",         required_argument, nullptr, 0 },
      {"minwordmatches",        required_argument, nullptr, 0 },
      {"mismatch",              required_argument, nullptr, 0 },
      {"mothur_shared_out",     required_argument, nullptr, 0 },
      {"msaout",                required_argument, nullptr, 0 },
      {"n_mismatch",            no_argument,       nullptr, 0 },
      {"no_progress",           no_argument,       nullptr, 0 },
      {"nonchimeras",           required_argument, nullptr, 0 },
      {"notmatched",            required_argument, nullptr, 0 },
      {"notmatchedfq",          required_argument, nullptr, 0 },
      {"notrunclabels",         no_argument,       nullptr, 0 },
      {"orient",                required_argument, nullptr, 0 },
      {"otutabout",             required_argument, nullptr, 0 },
      {"output",                required_argument, nullptr, 0 },
      {"output_no_hits",        no_argument,       nullptr, 0 },
      {"pattern",               required_argument, nullptr, 0 },
      {"profile",               required_argument, nullptr, 0 },
      {"qmask",                 required_argument, nullptr, 0 },
      {"qsegout",               required_argument, nullptr, 0 },
      {"query_cov",             required_argument, nullptr, 0 },
      {"quiet",                 no_argument,       nullptr, 0 },
      {"randseed",              required_argument, nullptr, 0 },
      {"relabel",               required_argument, nullptr, 0 },
      {"relabel_keep",          no_argument,       nullptr, 0 },
      {"relabel_md5",           no_argument,       nullptr, 0 },
      {"relabel_self",          no_argument,       nullptr, 0 },
      {"relabel_sha1",          no_argument,       nullptr, 0 },
      {"rereplicate",           required_argument, nullptr, 0 },
      {"reverse",               required_argument, nullptr, 0 },
      {"rightjust",             no_argument,       nullptr, 0 },
      {"rowlen",                required_argument, nullptr, 0 },
      {"samheader",             no_argument,       nullptr, 0 },
      {"samout",                required_argument, nullptr, 0 },
      {"sample",                required_argument, nullptr, 0 },
      {"sample_pct",            required_argument, nullptr, 0 },
      {"sample_size",           required_argument, nullptr, 0 },
      {"search_exact",          required_argument, nullptr, 0 },
      {"self",                  no_argument,       nullptr, 0 },
      {"selfid",                no_argument,       nullptr, 0 },
      {"sff_clip",              no_argument,       nullptr, 0 },
      {"sff_convert",           required_argument, nullptr, 0 },
      {"shuffle",               required_argument, nullptr, 0 },
      {"sintax",                required_argument, nullptr, 0 },
      {"sintax_cutoff",         required_argument, nullptr, 0 },
      {"sintax_random",         no_argument,       nullptr, 0 },
      {"sizein",                no_argument,       nullptr, 0 },
      {"sizeorder",             no_argument,       nullptr, 0 },
      {"sizeout",               no_argument,       nullptr, 0 },
      {"slots",                 required_argument, nullptr, 0 },
      {"sortbylength",          required_argument, nullptr, 0 },
      {"sortbysize",            required_argument, nullptr, 0 },
      {"strand",                required_argument, nullptr, 0 },
      {"subseq_end",            required_argument, nullptr, 0 },
      {"subseq_start",          required_argument, nullptr, 0 },
      {"tabbedout",             required_argument, nullptr, 0 },
      {"target_cov",            required_argument, nullptr, 0 },
      {"threads",               required_argument, nullptr, 0 },
      {"top_hits_only",         no_argument,       nullptr, 0 },
      {"topn",                  required_argument, nullptr, 0 },
      {"tsegout",               required_argument, nullptr, 0 },
      {"uc",                    required_argument, nullptr, 0 },
      {"uc_allhits",            no_argument,       nullptr, 0 },
      {"uchime2_denovo",        required_argument, nullptr, 0 },
      {"uchime3_denovo",        required_argument, nullptr, 0 },
      {"uchime_denovo",         required_argument, nullptr, 0 },
      {"uchime_ref",            required_argument, nullptr, 0 },
      {"uchimealns",            required_argument, nullptr, 0 },
      {"uchimeout",             required_argument, nullptr, 0 },
      {"uchimeout5",            no_argument,       nullptr, 0 },
      {"udb2fasta",             required_argument, nullptr, 0 },
      {"udbinfo",               required_argument, nullptr, 0 },
      {"udbstats",              required_argument, nullptr, 0 },
      {"unoise_alpha",          required_argument, nullptr, 0 },
      {"usearch_global",        required_argument, nullptr, 0 },
      {"userfields",            required_argument, nullptr, 0 },
      {"userout",               required_argument, nullptr, 0 },
      {"usersort",              no_argument,       nullptr, 0 },
      {"v",                     no_argument,       nullptr, 0 },
      {"version",               no_argument,       nullptr, 0 },
      {"weak_id",               required_argument, nullptr, 0 },
      {"wordlength",            required_argument, nullptr, 0 },
      {"xdrop_nw",              required_argument, nullptr, 0 },
      {"xee",                   no_argument,       nullptr, 0 },
      {"xlength",               no_argument,       nullptr, 0 },
      {"xn",                    required_argument, nullptr, 0 },
      {"xsize",                 no_argument,       nullptr, 0 },
      { nullptr,                0,                 nullptr, 0 }
      }};

  constexpr int options_count = (sizeof(long_options) / sizeof(struct option)) - 1;

  std::vector<bool> options_selected(options_count);

  int options_index = 0;
  int val = 0;  // recognized long option: return 'val' if 'flag' is nullptr

  while ((val = getopt_long_only(argc, argv, "", long_options.data(),
                               &options_index)) == 0)
    {
      if (options_index < options_count)
        {
          options_selected[options_index] = true;
        }

      switch(options_index)
        {
        case option_help:
          parameters.opt_help = true;
          break;

        case option_version:
          parameters.opt_version = true;
          break;

        case option_alnout:
          opt_alnout = optarg;
          break;

        case option_usearch_global:
          parameters.opt_usearch_global = optarg;
          break;

        case option_db:
          opt_db = optarg;
          parameters.opt_db = optarg;
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
            {
              opt_strand = 1;
              parameters.opt_strand = false;
            }
          else if (strcasecmp(optarg, "both") == 0)
            {
              opt_strand = 2;
              parameters.opt_strand = true;
            }
          else
            {
              fatal("The argument to --strand must be plus or both");
            }
          break;

        case option_threads:
          opt_threads = static_cast<int64_t>(args_getdouble(optarg));
          parameters.opt_threads = static_cast<int64_t>(args_getdouble(optarg));
          break;

        case option_gapopen:
          args_get_gap_penalty_string(optarg, true);
          break;

        case option_gapext:
          args_get_gap_penalty_string(optarg, false);
          break;

        case option_rowlen:
          opt_rowlen = args_getlong(optarg);
          break;

        case option_userfields:
          if (parse_userfields_arg(optarg) == 0)
            {
              fatal("Unrecognized userfield argument");
            }
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
          parameters.opt_uc = optarg;
          break;

        case option_weak_id:
          opt_weak_id = args_getdouble(optarg);
          break;

        case option_uc_allhits:
          opt_uc_allhits = 1;
          parameters.opt_uc_allhits = true;
          break;

        case option_notrunclabels:
          opt_notrunclabels = 1;
          parameters.opt_notrunclabels = true;
          break;

        case option_sortbysize:
          parameters.opt_sortbysize = optarg;
          break;

        case option_output:
          opt_output = optarg;
          parameters.opt_output = optarg;
          break;

        case option_minsize:
          opt_minsize = args_getlong(optarg);
          parameters.opt_minsize = args_getlong(optarg);
          if (parameters.opt_minsize <= 0)
            {
              fatal("The argument to --minsize must be at least 1");
            }
          break;

        case option_maxsize:
          opt_maxsize = args_getlong(optarg);
          parameters.opt_maxsize = args_getlong(optarg);
          break;

        case option_relabel:
          opt_relabel = optarg;
          parameters.opt_relabel = optarg;
          break;

        case option_sizeout:
          opt_sizeout = true;
          parameters.opt_sizeout = true;
          break;

        case option_derep_fulllength:
          parameters.opt_derep_fulllength = optarg;
          break;

        case option_minseqlength:
          opt_minseqlength = args_getlong(optarg);
          parameters.opt_minseqlength = args_getlong(optarg);
          if (parameters.opt_minseqlength < 0)
            {
              fatal("The argument to --minseqlength must not be negative");
            }
          break;

        case option_minuniquesize:
          opt_minuniquesize = args_getlong(optarg);
          parameters.opt_minuniquesize = args_getlong(optarg);
          break;

        case option_topn:
          opt_topn = args_getlong(optarg);
          parameters.opt_topn = args_getlong(optarg);
          if (parameters.opt_topn == 0)
            {
              fatal("The argument to --topn must be greater than zero");
            }
          break;

        case option_maxseqlength:
          opt_maxseqlength = args_getlong(optarg);
          parameters.opt_maxseqlength = args_getlong(optarg);
          break;

        case option_sizein:
          opt_sizein = true;
          parameters.opt_sizein = true;
          break;

        case option_sortbylength:
          parameters.opt_sortbylength = optarg;
          break;

        case option_matched:
          opt_matched = optarg;
          break;

        case option_notmatched:
          opt_notmatched = optarg;
          break;

        case option_dbmatched:
          opt_dbmatched = optarg;
          parameters.opt_dbmatched = optarg;
          break;

        case option_dbnotmatched:
          opt_dbnotmatched = optarg;
          parameters.opt_dbnotmatched = optarg;
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
          parameters.opt_fasta_width = args_getlong(optarg);
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
          parameters.opt_shuffle = optarg;
          break;

        case option_randseed:
          opt_randseed = args_getlong(optarg);
          parameters.opt_randseed = args_getlong(optarg);
          break;

        case option_maskfasta:
          parameters.opt_maskfasta = optarg;
          break;

        case option_hardmask:
          opt_hardmask = 1;
          parameters.opt_hardmask = true;
          break;

        case option_qmask:
          if (strcasecmp(optarg, "none") == 0)
            {
              opt_qmask = MASK_NONE;
              parameters.opt_qmask = MASK_NONE;
            }
          else if (strcasecmp(optarg, "dust") == 0)
            {
              opt_qmask = MASK_DUST;
              parameters.opt_qmask = MASK_DUST;
            }
          else if (strcasecmp(optarg, "soft") == 0)
            {
              opt_qmask = MASK_SOFT;
              parameters.opt_qmask = MASK_SOFT;
            }
          else
            {
              opt_qmask = MASK_ERROR;
              parameters.opt_qmask = MASK_ERROR;
            }
          break;

        case option_dbmask:
          if (strcasecmp(optarg, "none") == 0)
            {
              opt_dbmask = MASK_NONE;
            }
          else if (strcasecmp(optarg, "dust") == 0)
            {
              opt_dbmask = MASK_DUST;
            }
          else if (strcasecmp(optarg, "soft") == 0)
            {
              opt_dbmask = MASK_SOFT;
            }
          else
            {
              opt_dbmask = MASK_ERROR;
            }
          break;

        case option_cluster_smallmem:
          opt_cluster_smallmem = optarg;
          parameters.opt_cluster_smallmem = optarg;
          break;

        case option_cluster_fast:
          opt_cluster_fast = optarg;
          parameters.opt_cluster_fast = optarg;
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
          parameters.opt_maxuniquesize = args_getlong(optarg);
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
          parameters.opt_uchime_denovo = optarg;
          break;

        case option_uchime_ref:
          opt_uchime_ref = optarg;
          parameters.opt_uchime_ref = optarg;
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
          parameters.opt_allpairs_global = optarg;
          break;

        case option_acceptall:
          opt_acceptall = 1;
          break;

        case option_cluster_size:
          opt_cluster_size = optarg;
          parameters.opt_cluster_size = optarg;
          break;

        case option_samout:
          opt_samout = optarg;
          break;

        case option_log:
          opt_log = optarg;
          parameters.opt_log = optarg;
          break;

        case option_quiet:
          opt_quiet = true;
          parameters.opt_quiet = true;
          break;

        case option_fastx_subsample:
          parameters.opt_fastx_subsample = optarg;
          break;

        case option_sample_pct:
          opt_sample_pct = args_getdouble(optarg);
          parameters.opt_sample_pct = args_getdouble(optarg);
          break;

        case option_fastq_chars:
          parameters.opt_fastq_chars = optarg;
          break;

        case option_profile:
          opt_profile = optarg;
          break;

        case option_sample_size:
          opt_sample_size = args_getlong(optarg);
          parameters.opt_sample_size = args_getlong(optarg);
          break;

        case option_fastaout:
          opt_fastaout = optarg;
          parameters.opt_fastaout = optarg;
          break;

        case option_xsize:
          opt_xsize = true;
          parameters.opt_xsize = true;
          break;

        case option_clusterout_id:
          opt_clusterout_id = true;
          parameters.opt_clusterout_id = true;
          break;

        case option_clusterout_sort:
          opt_clusterout_sort = true;
          parameters.opt_clusterout_sort = true;
          break;

        case option_borderline:
          opt_borderline = optarg;
          break;

        case option_relabel_sha1:
          opt_relabel_sha1 = true;
          parameters.opt_relabel_sha1 = true;
          break;

        case option_relabel_md5:
          opt_relabel_md5 = true;
          parameters.opt_relabel_md5 = true;
          break;

        case option_derep_prefix:
          parameters.opt_derep_prefix = optarg;
          break;

        case option_fastq_filter:
          parameters.opt_fastq_filter = optarg;
          break;

        case option_fastqout:
          opt_fastqout = optarg;
          parameters.opt_fastqout = optarg;
          break;

        case option_fastaout_discarded:
          opt_fastaout_discarded = optarg;
          parameters.opt_fastaout_discarded = optarg;
          break;

        case option_fastqout_discarded:
          opt_fastqout_discarded = optarg;
          parameters.opt_fastqout_discarded = optarg;
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
          opt_eeout = true;
          parameters.opt_eeout = true;
          break;

        case option_fastq_ascii:
          opt_fastq_ascii = args_getlong(optarg);
          parameters.opt_fastq_ascii = args_getlong(optarg);
          break;

        case option_fastq_qmin:
          opt_fastq_qmin = args_getlong(optarg);
          parameters.opt_fastq_qmin = args_getlong(optarg);
          break;

        case option_fastq_qmax:
          opt_fastq_qmax = args_getlong(optarg);
          parameters.opt_fastq_qmax = args_getlong(optarg);
          break;

        case option_fastq_qmaxout:
          opt_fastq_qmaxout = args_getlong(optarg);
          parameters.opt_fastq_qmaxout = args_getlong(optarg);
          break;

        case option_fastq_stats:
          parameters.opt_fastq_stats = optarg;
          break;

        case option_fastq_tail:
          parameters.opt_fastq_tail = args_getlong(optarg);
          break;

        case option_fastx_revcomp:
          parameters.opt_fastx_revcomp = optarg;
          break;

        case option_label_suffix:
          opt_label_suffix = optarg;
          parameters.opt_label_suffix = optarg;
          break;

        case option_h:
          parameters.opt_help = true;
          break;

        case option_samheader:
          opt_samheader = true;
          parameters.opt_samheader = true;
          break;

        case option_sizeorder:
          opt_sizeorder = true;
          parameters.opt_sizeorder = true;
          break;

        case option_minwordmatches:
          opt_minwordmatches = args_getlong(optarg);
          if (opt_minwordmatches < 0)
            {
              fatal("The argument to --minwordmatches must not be negative");
            }
          break;

        case option_v:
          parameters.opt_version = true;
          break;

        case option_relabel_keep:
          opt_relabel_keep = true;
          parameters.opt_relabel_keep = true;
          break;

        case option_search_exact:
          parameters.opt_search_exact = optarg;
          break;

        case option_fastx_mask:
          parameters.opt_fastx_mask = optarg;
          break;

        case option_min_unmasked_pct:
          parameters.opt_min_unmasked_pct = args_getdouble(optarg);
          break;

        case option_max_unmasked_pct:
          parameters.opt_max_unmasked_pct = args_getdouble(optarg);
          break;

        case option_fastq_convert:
          opt_fastq_convert = optarg;
          parameters.opt_fastq_convert = optarg;
          break;

        case option_fastq_asciiout:
          opt_fastq_asciiout = args_getlong(optarg);
          parameters.opt_fastq_asciiout = args_getlong(optarg);
          break;

        case option_fastq_qminout:
          opt_fastq_qminout = args_getlong(optarg);
          parameters.opt_fastq_qminout = args_getlong(optarg);
          break;

        case option_fastq_mergepairs:
          parameters.opt_fastq_mergepairs = optarg;
          break;

        case option_fastq_eeout:
          opt_fastq_eeout = true;
          parameters.opt_fastq_eeout = true;
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
          opt_fastq_nostagger = true;
          parameters.opt_fastq_nostagger = true;
          break;

        case option_fastq_allowmergestagger:
          opt_fastq_allowmergestagger = true;
          parameters.opt_fastq_allowmergestagger = true;
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
          parameters.opt_reverse = optarg;
          break;

        case option_eetabbedout:
          opt_eetabbedout = optarg;
          break;

        case option_fasta_score:
          opt_fasta_score = true;
          parameters.opt_fasta_score = true;
          break;

        case option_fastq_eestats:
          parameters.opt_fastq_eestats = optarg;
          break;

        case option_rereplicate:
          parameters.opt_rereplicate = optarg;
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
          opt_gzip_decompress = true;
          parameters.opt_gzip_decompress = true;
          break;

        case option_bzip2_decompress:
          opt_bzip2_decompress = true;
          parameters.opt_bzip2_decompress = true;
          break;

        case option_fastq_maxlen:
          opt_fastq_maxlen = args_getlong(optarg);
          break;

        case option_fastq_truncee:
          opt_fastq_truncee = args_getdouble(optarg);
          break;

        case option_fastx_filter:
          parameters.opt_fastx_filter = optarg;
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
          opt_no_progress = true;
          parameters.opt_no_progress = true;
          break;

        case option_fastq_eestats2:
          parameters.opt_fastq_eestats2 = optarg;
          break;

        case option_ee_cutoffs:
          args_get_ee_cutoffs(optarg);
          break;

        case option_length_cutoffs:
          args_get_length_cutoffs(optarg);
          break;

        case option_makeudb_usearch:
          parameters.opt_makeudb_usearch = optarg;
          break;

        case option_udb2fasta:
          parameters.opt_udb2fasta = optarg;
          break;

        case option_udbinfo:
          parameters.opt_udbinfo = optarg;
          break;

        case option_udbstats:
          parameters.opt_udbstats = optarg;
          break;

        case option_cluster_unoise:
          opt_cluster_unoise = optarg;
          parameters.opt_cluster_unoise = optarg;
          break;

        case option_unoise_alpha:
          opt_unoise_alpha = args_getdouble(optarg);
          break;

        case option_uchime2_denovo:
          opt_uchime2_denovo = optarg;
          parameters.opt_uchime2_denovo = optarg;
          break;

        case option_uchime3_denovo:
          opt_uchime3_denovo = optarg;
          parameters.opt_uchime3_denovo = optarg;
          break;

        case option_sintax:
          parameters.opt_sintax = optarg;
          break;

        case option_sintax_cutoff:
          opt_sintax_cutoff = args_getdouble(optarg);
          break;

        case option_tabbedout:
          opt_tabbedout = optarg;
          parameters.opt_tabbedout = optarg;
          break;

        case option_fastq_maxdiffpct:
          opt_fastq_maxdiffpct = args_getdouble(optarg);
          break;

        case option_fastq_join:
          parameters.opt_fastq_join = optarg;
          break;

        case option_join_padgap:
          parameters.opt_join_padgap = optarg;
          break;

        case option_join_padgapq:
          parameters.opt_join_padgapq = optarg;
          parameters.opt_join_padgapq_set_by_user = true;
          break;

        case option_sff_convert:
          parameters.opt_sff_convert = optarg;
          break;

        case option_sff_clip:
          parameters.opt_sff_clip = true;
          break;

        case option_fastaout_rev:
          opt_fastaout_rev = optarg;
          parameters.opt_fastaout_rev = optarg;
          break;

        case option_fastaout_discarded_rev:
          opt_fastaout_discarded_rev = optarg;
          parameters.opt_fastaout_discarded_rev = optarg;
          break;

        case option_fastqout_rev:
          opt_fastqout_rev = optarg;
          parameters.opt_fastqout_rev = optarg;
          break;

        case option_fastqout_discarded_rev:
          opt_fastqout_discarded_rev = optarg;
          parameters.opt_fastqout_discarded_rev = optarg;
          break;

        case option_xee:
          opt_xee = true;
          parameters.opt_xee = true;
          break;

        case option_fastx_getseq:
          parameters.opt_fastx_getseq = optarg;
          break;

        case option_fastx_getseqs:
          parameters.opt_fastx_getseqs = optarg;
          break;

        case option_fastx_getsubseq:
          parameters.opt_fastx_getsubseq = optarg;
          break;

        case option_label_substr_match:
          opt_label_substr_match = true;
          parameters.opt_label_substr_match = true;
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
          parameters.opt_cut = optarg;
          break;

        case option_cut_pattern:
          parameters.opt_cut_pattern = optarg;
          break;

        case option_relabel_self:
          opt_relabel_self = true;
          parameters.opt_relabel_self = true;
          break;

        case option_derep_id:
          parameters.opt_derep_id = optarg;
          break;

        case option_orient:
          parameters.opt_orient = optarg;
          break;

        case option_fasta2fastq:
          parameters.opt_fasta2fastq = optarg;
          break;

        case option_lcaout:
          opt_lcaout = optarg;
          break;

        case option_lca_cutoff:
          opt_lca_cutoff = args_getdouble(optarg);
          break;

        case option_fastx_uniques:
          parameters.opt_fastx_uniques = optarg;
          break;

        case option_fastq_qout_max:
          parameters.opt_fastq_qout_max = true;
          break;

        case option_sample:
          opt_sample = optarg;
          parameters.opt_sample = optarg;
          break;

        case option_qsegout:
          opt_qsegout = optarg;
          break;

        case option_tsegout:
          opt_tsegout = optarg;
          break;

        case option_derep_smallmem:
          parameters.opt_derep_smallmem = optarg;
          break;

        case option_lengthout:
          opt_lengthout = true;
          parameters.opt_lengthout = true;
          break;

        case option_xlength:
          opt_xlength = true;
          parameters.opt_xlength = true;
          break;

        case option_chimeras_denovo:
          opt_chimeras_denovo = optarg;
          parameters.opt_chimeras_denovo = optarg;
          break;

        case option_chimeras_length_min:
          opt_chimeras_length_min = args_getlong(optarg);
          break;

        case option_chimeras_parts:
          opt_chimeras_parts = args_getlong(optarg);
          break;

        case option_chimeras_parents_max:
          opt_chimeras_parents_max = args_getlong(optarg);
          break;

        case option_chimeras_diff_pct:
          opt_chimeras_diff_pct = args_getdouble(optarg);
          break;

        case option_sintax_random:
          opt_sintax_random = true;
          parameters.opt_sintax_random = true;
          break;

        case option_n_mismatch:
          opt_n_mismatch = true;
          break;

        case option_fastq_minqual:
          opt_fastq_minqual = args_getlong(optarg);
          break;

        case option_fastq_truncee_rate:
          opt_fastq_truncee_rate = args_getdouble(optarg);
          break;

        default:
          fatal("Internal error in option parsing");
        }
    }

  /* Terminate if ambiguous or illegal options have been detected */
  if (val != -1)
    {
      exit(EXIT_FAILURE);
    }

  /* Terminate after reporting any extra non-option arguments */
  if (optind < argc)
    {
      fatal("Unrecognized string on command line (%s)", argv[optind]);
    }

  /* Below is a list of all command names, in alphabetical order. */

  static constexpr std::array<int, number_of_commands> command_options =
    {
      option_allpairs_global,
      option_chimeras_denovo,
      option_cluster_fast,
      option_cluster_size,
      option_cluster_smallmem,
      option_cluster_unoise,
      option_cut,
      option_derep_fulllength,
      option_derep_id,
      option_derep_prefix,
      option_derep_smallmem,
      option_fasta2fastq,
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
      option_fastx_uniques,
      option_h,
      option_help,
      option_makeudb_usearch,
      option_maskfasta,
      option_orient,
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

  static constexpr std::array<std::array<int, max_number_of_options>, number_of_commands> valid_options =
    {{
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
        option_label_suffix,
        option_leftjust,
        option_lengthout,
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
        option_n_mismatch,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_output_no_hits,
        option_pattern,
        option_qmask,
        option_qsegout,
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
        option_sample,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeout,
        option_slots,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_tsegout,
        option_uc,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_chimeras_denovo,
        option_abskew,
        option_alignwidth,
        option_alnout,
        option_chimeras,
        option_chimeras_diff_pct,
        option_chimeras_length_min,
        option_chimeras_parents_max,
        option_chimeras_parts,
        option_fasta_width,
        option_gapext,
        option_gapopen,
        option_hardmask,
        option_label_suffix,
        option_log,
        option_match,
        option_maxseqlength,
        option_minseqlength,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_tabbedout,
        option_threads,
        option_xee,
        option_xn,
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
        option_label_suffix,
        option_leftjust,
        option_lengthout,
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
        option_n_mismatch,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_profile,
        option_qmask,
        option_qsegout,
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
        option_sample,
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
        option_tsegout,
        option_uc,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_leftjust,
        option_lengthout,
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
        option_n_mismatch,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_profile,
        option_qmask,
        option_qsegout,
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
        option_sample,
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
        option_tsegout,
        option_uc,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_leftjust,
        option_lengthout,
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
        option_n_mismatch,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_profile,
        option_qmask,
        option_qsegout,
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
        option_sample,
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
        option_tsegout,
        option_uc,
        option_userfields,
        option_userout,
        option_usersort,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_leftjust,
        option_lengthout,
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
        option_n_mismatch,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_qsegout,
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
        option_sample,
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
        option_tsegout,
        option_uc,
        option_unoise_alpha,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
        option_log,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_derep_fulllength,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_topn,
        option_uc,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_derep_id,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_topn,
        option_uc,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_derep_prefix,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_topn,
        option_uc,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_derep_smallmem,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
        option_log,
        option_maxseqlength,
        option_maxuniquesize,
        option_minseqlength,
        option_minuniquesize,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_fasta2fastq,
        option_bzip2_decompress,
        option_fastq_asciiout,
        option_fastq_qmaxout,
        option_fastqout,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_fastq_minqual,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastq_stripleft,
        option_fastq_stripright,
        option_fastq_truncee,
        option_fastq_truncee_rate,
        option_fastq_trunclen,
        option_fastq_trunclen_keep,
        option_fastq_truncqual,
        option_fastqout,
        option_fastqout_discarded,
        option_fastqout_discarded_rev,
        option_fastqout_rev,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
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
        option_xlength,
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
        option_fastq_qminout,
        option_fastq_truncqual,
        option_fastqout,
        option_fastqout_notmerged_fwd,
        option_fastqout_notmerged_rev,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_reverse,
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_fastq_minqual,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastq_stripleft,
        option_fastq_stripright,
        option_fastq_truncee,
        option_fastq_truncee_rate,
        option_fastq_trunclen,
        option_fastq_trunclen_keep,
        option_fastq_truncqual,
        option_fastqout,
        option_fastqout_discarded,
        option_fastqout_discarded_rev,
        option_fastqout_rev,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_label_word,
        option_label_words,
        option_labels,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_subseq_end,
        option_subseq_start,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_lengthout,
        option_log,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sample_pct,
        option_sample_size,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_fastx_uniques,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastq_ascii,
        option_fastq_asciiout,
        option_fastq_qmax,
        option_fastq_qmaxout,
        option_fastq_qmin,
        option_fastq_qminout,
        option_fastq_qout_max,
        option_fastqout,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
        option_log,
        option_maxseqlength,
        option_maxuniquesize,
        option_minseqlength,
        option_minuniquesize,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_strand,
        option_tabbedout,
        option_threads,
        option_topn,
        option_uc,
        option_xee,
        option_xlength,
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
        option_maxseqlength,
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
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_orient,
        option_bzip2_decompress,
        option_db,
        option_dbmask,
        option_fasta_width,
        option_fastaout,
        option_fastqout,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
        option_log,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_qmask,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_tabbedout,
        option_threads,
        option_wordlength,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_rereplicate,
        option_bzip2_decompress,
        option_fasta_width,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lca_cutoff,
        option_lcaout,
        option_lengthout,
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
        option_qsegout,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_rowlen,
        option_samheader,
        option_samout,
        option_sample,
        option_self,
        option_sizein,
        option_sizeout,
        option_strand,
        option_threads,
        option_top_hits_only,
        option_tsegout,
        option_uc,
        option_uc_allhits,
        option_userfields,
        option_userout,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_sff_convert,
        option_fastq_asciiout,
        option_fastq_qmaxout,
        option_fastq_qminout,
        option_fastqout,
        option_label_suffix,
        option_lengthout,
        option_log,
        option_no_progress,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
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
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_topn,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_log,
        option_maxseqlength,
        option_minseqlength,
        option_no_progress,
        option_notrunclabels,
        option_quiet,
        option_randseed,
        option_sintax_cutoff,
        option_sintax_random,
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
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_topn,
        option_xee,
        option_xlength,
        option_xsize,
        -1 },

      { option_sortbysize,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_gzip_decompress,
        option_label_suffix,
        option_lengthout,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_topn,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
        option_log,
        option_match,
        option_maxseqlength,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_minseqlength,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_uchimealns,
        option_uchimeout,
        option_uchimeout5,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
        option_log,
        option_match,
        option_maxseqlength,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_minseqlength,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_uchimealns,
        option_uchimeout,
        option_uchimeout5,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
        option_log,
        option_match,
        option_maxseqlength,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_minseqlength,
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
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_uchimealns,
        option_uchimeout,
        option_uchimeout5,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lengthout,
        option_log,
        option_match,
        option_maxseqlength,
        option_mindiffs,
        option_mindiv,
        option_minh,
        option_minseqlength,
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
        option_sample,
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
        option_xlength,
        option_xn,
        option_xsize,
        -1 },

      { option_udb2fasta,
        option_fasta_width,
        option_label_suffix,
        option_lengthout,
        option_log,
        option_no_progress,
        option_output,
        option_quiet,
        option_relabel,
        option_relabel_keep,
        option_relabel_md5,
        option_relabel_self,
        option_relabel_sha1,
        option_sample,
        option_sizein,
        option_sizeout,
        option_threads,
        option_xee,
        option_xlength,
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
        option_label_suffix,
        option_lca_cutoff,
        option_lcaout,
        option_leftjust,
        option_lengthout,
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
        option_n_mismatch,
        option_no_progress,
        option_notmatched,
        option_notrunclabels,
        option_otutabout,
        option_output_no_hits,
        option_pattern,
        option_qmask,
        option_qsegout,
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
        option_sample,
        option_self,
        option_selfid,
        option_sizein,
        option_sizeout,
        option_slots,
        option_strand,
        option_target_cov,
        option_threads,
        option_top_hits_only,
        option_tsegout,
        option_uc,
        option_uc_allhits,
        option_userfields,
        option_userout,
        option_weak_id,
        option_wordlength,
        option_xdrop_nw,
        option_xee,
        option_xlength,
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
      }};

  /* check that only one commmand is specified */
  int commands = 0;
  int k = -1;
  for (int i = 0; i < commands_count; i++)
    {
      if (options_selected[command_options[i]])
        {
          ++commands;
          k = i;
        }
    }
  if (commands > 1)
    {
      fatal("More than one command specified");
    }

  /* check that only valid options are specified */
  int invalid_options = 0;

  if (commands == 0)
    {
      /* check if any options are specified */
      bool any_options = false;
      for (bool const i: options_selected)
        {
          if (i)
            {
              any_options = true;
            }
        }
      if (any_options)
        {
          fprintf(stderr, "WARNING: Options given, but no valid command specified.\n");
        }
    }
  else
    {
      for (int i = 0; i < options_count; i++)
        {
          if (options_selected[i])
            {
              int j = 0;
              bool option_is_valid = false;
              while (valid_options[k][j] >= 0)
                {
                  if (valid_options[k][j] == i)
                    {
                      option_is_valid = true;
                      break;
                    }
                  ++j;
                }
              if (not option_is_valid)
                {
                  ++invalid_options;

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
          for (int j = 1; valid_options[k][j] >= 0; j++)
            {
              fprintf(stderr, " --%s", long_options[valid_options[k][j]].name);
              ++count;
            }
          if (count == 0)
            {
              fprintf(stderr, " (none)");
            }
          fprintf(stderr, "\n");
          exit(EXIT_FAILURE);
        }
    }

  /* multi-threaded commands */

  if ((opt_threads < 0) or (opt_threads > n_threads_max))
    {
      fatal("The argument to --threads must be in the range 0 (default) to 1024");
    }

  if ((parameters.opt_allpairs_global != nullptr) or (parameters.opt_cluster_fast != nullptr) or (parameters.opt_cluster_size != nullptr) or
      (parameters.opt_cluster_smallmem != nullptr) or (parameters.opt_cluster_unoise != nullptr) or (parameters.opt_fastq_mergepairs != nullptr) or
      (parameters.opt_fastx_mask != nullptr) or (parameters.opt_maskfasta != nullptr) or (parameters.opt_search_exact != nullptr) or (parameters.opt_sintax != nullptr) or
      (parameters.opt_uchime_ref != nullptr) or (parameters.opt_usearch_global != nullptr))
    {
      if (parameters.opt_threads == 0)
        {
          opt_threads = arch_get_cores();
          parameters.opt_threads = arch_get_cores();
        }
    }
  else
    {
      if (parameters.opt_threads > 1)
        {
          fprintf(stderr, "WARNING: The %s command does not support multithreading.\nOnly 1 thread used.\n", long_options[command_options[k]].name);
        }
      opt_threads = 1;
      parameters.opt_threads = 1;
    }
  if ((parameters.opt_sintax != nullptr) and (parameters.opt_randseed != 0) and (parameters.opt_threads > 1))
    {
      fprintf(stderr, "WARNING: Using the --sintax command with the --randseed option may not work as intended with multiple threads. Use a single thread (--threads 1) to ensure reproducible results.\n");
    }

  if (parameters.opt_cluster_unoise != nullptr)
    {
      opt_weak_id = 0.90;
    }
  else
    if (opt_weak_id > opt_id)
      {
        opt_weak_id = opt_id;
      }

  if (opt_maxrejects == -1)
    {
      if (parameters.opt_cluster_fast != nullptr)
        {
          opt_maxrejects = 8;
        }
      else
        {
          opt_maxrejects = 32;
        }
    }

  if (opt_maxaccepts < 0)
    {
      fatal("The argument to --maxaccepts must not be negative");
    }

  if (opt_maxrejects < 0)
    {
      fatal("The argument to --maxrejects must not be negative");
    }

  if (opt_wordlength == 0)
    {
      /* set default word length */
      if (parameters.opt_orient != nullptr)
        {
          opt_wordlength = 12;
        }
      else
        {
          opt_wordlength = 8;
        }
    }

  if ((opt_wordlength < 3) or (opt_wordlength > 15))
    {
      fatal("The argument to --wordlength must be in the range 3 to 15");
    }

  if ((opt_iddef < 0) or (opt_iddef > 4))
    {
      fatal("The argument to --iddef must in the range 0 to 4");
    }

#if 0

  if (opt_match <= 0)
    fatal("The argument to --match must be positive");

  if (opt_mismatch >= 0)
    fatal("The argument to --mismatch must be negative");

#endif


  if (opt_alignwidth < 0)
    {
      fatal("The argument to --alignwidth must not be negative");
    }

  if (opt_rowlen < 0)
    {
      fatal("The argument to --rowlen must not be negative");
    }

  if (parameters.opt_qmask == MASK_ERROR)
    {
      fatal("The argument to --qmask must be none, dust or soft");
    }

  if (opt_dbmask == MASK_ERROR)
    {
      fatal("The argument to --dbmask must be none, dust or soft");
    }

  if ((opt_sample_pct < 0.0) or (opt_sample_pct > 100.0))
    {
      fatal("The argument to --sample_pct must be in the range 0.0 to 100.0");
    }

  if (opt_sample_size < 0)
    {
      fatal("The argument to --sample_size must not be negative");
    }

  if ((((parameters.opt_relabel != nullptr) ? 1 : 0) +
       static_cast<int>(parameters.opt_relabel_md5) +
       static_cast<int>(parameters.opt_relabel_self) +
       static_cast<int>(parameters.opt_relabel_sha1)) > 1)
    {
      fatal("Specify only one of --relabel, --relabel_self, --relabel_sha1, or --relabel_md5");
    }

  if (parameters.opt_fastq_tail < 1)
    {
      fatal("The argument to --fastq_tail must be positive");
    }

  if ((parameters.opt_min_unmasked_pct < 0.0) and (parameters.opt_min_unmasked_pct > 100.0))
    {
      fatal("The argument to --min_unmasked_pct must be between 0.0 and 100.0");
    }

  if ((parameters.opt_max_unmasked_pct < 0.0) and (parameters.opt_max_unmasked_pct > 100.0))
    {
      fatal("The argument to --max_unmasked_pct must be between 0.0 and 100.0");
    }

  if (parameters.opt_min_unmasked_pct > parameters.opt_max_unmasked_pct)
    {
      fatal("The argument to --min_unmasked_pct cannot be larger than --max_unmasked_pct");
    }

  if ((parameters.opt_fastq_ascii != 33) and (parameters.opt_fastq_ascii != 64))
    {
      fatal("The argument to --fastq_ascii must be 33 or 64");
    }

  if (opt_fastq_qmin > opt_fastq_qmax)
    {
      fatal("The argument to --fastq_qmin cannot be greater than --fastq_qmax");
    }

  if (parameters.opt_fastq_ascii + opt_fastq_qmin < 33)
    {
      fatal("Sum of arguments to --fastq_ascii and --fastq_qmin must be no less than 33");
    }

  if (parameters.opt_fastq_ascii + opt_fastq_qmax > 126)
    {
      fatal("Sum of arguments to --fastq_ascii and --fastq_qmax must be no more than 126");
    }

  if (parameters.opt_fastq_qminout > parameters.opt_fastq_qmaxout)
    {
      fatal("The argument to --fastq_qminout cannot be larger than --fastq_qmaxout");
    }

  if ((parameters.opt_fastq_asciiout != 33) and (parameters.opt_fastq_asciiout != 64))
    {
      fatal("The argument to --fastq_asciiout must be 33 or 64");
    }

  if (parameters.opt_fastq_asciiout + parameters.opt_fastq_qminout < 33)
    {
      fatal("Sum of arguments to --fastq_asciiout and --fastq_qminout must be no less than 33");
    }

  if (parameters.opt_fastq_asciiout + parameters.opt_fastq_qmaxout > 126)
    {
      fatal("Sum of arguments to --fastq_asciiout and --fastq_qmaxout must be no more than 126");
    }

  if (parameters.opt_gzip_decompress and parameters.opt_bzip2_decompress)
    {
      fatal("Specify either --gzip_decompress or --bzip2_decompress, not both");
    }

  if ((opt_sintax_cutoff < 0.0) or (opt_sintax_cutoff > 1.0))
    {
      fatal("The argument to sintax_cutoff must be in the range 0.0 to 1.0");
    }

  if ((opt_lca_cutoff <= 0.5) or (opt_lca_cutoff > 1.0))
    {
      fatal("The argument to lca_cutoff must be larger than 0.5, but not larger than 1.0");
    }

  if (parameters.opt_minuniquesize < 1)
    {
      fatal("The argument to minuniquesize must be at least 1");
    }

  if (parameters.opt_maxuniquesize < 1)
    {
      fatal("The argument to maxuniquesize must be at least 1");
    }

  if (parameters.opt_maxsize < 1)
    {
      fatal("The argument to maxsize must be at least 1");
    }

  if (opt_maxhits < 0)
    {
      fatal("The argument to maxhits cannot be negative");
    }

  if (opt_chimeras_length_min < 1)
    {
      fatal("The argument to chimeras_length_min must be at least 1");
    }

  if ((opt_chimeras_parents_max < 2) or (opt_chimeras_parents_max > maxparents))
    {
      std::array<char, 25> maxparents_string {{}};
      snprintf(maxparents_string.data(), maxparents_string.size(), "%d", maxparents);
      fatal("The argument to chimeras_parents_max must be in the range 2 to %s.\n", maxparents_string.data());
    }

  if ((opt_chimeras_diff_pct < 0.0) or (opt_chimeras_diff_pct > 50.0))
    {
      fatal("The argument to chimeras_diff_pct must be in the range 0.0 to 50.0");
    }

  if (options_selected[option_chimeras_parts] and
      ((opt_chimeras_parts < 2) or (opt_chimeras_parts > 100)))
    {
      fatal("The argument to chimeras_parts must be in the range 2 to 100");
    }

  if (parameters.opt_chimeras_denovo != nullptr)
    {
      if (not options_selected[option_alignwidth])
        {
          opt_alignwidth = 60;
        }
    }


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
    {
      opt_maxhits = int64_max;
    }

  if (opt_minwordmatches < 0)
    {
      opt_minwordmatches = minwordmatches_defaults[opt_wordlength];
    }

  /* set default opt_minsize depending on command */
  if (parameters.opt_minsize == 0)
    {
      if (parameters.opt_cluster_unoise != nullptr)
        {
          opt_minsize = 8;
          parameters.opt_minsize = 8;
        }
      else
        {
          opt_minsize = 1;
          parameters.opt_minsize = 1;
        }
    }

  /* set default opt_abskew depending on command */
  if (not options_selected[option_abskew])
    {
      if (parameters.opt_chimeras_denovo != nullptr)
        {
          opt_abskew = 1.0;
        }
      else if (opt_uchime3_denovo != nullptr)
        {
          opt_abskew = 16.0;
        }
      else
        {
          opt_abskew = 2.0;
        }
    }

  /* set default opt_minseqlength depending on command */

  if (parameters.opt_minseqlength < 0)
    {
      if ((parameters.opt_cluster_fast != nullptr) or
          (parameters.opt_cluster_size != nullptr) or
          (parameters.opt_cluster_smallmem != nullptr) or
          (parameters.opt_cluster_unoise != nullptr) or
          (parameters.opt_derep_fulllength != nullptr) or
          (parameters.opt_derep_id != nullptr) or
          (parameters.opt_derep_prefix != nullptr) or
          (parameters.opt_makeudb_usearch != nullptr) or
          (parameters.opt_sintax != nullptr) or
          (parameters.opt_usearch_global != nullptr))
        {
          opt_minseqlength = 32;
          parameters.opt_minseqlength = 32;
        }
      else
        {
          opt_minseqlength = 1;
          parameters.opt_minseqlength = 1;
        }
    }

  if (parameters.opt_sintax != nullptr)
    {
    opt_notrunclabels = 1;
    parameters.opt_notrunclabels = true;
    }
}


auto show_publication() -> void
{
  fprintf(stdout,
          "Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)\n"
          "VSEARCH: a versatile open source tool for metagenomics\n"
          "PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584\n"
          "\n");
}


auto cmd_version(struct Parameters const & parameters) -> void
{
  if (parameters.opt_quiet) { return ; }

  show_publication();

#ifdef HAVE_ZLIB_H
  printf("Compiled with support for gzip-compressed files,");
  if (gz_lib != nullptr)
    {
      printf(" and the library is loaded.\n");

      char * (*zlibVersion_p)();
      zlibVersion_p = (char * (*)()) arch_dlsym(gz_lib,
                                                "zlibVersion");
      char * gz_version = (*zlibVersion_p)();
      uLong (*zlibCompileFlags_p)();
      zlibCompileFlags_p = (uLong (*)()) arch_dlsym(gz_lib,
                                                    "zlibCompileFlags");
      uLong const flags = (*zlibCompileFlags_p)();

      printf("zlib version %s, compile flags %lx", gz_version, flags);
      static constexpr auto check_10th_bit = 1024U; // 0x0400
      if ((flags & check_10th_bit) != 0U)
        {
          printf(" (ZLIB_WINAPI)");
        }
      printf("\n");
    }
  else
    {
      printf(" but the library was not found.\n");
    }
#else
  printf("Compiled without support for gzip-compressed files.\n");
#endif

#ifdef HAVE_BZLIB_H
  printf("Compiled with support for bzip2-compressed files,");
  if (bz2_lib != nullptr)
    {
      printf(" and the library is loaded.\n");
    }
  else
    {
      printf(" but the library was not found.\n");
    }
#else
  printf("Compiled without support for bzip2-compressed files.\n");
#endif
}


auto cmd_help(struct Parameters const & parameters) -> void {
  if (parameters.opt_quiet) { return ; }

  show_publication();

  /*       0         1         2         3         4         5         6         7          */
  /*       01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
  fprintf(stdout,
          "Usage: %s [OPTIONS]\n", parameters.progname);

  fprintf(stdout,
          "\n"
          "For further details, please consult the manual by entering: man vsearch\n"
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
          "Chimera detection with new algorithm\n"
          "  --chimeras_denovo FILENAME  detect chimeras de novo in long exact sequences\n"
          " Parameters\n"
          "  --abskew REAL               minimum abundance ratio (1.0)\n"
          "  --chimeras_diff_pct REAL    mismatch %% allowed in each chimeric region (0.0)\n"
          "  --chimeras_length_min INT   minimum length of each chimeric region (10)\n"
          "  --chimeras_parents_max INT  maximum number of parent sequences (3)\n"
          "  --chimeras_parts INT        number of parts to divide sequences (length/100)\n"
          "  --sizein                    propagate abundance annotation from input\n"
          " Output\n"
          "  --alignwidth INT            width of alignments in alignment output file (60)\n"
          "  --alnout FILENAME           output chimera alignments to file\n"
          "  --chimeras FILENAME         output chimeric sequences to file\n"
          "  --nonchimeras FILENAME      output non-chimeric sequences to file\n"
          "  --relabel STRING            relabel nonchimeras with this prefix string\n"
          "  --relabel_keep              keep the old label after the new when relabelling\n"
          "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
          "  --relabel_self              relabel with the sequence itself as label\n"
          "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
          "  --sizeout                   include abundance information when relabelling\n"
          "  --tabbedout FILENAME        output chimera info to tab-separated file\n"
          "  --xsize                     strip abundance information in output\n"
          "\n"
          "Chimera detection with UCHIME algorithms\n"
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
          "  --fasta_score               include chimera score in FASTA output\n"
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
          "  --derep_id FILENAME         dereplicate using both identifiers and sequences\n"
          "  --derep_prefix FILENAME     dereplicate sequences in file based on prefixes\n"
          "  --derep_smallmem FILENAME   dereplicate sequences in file using less memory\n"
          "  --fastx_uniques FILENAME    dereplicate sequences in the FASTA/FASTQ file\n"
          "  --rereplicate FILENAME      rereplicate sequences in the given FASTA file\n"
          " Parameters\n"
          "  --maxuniquesize INT         maximum abundance for output from dereplication\n"
          "  --minuniquesize INT         minimum abundance for output from dereplication\n"
          "  --sizein                    propagate abundance annotation from input\n"
          "  --strand plus|both          dereplicate plus or both strands (plus)\n"
          " Output\n"
          "  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)\n"
          "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
          "  --fastq_qmaxout INT         maximum base quality value for FASTQ output (41)\n"
          "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
          "  --fastq_qminout INT         minimum base quality value for FASTQ output (0)\n"
          "  --fastaout FILENAME         output FASTA file (for fastx_uniques)\n"
          "  --fastqout FILENAME         output FASTQ file (for fastx_uniques)\n"
          "  --output FILENAME           output FASTA file (not for fastx_uniques)\n"
          "  --relabel STRING            relabel with this prefix string\n"
          "  --relabel_keep              keep the old label after the new when relabelling\n"
          "  --relabel_md5               relabel with md5 digest of normalized sequence\n"
          "  --relabel_self              relabel with the sequence itself as label\n"
          "  --relabel_sha1              relabel with sha1 digest of normalized sequence\n"
          "  --sizeout                   write abundance annotation to output\n"
          "  --tabbedout FILENAME        write cluster info to tsv file for fastx_uniques\n"
          "  --topn INT                  output only n most abundant sequences after derep\n"
          "  --uc FILENAME               filename for UCLUST-like dereplication output\n"
          "  --xsize                     strip abundance information in derep output\n"
          "\n"
          "FASTA to FASTQ conversion\n"
          "  --fasta2fastq FILENAME      convert from FASTA to FASTQ, fake quality scores\n"
          " Parameters\n"
          "  --fastq_asciiout INT        FASTQ output quality score ASCII base char (33)\n"
          "  --fastq_qmaxout INT         fake quality score for FASTQ output (41)\n"
          " Output\n"
          "  --fastqout FILENAME         FASTQ output filename for converted sequences\n"
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
          "Orient sequences in forward or reverse direction\n"
          "  --orient FILENAME           orient sequences in given FASTA/FASTQ file\n"
          " Data\n"
          "  --db FILENAME               database of sequences in correct orientation\n"
          "  --dbmask none|dust|soft     mask db seqs with dust, soft or no method (dust)\n"
          "  --qmask none|dust|soft      mask query with dust, soft or no method (dust)\n"
          "  --wordlength INT            length of words used for matching 3-15 (12)\n"
          " Output\n"
          "  --fastaout FILENAME         FASTA output filename for oriented sequences\n"
          "  --fastqout FILENAME         FASTQ output filenamr for oriented sequences\n"
          "  --notmatched FILENAME       output filename for undetermined sequences\n"
          "  --tabbedout FILENAME        output filename for result information\n"
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
          "  --label_suffix STRING       suffix to append to label of merged sequences\n"
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
          "  --lca_cutoff REAL           fraction of matching hits required for LCA (1.0)\n"
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
          "  --n_mismatch                consider aligning with N's as mismatches\n"
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
          "  --lcaout FILENAME           output LCA of matching sequences to file\n"
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
          "  --randseed INT              seed for PRNG, zero to use random data source (0)\n"
          "  --sintax_cutoff REAL        confidence value cutoff level (0.0)\n"
          "  --sintax_random             use random sequence, not shortest, if equal match\n"
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
          "  --fastq_minqual INT         discard if any base quality value lower (0)\n"
          "  --fastq_qmax INT            maximum base quality value for FASTQ input (41)\n"
          "  --fastq_qmin INT            minimum base quality value for FASTQ input (0)\n"
          "  --fastq_stripleft INT       delete given number of bases from the 5' end\n"
          "  --fastq_stripright INT      delete given number of bases from the 3' end\n"
          "  --fastq_truncee REAL        truncate to given maximum expected error\n"
          "  --fastq_truncee_rate REAL   truncate to given maximum expected error rate\n"
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


auto cmd_allpairs_global(struct Parameters const & parameters) -> void
{
  /* check options */

  if ((opt_alnout == nullptr) and (opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (opt_blast6out == nullptr) and
      (opt_matched == nullptr) and (opt_notmatched == nullptr) and
      (opt_samout == nullptr) and (opt_fastapairs == nullptr))
    {
      fatal("No output files specified");
    }

  if (not ((opt_acceptall != 0) or ((opt_id >= 0.0) and (opt_id <= 1.0))))
    {
      fatal("Specify either --acceptall or --id with an identity from 0.0 to 1.0");
    }

  allpairs_global(parameters, cmdline, prog_header.data());
}


auto cmd_usearch_global(struct Parameters const & parameters) -> void
{
  /* check options */

  if ((opt_alnout == nullptr) and (opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (opt_blast6out == nullptr) and
      (opt_matched == nullptr) and (opt_notmatched == nullptr) and
      (parameters.opt_dbmatched == nullptr) and (parameters.opt_dbnotmatched == nullptr) and
      (opt_samout == nullptr) and (opt_otutabout == nullptr) and
      (opt_biomout == nullptr) and (opt_mothur_shared_out == nullptr) and
      (opt_fastapairs == nullptr) and (opt_lcaout == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_db == nullptr)
    {
      fatal("Database filename not specified with --db");
    }

  if ((opt_id < 0.0) or (opt_id > 1.0))
    {
      fatal("Identity between 0.0 and 1.0 must be specified with --id");
    }

  usearch_global(parameters, cmdline, prog_header.data());
}


auto cmd_search_exact(struct Parameters const & parameters) -> void
{
  /* check options */

  if ((opt_alnout == nullptr) and (opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (opt_blast6out == nullptr) and
      (opt_matched == nullptr) and (opt_notmatched == nullptr) and
      (parameters.opt_dbmatched == nullptr) and (parameters.opt_dbnotmatched == nullptr) and
      (opt_samout == nullptr) and (opt_otutabout == nullptr) and
      (opt_biomout == nullptr) and (opt_mothur_shared_out == nullptr) and
      (opt_fastapairs == nullptr) and (opt_lcaout == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_db == nullptr)
    {
      fatal("Database filename not specified with --db");
    }

  search_exact(parameters, cmdline, prog_header.data());
}


auto cmd_subsample(struct Parameters const & parameters) -> void
{
  if ((opt_fastaout == nullptr) and (opt_fastqout == nullptr))
    {
      fatal("Specify output files for subsampling with --fastaout and/or --fastqout");
    }

  if ((opt_sample_pct > 0) == (opt_sample_size > 0))
    {
      fatal("Specify either --sample_pct or --sample_size");
    }

  subsample(parameters);
}


auto cmd_none(struct Parameters const & parameters) -> void {
  if (parameters.opt_quiet) { return ; }
  fprintf(stderr,
          "For more help, please enter: %s --help\n"
          "For further details, please consult the manual by entering: man vsearch\n"
          "\n"
          "Selected command examples:\n"
          "\n"
          "vsearch --allpairs_global FILENAME --id 0.5 --alnout FILENAME\n"
          "vsearch --cluster_size FILENAME --id 0.97 --centroids FILENAME\n"
          "vsearch --cut FILENAME --cut_pattern G^AATT_C --fastaout FILENAME\n"
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
          "vsearch --fastx_uniques FILENAME --fastaout FILENAME\n"
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
          "Other commands: cluster_fast, cluster_smallmem, cluster_unoise, cut,\n"
          "                derep_id, derep_fulllength, derep_prefix, derep_smallmem,\n"
          "                fasta2fastq, fastq_filter, fastq_join, fastx_getseqs,\n"
          "                fastx_getsubseq, maskfasta, orient, rereplicate, uchime2_denovo,\n"
          "                uchime3_denovo, udb2fasta, udbinfo, udbstats, version\n"
          "\n",
          parameters.progname);
}


auto cmd_cluster(struct Parameters const & parameters) -> void
{
  if ((opt_alnout == nullptr) and (opt_userout == nullptr) and
      (parameters.opt_uc == nullptr) and (opt_blast6out == nullptr) and
      (opt_matched == nullptr) and (opt_notmatched == nullptr) and
      (opt_centroids == nullptr) and (opt_clusters == nullptr) and
      (opt_consout == nullptr) and (opt_msaout == nullptr) and
      (opt_samout == nullptr) and (opt_profile == nullptr) and
      (opt_otutabout == nullptr) and (opt_biomout == nullptr) and
      (opt_mothur_shared_out == nullptr))
    {
      fatal("No output files specified");
    }

  if (parameters.opt_cluster_unoise == nullptr)
    {
      if ((opt_id < 0.0) or (opt_id > 1.0))
        {
          fatal("Identity between 0.0 and 1.0 must be specified with --id");
        }
    }

  if (parameters.opt_cluster_fast != nullptr)
    {
      cluster_fast(cmdline, prog_header.data());
    }
  else if (parameters.opt_cluster_smallmem != nullptr)
    {
      cluster_smallmem(cmdline, prog_header.data());
    }
  else if (parameters.opt_cluster_size != nullptr)
    {
      cluster_size(cmdline, prog_header.data());
    }
  else if (parameters.opt_cluster_unoise != nullptr)
    {
      cluster_unoise(cmdline, prog_header.data());
    }
}


auto cmd_chimera(struct Parameters const & parameters) -> void
{
  if ((opt_chimeras == nullptr)  and (opt_nonchimeras == nullptr) and
      (opt_uchimeout == nullptr) and (opt_uchimealns == nullptr))
    {
      fatal("No output files specified");
    }

  if ((parameters.opt_uchime_ref != nullptr) and (parameters.opt_db == nullptr))
    {
      fatal("Database filename not specified with --db");
    }

  if (opt_abskew < 1.0)
    {
      fatal("Argument to --abskew must be >= 1.0");
    }

  if (opt_xn <= 1.0)
    {
      fatal("Argument to --xn must be > 1");
    }

  if (opt_dn <= 0.0)
    {
      fatal("Argument to --dn must be > 0");
    }

  if ((parameters.opt_uchime2_denovo == nullptr) and (parameters.opt_uchime3_denovo == nullptr))
    {
      if (opt_mindiffs <= 0)
        {
          fatal("Argument to --mindiffs must be > 0");
        }

      if (opt_mindiv <= 0.0)
        {
          fatal("Argument to --mindiv must be > 0");
        }

      if (opt_minh <= 0.0)
        {
          fatal("Argument to --minh must be > 0");
        }
    }

  chimera(parameters);
}


auto cmd_fastq_mergepairs(struct Parameters const & parameters) -> void
{
  if (parameters.opt_reverse == nullptr)
    {
      fatal("No reverse reads file specified with --reverse");
    }
  if ((parameters.opt_fastqout == nullptr) and
      (parameters.opt_fastaout == nullptr) and
      (opt_fastqout_notmerged_fwd == nullptr) and
      (opt_fastqout_notmerged_rev == nullptr) and
      (opt_fastaout_notmerged_fwd == nullptr) and
      (opt_fastaout_notmerged_rev == nullptr) and
      (opt_eetabbedout == nullptr))
    {
      fatal("No output files specified");
    }
  if (opt_fastq_maxdiffs < 0) {
    fatal("Argument to --fastq_maxdiffs must be positive");
  }
  fastq_mergepairs(parameters);
}


auto fill_prog_header() -> void
{
  static constexpr auto one_gigabyte = double{1024 * 1024 * 1024};
  auto const * const format = "%s v%s_%s, %.1fGB RAM, %ld cores";
  static_cast<void>(snprintf(
      prog_header.data(), max_line_length, format, PROG_NAME, PROG_VERSION,
      PROG_ARCH, static_cast<double>(arch_get_memtotal()) / one_gigabyte,
      arch_get_cores()));
}


auto getentirecommandline(int argc, char** argv) -> void
{
  int len = 0;
  for (int i = 0; i < argc; i++)
    {
      len += std::strlen(argv[i]);
    }

  cmdline = (char *) xmalloc(len + argc);
  cmdline[0] = 0;

  for (int i = 0; i < argc; i++)
    {
      if (i > 0)
        {
          std::strcat(cmdline, " ");
        }
      std::strcat(cmdline, argv[i]);
    }
}


auto show_header(struct Parameters const & parameters) -> void {
  if (parameters.opt_quiet) { return ; }
  fprintf(stderr, "%s\n", prog_header.data());
  fprintf(stderr, "https://github.com/torognes/vsearch\n");
  fprintf(stderr, "\n");
}


auto main(int argc, char** argv) -> int
{
  struct Parameters parameters;

  fill_prog_header();

  getentirecommandline(argc, argv);
  parameters.command_line = std::string{cmdline, std::strlen(cmdline)};

  cpu_features_detect();

  args_init(argc, argv, parameters);

  if (parameters.opt_log != nullptr)
    {
      fp_log = fopen_output(opt_log);
      parameters.fp_log = fp_log;
      if (fp_log == nullptr)
        {
          fatal("Unable to open log file for writing");
        }
      fprintf(fp_log, "%s\n", prog_header.data());
      fprintf(fp_log, "%s\n", cmdline);

      std::array<char, 26> time_string {{}};
      time_start = time(nullptr);
      struct tm * tm_start = localtime(& time_start);
      strftime(time_string.data(), time_string.size(), "%Y-%m-%dT%H:%M:%S", tm_start);
      fprintf(fp_log, "Started  %s\n", time_string.data());
    }

  random_init();

  show_header(parameters);

  dynlibs_open();

#ifdef __x86_64__
  if (sse2_present == 0)
    {
      fatal("Sorry, this program requires a cpu with SSE2.");
    }
#endif

  if (parameters.opt_help)
    {
      cmd_help(parameters);
    }
  else if (parameters.opt_allpairs_global != nullptr)
    {
      opt_strand = 1;
      parameters.opt_strand = false;
      opt_uc_allhits = 1;
      parameters.opt_uc_allhits = true;
      cmd_allpairs_global(parameters);
    }
  else if (parameters.opt_usearch_global != nullptr)
    {
      cmd_usearch_global(parameters);
    }
  else if (parameters.opt_sortbysize != nullptr)
    {
      sortbysize(parameters);
    }
  else if (parameters.opt_sortbylength != nullptr)
    {
      sortbylength(parameters);
    }
  else if (parameters.opt_derep_fulllength != nullptr)
    {
      derep(parameters, parameters.opt_derep_fulllength, false);
    }
  else if (parameters.opt_derep_prefix != nullptr)
    {
      derep_prefix(parameters);
    }
  else if (parameters.opt_derep_smallmem != nullptr)
    {
      derep_smallmem(parameters);
    }
  else if (parameters.opt_derep_id != nullptr)
    {
      derep(parameters, parameters.opt_derep_id, true);
    }
  else if (parameters.opt_shuffle != nullptr)
    {
      shuffle(parameters);
    }
  else if (parameters.opt_fastx_subsample != nullptr)
    {
      cmd_subsample(parameters);
    }
  else if (parameters.opt_maskfasta != nullptr)
    {
      maskfasta(parameters);
    }
  else if ((parameters.opt_cluster_smallmem != nullptr) or (parameters.opt_cluster_fast != nullptr) or (parameters.opt_cluster_size != nullptr) or (parameters.opt_cluster_unoise != nullptr))
    {
      cmd_cluster(parameters);
    }
  else if ((parameters.opt_uchime_denovo != nullptr) or (parameters.opt_uchime_ref != nullptr) or (parameters.opt_uchime2_denovo != nullptr) or (parameters.opt_uchime3_denovo != nullptr) or (parameters.opt_chimeras_denovo != nullptr))
    {
      cmd_chimera(parameters);
    }
  else if (parameters.opt_fastq_chars != nullptr)
    {
      fastq_chars(parameters);
    }
  else if (parameters.opt_fastq_stats != nullptr)
    {
      fastq_stats(parameters);
    }
  else if (parameters.opt_fastq_filter != nullptr)
    {
      fastq_filter(parameters);
    }
  else if (parameters.opt_fastx_filter != nullptr)
    {
      fastx_filter(parameters);
    }
  else if (parameters.opt_fastx_revcomp != nullptr)
    {
      fastx_revcomp(parameters);
    }
  else if (parameters.opt_search_exact != nullptr)
    {
      opt_id = 1.0;
      cmd_search_exact(parameters);
    }
  else if (parameters.opt_fastx_mask != nullptr)
    {
      fastx_mask(parameters);
    }
  else if (parameters.opt_fastq_convert != nullptr)
    {
      fastq_convert(parameters);
    }
  else if (parameters.opt_fastq_mergepairs != nullptr)
    {
      cmd_fastq_mergepairs(parameters);
    }
  else if (parameters.opt_fastq_eestats != nullptr)
    {
      fastq_eestats(parameters);
    }
  else if (parameters.opt_fastq_eestats2 != nullptr)
    {
      fastq_eestats2(parameters);
    }
  else if (parameters.opt_fastq_join != nullptr)
    {
      if ((not parameters.opt_join_padgapq_set_by_user) and
          (parameters.opt_fastq_ascii != default_ascii_offset)) {
        parameters.opt_join_padgapq = alternative_quality_padding;
      }
      fastq_join(parameters);
    }
  else if (parameters.opt_rereplicate != nullptr)
    {
      opt_xsize = true;
      parameters.opt_xsize = true;
      rereplicate(parameters);
    }
  else if (parameters.opt_version)
    {
      cmd_version(parameters);
    }
  else if (parameters.opt_makeudb_usearch != nullptr)
    {
      udb_make(parameters);
    }
  else if (parameters.opt_udb2fasta != nullptr)
    {
      udb_fasta(parameters);
    }
  else if (parameters.opt_udbinfo != nullptr)
    {
      udb_info(parameters);
    }
  else if (parameters.opt_udbstats != nullptr)
    {
      udb_stats(parameters);
    }
  else if (parameters.opt_sintax != nullptr)
    {
      sintax(parameters);
    }
  else if (parameters.opt_sff_convert != nullptr)
    {
      sff_convert(parameters);
    }
  else if (parameters.opt_fastx_getseq != nullptr)
    {
      fastx_getseq(parameters);
    }
  else if (parameters.opt_fastx_getseqs != nullptr)
    {
      fastx_getseqs(parameters);
    }
  else if (parameters.opt_fastx_getsubseq != nullptr)
    {
      fastx_getsubseq(parameters);
    }
  else if (parameters.opt_cut != nullptr)
    {
      cut(parameters);
    }
  else if (parameters.opt_orient != nullptr)
    {
      orient(parameters);
    }
  else if (parameters.opt_fasta2fastq != nullptr)
    {
      fasta2fastq(parameters);
    }
  else if (parameters.opt_fastx_uniques != nullptr)
    {
      derep(parameters, parameters.opt_fastx_uniques, false);
    }
  else
    {
      cmd_none(parameters);
    }

  if (parameters.opt_log != nullptr)
    {
      time_finish = time(nullptr);
      struct tm * tm_finish = localtime(& time_finish);
      std::array<char, 26> time_string {{}};
      strftime(time_string.data(), time_string.size(), "%Y-%m-%dT%H:%M:%S", tm_finish);
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Finished %s", time_string.data());

      double const time_diff = difftime(time_finish, time_start);
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Elapsed time %02.0lf:%02.0lf\n",
              floor(time_diff / 60.0),
              floor(time_diff - (60.0 * floor(time_diff / 60.0))));
      double const maxmem = arch_get_memused() / 1048576.0;
      if (maxmem < 1024.0)
        {
          fprintf(fp_log, "Max memory %.1lfMB\n", maxmem);
        }
      else
        {
          fprintf(fp_log, "Max memory %.1lfGB\n", maxmem / 1024.0);
        }
      fclose(fp_log);
    }

  if (opt_ee_cutoffs_values != nullptr)
    {
      xfree(opt_ee_cutoffs_values);
    }
  opt_ee_cutoffs_values = nullptr;

  xfree(cmdline);
  dynlibs_close();
}
