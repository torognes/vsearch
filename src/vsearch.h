/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#define __STDC_FORMAT_MACROS 1
#define __restrict

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <regex.h>

#include <string>
#include <set>
#include <map>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <float.h>
#include <inttypes.h>
#include <stdarg.h>

#define PROG_NAME PACKAGE
#define PROG_VERSION PACKAGE_VERSION

#ifdef __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#else

#define PROG_CPU "x86_64"

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

#endif


#ifdef _WIN32

#define PROG_OS "win"
#include <windows.h>
#include <psapi.h>

#else

#ifdef __APPLE__

#define PROG_OS "macos"
#include <sys/sysctl.h>

#else

#ifdef __linux__
#define PROG_OS "linux"
#else
#define PROG_OS "unknown"
#endif

#include <sys/sysinfo.h>

#endif

#include <sys/resource.h>

#endif


#define PROG_ARCH PROG_OS "_" PROG_CPU

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB_H
#include <bzlib.h>
#endif

#include "city.h"
#include "md5.h"
#include "sha1.h"

#include "arch.h"
#include "dynlibs.h"
#include "util.h"
#include "xstring.h"
#include "align_simd.h"
#include "maps.h"
#include "abundance.h"
#include "db.h"
#include "align.h"
#include "unique.h"
#include "bitmap.h"
#include "dbindex.h"
#include "minheap.h"
#include "search.h"
#include "linmemalign.h"
#include "searchcore.h"
#include "showalign.h"
#include "userfields.h"
#include "results.h"
#include "sortbysize.h"
#include "sortbylength.h"
#include "derep.h"
#include "shuffle.h"
#include "mask.h"
#include "cluster.h"
#include "msa.h"
#include "chimera.h"
#include "cpu.h"
#include "allpairs.h"
#include "subsample.h"
#include "fastx.h"
#include "fasta.h"
#include "fastq.h"
#include "fastqops.h"
#include "dbhash.h"
#include "searchexact.h"
#include "mergepairs.h"
#include "eestats.h"
#include "rerep.h"
#include "otutable.h"
#include "udb.h"
#include "kmerhash.h"

/* options */

extern bool opt_bzip2_decompress;
extern bool opt_clusterout_id;
extern bool opt_clusterout_sort;
extern bool opt_eeout;
extern bool opt_fasta_score;
extern bool opt_fastq_allowmergestagger;
extern bool opt_fastq_eeout;
extern bool opt_fastq_nostagger;
extern bool opt_gzip_decompress;
extern bool opt_no_progress;
extern bool opt_quiet;
extern bool opt_relabel_keep;
extern bool opt_relabel_md5;
extern bool opt_relabel_sha1;
extern bool opt_samheader;
extern bool opt_sizeorder;
extern bool opt_xsize;
extern char * opt_allpairs_global;
extern char * opt_alnout;
extern char * opt_blast6out;
extern char * opt_biomout;
extern char * opt_borderline;
extern char * opt_centroids;
extern char * opt_chimeras;
extern char * opt_cluster_fast;
extern char * opt_cluster_size;
extern char * opt_cluster_smallmem;
extern char * opt_cluster_unoise;
extern char * opt_clusters;
extern char * opt_consout;
extern char * opt_db;
extern char * opt_dbmatched;
extern char * opt_dbnotmatched;
extern char * opt_derep_fulllength;
extern char * opt_derep_prefix;
extern char * opt_eetabbedout;
extern char * opt_fastaout;
extern char * opt_fastaout_discarded;
extern char * opt_fastaout_notmerged_fwd;
extern char * opt_fastaout_notmerged_rev;
extern char * opt_fastapairs;
extern char * opt_fastq_chars;
extern char * opt_fastq_convert;
extern char * opt_fastq_eestats;
extern char * opt_fastq_eestats2;
extern char * opt_fastq_filter;
extern char * opt_fastq_mergepairs;
extern char * opt_fastq_stats;
extern char * opt_fastqout;
extern char * opt_fastqout_discarded;
extern char * opt_fastqout_notmerged_fwd;
extern char * opt_fastqout_notmerged_rev;
extern char * opt_fastx_filter;
extern char * opt_fastx_mask;
extern char * opt_fastx_revcomp;
extern char * opt_fastx_subsample;
extern char * opt_label_suffix;
extern char * opt_log;
extern char * opt_makeudb_usearch;
extern char * opt_maskfasta;
extern char * opt_matched;
extern char * opt_mothur_shared_out;
extern char * opt_msaout;
extern char * opt_nonchimeras;
extern char * opt_notmatched;
extern char * opt_otutabout;
extern char * opt_output;
extern char * opt_pattern;
extern char * opt_profile;
extern char * opt_relabel;
extern char * opt_rereplicate;
extern char * opt_reverse;
extern char * opt_samout;
extern char * opt_search_exact;
extern char * opt_shuffle;
extern char * opt_sortbylength;
extern char * opt_sortbysize;
extern char * opt_udb2fasta;
extern char * opt_udbinfo;
extern char * opt_udbstats;
extern char * opt_uc;
extern char * opt_uchime_denovo;
extern char * opt_uchime2_denovo;
extern char * opt_uchime3_denovo;
extern char * opt_uchime_ref;
extern char * opt_uchimealns;
extern char * opt_uchimeout;
extern char * opt_usearch_global;
extern char * opt_userout;
extern double * opt_ee_cutoffs_values;
extern double opt_abskew;
extern double opt_dn;
extern double opt_fastq_maxee;
extern double opt_fastq_maxee_rate;
extern double opt_fastq_truncee;
extern double opt_id;
extern double opt_max_unmasked_pct;
extern double opt_maxid;
extern double opt_maxqt;
extern double opt_maxsizeratio;
extern double opt_maxsl;
extern double opt_mid;
extern double opt_min_unmasked_pct;
extern double opt_mindiv;
extern double opt_minh;
extern double opt_minqt;
extern double opt_minsizeratio;
extern double opt_minsl;
extern double opt_query_cov;
extern double opt_sample_pct;
extern double opt_target_cov;
extern double opt_unoise_alpha;
extern double opt_weak_id;
extern double opt_xn;
extern int opt_acceptall;
extern int opt_alignwidth;
extern int opt_cons_truncate;
extern int opt_ee_cutoffs_count;
extern int opt_length_cutoffs_increment;
extern int opt_length_cutoffs_longest;
extern int opt_length_cutoffs_shortest;
extern int opt_gap_extension_query_interior;
extern int opt_gap_extension_query_left;
extern int opt_gap_extension_query_right;
extern int opt_gap_extension_target_interior;
extern int opt_gap_extension_target_left;
extern int opt_gap_extension_target_right;
extern int opt_gap_open_query_interior;
extern int opt_gap_open_query_left;
extern int opt_gap_open_query_right;
extern int opt_gap_open_target_interior;
extern int opt_gap_open_target_left;
extern int opt_gap_open_target_right;
extern int opt_help;
extern int opt_mindiffs;
extern int opt_slots;
extern int opt_uchimeout5;
extern int opt_usersort;
extern int opt_version;
extern int64_t opt_dbmask;
extern int64_t opt_fasta_width;
extern int64_t opt_fastq_ascii;
extern int64_t opt_fastq_asciiout;
extern int64_t opt_fastq_maxdiffs;
extern int64_t opt_fastq_maxlen;
extern int64_t opt_fastq_maxmergelen;
extern int64_t opt_fastq_maxns;
extern int64_t opt_fastq_minlen;
extern int64_t opt_fastq_minmergelen;
extern int64_t opt_fastq_minovlen;
extern int64_t opt_fastq_qmax;
extern int64_t opt_fastq_qmaxout;
extern int64_t opt_fastq_qmin;
extern int64_t opt_fastq_qminout;
extern int64_t opt_fastq_stripleft;
extern int64_t opt_fastq_stripright;
extern int64_t opt_fastq_tail;
extern int64_t opt_fastq_trunclen;
extern int64_t opt_fastq_trunclen_keep;
extern int64_t opt_fastq_truncqual;
extern int64_t opt_fulldp;
extern int64_t opt_hardmask;
extern int64_t opt_iddef;
extern int64_t opt_idprefix;
extern int64_t opt_idsuffix;
extern int64_t opt_leftjust;
extern int64_t opt_match;
extern int64_t opt_maxaccepts;
extern int64_t opt_maxdiffs;
extern int64_t opt_maxgaps;
extern int64_t opt_maxhits;
extern int64_t opt_maxqsize;
extern int64_t opt_maxrejects;
extern int64_t opt_maxseqlength;
extern int64_t opt_maxsize;
extern int64_t opt_maxsubs;
extern int64_t opt_maxuniquesize;
extern int64_t opt_mincols;
extern int64_t opt_minseqlength;
extern int64_t opt_minsize;
extern int64_t opt_mintsize;
extern int64_t opt_minuniquesize;
extern int64_t opt_minwordmatches;
extern int64_t opt_mismatch;
extern int64_t opt_notrunclabels;
extern int64_t opt_output_no_hits;
extern int64_t opt_qmask;
extern int64_t opt_randseed;
extern int64_t opt_rightjust;
extern int64_t opt_rowlen;
extern int64_t opt_sample_size;
extern int64_t opt_self;
extern int64_t opt_selfid;
extern int64_t opt_sizein;
extern int64_t opt_sizeout;
extern int64_t opt_strand;
extern int64_t opt_threads;
extern int64_t opt_top_hits_only;
extern int64_t opt_topn;
extern int64_t opt_uc_allhits;
extern int64_t opt_wordlength;

extern int64_t altivec_present;
extern int64_t mmx_present;
extern int64_t sse_present;
extern int64_t sse2_present;
extern int64_t sse3_present;
extern int64_t ssse3_present;
extern int64_t sse41_present;
extern int64_t sse42_present;
extern int64_t popcnt_present;
extern int64_t avx_present;
extern int64_t avx2_present;

extern FILE * fp_log;
