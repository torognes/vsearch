/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#define _GNU_SOURCE 1
#define __STDC_CONSTANT_MACROS 1
#define __STDC_FORMAT_MACROS 1
#define __STDC_LIMIT_MACROS 1
#define __restrict

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cctype>
#include <cfloat>
#include <cinttypes>
#include <climits>
#include <clocale>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <set>
#include <string>
#include <cassert>
#include <limits>

/* include appropriate regex library */

#ifdef HAVE_REGEX_H
#include <regex.h>
#else
#include <regex>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <pthread.h>
#include <getopt.h>
#include <fcntl.h>
#include <unistd.h>

#define PROG_NAME PACKAGE
#define PROG_VERSION PACKAGE_VERSION

#ifdef __x86_64__

#define PROG_CPU "x86_64"
#include <x86intrin.h>

#elif __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#include <altivec.h>
#undef bool
#else
#error Big endian ppc64 CPUs not supported
#endif

#elif __aarch64__

#define PROG_CPU "aarch64"
#include <arm_neon.h>

#else

#define PROG_CPU "simde"
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/avx512.h>
#endif


#ifdef _WIN32

#define PROG_OS "win"
#include <windows.h>
#include <psapi.h>
#include <shlwapi.h>
#define bswap_16(x) _byteswap_ushort(x)
#define bswap_32(x) _byteswap_ulong(x)
#define bswap_64(x) _byteswap_uint64(x)

#elif __APPLE__

#define PROG_OS "macos"
#include <sys/sysctl.h>
#include <libkern/OSByteOrder.h>
#include <sys/resource.h>
#define bswap_16(x) OSSwapInt16(x)
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#elif __linux__

#define PROG_OS "linux"
#include <sys/sysinfo.h>
#include <byteswap.h>
#include <sys/resource.h>

#elif __FreeBSD__

#define PROG_OS "freebsd"
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <sys/endian.h>
#define bswap_16(x) bswap16(x)
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)

#elif __NetBSD__

#define PROG_OS "netbsd"
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/bswap.h>
#define bswap_16(x) bswap16(x)
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)
/* Alters behavior, but NetBSD 7 does not have getopt_long_only() */
#define getopt_long_only getopt_long

#else

#define PROG_OS "unknown"
#include <sys/sysinfo.h>
#include <byteswap.h>
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
#include "sha1.h"

#include "arch.h"
#include "util.h"
#include "xstring.h"
#include "db.h"
#include "linmemalign.h"
#include "searchcore.h"
#include "results.h"
#include "cpu.h"
#include "fastx.h"
#include "fasta.h"
#include "fastq.h"
#include "dbhash.h"
#include "kmerhash.h"
#include "pcr.h"

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
extern bool opt_label_substr_match;
extern bool opt_lengthout;
extern bool opt_n_mismatch;
extern bool opt_no_progress;
extern bool opt_quiet;
extern bool opt_relabel_keep;
extern bool opt_relabel_md5;
extern bool opt_relabel_self;
extern bool opt_relabel_sha1;
extern bool opt_samheader;
extern bool opt_sff_clip;
extern bool opt_sintax_random;
extern bool opt_sizein;
extern bool opt_sizeorder;
extern bool opt_sizeout;
extern bool opt_xee;
extern bool opt_xlength;
extern bool opt_xsize;
extern char * opt_allpairs_global;
extern char * opt_alnout;
extern char * opt_biomout;
extern char * opt_blast6out;
extern char * opt_borderline;
extern char * opt_centroids;
extern char * opt_chimeras;
extern char * opt_chimeras_denovo;
extern char * opt_cluster_fast;
extern char * opt_cluster_size;
extern char * opt_cluster_smallmem;
extern char * opt_cluster_unoise;
extern char * opt_clusters;
extern char * opt_consout;
extern char * opt_db;
extern char * opt_dbmatched;
extern char * opt_dbnotmatched;
extern char * opt_eetabbedout;
extern char * opt_fastaout;
extern char * opt_fastaout_discarded;
extern char * opt_fastaout_discarded_rev;
extern char * opt_fastaout_notmerged_fwd;
extern char * opt_fastaout_notmerged_rev;
extern char * opt_fastaout_rev;
extern char * opt_fastapairs;
extern char * opt_fastq_convert;
extern char * opt_fastq_eestats2;
extern char * opt_fastq_eestats;
extern char * opt_fastq_filter;
extern char * opt_fastq_mergepairs;
extern char * opt_fastq_stats;
extern char * opt_fastqout;
extern char * opt_fastqout_discarded;
extern char * opt_fastqout_discarded_rev;
extern char * opt_fastqout_notmerged_fwd;
extern char * opt_fastqout_notmerged_rev;
extern char * opt_fastqout_rev;
extern char * opt_fastx_filter;
extern char * opt_fastx_getseq;
extern char * opt_fastx_getseqs;
extern char * opt_fastx_getsubseq;
extern char * opt_fastx_mask;
extern char * opt_fastx_revcomp;
extern char * opt_label;
extern char * opt_label_field;
extern char * opt_label_suffix;
extern char * opt_label_word;
extern char * opt_label_words;
extern char * opt_labels;
extern char * opt_lcaout;
extern char * opt_log;
extern char * opt_makeudb_usearch;
extern char * opt_maskfasta;
extern char * opt_matched;
extern char * opt_mothur_shared_out;
extern char * opt_msaout;
extern char * opt_nonchimeras;
extern char * opt_notmatched;
extern char * opt_notmatchedfq;
extern char * opt_orient;
extern char * opt_otutabout;
extern char * opt_output;
extern char * opt_pattern;
extern char * opt_pcr_sim;
extern char * opt_profile;
extern char * opt_qsegout;
extern char * opt_relabel;
extern char * opt_reverse;
extern char * opt_samout;
extern char * opt_sample;
extern char * opt_search_exact;
extern char * opt_sff_convert;
extern char * opt_sintax;
extern char * opt_tabbedout;
extern char * opt_tsegout;
extern char * opt_uc;
extern char * opt_uchime2_denovo;
extern char * opt_uchime3_denovo;
extern char * opt_uchime_denovo;
extern char * opt_uchime_ref;
extern char * opt_uchimealns;
extern char * opt_uchimeout;
extern char * opt_udb2fasta;
extern char * opt_udbinfo;
extern char * opt_udbstats;
extern char * opt_usearch_global;
extern char * opt_userout;
extern double * opt_ee_cutoffs_values;
extern double opt_abskew;
extern double opt_chimeras_diff_pct;
extern double opt_dn;
extern double opt_fastq_maxdiffpct;
extern double opt_fastq_maxee;
extern double opt_fastq_maxee_rate;
extern double opt_fastq_truncee;
extern double opt_id;
extern double opt_lca_cutoff;
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
extern double opt_pcr_chimera_p;
extern double opt_pcr_subst_p;
extern double opt_query_cov;
extern double opt_sample_pct;
extern double opt_sintax_cutoff;
extern double opt_target_cov;
extern double opt_unoise_alpha;
extern double opt_weak_id;
extern double opt_xn;
extern int opt_acceptall;
extern int opt_alignwidth;
extern int opt_chimeras_length_min;
extern int opt_chimeras_parents_max;
extern int opt_chimeras_parts;
extern int opt_cons_truncate;
extern int opt_ee_cutoffs_count;
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
extern int opt_length_cutoffs_increment;
extern int opt_length_cutoffs_longest;
extern int opt_length_cutoffs_shortest;
extern int opt_mindiffs;
extern int opt_slots;
extern int opt_uchimeout5;
extern int opt_usersort;
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
extern int64_t opt_pcr_cycles;
extern int64_t opt_qmask;
extern int64_t opt_randseed;
extern int64_t opt_rightjust;
extern int64_t opt_rowlen;
extern int64_t opt_sample_size;
extern int64_t opt_self;
extern int64_t opt_selfid;
extern int64_t opt_strand;
extern int64_t opt_subseq_end;
extern int64_t opt_subseq_start;
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

extern std::FILE * fp_log;

constexpr auto tax_levels = 9;
constexpr int64_t default_maxseqlength = 50000;
constexpr int64_t default_ascii_offset = 33;
constexpr char alternative_ascii_offset = 64;
constexpr int64_t default_max_quality = 41;
constexpr auto int64_max = std::numeric_limits<int64_t>::max();
std::string const default_quality_padding = "IIIIIIII";  // Q40 with an offset of 33
std::string const alternative_quality_padding = "hhhhhhhh";  // Q40 with an offset of 64
std::string const default_sequence_padding = "NNNNNNNN";

struct Parameters {
  char * opt_cut = nullptr;
  std::string opt_cut_pattern;
  char * opt_derep_fulllength = nullptr;
  char * opt_derep_id = nullptr;
  char * opt_derep_prefix = nullptr;
  char * opt_derep_smallmem = nullptr;
  char * opt_fasta2fastq = nullptr;
  char * opt_fastaout = nullptr;
  char * opt_fastaout_rev = nullptr;
  char * opt_fastaout_discarded = nullptr;
  char * opt_fastaout_discarded_rev = nullptr;
  char * opt_fastq_chars = nullptr;
  char * opt_fastq_join = nullptr;
  char * opt_fastqout = nullptr;
  char * opt_fastqout_rev = nullptr;
  char * opt_fastqout_discarded = nullptr;
  char * opt_fastqout_discarded_rev = nullptr;
  char * opt_fastx_subsample = nullptr;
  char * opt_fastx_uniques = nullptr;
  std::string opt_join_padgap = default_sequence_padding;
  std::string opt_join_padgapq = default_quality_padding;
  char * opt_log = nullptr;
  char * opt_output = nullptr;
  char * opt_relabel = nullptr;
  char * opt_rereplicate = nullptr;
  char * opt_reverse = nullptr;
  char * opt_shuffle = nullptr;
  char * opt_sortbylength = nullptr;
  char * opt_sortbysize = nullptr;
  char * opt_tabbedout = nullptr;
  char * opt_uc = nullptr;
  char * progname = nullptr;
  std::FILE * fp_log = nullptr;
  double opt_sample_pct = 0;
  int64_t opt_fastq_ascii = default_ascii_offset;
  int64_t opt_fastq_asciiout = default_ascii_offset;
  int64_t opt_fastq_qmaxout = default_max_quality;
  int64_t opt_fastq_qminout = 0;
  int64_t opt_fastq_tail = 4;
  int64_t opt_maxseqlength = default_maxseqlength;
  int64_t opt_maxsize = int64_max;
  int64_t opt_maxuniquesize = int64_max;
  int64_t opt_minseqlength = -1;
  int64_t opt_minsize = 0;
  int64_t opt_minuniquesize = 1;
  int64_t opt_randseed = 0;
  int64_t opt_sample_size = 0;
  int64_t opt_threads = 0;
  int64_t opt_topn = int64_max;
  bool opt_fastq_qout_max = false;
  bool opt_help = false;
  bool opt_join_padgapq_set_by_user = false;
  bool opt_notrunclabels = false;
  bool opt_quiet = false;
  bool opt_sizein = false;
  bool opt_strand = false;
  bool opt_version = false;
  bool opt_xsize = false;
};
