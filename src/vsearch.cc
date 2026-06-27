/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "vsearch_api.h"
#include "allpairs.h"
#include "chimera.h"
#include "cli.h"
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
#include "fastq_mergepairs.h"
#include "fastq_stats.h"
#include "fastqops.h"
#include "fastx_syncpairs.h"
#include "filter.h"
#include "getseq.h"
#include "mask.h"
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
#include "utils/compare_strings_nocase.hpp"
#include "utils/cpu_features.hpp"
#include "utils/fatal.hpp"
#include <algorithm>  // std::count, std::any_of
#include <array>
#include <cerrno>  // errno, ERANGE
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::floor
#include <ctime>  // std::strftime, std::localtime, std::time, std::time_t, std::tm, std::difftime
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::size_t, std::sscanf, std::fclose, std::snprintf, std::printf, std::strcat
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <cstring>  // std::strlen, std::memset, std::strcspn
#include <getopt.h>  // getopt_long_only, optarg, optind, opterr, struct
                     // option (no_argument, required_argument)
#include <limits>
#include <mutex>  // std::mutex
#include <new>  // std::set_new_handler
#include <string>
#include <unistd.h>  // write, _exit, STDERR_FILENO
#include <vector>


constexpr auto max_line_length = std::size_t{80};

/* options */

bool opt_bzip2_decompress = false;
bool opt_centroid_sizeout = false;
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

static std::array<char, max_line_length> prog_header {{}};
static char * cmdline;

#ifndef VSEARCH_NO_MAIN
static time_t time_start;
static time_t time_finish;
#endif

std::FILE * fp_log = nullptr;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // performance: lambda compiles to a single x86-64 instruction (shr al, 7),
  // equivalent to (c & ~0x7F) == 0
  auto is_not_ASCII(std::string const & user_string) -> bool {
    static constexpr auto ascii_max = static_cast<unsigned char>(std::numeric_limits<signed char>::max());
    auto is_not_in_range = [](char const user_char) -> bool {
      auto const unsigned_user_char = static_cast<unsigned char>(user_char);
      return (unsigned_user_char > ascii_max);
    };
    return std::any_of(user_string.begin(), user_string.end(), is_not_in_range);
  }

}  // end of anonymous namespace


/* Initialize all global opt_* variables to their default values.
   This is the library-API equivalent of the defaults block in args_init().
   Must be called before using any vsearch library function.

   After calling this, the caller should override any options needed
   for their specific use case (e.g., chimera detection parameters).

   Note: allocates opt_ee_cutoffs_values via xmalloc. Caller is responsible
   for freeing it (or calling vsearch_init_defaults again, which leaks the
   old allocation — acceptable for single-init library use). */

static std::mutex session_mutex;

auto vsearch_api_version() -> int
{
  return VSEARCH_API_VERSION;
}

auto vsearch_api_version_string() -> const char *
{
  return VSEARCH_API_VERSION_STRING;
}

auto vsearch_init_defaults() -> void
{
  session_mutex.lock();
  static constexpr auto int_max = std::numeric_limits<int>::max();
  static constexpr auto long_min = std::numeric_limits<long>::min();

  opt_abskew = 0.0;
  opt_acceptall = 0;
  opt_alignwidth = 80;
  opt_alnout = nullptr;
  opt_biomout = nullptr;
  opt_blast6out = nullptr;
  opt_borderline = nullptr;
  opt_centroid_sizeout = false;
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
  if (opt_ee_cutoffs_values != nullptr)
    {
      xfree(opt_ee_cutoffs_values);
    }
  opt_ee_cutoffs_values = static_cast<double *>(xmalloc(static_cast<size_t>(opt_ee_cutoffs_count) * sizeof(double)));
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
  opt_maxqsize = int64_max;
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
  opt_no_progress = true;
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
  opt_quiet = true;
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
}


auto vsearch_session_end() -> void
{
  session_mutex.unlock();
}


/* Apply post-parsing fixups to global options.
   Resolves sentinel values (e.g., opt_minwordmatches=-1, opt_maxhits=0)
   to their computed defaults. Call after vsearch_init_defaults() and
   after setting any option overrides (e.g., opt_wordlength).

   In the CLI path, args_init() calls this implicitly at the end of
   option parsing. Library users must call it explicitly after setting
   their overrides. */
auto vsearch_apply_defaults_fixups() -> void
{
  if (opt_maxhits == 0)
    {
      opt_maxhits = int64_max;
    }

  if (opt_minwordmatches < 0)
    {
      if (opt_wordlength >= 0 and
          opt_wordlength < static_cast<int64_t>(minwordmatches_defaults.size()))
        {
          opt_minwordmatches = minwordmatches_defaults[static_cast<size_t>(opt_wordlength)];
        }
      else
        {
          opt_minwordmatches = 0;
        }
    }


  /* opt_weak_id default (10.0) is a sentinel meaning "not set by user".
     Clamp to opt_id so the acceptance check in search_acceptable_aligned
     doesn't reject everything. Only clamp when opt_id has been set to
     a real value (>= 0). For chimera detection, opt_id is set later by
     chimera_detect_init, which re-checks opt_weak_id. */
  if (opt_id >= 0.0 and opt_weak_id > opt_id)
    {
      opt_weak_id = opt_id;
    }

  /* Resolve opt_threads sentinel: 0 means "auto-detect".
     The CLI resolves this in command dispatch; library callers need it
     resolved here so that dust_all() and other threaded functions work. */
  if (opt_threads == 0)
    {
      opt_threads = arch_get_cores();
    }

  /* Resolve opt_maxrejects sentinel: -1 means "not set by user".
     CLI resolves to 8 (cluster_fast) or 32 (others). Default to 32. */
  if (opt_maxrejects < 0)
    {
      opt_maxrejects = 32;
    }

  /* Validate opt_wordlength — 0 is a sentinel meaning "not set".
     A wordlength of 0 destroys the k-mer index (hash size = 1).
     Library users MUST set this before calling fixups. */
  if (opt_wordlength == 0)
    {
      opt_wordlength = 8;
    }

  /* Adjust gap-open penalties: subtract the first-nucleotide extension cost.
     The rest of the code assumes gap_open does NOT include the first
     extension. Safe to call repeatedly — vsearch_init_defaults() resets
     opt_gap_open_* to their raw (pre-adjustment) values each time. */
  opt_gap_open_query_left -= opt_gap_extension_query_left;
  opt_gap_open_target_left -= opt_gap_extension_target_left;
  opt_gap_open_query_interior -= opt_gap_extension_query_interior;
  opt_gap_open_target_interior -= opt_gap_extension_target_interior;
  opt_gap_open_query_right -= opt_gap_extension_query_right;
  opt_gap_open_target_right -= opt_gap_extension_target_right;

  /* Note: opt_minsize is NOT resolved here — it has command-specific
     defaults (1 for most commands, 8 for cluster_unoise).
     Library users should set opt_minsize explicitly if needed. */
}


auto show_publication() -> void
{
  std::fprintf(stdout,
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
  std::printf("Compiled with support for gzip-compressed files,");
  if (gz_lib != nullptr)
    {
      std::printf(" and the library is loaded.\n");

      char * (*zlibVersion_p)();
      zlibVersion_p = reinterpret_cast<char * (*)()>(arch_dlsym(gz_lib, "zlibVersion"));
      char const * gz_version = (*zlibVersion_p)();

      long unsigned int (*zlibCompileFlags_p)();
      zlibCompileFlags_p = reinterpret_cast<long unsigned int (*)()>(arch_dlsym(gz_lib, "zlibCompileFlags"));
      long unsigned int const flags = (*zlibCompileFlags_p)();

      std::printf("zlib version %s, compile flags %lx", gz_version, flags);
      static constexpr auto check_10th_bit = 1024U; // 0x0400
      if ((flags & check_10th_bit) != 0U)
        {
          std::printf(" (ZLIB_WINAPI)");
        }
      std::printf("\n");
    }
  else
    {
      std::printf(" but the library was not found.\n");
    }
#else
  std::printf("Compiled without support for gzip-compressed files.\n");
#endif

#ifdef HAVE_BZLIB_H
  std::printf("Compiled with support for bzip2-compressed files,");
  if (bz2_lib != nullptr)
    {
      std::printf(" and the library is loaded.\n");
    }
  else
    {
      std::printf(" but the library was not found.\n");
    }
#else
  std::printf("Compiled without support for bzip2-compressed files.\n");
#endif
}


auto cmd_help(struct Parameters const & parameters) -> void {
  if (parameters.opt_quiet) { return ; }

  show_publication();

  /*       0         1         2         3         4         5         6         7          */
  /*       01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
  std::fprintf(stdout,
          "Usage: %s [OPTIONS]\n", parameters.progname);

  std::fprintf(stdout,
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
          "  --centroid_sizeout          write centroid abundances to centroid file\n"
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
          "Paired-end reads synchronizing\n"
          "  --fastx_syncpairs FILENAME  reorder paired reads to match across both files\n"
          " Data\n"
          "  --reverse FILENAME          specify FASTA/FASTQ file with reverse reads\n"
          " Parameters\n"
          "  --read_separators STRING    characters preceding the mate number 1 or 2 (/)\n"
          " Output\n"
          "  --fastaout FILENAME         FASTA filename for synced forward reads\n"
          "  --fastaout_rev FILENAME     FASTA filename for synced reverse reads\n"
          "  --fastaout_orphans FN       FASTA filename for unpaired forward reads\n"
          "  --fastaout_orphans_rev FN   FASTA filename for unpaired reverse reads\n"
          "  --fastqout FILENAME         FASTQ filename for synced forward reads\n"
          "  --fastqout_rev FILENAME     FASTQ filename for synced reverse reads\n"
          "  --fastqout_orphans FN       FASTQ filename for unpaired forward reads\n"
          "  --fastqout_orphans_rev FN   FASTQ filename for unpaired reverse reads\n"
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
          "  --db FILENAME               FASTA or UDB database (only FASTA for search_exact)\n"
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
  std::fprintf(stderr,
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
          "                fastx_getsubseq, fastx_syncpairs, maskfasta, orient, rereplicate,\n"
          "                uchime2_denovo, uchime3_denovo, udb2fasta, udbinfo, udbstats,\n"
          "                version\n"
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
      (opt_uchimeout == nullptr) and (opt_uchimealns == nullptr) and
      (opt_tabbedout == nullptr) and (opt_alnout == nullptr))
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


auto cmd_fastq_join(struct Parameters & parameters) -> void
{
  if (is_not_ASCII(parameters.opt_join_padgap)) {
    fatal("Option --join_padgap contains non-ASCII characters");
  }
  if (is_not_ASCII(parameters.opt_join_padgapq)) {
    fatal("Option --join_padgapq contains non-ASCII characters");
  }
  if ((not parameters.opt_join_padgapq_set_by_user) and
      (parameters.opt_fastq_ascii != default_ascii_offset)) {
    parameters.opt_join_padgapq = alternative_quality_padding;
  }
  fastq_join(parameters);
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
  if (opt_fastq_maxee <= 0.0) {
    /* expected error is the sum of per-base error probabilities;
       probabilities are strictly positive (min quality score is 93,
       corresponding to probability ~1e-9.3), so the sum cannot be
       zero or negative. A null or negative threshold would always
       reject every read and is almost certainly a user mistake. */
    fatal("Argument to --fastq_maxee must be a strictly positive number");
  }
  if (opt_fastq_maxlen < 1) {
    fatal("Argument to --fastq_maxlen must be a positive integer");
  }
  if (opt_fastq_minlen < 1) {
    fatal("Argument to --fastq_minlen must be a positive integer");
  }
  if (opt_fastq_maxns < 0) {
    fatal("Argument to --fastq_maxns must be a non-negative integer");
  }
  if (opt_fastq_maxmergelen < 1) {
    fatal("Argument to --fastq_maxmergelen must be a positive integer");
  }
  if (opt_fastq_minmergelen < 0) {
    fatal("Argument to --fastq_minmergelen must be a non-negative integer");
  }
  {
    /* Quality score range: 0..93 (Phred scores encoded in fastq).
       The default value is std::numeric_limits<long>::min(), meaning
       "no truncation"; skip the range check in that case so the default
       is preserved. */
    static constexpr auto long_min = std::numeric_limits<long>::min();
    if ((opt_fastq_truncqual != long_min) and
        ((opt_fastq_truncqual < 0) or (opt_fastq_truncqual > 93))) {
      fatal("Argument to --fastq_truncqual must be in range 0..93");
    }
  }
  fastq_mergepairs(parameters);
}


auto fill_prog_header() -> void
{
  static constexpr auto one_gigabyte = double{1024 * 1024 * 1024};
  auto const * const format = "%s v%s_%s, %.1fGB RAM, %ld cores";
  static_cast<void>(std::snprintf(
      prog_header.data(), max_line_length, format, PROG_NAME, PROG_VERSION,
      PROG_ARCH, static_cast<double>(arch_get_memtotal()) / one_gigabyte,
      arch_get_cores()));
}


auto getentirecommandline(int argc, char** argv) -> void
{
  size_t len = 0;
  for (int i = 0; i < argc; i++)
    {
      len += std::strlen(argv[i]);
    }

  cmdline = static_cast<char *>(xmalloc(len + static_cast<size_t>(argc)));
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
  std::fprintf(stderr, "%s\n", prog_header.data());
  std::fprintf(stderr, "https://github.com/torognes/vsearch\n");
  std::fprintf(stderr, "\n");
}


/* When building as a library (VSEARCH_NO_MAIN), exclude main() to avoid
   symbol conflicts with the embedding application's entry point.
   All global variable definitions and helper functions above remain
   available — only the CLI entry point is excluded. */
#ifndef VSEARCH_NO_MAIN
/* Installed via std::set_new_handler so that a failed C++ allocation
   (std::vector, std::string, ...) reports the same friendly message as the
   xmalloc path instead of throwing std::bad_alloc, which under -fno-exceptions
   would call std::terminate()/abort(). The handler runs while memory is
   exhausted, so it must not allocate: it writes a fixed string with write(2)
   and leaves via _exit() without flushing stdio. This catches a refused
   allocation (operator new gets nullptr); it cannot catch a kernel OOM kill
   of an overcommitted allocation. */
static auto vsearch_new_handler() -> void
{
  static char const message[] = "\n\nFatal error: Unable to allocate enough memory.\n";
  ssize_t const written = write(STDERR_FILENO, message, sizeof(message) - 1);
  (void) written;
  _exit(EXIT_FAILURE);
}

auto main(int argc, char** argv) -> int
{
  std::set_new_handler(vsearch_new_handler);

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
      std::fprintf(fp_log, "%s\n", prog_header.data());
      std::fprintf(fp_log, "%s\n", cmdline);

      std::array<char, 26> time_string {{}};
      time_start = std::time(nullptr);
      struct tm const * tm_start = localtime(& time_start);
      std::strftime(time_string.data(), time_string.size(), "%Y-%m-%dT%H:%M:%S", tm_start);
      std::fprintf(fp_log, "Started  %s\n", time_string.data());
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
  else if (parameters.opt_fastx_syncpairs != nullptr)
    {
      fastx_syncpairs(parameters);
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
      cmd_fastq_join(parameters);
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
      time_finish = std::time(nullptr);
      struct tm const * tm_finish = localtime(& time_finish);
      std::array<char, 26> time_string {{}};
      std::strftime(time_string.data(), time_string.size(), "%Y-%m-%dT%H:%M:%S", tm_finish);
      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "Finished %s", time_string.data());

      double const time_diff = std::difftime(time_finish, time_start);
      std::fprintf(fp_log, "\n");
      std::fprintf(fp_log, "Elapsed time %02.0lf:%02.0lf\n",
              std::floor(time_diff / 60.0),
              std::floor(time_diff - (60.0 * std::floor(time_diff / 60.0))));
      double const maxmem = static_cast<double>(arch_get_memused()) / 1048576.0;
      if (maxmem < 1024.0)
        {
          std::fprintf(fp_log, "Max memory %.1lfMB\n", maxmem);
        }
      else
        {
          std::fprintf(fp_log, "Max memory %.1lfGB\n", maxmem / 1024.0);
        }
      std::fclose(fp_log);
    }

  if (opt_ee_cutoffs_values != nullptr)
    {
      xfree(opt_ee_cutoffs_values);
    }
  opt_ee_cutoffs_values = nullptr;

  xfree(cmdline);
  dynlibs_close();
}
#endif /* VSEARCH_NO_MAIN */
