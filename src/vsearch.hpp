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

#pragma once

/* _GNU_SOURCE and the __STDC_*_MACROS feature-test macros are defined in
   config.h (see configure.ac), which is included first below so they precede
   the standard library headers. */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdint>
#include <cstdio>  // std::size_t
#include <string>
#include <vector>
#include <limits>

#include "core/mask.hpp"  // Masking
#include "utils/quality_encoding.hpp"  // sanger_ascii_offset (opt_fastq_ascii default)

// C++20 refactoring: constexpr
std::string const default_quality_padding = "IIIIIIII";  // Q40 with an offset of 33
std::string const default_sequence_padding = "NNNNNNNN";

class DynamicLibraries;  // owned by main(), referenced through parameters.dyn_libs

struct Parameters {
private:
  /* Canonical option defaults referenced only by this struct's own member
     initializers (formerly namespace-scope constants in this header). Private
     so they do not widen the public Parameters surface; the initializers,
     being members, may still name them. */
  static constexpr int64_t default_fasta_width = 80;
  static constexpr int64_t default_fastq_tail = 4;
  static constexpr int64_t default_maxseqlength = 50000;
  static constexpr int64_t default_max_quality = 41;

public:
  std::string prog_header;
  std::string command_line;

  /* CPU features detected at startup (see utils/cpu_features) */
  int64_t altivec_present {0};  // unused
  int64_t neon_present {0};     // unused
  int64_t mmx_present {0};      // unused
  int64_t sse_present {0};      // unused
  int64_t sse2_present {0};
  int64_t sse3_present {0};     // unused
  int64_t ssse3_present {0};
  int64_t sse41_present {0};    // unused
  int64_t sse42_present {0};    // unused
  int64_t popcnt_present {0};   // unused
  int64_t avx_present {0};      // unused
  int64_t avx2_present {0};     // unused
  char * opt_allpairs_global = nullptr;
  char * opt_chimeras_denovo = nullptr;
  char * opt_cluster_fast = nullptr;
  char * opt_cluster_size = nullptr;
  char * opt_cluster_smallmem = nullptr;
  char * opt_cluster_unoise = nullptr;
  char * opt_cut = nullptr;
  std::string opt_cut_pattern;
  char * opt_db = nullptr;
  char * opt_dbmatched = nullptr;
  char * opt_dbnotmatched = nullptr;
  char * opt_derep_fulllength = nullptr;
  char * opt_derep_id = nullptr;
  char * opt_derep_prefix = nullptr;
  char * opt_derep_smallmem = nullptr;
  char * opt_fasta2fastq = nullptr;
  char * opt_fastaout = nullptr;
  char * opt_fastaout_rev = nullptr;
  char * opt_fastaout_discarded = nullptr;
  char * opt_fastaout_discarded_rev = nullptr;
  char * opt_fastaout_orphans = nullptr;
  char * opt_fastaout_orphans_rev = nullptr;
  char * opt_fastq_chars = nullptr;
  char * opt_fastq_convert = nullptr;
  char * opt_fastq_eestats2 = nullptr;
  char * opt_fastq_eestats = nullptr;
  char * opt_fastq_filter = nullptr;
  char * opt_fastq_join = nullptr;
  char * opt_fastq_mergepairs = nullptr;
  char * opt_fastq_stats = nullptr;
  char * opt_fastqout = nullptr;
  char * opt_fastqout_rev = nullptr;
  char * opt_fastqout_discarded = nullptr;
  char * opt_fastqout_discarded_rev = nullptr;
  char * opt_fastqout_orphans = nullptr;
  char * opt_fastqout_orphans_rev = nullptr;
  char * opt_fastx_filter = nullptr;
  char * opt_fastx_getseq = nullptr;
  char * opt_fastx_getseqs = nullptr;
  char * opt_fastx_getsubseq = nullptr;
  char * opt_fastx_mask = nullptr;
  char * opt_fastx_revcomp = nullptr;
  char * opt_fastx_subsample = nullptr;
  char * opt_fastx_syncpairs = nullptr;
  char * opt_fastx_uniques = nullptr;
  std::string opt_join_padgap = default_sequence_padding;
  std::string opt_join_padgapq = default_quality_padding;
  char * opt_label_suffix = nullptr;
  char * opt_log = nullptr;
  char * opt_makeudb_usearch = nullptr;
  char * opt_maskfasta = nullptr;
  char * opt_orient = nullptr;
  char * opt_output = nullptr;
  char * opt_relabel = nullptr;
  char * opt_read_separators = nullptr;
  char * opt_rereplicate = nullptr;
  char * opt_reverse = nullptr;
  char * opt_sample = nullptr;
  char * opt_search_exact = nullptr;
  char * opt_sff_convert = nullptr;
  char * opt_shuffle = nullptr;
  char * opt_sintax = nullptr;
  char * opt_sortbylength = nullptr;
  char * opt_sortbysize = nullptr;
  char * opt_tabbedout = nullptr;
  char * opt_uc = nullptr;
  char * opt_uchime2_denovo = nullptr;
  char * opt_uchime3_denovo = nullptr;
  char * opt_uchime_denovo = nullptr;
  char * opt_uchime_ref = nullptr;
  char * opt_udb2fasta = nullptr;
  char * opt_udbinfo = nullptr;
  char * opt_udbstats = nullptr;
  char * opt_usearch_global = nullptr;
  char * progname = nullptr;  // refactoring: unused?
  std::FILE * fp_log = nullptr;
  DynamicLibraries const * dyn_libs = nullptr;
  double opt_fastq_truncee_rate = std::numeric_limits<double>::max();
  double opt_max_unmasked_pct = 100.0;
  double opt_min_unmasked_pct = 0;
  double opt_sample_pct = 0;
  int64_t opt_fasta_width = default_fasta_width;
  int64_t opt_fastq_ascii = sanger_ascii_offset;
  int64_t opt_fastq_asciiout = sanger_ascii_offset;
  int64_t opt_fastq_qmax = default_max_quality;
  int64_t opt_fastq_qmaxout = default_max_quality;
  int64_t opt_fastq_qmin = 0;
  int64_t opt_fastq_qminout = 0;
  int64_t opt_fastq_minqual = 0;
  int64_t opt_fastq_tail = default_fastq_tail;
  int64_t opt_maxseqlength = default_maxseqlength;
  int64_t opt_maxsize = std::numeric_limits<int64_t>::max();
  int64_t opt_maxuniquesize = std::numeric_limits<int64_t>::max();
  int64_t opt_minseqlength = -1;
  int64_t opt_minsize = 0;
  int64_t opt_minuniquesize = 1;
  Masking opt_qmask = Masking::dust;
  int64_t opt_randseed = 0;
  int64_t opt_sample_size = 0;
  int64_t opt_threads = 0;
  int64_t opt_topn = std::numeric_limits<int64_t>::max();
  bool opt_bzip2_decompress = false;
  bool opt_clusterout_id = false;
  bool opt_clusterout_sort = false;
  bool opt_eeout = false;
  bool opt_fasta_score = false;
  bool opt_fastq_allowmergestagger = false;
  bool opt_fastq_eeout = false;
  bool opt_fastq_nostagger = true;
  bool opt_fastq_qout_max = false;
  bool opt_gzip_decompress = false;
  bool opt_hardmask = false;
  bool opt_label_substr_match = false;
  bool opt_lengthout = false;
  bool opt_no_progress = true;   // library default (quiet); the CLI overrides in args_init
  bool opt_help = false;
  bool opt_join_padgapq_set_by_user = false;
  bool opt_notrunclabels = false;
  bool opt_quiet = true;   // library default (quiet); the CLI overrides in args_init
  bool opt_relabel_keep = false;
  bool opt_relabel_md5 = false;
  bool opt_relabel_self = false;
  bool opt_relabel_sha1 = false;
  bool opt_samheader = false;
  bool opt_sff_clip = false;
  bool opt_sintax_random = false;
  bool opt_sizein = false;
  bool opt_sizeorder = false;
  bool opt_sizeout = false;
  bool opt_stderr_is_tty = false;
  bool opt_strand = false;
  bool opt_uc_allhits = false;
  bool opt_version = false;
  bool opt_xee = false;
  bool opt_xlength = false;
  bool opt_xsize = false;

  /* Options migrated from the opt_* globals (E1/F1). Types match the globals
     and the default values are the canonical library defaults, so a
     default-constructed Parameters is a fully-defaulted configuration that
     vsearch_session_begin() then derives the globals from. */
  char *    opt_alnout                       = nullptr;
  char *    opt_biomout                      = nullptr;
  char *    opt_blast6out                    = nullptr;
  char *    opt_borderline                   = nullptr;
  char *    opt_centroids                    = nullptr;
  char *    opt_chimeras                     = nullptr;
  char *    opt_clusters                     = nullptr;
  char *    opt_consout                      = nullptr;
  char *    opt_eetabbedout                  = nullptr;
  char *    opt_fastaout_notmerged_fwd       = nullptr;
  char *    opt_fastaout_notmerged_rev       = nullptr;
  char *    opt_fastapairs                   = nullptr;
  char *    opt_fastqout_notmerged_fwd       = nullptr;
  char *    opt_fastqout_notmerged_rev       = nullptr;
  char *    opt_label_field                  = nullptr;
  char *    opt_label                        = nullptr;
  char *    opt_labels                       = nullptr;
  char *    opt_label_word                   = nullptr;
  char *    opt_label_words                  = nullptr;
  char *    opt_lcaout                       = nullptr;
  char *    opt_matched                      = nullptr;
  char *    opt_mothur_shared_out            = nullptr;
  char *    opt_msaout                       = nullptr;
  char *    opt_nonchimeras                  = nullptr;
  char *    opt_notmatchedfq                 = nullptr;
  char *    opt_notmatched                   = nullptr;
  char *    opt_otutabout                    = nullptr;
  char *    opt_pattern                      = nullptr;
  char *    opt_profile                      = nullptr;
  char *    opt_qsegout                      = nullptr;
  char *    opt_samout                       = nullptr;
  char *    opt_tsegout                      = nullptr;
  char *    opt_uchimealns                   = nullptr;
  char *    opt_uchimeout                    = nullptr;
  char *    opt_userout                      = nullptr;

  double    opt_abskew                       = 0.0;
  double    opt_chimeras_diff_pct            = 0.0;
  double    opt_dn                           = 1.4;
  double    opt_fastq_maxdiffpct             = 100.0;
  double    opt_fastq_maxee                  = std::numeric_limits<double>::max();
  double    opt_fastq_maxee_rate             = std::numeric_limits<double>::max();
  double    opt_fastq_truncee                = std::numeric_limits<double>::max();
  double    opt_id                           = -1.0;
  double    opt_lca_cutoff                   = 1.0;
  double    opt_maxid                        = 1.0;
  double    opt_maxqt                        = std::numeric_limits<double>::max();
  double    opt_maxsizeratio                 = std::numeric_limits<double>::max();
  double    opt_maxsl                        = std::numeric_limits<double>::max();
  double    opt_mid                          = 0.0;
  double    opt_mindiv                       = 0.8;
  double    opt_minh                         = 0.28;
  double    opt_minqt                        = 0.0;
  double    opt_minsizeratio                 = 0.0;
  double    opt_minsl                        = 0.0;
  double    opt_query_cov                    = 0.0;
  double    opt_sintax_cutoff                = 0.0;
  double    opt_target_cov                   = 0.0;
  double    opt_unoise_alpha                 = 2.0;
  double    opt_weak_id                      = 10.0;
  double    opt_xn                           = 8.0;

  int       opt_acceptall                    = 0;
  int       opt_alignwidth                   = 80;
  int       opt_chimeras_length_min          = 10;
  int       opt_chimeras_parents_max         = 3;
  int       opt_chimeras_parts               = 0;
  int       opt_cons_truncate                = 0;
  int       opt_gap_extension_query_interior = 2;
  int       opt_gap_extension_query_left     = 1;
  int       opt_gap_extension_query_right    = 1;
  int       opt_gap_extension_target_interior = 2;
  int       opt_gap_extension_target_left    = 1;
  int       opt_gap_extension_target_right   = 1;
  int       opt_gap_open_query_interior      = 20;
  int       opt_gap_open_query_left          = 2;
  int       opt_gap_open_query_right         = 2;
  int       opt_gap_open_target_interior     = 20;
  int       opt_gap_open_target_left         = 2;
  int       opt_gap_open_target_right        = 2;
  /* '*' (infinite) gap-penalty sentinels, recorded at parse time. The numeric
     penalties above are kept for scoring; these enforce the hard forbid at
     accept time: an infinite gap-open forbids that gap class outright, an
     infinite gap-extension forbids gaps of that class longer than one.
     opt_gap_penalty_has_infinite is the OR of the twelve, computed in
     vsearch_apply_defaults_fixups so the hot accept path can skip the scan. */
  bool      opt_gap_open_query_left_infinite          = false;
  bool      opt_gap_open_query_interior_infinite      = false;
  bool      opt_gap_open_query_right_infinite         = false;
  bool      opt_gap_open_target_left_infinite         = false;
  bool      opt_gap_open_target_interior_infinite     = false;
  bool      opt_gap_open_target_right_infinite        = false;
  bool      opt_gap_extension_query_left_infinite     = false;
  bool      opt_gap_extension_query_interior_infinite = false;
  bool      opt_gap_extension_query_right_infinite    = false;
  bool      opt_gap_extension_target_left_infinite    = false;
  bool      opt_gap_extension_target_interior_infinite = false;
  bool      opt_gap_extension_target_right_infinite   = false;
  bool      opt_gap_penalty_has_infinite              = false;
  int       opt_length_cutoffs_increment     = 50;
  int       opt_length_cutoffs_longest       = std::numeric_limits<int>::max();
  int       opt_length_cutoffs_shortest      = 50;
  int       opt_mindiffs                     = 3;
  int       opt_slots                        = 0;
  int       opt_uchimeout5                   = 0;
  int       opt_usersort                     = 0;

  Masking   opt_dbmask                       = Masking::dust;
  int64_t   opt_fastq_maxdiffs               = 10;
  int64_t   opt_fastq_maxlen                 = std::numeric_limits<int64_t>::max();
  int64_t   opt_fastq_maxmergelen            = 1000000;
  int64_t   opt_fastq_maxns                  = std::numeric_limits<int64_t>::max();
  int64_t   opt_fastq_minlen                 = 1;
  int64_t   opt_fastq_minmergelen            = 0;
  int64_t   opt_fastq_minovlen               = 10;
  int64_t   opt_fastq_stripleft              = 0;
  int64_t   opt_fastq_stripright             = 0;
  int64_t   opt_fastq_trunclen               = -1;
  int64_t   opt_fastq_trunclen_keep          = -1;
  int64_t   opt_fastq_truncqual              = std::numeric_limits<long>::min();
  int64_t   opt_fulldp                       = 0;
  int64_t   opt_iddef                        = 2;
  int64_t   opt_idprefix                     = 0;
  int64_t   opt_idsuffix                     = 0;
  int64_t   opt_leftjust                     = 0;
  int64_t   opt_match                        = 2;
  int64_t   opt_maxaccepts                   = 1;
  int64_t   opt_maxdiffs                     = std::numeric_limits<int>::max();
  int64_t   opt_maxgaps                      = std::numeric_limits<int>::max();
  int64_t   opt_maxhits                      = 0;
  int64_t   opt_maxqsize                     = std::numeric_limits<int64_t>::max();
  int64_t   opt_maxrejects                   = -1;
  int64_t   opt_maxsubs                      = std::numeric_limits<int>::max();
  int64_t   opt_mincols                      = 0;
  int64_t   opt_mintsize                     = 0;
  int64_t   opt_minwordmatches               = -1;
  int64_t   opt_mismatch                     = -4;
  int64_t   opt_output_no_hits               = 0;
  int64_t   opt_rightjust                    = 0;
  int64_t   opt_rowlen                       = 64;
  int64_t   opt_self                         = 0;
  int64_t   opt_selfid                       = 0;
  int64_t   opt_subseq_end                   = std::numeric_limits<int64_t>::max();
  int64_t   opt_subseq_start                 = 1;
  int64_t   opt_top_hits_only                = 0;
  int64_t   opt_wordlength                   = 0;

  bool      opt_centroid_sizeout             = false;
  bool      opt_n_mismatch                   = false;

  std::vector<double> opt_ee_cutoffs = {0.5, 1.0, 2.0};  // was opt_ee_cutoffs_values/_count
  std::vector<int> opt_userfields;  // was userfields_requested/_count (globals)

  /* Internal state (not an option): guards the once-only gap-open penalty
     adjustment in vsearch_apply_defaults_fixups() so a repeated call on the
     same struct is idempotent rather than double-subtracting. */
  bool gap_penalties_adjusted = false;
};

/* The shared parameter-resolution / session-lifecycle declarations that used
   to sit here now live in parameters.hpp (paired with parameters.cpp); it is
   included so the many translation units that reach these functions through
   vsearch.hpp keep compiling unchanged. */
#include "parameters.hpp"
