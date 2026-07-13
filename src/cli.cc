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

#include "cli.h"
#include "vsearch.h"
#include "vsearch_api.h"
#include "os/system.hpp"  // system_get_cores
#include "core/chimera.hpp"  // maxparents
#include "core/mask.hpp"  // Masking
#include "utils/userfields.hpp"  // parse_userfields_arg
#include "utils/compare_strings_nocase.hpp"  // are_same_string
#include "utils/fatal.hpp"  // fatal
#include <algorithm>  // std::count
#include <array>
#include <cerrno>  // errno, ERANGE
#include <cinttypes>  // macro SCNd64
#include <cmath>  // std::isfinite
#include <cstdint>  // int64_t
#include <cstdio>  // std::sscanf, std::fprintf, fprintf, stderr, stdout
#include <cstdlib>  // exit, EXIT_FAILURE
#include <cstring>  // std::strlen
#include <getopt.h>  // getopt_long_only, optarg, optind, opterr, struct
                     // option (no_argument, required_argument)
#include <limits>  // std::numeric_limits
#include <unistd.h>  // isatty, fileno


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto args_get_ee_cutoffs(char * arg, struct Parameters & parameters) -> void
  {
    /* get comma-separated list of floating point numbers */
    /* save in parameters.opt_ee_cutoffs */

    parameters.opt_ee_cutoffs.clear();

    char const * cursor = arg;
    while (true)
      {
        double val = 0;
        int skip = 0;

        if ((std::sscanf(cursor, "%lf%n", &val, &skip) != 1) or (val <= 0.0))
          {
            fatal("Invalid arguments to ee_cutoffs");
          }

        parameters.opt_ee_cutoffs.push_back(val);

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


  auto args_get_length_cutoffs(char const * arg, struct Parameters & parameters) -> void
  {
    /* get comma-separated list of 3 integers: */
    /* smallest, largest and increment. */
    /* second value may be * indicating no limit */
    /* save in length_cutoffs_{smallest,largest,increment} */

    // refactoring: std::stoi(), faster than sscanf()
    static constexpr auto n_of_expected_assignments= 3;
    int skip = 0;  // receives the number of characters read so far ('%n')
    if (std::sscanf(arg, "%d,%d,%d%n", &parameters.opt_length_cutoffs_shortest, &parameters.opt_length_cutoffs_longest, &parameters.opt_length_cutoffs_increment, & skip) == n_of_expected_assignments)
      {
        if (static_cast<size_t>(skip) < std::strlen(arg))
          {
            fatal("Invalid arguments to length_cutoffs");
          }
      }
    else if (std::sscanf(arg, "%d,*,%d%n", &parameters.opt_length_cutoffs_shortest, &parameters.opt_length_cutoffs_increment, &skip) == 2)
      {
        if (static_cast<size_t>(skip) < std::strlen(arg))
          {
            fatal("Invalid arguments to length_cutoffs");
          }
        parameters.opt_length_cutoffs_longest = std::numeric_limits<int>::max();
      }
    else
      {
        fatal("Invalid arguments to length_cutoffs");
      }

    if ((parameters.opt_length_cutoffs_shortest < 1) or
        (parameters.opt_length_cutoffs_shortest > parameters.opt_length_cutoffs_longest) or
        (parameters.opt_length_cutoffs_increment < 1))
      {
        fatal("Invalid arguments to length_cutoffs");
      }
  }


  /* The SIMD aligner (search16, align_simd.cc) represents alignment scores in
     signed 16-bit cells (CELL == signed short). --match and --mismatch are
     stored verbatim in the score matrix, so they must fit a CELL. Gap penalties
     feed the block-seed value -(open + CDEPTH * ext) (CDEPTH == 4 in
     align_simd.cc) through a non-saturating narrowing cast, so bounding each
     penalty by SHRT_MAX / (1 + CDEPTH) keeps open + CDEPTH * ext within CELL
     range for any open/ext combination. Values outside these ranges would wrap
     silently and corrupt alignments, so they are rejected at parse time. The
     '*' (infinite) gap penalty is a documented sentinel and is left untouched
     here; its handling in the SIMD path is addressed separately. */
  static constexpr auto max_alignment_score = int64_t{32767};  // SHRT_MAX
  static constexpr auto max_gap_penalty = int64_t{32767 / (1 + 4)};  // SHRT_MAX / (1 + CDEPTH)


  auto args_get_gap_penalty_string(char * arg, bool const is_open, struct Parameters & parameters) -> void
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
    static constexpr auto max_penality = std::numeric_limits<int>::max();

    char * cursor = arg;

    while (*cursor != '\0')
      {
        int skip = 0;
        int pen = 0;
        bool is_infinite = false;

        if (std::sscanf(cursor, "%d%n", &pen, &skip) == 1)
          {
            cursor += skip;
            if ((pen < 0) or (pen > max_gap_penalty))
              {
                std::string const message =
                  "A finite gap penalty must be in the range 0 to "
                  + std::to_string(max_gap_penalty)
                  + "; use '*' to declare an infinite penalty";
                fatal(message.c_str());
              }
          }
        else if (*cursor == '*')
          {
            pen = max_penality;
            is_infinite = true;
            ++cursor;
          }
        else
          {
            fatal("Invalid gap penalty argument (%s)", cursor);
          }

        char const * q = cursor;

        int set_E = 0;
        int set_I = 0;
        int set_L = 0;
        int set_R = 0;
        int set_Q = 0;
        int set_T = 0;

        while ((*cursor != '\0') and (*cursor != '/'))
          {
            switch (*cursor)
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

        /* Each assignment also records whether the penalty is the '*'
           sentinel, so a later finite token in the same string (e.g. '*I/20I')
           clears it -- last token wins, matching the numeric penalty. */
        if (is_open)
          {
            if (set_Q != 0)
              {
                if (set_L != 0)
                  {
                    parameters.opt_gap_open_query_left = pen;
                    parameters.opt_gap_open_query_left_infinite = is_infinite;
                  }
                if (set_I != 0)
                  {
                    parameters.opt_gap_open_query_interior = pen;
                    parameters.opt_gap_open_query_interior_infinite = is_infinite;
                  }
                if (set_R != 0)
                  {
                    parameters.opt_gap_open_query_right = pen;
                    parameters.opt_gap_open_query_right_infinite = is_infinite;
                  }
              }
            if (set_T != 0)
              {
                if (set_L != 0)
                  {
                    parameters.opt_gap_open_target_left = pen;
                    parameters.opt_gap_open_target_left_infinite = is_infinite;
                  }
                if (set_I != 0)
                  {
                    parameters.opt_gap_open_target_interior = pen;
                    parameters.opt_gap_open_target_interior_infinite = is_infinite;
                  }
                if (set_R != 0)
                  {
                    parameters.opt_gap_open_target_right = pen;
                    parameters.opt_gap_open_target_right_infinite = is_infinite;
                  }
              }
          }
        else
          {
            if (set_Q != 0)
              {
                if (set_L != 0)
                  {
                    parameters.opt_gap_extension_query_left = pen;
                    parameters.opt_gap_extension_query_left_infinite = is_infinite;
                  }
                if (set_I != 0)
                  {
                    parameters.opt_gap_extension_query_interior = pen;
                    parameters.opt_gap_extension_query_interior_infinite = is_infinite;
                  }
                if (set_R != 0)
                  {
                    parameters.opt_gap_extension_query_right = pen;
                    parameters.opt_gap_extension_query_right_infinite = is_infinite;
                  }
              }
            if (set_T != 0)
              {
                if (set_L != 0)
                  {
                    parameters.opt_gap_extension_target_left = pen;
                    parameters.opt_gap_extension_target_left_infinite = is_infinite;
                  }
                if (set_I != 0)
                  {
                    parameters.opt_gap_extension_target_interior = pen;
                    parameters.opt_gap_extension_target_interior_infinite = is_infinite;
                  }
                if (set_R != 0)
                  {
                    parameters.opt_gap_extension_target_right = pen;
                    parameters.opt_gap_extension_target_right_infinite = is_infinite;
                  }
              }
          }
      }
  }


  auto args_getlong(char const * arg) -> int64_t
  {
    int len = 0;
    int64_t temp = 0;
    errno = 0;
    auto const ret = std::sscanf(arg, "%" SCNd64 "%n", &temp, &len);
    if ((ret == 0) or ((static_cast<unsigned int>(len)) < std::strlen(arg)) or (errno == ERANGE))
      {
        fatal("Illegal option argument");
      }
    return temp;
  }


  auto args_getdouble(char const * arg) -> double
  {
    int len = 0;
    double temp = 0;
    errno = 0;
    auto const ret = std::sscanf(arg, "%lf%n", &temp, &len);

    if ((ret == 0) or ((static_cast<unsigned int>(len)) < std::strlen(arg)) or (errno == ERANGE) or (not std::isfinite(temp)))
      {
        fatal("Illegal option argument");
      }
    return temp;
  }


  static constexpr auto number_of_commands = std::size_t{51};
  static constexpr auto number_of_options = std::size_t{254};
  static constexpr auto max_number_of_options_per_command = std::size_t{100};

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
      option_centroid_sizeout,
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
      option_fastaout_orphans,
      option_fastaout_orphans_rev,
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
      option_fastqout_orphans,
      option_fastqout_orphans_rev,
      option_fastqout_rev,
      option_fastx_filter,
      option_fastx_getseq,
      option_fastx_getseqs,
      option_fastx_getsubseq,
      option_fastx_mask,
      option_fastx_revcomp,
      option_fastx_subsample,
      option_fastx_syncpairs,
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
      option_read_separators,
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
      option_xsize,
      option_count  // number of options; keep last (sizes option_specs below)
    };

  /*
    Single source of truth for the command-line option metadata. The getopt
    long_options array is derived from this table (see below), so a new
    option's name and argument requirement are declared in exactly one place.
    The entries mirror the option_* enum one-for-one and MUST stay in the same
    (alphabetical) order: getopt reports a match by its index into
    long_options, which is used directly as the enum value in the switch.
    needs_arg selects required_argument (true) versus no_argument (false).
    vsearch has no short options, so getopt_long_only is given an empty short
    string and every entry's flag/val are nullptr/0.
  */
  struct OptionSpec
  {
    char const * name;
    bool needs_arg;
  };

  static constexpr std::array<OptionSpec, option_count> option_specs =
    {{
      {"abskew",                     true },
      {"acceptall",                  false },
      {"alignwidth",                 true },
      {"allpairs_global",            true },
      {"alnout",                     true },
      {"band",                       true },
      {"biomout",                    true },
      {"blast6out",                  true },
      {"borderline",                 true },
      {"bzip2_decompress",           false },
      {"centroid_sizeout",           false },
      {"centroids",                  true },
      {"chimeras",                   true },
      {"chimeras_denovo",            true },
      {"chimeras_diff_pct",          true },
      {"chimeras_length_min",        true },
      {"chimeras_parents_max",       true },
      {"chimeras_parts",             true },
      {"cluster_fast",               true },
      {"cluster_size",               true },
      {"cluster_smallmem",           true },
      {"cluster_unoise",             true },
      {"clusterout_id",              false },
      {"clusterout_sort",            false },
      {"clusters",                   true },
      {"cons_truncate",              false },
      {"consout",                    true },
      {"cut",                        true },
      {"cut_pattern",                true },
      {"db",                         true },
      {"dbmask",                     true },
      {"dbmatched",                  true },
      {"dbnotmatched",               true },
      {"derep_fulllength",           true },
      {"derep_id",                   true },
      {"derep_prefix",               true },
      {"derep_smallmem",             true },
      {"dn",                         true },
      {"ee_cutoffs",                 true },
      {"eeout",                      false },
      {"eetabbedout",                true },
      {"fasta2fastq",                true },
      {"fasta_score",                false },
      {"fasta_width",                true },
      {"fastaout",                   true },
      {"fastaout_discarded",         true },
      {"fastaout_discarded_rev",     true },
      {"fastaout_notmerged_fwd",     true },
      {"fastaout_notmerged_rev",     true },
      {"fastaout_orphans",           true },
      {"fastaout_orphans_rev",       true },
      {"fastaout_rev",               true },
      {"fastapairs",                 true },
      {"fastq_allowmergestagger",    false },
      {"fastq_ascii",                true },
      {"fastq_asciiout",             true },
      {"fastq_chars",                true },
      {"fastq_convert",              true },
      {"fastq_eeout",                false },
      {"fastq_eestats",              true },
      {"fastq_eestats2",             true },
      {"fastq_filter",               true },
      {"fastq_join",                 true },
      {"fastq_maxdiffpct",           true },
      {"fastq_maxdiffs",             true },
      {"fastq_maxee",                true },
      {"fastq_maxee_rate",           true },
      {"fastq_maxlen",               true },
      {"fastq_maxmergelen",          true },
      {"fastq_maxns",                true },
      {"fastq_mergepairs",           true },
      {"fastq_minlen",               true },
      {"fastq_minmergelen",          true },
      {"fastq_minovlen",             true },
      {"fastq_minqual",              true },
      {"fastq_nostagger",            false },
      {"fastq_qmax",                 true },
      {"fastq_qmaxout",              true },
      {"fastq_qmin",                 true },
      {"fastq_qminout",              true },
      {"fastq_qout_max",             false },
      {"fastq_stats",                true },
      {"fastq_stripleft",            true },
      {"fastq_stripright",           true },
      {"fastq_tail",                 true },
      {"fastq_truncee",              true },
      {"fastq_truncee_rate",         true },
      {"fastq_trunclen",             true },
      {"fastq_trunclen_keep",        true },
      {"fastq_truncqual",            true },
      {"fastqout",                   true },
      {"fastqout_discarded",         true },
      {"fastqout_discarded_rev",     true },
      {"fastqout_notmerged_fwd",     true },
      {"fastqout_notmerged_rev",     true },
      {"fastqout_orphans",           true },
      {"fastqout_orphans_rev",       true },
      {"fastqout_rev",               true },
      {"fastx_filter",               true },
      {"fastx_getseq",               true },
      {"fastx_getseqs",              true },
      {"fastx_getsubseq",            true },
      {"fastx_mask",                 true },
      {"fastx_revcomp",              true },
      {"fastx_subsample",            true },
      {"fastx_syncpairs",            true },
      {"fastx_uniques",              true },
      {"fulldp",                     false },
      {"gapext",                     true },
      {"gapopen",                    true },
      {"gzip_decompress",            false },
      {"h",                          false },
      {"hardmask",                   false },
      {"help",                       false },
      {"hspw",                       true },
      {"id",                         true },
      {"iddef",                      true },
      {"idprefix",                   true },
      {"idsuffix",                   true },
      {"join_padgap",                true },
      {"join_padgapq",               true },
      {"label",                      true },
      {"label_field",                true },
      {"label_substr_match",         false },
      {"label_suffix",               true },
      {"label_word",                 true },
      {"label_words",                true },
      {"labels",                     true },
      {"lca_cutoff",                 true },
      {"lcaout",                     true },
      {"leftjust",                   false },
      {"length_cutoffs",             true },
      {"lengthout",                  false },
      {"log",                        true },
      {"makeudb_usearch",            true },
      {"maskfasta",                  true },
      {"match",                      true },
      {"matched",                    true },
      {"max_unmasked_pct",           true },
      {"maxaccepts",                 true },
      {"maxdiffs",                   true },
      {"maxgaps",                    true },
      {"maxhits",                    true },
      {"maxid",                      true },
      {"maxqsize",                   true },
      {"maxqt",                      true },
      {"maxrejects",                 true },
      {"maxseqlength",               true },
      {"maxsize",                    true },
      {"maxsizeratio",               true },
      {"maxsl",                      true },
      {"maxsubs",                    true },
      {"maxuniquesize",              true },
      {"mid",                        true },
      {"min_unmasked_pct",           true },
      {"mincols",                    true },
      {"mindiffs",                   true },
      {"mindiv",                     true },
      {"minh",                       true },
      {"minhsp",                     true },
      {"minqt",                      true },
      {"minseqlength",               true },
      {"minsize",                    true },
      {"minsizeratio",               true },
      {"minsl",                      true },
      {"mintsize",                   true },
      {"minuniquesize",              true },
      {"minwordmatches",             true },
      {"mismatch",                   true },
      {"mothur_shared_out",          true },
      {"msaout",                     true },
      {"n_mismatch",                 false },
      {"no_progress",                false },
      {"nonchimeras",                true },
      {"notmatched",                 true },
      {"notmatchedfq",               true },
      {"notrunclabels",              false },
      {"orient",                     true },
      {"otutabout",                  true },
      {"output",                     true },
      {"output_no_hits",             false },
      {"pattern",                    true },
      {"profile",                    true },
      {"qmask",                      true },
      {"qsegout",                    true },
      {"query_cov",                  true },
      {"quiet",                      false },
      {"randseed",                   true },
      {"read_separators",            true },
      {"relabel",                    true },
      {"relabel_keep",               false },
      {"relabel_md5",                false },
      {"relabel_self",               false },
      {"relabel_sha1",               false },
      {"rereplicate",                true },
      {"reverse",                    true },
      {"rightjust",                  false },
      {"rowlen",                     true },
      {"samheader",                  false },
      {"samout",                     true },
      {"sample",                     true },
      {"sample_pct",                 true },
      {"sample_size",                true },
      {"search_exact",               true },
      {"self",                       false },
      {"selfid",                     false },
      {"sff_clip",                   false },
      {"sff_convert",                true },
      {"shuffle",                    true },
      {"sintax",                     true },
      {"sintax_cutoff",              true },
      {"sintax_random",              false },
      {"sizein",                     false },
      {"sizeorder",                  false },
      {"sizeout",                    false },
      {"slots",                      true },
      {"sortbylength",               true },
      {"sortbysize",                 true },
      {"strand",                     true },
      {"subseq_end",                 true },
      {"subseq_start",               true },
      {"tabbedout",                  true },
      {"target_cov",                 true },
      {"threads",                    true },
      {"top_hits_only",              false },
      {"topn",                       true },
      {"tsegout",                    true },
      {"uc",                         true },
      {"uc_allhits",                 false },
      {"uchime2_denovo",             true },
      {"uchime3_denovo",             true },
      {"uchime_denovo",              true },
      {"uchime_ref",                 true },
      {"uchimealns",                 true },
      {"uchimeout",                  true },
      {"uchimeout5",                 false },
      {"udb2fasta",                  true },
      {"udbinfo",                    true },
      {"udbstats",                   true },
      {"unoise_alpha",               true },
      {"usearch_global",             true },
      {"userfields",                 true },
      {"userout",                    true },
      {"usersort",                   false },
      {"v",                          false },
      {"version",                    false },
      {"weak_id",                    true },
      {"wordlength",                 true },
      {"xdrop_nw",                   true },
      {"xee",                        false },
      {"xlength",                    false },
      {"xn",                         true },
      {"xsize",                      false },
      }};

  static_assert(option_specs.size() + 1 == number_of_options,
                "number_of_options must be the option count plus the getopt sentinel");

  auto build_long_options() -> std::array<struct option, number_of_options>
  {
    /* Derive the getopt long_options array from the single-source table above.
       The array is value-initialised, so the trailing element is left as the
       {nullptr, 0, nullptr, 0} sentinel getopt_long_only expects. */
    std::array<struct option, number_of_options> long_options {};
    for (std::size_t idx = 0; idx < option_specs.size(); ++idx)
      {
        long_options[idx].name    = option_specs[idx].name;
        long_options[idx].has_arg = option_specs[idx].needs_arg ? required_argument : no_argument;
        long_options[idx].flag    = nullptr;
        long_options[idx].val     = 0;
      }
    return long_options;
  }

  /*
    Below is a list of all the options that are valid for each command.
    The first line is the command and the lines below are the valid options.
    There is one row per command, with the commands in alphabetical order by
    name; each row is terminated by -1. A command is recognised, and its name
    reported, through the first entry of its row (valid_options[k][0]) -- this
    is the single source for the command list.
  */

  static constexpr std::array<std::array<int, max_number_of_options_per_command>, number_of_commands> valid_options =
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
        option_lengthout,
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
        option_xlength,
        option_xn,
        option_xsize,
        -1 },

      { option_cluster_fast,
        option_alnout,
        option_band,
        option_biomout,
        option_blast6out,
        option_bzip2_decompress,
        option_centroid_sizeout,
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
        option_centroid_sizeout,
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
        option_centroid_sizeout,
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
        option_centroid_sizeout,
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

      { option_fastx_syncpairs,
        option_bzip2_decompress,
        option_fasta_width,
        option_fastaout,
        option_fastaout_orphans,
        option_fastaout_orphans_rev,
        option_fastaout_rev,
        option_fastq_ascii,
        option_fastq_qmax,
        option_fastq_qmin,
        option_fastqout,
        option_fastqout_orphans,
        option_fastqout_orphans_rev,
        option_fastqout_rev,
        option_gzip_decompress,
        option_log,
        option_no_progress,
        option_quiet,
        option_read_separators,
        option_reverse,
        option_threads,
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
        option_hardmask,
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

  /* Parse the command line with getopt: set each recognised option's
     Parameters field, record which options were seen in the returned
     vector (indexed by the option_* enum), and terminate on an ambiguous
     or illegal option or a stray non-option argument. */
  auto parse_options(int argc, char ** argv,
                     std::array<struct option, number_of_options> const & long_options,
                     struct Parameters & parameters) -> std::vector<bool>
  {
    std::vector<bool> options_selected(static_cast<std::size_t>(option_count));

    int options_index = 0;
    int val = 0;  // recognized long option: return 'val' if 'flag' is nullptr

    while ((val = getopt_long_only(argc, argv, "", long_options.data(),
                                 &options_index)) == 0)
      {
        if (options_index < option_count)
          {
            options_selected[static_cast<size_t>(options_index)] = true;
          }

        switch (options_index)
          {
          case option_help:
            parameters.opt_help = true;
            break;

          case option_version:
            parameters.opt_version = true;
            break;

          case option_alnout:
            parameters.opt_alnout = optarg;
            break;

          case option_usearch_global:
            parameters.opt_usearch_global = optarg;
            break;

          case option_db:
            parameters.opt_db = optarg;
            break;

          case option_id:
            parameters.opt_id = args_getdouble(optarg);
            break;

          case option_maxaccepts:
            parameters.opt_maxaccepts = args_getlong(optarg);
            break;

          case option_maxrejects:
            parameters.opt_maxrejects = args_getlong(optarg);
            break;

          case option_wordlength:
            parameters.opt_wordlength = args_getlong(optarg);
            /* 0 is reserved as the "not set" sentinel (resolved later
               to a command-specific default); reject an explicit 0 so
               the user-facing range matches the manpage (3..15). */
            if ((parameters.opt_wordlength < 3) or (parameters.opt_wordlength > 15))
              {
                fatal("Argument to --wordlength must be in the range 3 to 15");
              }
            break;

          case option_match:
            parameters.opt_match = args_getlong(optarg);
            break;

          case option_mismatch:
            parameters.opt_mismatch = args_getlong(optarg);
            break;

          case option_fulldp:
            parameters.opt_fulldp = 1;
            std::fprintf(stderr, "WARNING: Option --fulldp is ignored\n");
            break;

          case option_strand:
            if (are_same_string(optarg, "plus"))
              {
                parameters.opt_strand = false;
              }
            else if (are_same_string(optarg, "both"))
              {
                parameters.opt_strand = true;
              }
            else
              {
                fatal("The argument to --strand must be plus or both");
              }
            break;

          case option_threads:
            parameters.opt_threads = static_cast<int64_t>(args_getdouble(optarg));
            break;

          case option_gapopen:
            args_get_gap_penalty_string(optarg, true, parameters);
            break;

          case option_gapext:
            args_get_gap_penalty_string(optarg, false, parameters);
            break;

          case option_rowlen:
            parameters.opt_rowlen = args_getlong(optarg);
            break;

          case option_userfields:
            if (not parse_userfields_arg(optarg, parameters))
              {
                fatal("Unrecognized userfield argument");
              }
            break;

          case option_userout:
            parameters.opt_userout = optarg;
            break;

          case option_self:
            parameters.opt_self = 1;
            break;

          case option_blast6out:
            parameters.opt_blast6out = optarg;
            break;

          case option_uc:
            parameters.opt_uc = optarg;
            break;

          case option_weak_id:
            parameters.opt_weak_id = args_getdouble(optarg);
            break;

          case option_uc_allhits:
            parameters.opt_uc_allhits = true;
            break;

          case option_notrunclabels:
            parameters.opt_notrunclabels = true;
            break;

          case option_sortbysize:
            parameters.opt_sortbysize = optarg;
            break;

          case option_output:
            parameters.opt_output = optarg;
            break;

          case option_minsize:
            parameters.opt_minsize = args_getlong(optarg);
            if (parameters.opt_minsize <= 0)
              {
                fatal("The argument to --minsize must be at least 1");
              }
            break;

          case option_maxsize:
            parameters.opt_maxsize = args_getlong(optarg);
            break;

          case option_relabel:
            parameters.opt_relabel = optarg;
            break;

          case option_sizeout:
            parameters.opt_sizeout = true;
            break;

          case option_derep_fulllength:
            parameters.opt_derep_fulllength = optarg;
            break;

          case option_minseqlength:
            parameters.opt_minseqlength = args_getlong(optarg);
            if (parameters.opt_minseqlength < 0)
              {
                fatal("The argument to --minseqlength must not be negative");
              }
            break;

          case option_minuniquesize:
            parameters.opt_minuniquesize = args_getlong(optarg);
            break;

          case option_topn:
            parameters.opt_topn = args_getlong(optarg);
            if (parameters.opt_topn == 0)
              {
                fatal("The argument to --topn must be greater than zero");
              }
            break;

          case option_maxseqlength:
            parameters.opt_maxseqlength = args_getlong(optarg);
            break;

          case option_sizein:
            parameters.opt_sizein = true;
            break;

          case option_sortbylength:
            parameters.opt_sortbylength = optarg;
            break;

          case option_matched:
            parameters.opt_matched = optarg;
            break;

          case option_notmatched:
            parameters.opt_notmatched = optarg;
            break;

          case option_dbmatched:
            parameters.opt_dbmatched = optarg;
            break;

          case option_dbnotmatched:
            parameters.opt_dbnotmatched = optarg;
            break;

          case option_fastapairs:
            parameters.opt_fastapairs = optarg;
            break;

          case option_output_no_hits:
            parameters.opt_output_no_hits = 1;
            break;

          case option_maxhits:
            parameters.opt_maxhits = args_getlong(optarg);
            break;

          case option_top_hits_only:
            parameters.opt_top_hits_only = 1;
            break;

          case option_fasta_width:
            parameters.opt_fasta_width = args_getlong(optarg);
            break;

          case option_query_cov:
            parameters.opt_query_cov = args_getdouble(optarg);
            break;

          case option_target_cov:
            parameters.opt_target_cov = args_getdouble(optarg);
            break;

          case option_idprefix:
            parameters.opt_idprefix = args_getlong(optarg);
            break;

          case option_idsuffix:
            parameters.opt_idsuffix = args_getlong(optarg);
            break;

          case option_minqt:
            parameters.opt_minqt = args_getdouble(optarg);
            break;

          case option_maxqt:
            parameters.opt_maxqt = args_getdouble(optarg);
            break;

          case option_minsl:
            parameters.opt_minsl = args_getdouble(optarg);
            break;

          case option_maxsl:
            parameters.opt_maxsl = args_getdouble(optarg);
            break;

          case option_leftjust:
            parameters.opt_leftjust = 1;
            break;

          case option_rightjust:
            parameters.opt_rightjust = 1;
            break;

          case option_selfid:
            parameters.opt_selfid = 1;
            break;

          case option_maxid:
            parameters.opt_maxid = args_getdouble(optarg);
            break;

          case option_minsizeratio:
            parameters.opt_minsizeratio = args_getdouble(optarg);
            break;

          case option_maxsizeratio:
            parameters.opt_maxsizeratio = args_getdouble(optarg);
            break;

          case option_maxdiffs:
            parameters.opt_maxdiffs = args_getlong(optarg);
            break;

          case option_maxsubs:
            parameters.opt_maxsubs = args_getlong(optarg);
            break;

          case option_maxgaps:
            parameters.opt_maxgaps = args_getlong(optarg);
            break;

          case option_mincols:
            parameters.opt_mincols = args_getlong(optarg);
            break;

          case option_maxqsize:
            parameters.opt_maxqsize = args_getlong(optarg);
            break;

          case option_mintsize:
            parameters.opt_mintsize = args_getlong(optarg);
            break;

          case option_mid:
            parameters.opt_mid = args_getdouble(optarg);
            break;

          case option_shuffle:
            parameters.opt_shuffle = optarg;
            break;

          case option_randseed:
            parameters.opt_randseed = args_getlong(optarg);
            break;

          case option_maskfasta:
            parameters.opt_maskfasta = optarg;
            break;

          case option_hardmask:
            parameters.opt_hardmask = true;
            break;

          case option_qmask:
            if (are_same_string(optarg, "none"))
              {
                parameters.opt_qmask = Masking::none;
              }
            else if (are_same_string(optarg, "dust"))
              {
                parameters.opt_qmask = Masking::dust;
              }
            else if (are_same_string(optarg, "soft"))
              {
                parameters.opt_qmask = Masking::soft;
              }
            else
              {
                parameters.opt_qmask = Masking::error;
              }
            break;

          case option_dbmask:
            if (are_same_string(optarg, "none"))
              {
                parameters.opt_dbmask = Masking::none;
              }
            else if (are_same_string(optarg, "dust"))
              {
                parameters.opt_dbmask = Masking::dust;
              }
            else if (are_same_string(optarg, "soft"))
              {
                parameters.opt_dbmask = Masking::soft;
              }
            else
              {
                parameters.opt_dbmask = Masking::error;
              }
            break;

          case option_cluster_smallmem:
            parameters.opt_cluster_smallmem = optarg;
            break;

          case option_cluster_fast:
            parameters.opt_cluster_fast = optarg;
            break;

          case option_centroids:
            parameters.opt_centroids = optarg;
            break;

          case option_clusters:
            parameters.opt_clusters = optarg;
            break;

          case option_consout:
            parameters.opt_consout = optarg;
            break;

          case option_cons_truncate:
            std::fprintf(stderr, "WARNING: Option --cons_truncate is ignored\n");
            parameters.opt_cons_truncate = 1;
            break;

          case option_msaout:
            parameters.opt_msaout = optarg;
            break;

          case option_usersort:
            parameters.opt_usersort = 1;
            break;

          case option_xn:
            parameters.opt_xn = args_getdouble(optarg);
            break;

          case option_iddef:
            parameters.opt_iddef = args_getlong(optarg);
            break;

          case option_slots:
            std::fprintf(stderr, "WARNING: Option --slots is ignored\n");
            parameters.opt_slots = static_cast<int>(args_getlong(optarg));
            break;

          case option_pattern:
            std::fprintf(stderr, "WARNING: Option --pattern is ignored\n");
            parameters.opt_pattern = optarg;
            break;

          case option_maxuniquesize:
            parameters.opt_maxuniquesize = args_getlong(optarg);
            break;

          case option_abskew:
            parameters.opt_abskew = args_getdouble(optarg);
            break;

          case option_chimeras:
            parameters.opt_chimeras = optarg;
            break;

          case option_dn:
            parameters.opt_dn = args_getdouble(optarg);
            break;

          case option_mindiffs:
            parameters.opt_mindiffs = static_cast<int>(args_getlong(optarg));
            break;

          case option_mindiv:
            parameters.opt_mindiv = args_getdouble(optarg);
            break;

          case option_minh:
            parameters.opt_minh = args_getdouble(optarg);
            break;

          case option_nonchimeras:
            parameters.opt_nonchimeras = optarg;
            break;

          case option_uchime_denovo:
            parameters.opt_uchime_denovo = optarg;
            break;

          case option_uchime_ref:
            parameters.opt_uchime_ref = optarg;
            break;

          case option_uchimealns:
            parameters.opt_uchimealns = optarg;
            break;

          case option_uchimeout:
            parameters.opt_uchimeout = optarg;
            break;

          case option_uchimeout5:
            parameters.opt_uchimeout5 = 1;
            break;

          case option_alignwidth:
            parameters.opt_alignwidth = static_cast<int>(args_getlong(optarg));
            break;

          case option_allpairs_global:
            parameters.opt_allpairs_global = optarg;
            break;

          case option_acceptall:
            parameters.opt_acceptall = 1;
            break;

          case option_cluster_size:
            parameters.opt_cluster_size = optarg;
            break;

          case option_samout:
            parameters.opt_samout = optarg;
            break;

          case option_log:
            parameters.opt_log = optarg;
            break;

          case option_quiet:
            parameters.opt_quiet = true;
            break;

          case option_fastx_subsample:
            parameters.opt_fastx_subsample = optarg;
            break;

          case option_sample_pct:
            parameters.opt_sample_pct = args_getdouble(optarg);
            break;

          case option_fastq_chars:
            parameters.opt_fastq_chars = optarg;
            break;

          case option_profile:
            parameters.opt_profile = optarg;
            break;

          case option_sample_size:
            parameters.opt_sample_size = args_getlong(optarg);
            break;

          case option_fastaout:
            parameters.opt_fastaout = optarg;
            break;

          case option_xsize:
            parameters.opt_xsize = true;
            break;

          case option_clusterout_id:
            parameters.opt_clusterout_id = true;
            break;

          case option_clusterout_sort:
            parameters.opt_clusterout_sort = true;
            break;

          case option_borderline:
            parameters.opt_borderline = optarg;
            break;

          case option_relabel_sha1:
            parameters.opt_relabel_sha1 = true;
            break;

          case option_relabel_md5:
            parameters.opt_relabel_md5 = true;
            break;

          case option_derep_prefix:
            parameters.opt_derep_prefix = optarg;
            break;

          case option_fastq_filter:
            parameters.opt_fastq_filter = optarg;
            break;

          case option_fastqout:
            parameters.opt_fastqout = optarg;
            break;

          case option_fastaout_discarded:
            parameters.opt_fastaout_discarded = optarg;
            break;

          case option_fastqout_discarded:
            parameters.opt_fastqout_discarded = optarg;
            break;

          case option_fastx_syncpairs:
            parameters.opt_fastx_syncpairs = optarg;
            break;

          case option_fastaout_orphans:
            parameters.opt_fastaout_orphans = optarg;
            break;

          case option_fastaout_orphans_rev:
            parameters.opt_fastaout_orphans_rev = optarg;
            break;

          case option_fastqout_orphans:
            parameters.opt_fastqout_orphans = optarg;
            break;

          case option_fastqout_orphans_rev:
            parameters.opt_fastqout_orphans_rev = optarg;
            break;

          case option_read_separators:
            parameters.opt_read_separators = optarg;
            break;

          case option_fastq_truncqual:
            parameters.opt_fastq_truncqual = args_getlong(optarg);
            break;

          case option_fastq_maxee:
            parameters.opt_fastq_maxee = args_getdouble(optarg);
            break;

          case option_fastq_trunclen:
            parameters.opt_fastq_trunclen = args_getlong(optarg);
            break;

          case option_fastq_minlen:
            parameters.opt_fastq_minlen = args_getlong(optarg);
            break;

          case option_fastq_stripleft:
            parameters.opt_fastq_stripleft = args_getlong(optarg);
            break;

          case option_fastq_maxee_rate:
            parameters.opt_fastq_maxee_rate = args_getdouble(optarg);
            break;

          case option_fastq_maxns:
            parameters.opt_fastq_maxns = args_getlong(optarg);
            break;

          case option_eeout:
            parameters.opt_eeout = true;
            break;

          case option_fastq_ascii:
            parameters.opt_fastq_ascii = args_getlong(optarg);
            break;

          case option_fastq_qmin:
            parameters.opt_fastq_qmin = args_getlong(optarg);
            break;

          case option_fastq_qmax:
            parameters.opt_fastq_qmax = args_getlong(optarg);
            break;

          case option_fastq_qmaxout:
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
            parameters.opt_label_suffix = optarg;
            break;

          case option_h:
            parameters.opt_help = true;
            break;

          case option_samheader:
            parameters.opt_samheader = true;
            break;

          case option_sizeorder:
            parameters.opt_sizeorder = true;
            break;

          case option_minwordmatches:
            parameters.opt_minwordmatches = args_getlong(optarg);
            if (parameters.opt_minwordmatches < 0)
              {
                fatal("The argument to --minwordmatches must not be negative");
              }
            break;

          case option_v:
            parameters.opt_version = true;
            break;

          case option_relabel_keep:
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
            parameters.opt_fastq_convert = optarg;
            break;

          case option_fastq_asciiout:
            parameters.opt_fastq_asciiout = args_getlong(optarg);
            break;

          case option_fastq_qminout:
            parameters.opt_fastq_qminout = args_getlong(optarg);
            break;

          case option_fastq_mergepairs:
            parameters.opt_fastq_mergepairs = optarg;
            break;

          case option_fastq_eeout:
            parameters.opt_fastq_eeout = true;
            break;

          case option_fastqout_notmerged_fwd:
            parameters.opt_fastqout_notmerged_fwd = optarg;
            break;

          case option_fastqout_notmerged_rev:
            parameters.opt_fastqout_notmerged_rev = optarg;
            break;

          case option_fastq_minovlen:
            parameters.opt_fastq_minovlen = args_getlong(optarg);
            break;

          case option_fastq_minmergelen:
            parameters.opt_fastq_minmergelen = args_getlong(optarg);
            break;

          case option_fastq_maxmergelen:
            parameters.opt_fastq_maxmergelen = args_getlong(optarg);
            break;

          case option_fastq_nostagger:
            parameters.opt_fastq_nostagger = true;
            break;

          case option_fastq_allowmergestagger:
            parameters.opt_fastq_allowmergestagger = true;
            break;

          case option_fastq_maxdiffs:
            parameters.opt_fastq_maxdiffs = args_getlong(optarg);
            break;

          case option_fastaout_notmerged_fwd:
            parameters.opt_fastaout_notmerged_fwd = optarg;
            break;

          case option_fastaout_notmerged_rev:
            parameters.opt_fastaout_notmerged_rev = optarg;
            break;

          case option_reverse:
            parameters.opt_reverse = optarg;
            break;

          case option_eetabbedout:
            parameters.opt_eetabbedout = optarg;
            break;

          case option_fasta_score:
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
            std::fprintf(stderr, "WARNING: Option --xdrop_nw is ignored\n");
            break;

          case option_minhsp:
            /* minhsp ignored */
            std::fprintf(stderr, "WARNING: Option --minhsp is ignored\n");
            break;

          case option_band:
            /* band ignored */
            std::fprintf(stderr, "WARNING: Option --band is ignored\n");
            break;

          case option_hspw:
            /* hspw ignored */
            std::fprintf(stderr, "WARNING: Option --hspw is ignored\n");
            break;

          case option_gzip_decompress:
            parameters.opt_gzip_decompress = true;
            break;

          case option_bzip2_decompress:
            parameters.opt_bzip2_decompress = true;
            break;

          case option_fastq_maxlen:
            parameters.opt_fastq_maxlen = args_getlong(optarg);
            break;

          case option_fastq_truncee:
            parameters.opt_fastq_truncee = args_getdouble(optarg);
            break;

          case option_fastx_filter:
            parameters.opt_fastx_filter = optarg;
            break;

          case option_otutabout:
            parameters.opt_otutabout = optarg;
            break;

          case option_mothur_shared_out:
            parameters.opt_mothur_shared_out = optarg;
            break;

          case option_biomout:
            parameters.opt_biomout = optarg;
            break;

          case option_fastq_trunclen_keep:
            parameters.opt_fastq_trunclen_keep = args_getlong(optarg);
            break;

          case option_fastq_stripright:
            parameters.opt_fastq_stripright = args_getlong(optarg);
            break;

          case option_no_progress:
            parameters.opt_no_progress = true;
            break;

          case option_fastq_eestats2:
            parameters.opt_fastq_eestats2 = optarg;
            break;

          case option_ee_cutoffs:
            args_get_ee_cutoffs(optarg, parameters);
            break;

          case option_length_cutoffs:
            args_get_length_cutoffs(optarg, parameters);
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
            parameters.opt_cluster_unoise = optarg;
            break;

          case option_unoise_alpha:
            parameters.opt_unoise_alpha = args_getdouble(optarg);
            break;

          case option_uchime2_denovo:
            parameters.opt_uchime2_denovo = optarg;
            break;

          case option_uchime3_denovo:
            parameters.opt_uchime3_denovo = optarg;
            break;

          case option_sintax:
            parameters.opt_sintax = optarg;
            break;

          case option_sintax_cutoff:
            parameters.opt_sintax_cutoff = args_getdouble(optarg);
            break;

          case option_tabbedout:
            parameters.opt_tabbedout = optarg;
            break;

          case option_fastq_maxdiffpct:
            parameters.opt_fastq_maxdiffpct = args_getdouble(optarg);
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
            parameters.opt_fastaout_rev = optarg;
            break;

          case option_fastaout_discarded_rev:
            parameters.opt_fastaout_discarded_rev = optarg;
            break;

          case option_fastqout_rev:
            parameters.opt_fastqout_rev = optarg;
            break;

          case option_fastqout_discarded_rev:
            parameters.opt_fastqout_discarded_rev = optarg;
            break;

          case option_xee:
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
            parameters.opt_label_substr_match = true;
            break;

          case option_label:
            parameters.opt_label = optarg;
            break;

          case option_subseq_start:
            parameters.opt_subseq_start = args_getlong(optarg);
            break;

          case option_subseq_end:
            parameters.opt_subseq_end = args_getlong(optarg);
            break;

          case option_notmatchedfq:
            parameters.opt_notmatchedfq = optarg;
            break;

          case option_label_field:
            parameters.opt_label_field = optarg;
            break;

          case option_label_word:
            parameters.opt_label_word = optarg;
            break;

          case option_label_words:
            parameters.opt_label_words = optarg;
            break;

          case option_labels:
            parameters.opt_labels = optarg;
            break;

          case option_cut:
            parameters.opt_cut = optarg;
            break;

          case option_cut_pattern:
            parameters.opt_cut_pattern = optarg;
            break;

          case option_relabel_self:
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
            parameters.opt_lcaout = optarg;
            break;

          case option_lca_cutoff:
            parameters.opt_lca_cutoff = args_getdouble(optarg);
            break;

          case option_fastx_uniques:
            parameters.opt_fastx_uniques = optarg;
            break;

          case option_fastq_qout_max:
            parameters.opt_fastq_qout_max = true;
            break;

          case option_sample:
            /* truncate at first ';' or blank, per manpage: these
               characters are header separators in fasta/fastq labels
               and would break downstream parsing if left in. */
            optarg[std::strcspn(optarg, "; \t\r\n\v\f")] = '\0';
            parameters.opt_sample = optarg;
            break;

          case option_qsegout:
            parameters.opt_qsegout = optarg;
            break;

          case option_tsegout:
            parameters.opt_tsegout = optarg;
            break;

          case option_derep_smallmem:
            parameters.opt_derep_smallmem = optarg;
            break;

          case option_lengthout:
            parameters.opt_lengthout = true;
            break;

          case option_xlength:
            parameters.opt_xlength = true;
            break;

          case option_chimeras_denovo:
            parameters.opt_chimeras_denovo = optarg;
            break;

          case option_chimeras_length_min:
            parameters.opt_chimeras_length_min = static_cast<int>(args_getlong(optarg));
            break;

          case option_chimeras_parts:
            parameters.opt_chimeras_parts = static_cast<int>(args_getlong(optarg));
            break;

          case option_chimeras_parents_max:
            parameters.opt_chimeras_parents_max = static_cast<int>(args_getlong(optarg));
            break;

          case option_chimeras_diff_pct:
            parameters.opt_chimeras_diff_pct = args_getdouble(optarg);
            break;

          case option_sintax_random:
            parameters.opt_sintax_random = true;
            break;

          case option_n_mismatch:
            parameters.opt_n_mismatch = true;
            break;

          case option_fastq_minqual:
            parameters.opt_fastq_minqual = args_getlong(optarg);
            break;

          case option_fastq_truncee_rate:
            parameters.opt_fastq_truncee_rate = args_getdouble(optarg);
            break;

          case option_centroid_sizeout:
            parameters.opt_centroid_sizeout = true;
            break;

          default:
            fatal("Internal error in option parsing");
          }
      }

    /* Terminate if ambiguous or illegal options have been detected */
    if (val != -1)
      {
        std::exit(EXIT_FAILURE);
      }

    /* Terminate after reporting any extra non-option arguments */
    if (optind < argc)
      {
        fatal("Unrecognized string on command line (%s)", argv[optind]);
      }

    return options_selected;
  }

  /* Identify the single requested command (its row in valid_options), warn if
     options were given with no command, and fatally reject any option that is
     not valid for the chosen command. Returns the command's row index. */
  auto resolve_command(std::vector<bool> const & options_selected,
                       std::array<struct option, number_of_options> const & long_options) -> int
  {
    /* check that only one commmand is specified */
    int commands = 0;
    int k = -1;
    for (int i = 0; i < static_cast<int>(number_of_commands); i++)
      {
        if (options_selected[static_cast<size_t>(valid_options[static_cast<size_t>(i)][0])])
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
            std::fprintf(stderr, "WARNING: Options given, but no valid command specified.\n");
          }
      }
    else
      {
        int invalid_options = 0;
        for (int i = 0; i < option_count; i++)
          {
            if (options_selected[static_cast<size_t>(i)])
              {
                int j = 0;
                bool option_is_valid = false;
                while ((static_cast<size_t>(j) < max_number_of_options_per_command) and
                       (valid_options[static_cast<size_t>(k)][static_cast<size_t>(j)] >= 0))
                  {
                    if (valid_options[static_cast<size_t>(k)][static_cast<size_t>(j)] == i)
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
                        std::fprintf(stderr,
                                "Fatal error: Invalid options to command %s\n",
                                long_options[static_cast<size_t>(valid_options[static_cast<size_t>(k)][0])].name);
                        std::fprintf(stderr,
                                "Invalid option(s):");
                      }
                    std::fprintf(stderr, " --%s",
                            long_options[static_cast<size_t>(i)].name);
                  }
              }
          }

        if (invalid_options > 0)
          {
            std::fprintf(stderr, "\nThe valid options for the %s command are:",
                    long_options[static_cast<size_t>(valid_options[static_cast<size_t>(k)][0])].name);
            int count = 0;
            for (int j = 1;
                 (static_cast<size_t>(j) < max_number_of_options_per_command) and
                 (valid_options[static_cast<size_t>(k)][static_cast<size_t>(j)] >= 0);
                 j++)
              {
                std::fprintf(stderr, " --%s", long_options[static_cast<size_t>(valid_options[static_cast<size_t>(k)][static_cast<size_t>(j)])].name);
                ++count;
              }
            if (count == 0)
              {
                std::fprintf(stderr, " (none)");
              }
            std::fprintf(stderr, "\n");
            std::exit(EXIT_FAILURE);
          }
      }
    return k;
  }

  /* Resolve the thread count: validate the --threads range, use all cores for
     the multithreaded commands (otherwise force a single thread, warning if the
     user asked for more), and warn about --sintax --randseed across threads. */
  auto configure_threads(int k,
                         std::array<struct option, number_of_options> const & long_options,
                         struct Parameters & parameters) -> void
  {
    /* multi-threaded commands */

    if ((parameters.opt_threads < 0) or (parameters.opt_threads > n_threads_max))
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
            parameters.opt_threads = system_get_cores();
          }
      }
    else
      {
        if (parameters.opt_threads > 1)
          {
            std::fprintf(stderr, "WARNING: The %s command does not support multithreading.\nOnly 1 thread used.\n", long_options[static_cast<size_t>(valid_options[static_cast<size_t>(k)][0])].name);
          }
        parameters.opt_threads = 1;
      }
    if ((parameters.opt_sintax != nullptr) and (parameters.opt_randseed != 0) and (parameters.opt_threads > 1))
      {
        std::fprintf(stderr, "WARNING: Using the --sintax command with the --randseed option may not work as intended with multiple threads. Use a single thread (--threads 1) to ensure reproducible results.\n");
      }
  }

  /* Validate option value ranges (fatal on out-of-range input) and resolve the
     co-dependent defaults that must be settled alongside the range checks
     (weak_id, maxrejects, wordlength, ...). */
  auto validate_option_values(std::vector<bool> const & options_selected,
                              struct Parameters & parameters) -> void
  {
    if (parameters.opt_cluster_unoise != nullptr)
      {
        parameters.opt_weak_id = 0.90;
      }
    else
      if (parameters.opt_weak_id > parameters.opt_id)
        {
          parameters.opt_weak_id = parameters.opt_id;
        }

    if (parameters.opt_maxrejects == -1)
      {
        if (parameters.opt_cluster_fast != nullptr)
          {
            parameters.opt_maxrejects = 8;
          }
        else
          {
            parameters.opt_maxrejects = 32;
          }
      }

    if (parameters.opt_maxaccepts < 0)
      {
        fatal("The argument to --maxaccepts must not be negative");
      }

    if (parameters.opt_maxrejects < 0)
      {
        fatal("The argument to --maxrejects must not be negative");
      }

    if (parameters.opt_wordlength == 0)
      {
        /* set default word length */
        if (parameters.opt_orient != nullptr)
          {
            parameters.opt_wordlength = 12;
          }
        else
          {
            parameters.opt_wordlength = 8;
          }
      }

    if ((parameters.opt_wordlength < 3) or (parameters.opt_wordlength > 15))
      {
        fatal("The argument to --wordlength must be in the range 3 to 15");
      }

    if ((parameters.opt_iddef < 0) or (parameters.opt_iddef > 4))
      {
        fatal("The argument to --iddef must in the range 0 to 4");
      }

  #if 0

    if (parameters.opt_match <= 0)
      fatal("The argument to --match must be positive");

    if (parameters.opt_mismatch >= 0)
      fatal("The argument to --mismatch must be negative");

  #endif

    if ((parameters.opt_match > max_alignment_score) or
        (parameters.opt_match < -max_alignment_score))
      {
        fatal("The argument to --match must be in the range -32767 to 32767");
      }

    if ((parameters.opt_mismatch > max_alignment_score) or
        (parameters.opt_mismatch < -max_alignment_score))
      {
        fatal("The argument to --mismatch must be in the range -32767 to 32767");
      }


    if (parameters.opt_alignwidth < 0)
      {
        fatal("The argument to --alignwidth must not be negative");
      }

    if (parameters.opt_rowlen < 0)
      {
        fatal("The argument to --rowlen must not be negative");
      }

    if (parameters.opt_qmask == Masking::error)
      {
        fatal("The argument to --qmask must be none, dust or soft");
      }

    if (parameters.opt_dbmask == Masking::error)
      {
        fatal("The argument to --dbmask must be none, dust or soft");
      }

    if ((parameters.opt_sample_pct < 0.0) or (parameters.opt_sample_pct > 100.0))
      {
        fatal("The argument to --sample_pct must be in the range 0.0 to 100.0");
      }

    if (parameters.opt_sample_size < 0)
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
        fatal("The argument to --fastq_tail must be greater than zero");
      }

    if ((parameters.opt_min_unmasked_pct < 0.0) or (parameters.opt_min_unmasked_pct > 100.0))
      {
        fatal("The argument to --min_unmasked_pct must be between 0.0 and 100.0");
      }

    if ((parameters.opt_max_unmasked_pct < 0.0) or (parameters.opt_max_unmasked_pct > 100.0))
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

    if (parameters.opt_fastq_qmin > parameters.opt_fastq_qmax)
      {
        fatal("The argument to --fastq_qmin cannot be greater than --fastq_qmax");
      }

    if (parameters.opt_fastq_ascii + parameters.opt_fastq_qmin < 33)
      {
        fatal("Sum of arguments to --fastq_ascii and --fastq_qmin must be no less than 33");
      }

    if (parameters.opt_fastq_ascii + parameters.opt_fastq_qmax > 126)
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

    /* --db and the query both map "-" to a duplicate of standard input, so
       giving "-" to both makes the two readers race over the same stream and
       silently return no hits. Reject that ambiguous spelling; feed one of
       them through an explicit stream path (/dev/stdin, a named pipe, or a
       process substitution) instead. --db is valid only for the five search
       commands below, and only the active one's query field is set, so at
       most one branch is non-null. */
    {
      auto const reads_stdin = [](char const * name) -> bool {
        return (name != nullptr) and (std::strcmp(name, "-") == 0);
      };
      char const * const query =
        (parameters.opt_usearch_global != nullptr) ? parameters.opt_usearch_global :
        (parameters.opt_search_exact   != nullptr) ? parameters.opt_search_exact   :
        (parameters.opt_sintax         != nullptr) ? parameters.opt_sintax         :
        (parameters.opt_orient         != nullptr) ? parameters.opt_orient         :
        (parameters.opt_uchime_ref     != nullptr) ? parameters.opt_uchime_ref     :
        nullptr;
      if (reads_stdin(parameters.opt_db) and reads_stdin(query))
        {
          fatal("Cannot read both the query and the database from standard "
                "input; give one of them an explicit path such as /dev/stdin, "
                "a named pipe, or a process substitution");
        }
    }

    if ((parameters.opt_sintax_cutoff < 0.0) or (parameters.opt_sintax_cutoff > 1.0))
      {
        fatal("The argument to sintax_cutoff must be in the range 0.0 to 1.0");
      }

    if ((parameters.opt_lca_cutoff <= 0.5) or (parameters.opt_lca_cutoff > 1.0))
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

    if (parameters.opt_maxhits < 0)
      {
        fatal("The argument to maxhits cannot be negative");
      }

    if (parameters.opt_chimeras_length_min < 1)
      {
        fatal("The argument to chimeras_length_min must be at least 1");
      }

    if ((parameters.opt_chimeras_parents_max < 2) or (parameters.opt_chimeras_parents_max > maxparents))
      {
        std::array<char, 25> maxparents_string {{}};
        std::snprintf(maxparents_string.data(), maxparents_string.size(), "%d", maxparents);
        fatal("The argument to chimeras_parents_max must be in the range 2 to %s.\n", maxparents_string.data());
      }

    if ((parameters.opt_chimeras_diff_pct < 0.0) or (parameters.opt_chimeras_diff_pct > 50.0))
      {
        fatal("The argument to chimeras_diff_pct must be in the range 0.0 to 50.0");
      }

    if (options_selected[option_chimeras_parts] and
        ((parameters.opt_chimeras_parts < 2) or (parameters.opt_chimeras_parts > 100)))
      {
        fatal("The argument to chimeras_parts must be in the range 2 to 100");
      }

    /* --fasta_width accepts 0 to disable line wrapping (documented);
       reject only negative values. */
    if (parameters.opt_fasta_width < 0)
      {
        fatal("The argument to --fasta_width cannot be negative");
      }

    if (parameters.opt_maxseqlength < 1)
      {
        fatal("The argument to --maxseqlength must be a positive integer");
      }

    // The sequence length is narrowed to int in the search/cluster/chimera
    // engine (searchinfo_s::qseqlen, sized as qseqlen + buffer_headroom), so a
    // longer sequence would wrap negative and overflow the per-query buffer. Cap
    // the option at INT_MAX - buffer_headroom so that cannot happen; this mirrors
    // the header length limit enforced in fastx_filter_header.
    static constexpr int maxseqlength_limit =
      std::numeric_limits<int>::max() - buffer_headroom;
    if (parameters.opt_maxseqlength > maxseqlength_limit)
      {
        std::array<char, 128> message {{}};
        std::snprintf(message.data(), message.size(),
                      "The argument to --maxseqlength cannot exceed %d (INT_MAX - %d)",
                      maxseqlength_limit, buffer_headroom);
        fatal(message.data());
      }

    if (parameters.opt_chimeras_denovo != nullptr)
      {
        if (not options_selected[option_alignwidth])
          {
            parameters.opt_alignwidth = 60;
          }
      }
  }

  /* Apply the generic sentinel fixups and the command-specific defaults
     (minsize, abskew, minseqlength, sintax label handling). */
  auto apply_command_defaults(std::vector<bool> const & options_selected,
                              struct Parameters & parameters) -> void
  {
    /* TODO: check valid range of gap penalties */

    /* adapt/adjust parameters */

    /* Resolve sentinel defaults and adjust gap penalties on `parameters`.
       Generic fixups (including gap-open adjustment) are in
       vsearch_apply_defaults_fixups(Parameters&); command-specific overrides
       (abskew, minsize for unoise) follow below. */
    vsearch_apply_defaults_fixups(parameters);

    /* set default opt_minsize depending on command */
    if (parameters.opt_minsize == 0)
      {
        if (parameters.opt_cluster_unoise != nullptr)
          {
            parameters.opt_minsize = 8;
          }
        else
          {
            parameters.opt_minsize = 1;
          }
      }

    /* set default opt_abskew depending on command */
    if (not options_selected[option_abskew])
      {
        if (parameters.opt_chimeras_denovo != nullptr)
          {
            parameters.opt_abskew = 1.0;
          }
        else if (parameters.opt_uchime3_denovo != nullptr)
          {
            parameters.opt_abskew = 16.0;
          }
        else
          {
            parameters.opt_abskew = 2.0;
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
            parameters.opt_minseqlength = 32;
          }
        else
          {
            parameters.opt_minseqlength = 1;
          }
      }

    if (parameters.opt_sintax != nullptr)
      {
      parameters.opt_notrunclabels = true;
      }
  }


  /* Validate the chimera-detection options. Relocated verbatim from the former
     cmd_chimera(); runs only when a chimera command is active and after
     apply_command_defaults() has resolved opt_abskew's command-specific default,
     so it sees exactly the configuration the command dispatch used to. */
  auto validate_chimera_options(struct Parameters const & parameters) -> void
  {
    bool const is_chimera_command =
      (parameters.opt_uchime_denovo != nullptr) or
      (parameters.opt_uchime_ref != nullptr) or
      (parameters.opt_uchime2_denovo != nullptr) or
      (parameters.opt_uchime3_denovo != nullptr) or
      (parameters.opt_chimeras_denovo != nullptr);
    if (not is_chimera_command)
      {
        return;
      }

    if ((parameters.opt_chimeras == nullptr)  and (parameters.opt_nonchimeras == nullptr) and
        (parameters.opt_uchimeout == nullptr) and (parameters.opt_uchimealns == nullptr) and
        (parameters.opt_tabbedout == nullptr) and (parameters.opt_alnout == nullptr))
      {
        fatal("No output files specified");
      }

    if ((parameters.opt_uchime_ref != nullptr) and (parameters.opt_db == nullptr))
      {
        fatal("Database filename not specified with --db");
      }

    if (parameters.opt_abskew < 1.0)
      {
        fatal("Argument to --abskew must be >= 1.0");
      }

    if (parameters.opt_xn <= 1.0)
      {
        fatal("Argument to --xn must be > 1");
      }

    if (parameters.opt_dn <= 0.0)
      {
        fatal("Argument to --dn must be > 0");
      }

    if ((parameters.opt_uchime2_denovo == nullptr) and (parameters.opt_uchime3_denovo == nullptr))
      {
        if (parameters.opt_mindiffs <= 0)
          {
            fatal("Argument to --mindiffs must be > 0");
          }

        if (parameters.opt_mindiv <= 0.0)
          {
            fatal("Argument to --mindiv must be > 0");
          }

        if (parameters.opt_minh <= 0.0)
          {
            fatal("Argument to --minh must be > 0");
          }
      }
  }
}  // end of anonymous namespace


auto args_init(int argc, char ** argv, struct Parameters & parameters) -> void
{
  parameters.progname = argv[0];

  /* The option parser (parse_options) populates this Parameters
     (self-defaulting), which is the single configuration source the run reads
     from. */

  /* Library defaults to quiet; the CLI needs output. */
  parameters.opt_quiet = false;
  parameters.opt_no_progress = false;

  opterr = 1;

  auto const long_options = build_long_options();

  auto const options_selected = parse_options(argc, argv, long_options, parameters);

  int const k = resolve_command(options_selected, long_options);

  configure_threads(k, long_options, parameters);

  validate_option_values(options_selected, parameters);

  apply_command_defaults(options_selected, parameters);

  validate_chimera_options(parameters);

  // refactoring: C++17 <filesystem> std::filesystem::is_regular_file
  // check if stderr is referring to a terminal
  //  - fileno() returns a file descriptor (fd)
  //  - isatty() returns 1 if a file descriptor refers to a terminal
  parameters.opt_stderr_is_tty = (isatty(fileno(stderr)) == 1);

  /* Configuration is now fully resolved in `parameters` (parsed options,
     command-specific defaults and the sentinel fixups above); the compute
     engines and command dispatch read it directly. */
}
