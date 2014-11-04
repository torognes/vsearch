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

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include <bzlib.h>
#endif

#include "align_simd.h"
#include "maps.h"
#include "arch.h"
#include "util.h"
#include "db.h"
#include "query.h"
#include "align.h"
#include "unique.h"
#include "bitmap.h"
#include "dbindex.h"
#include "minheap.h"
#include "search.h"
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

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
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

#define PROG_NAME "vsearch"
#define PROG_VERSION "v0.2.1"


/* options */

extern char * ucfilename;
extern char * useroutfilename;
extern long match_score;
extern long maxaccepts;
extern long maxrejects;
extern long mismatch_score;

extern char * opt_alnout;
extern char * opt_blast6out;
extern char * opt_centroids;
extern char * opt_cluster_fast;
extern char * opt_cluster_smallmem;
extern char * opt_clusters;
extern char * opt_consout;
extern char * opt_db;
extern char * opt_dbmatched;
extern char * opt_dbnotmatched;
extern char * opt_derep_fulllength;
extern char * opt_fastapairs;
extern char * opt_maskfasta;
extern char * opt_matched;
extern char * opt_msaout;
extern char * opt_notmatched;
extern char * opt_output;
extern char * opt_relabel;
extern char * opt_shuffle;
extern char * opt_sortbylength;
extern char * opt_sortbysize;
extern char * opt_vsearch_global;
extern double opt_id;
extern double opt_maxid;
extern double opt_maxqt;
extern double opt_maxsizeratio;
extern double opt_maxsl;
extern double opt_mid;
extern double opt_minqt;
extern double opt_minsizeratio;
extern double opt_minsl;
extern double opt_query_cov;
extern double opt_target_cov;
extern double opt_weak_id;
extern int opt_cons_truncate;
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
extern int opt_usersort;
extern long opt_dbmask;
extern long opt_fasta_width;
extern long opt_fulldp;
extern long opt_hardmask;
extern long opt_iddef;
extern long opt_idprefix;
extern long opt_idsuffix;
extern long opt_leftjust;
extern long opt_maxdiffs;
extern long opt_maxgaps;
extern long opt_maxhits;
extern long opt_maxqsize;
extern long opt_maxseqlength;
extern long opt_maxsize;
extern long opt_maxsubs;
extern long opt_mincols;
extern long opt_minseqlength;
extern long opt_minsize;
extern long opt_mintsize;
extern long opt_minuniquesize;
extern long opt_notrunclabels;
extern long opt_output_no_hits;
extern long opt_qmask;
extern long opt_rightjust;
extern long opt_seed;
extern long opt_self;
extern long opt_selfid;
extern long opt_sizein;
extern long opt_sizeout;
extern long opt_strand;
extern long opt_threads;
extern long opt_top_hits_only;
extern long opt_topn;
extern long opt_uc_allhits;
extern long opt_rowlen;
extern long opt_wordlength;
