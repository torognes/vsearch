/*
    Copyright (C) 2014-2015 Torbjorn Rognes & Frederic Mahe

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

#include "vsearch.h"

/* options */

bool opt_cluster_id;
bool opt_cluster_sort;
bool opt_quiet;
bool opt_xsize;
char * opt_alnout;
char * opt_allpairs_global;
char * opt_blast6out;
char * opt_centroids;
char * opt_chimeras;
char * opt_cluster_fast;
char * opt_cluster_smallmem;
char * opt_cluster_size;
char * opt_clusters;
char * opt_consout;
char * opt_db;
char * opt_dbmatched;
char * opt_dbnotmatched;
char * opt_derep_fulllength;
char * opt_fastaout;
char * opt_fastapairs;
char * opt_fastq_chars;
char * opt_fastx_subsample;
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
double opt_abskew;
double opt_dn;
double opt_id;
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
long opt_dbmask;
long opt_fasta_width;
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

FILE * fp_log = 0;
time_t time_start;
time_t time_finish;

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
          pen = 10000;
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

void args_init(int argc, char **argv)
{
  /* Set defaults */

  progname = argv[0];

  opt_abskew = 2.0;
  opt_acceptall = 0;
  opt_alignwidth = 80;
  opt_allpairs_global = 0;
  opt_alnout = 0;
  opt_centroids = 0;
  opt_chimeras = 0;
  opt_cluster_fast = 0;
  opt_cluster_id = 0;
  opt_cluster_size = 0;
  opt_cluster_smallmem = 0;
  opt_cluster_sort = 0;
  opt_clusters = 0;
  opt_cons_truncate = 0;
  opt_consout = 0;
  opt_db = 0;
  opt_dbmask = MASK_DUST;
  opt_dbmatched = 0;
  opt_dbnotmatched = 0;
  opt_derep_fulllength = 0;
  opt_dn = 1.4;
  opt_fasta_width = 80;
  opt_fastapairs = 0;
  opt_fastaout = 0;
  opt_fastq_chars = 0;
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
  opt_idprefix = 0;
  opt_idsuffix = 0;
  opt_leftjust = 0;
  opt_log = 0;
  opt_match = 2;
  opt_matched = 0;
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
  opt_rightjust = 0;
  opt_rowlen = 64;
  opt_samout = 0;
  opt_sample_pct = 0;
  opt_sample_size = 0;
  opt_self = 0;
  opt_selfid = 0;
  opt_shuffle = 0;
  opt_sizein = 0;
  opt_sizeout = 0;
  opt_slots = 0;
  opt_sortbylength = 0;
  opt_sortbysize = 0;
  opt_strand = 1;
  opt_target_cov = 0.0;
  opt_threads = 0;
  opt_top_hits_only = 0;
  opt_topn = LONG_MAX;
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
    {"cluster_id",            no_argument,       0, 0 },
    {"cluster_sort",          no_argument,       0, 0 },
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
          /* help */
          opt_help = 1;
          break;
              
        case 1:
          /* version */
          opt_version = 1;
          break;

        case 2:
          /* alnout */
          opt_alnout = optarg;
          break;
          
        case 3:
          /* usearch_global */
          opt_usearch_global = optarg;
          break;

        case 4:
          /* db */
          opt_db = optarg;
          break;

        case 5:
          /* id */
          opt_id = atof(optarg);
          break;

        case 6:
          /* maxaccepts */
          opt_maxaccepts = args_getlong(optarg);
          break;

        case 7:
          /* maxrejects */
          opt_maxrejects = args_getlong(optarg);
          break;

        case 8:
          /* wordlength */
          opt_wordlength = args_getlong(optarg);
          break;

        case 9:
          /* match */
          opt_match = args_getlong(optarg);
          break;

        case 10:
          /* mismatch */
          opt_mismatch = args_getlong(optarg);
          break;

        case 11:
          /* fulldp */
          opt_fulldp = 1;
          break;

        case 12:
          /* strand */
          if (strcasecmp(optarg, "plus") == 0)
            opt_strand = 1;
          else if (strcasecmp(optarg, "both") == 0)
            opt_strand = 2;
          else
            opt_strand = 0;
          break;

        case 13:
          /* threads */
          opt_threads = args_getlong(optarg);
          break;

        case 14:
          /* gapopen */
          args_get_gap_penalty_string(optarg, 1);
          break;

        case 15:
          /* gapext */
          args_get_gap_penalty_string(optarg, 0);
          break;

        case 16:
          /* rowlen */
          opt_rowlen = args_getlong(optarg);
          break;

        case 17:
          /* userfields */
          if (!parse_userfields_arg(optarg))
            fatal("Unrecognized userfield argument");
          break;

        case 18:
          /* userout */
          opt_userout = optarg;
          break;
      
        case 19:
          /* self */
          opt_self = 1;
          break;
      
        case 20:
          /* blast6out */
          opt_blast6out = optarg;
          break;
      
        case 21:
          /* uc */
          opt_uc = optarg;
          break;
      
        case 22:
          /* weak_id */
          opt_weak_id = atof(optarg);
          break;

        case 23:
          /* uc_allhits */
          opt_uc_allhits = 1;
          break;

        case 24:
          /* notrunclabels */
          opt_notrunclabels = 1;
          break;

        case 25:
          /* sortbysize */
          opt_sortbysize = optarg;
          break;

        case 26:
          /* output */
          opt_output = optarg;
          break;

        case 27:
          /* minsize */
          opt_minsize = args_getlong(optarg);
          break;

        case 28:
          /* maxsize */
          opt_maxsize = args_getlong(optarg);
          break;

        case 29:
          /* relabel */
          opt_relabel = optarg;
          break;

        case 30:
          /* sizeout */
          opt_sizeout = 1;
          break;

        case 31:
          /* derep_fulllength */
          opt_derep_fulllength = optarg;
          break;

        case 32:
          /* minseqlength */
          opt_minseqlength = args_getlong(optarg);
          break;

        case 33:
          /* minuniquesize */
          opt_minuniquesize = args_getlong(optarg);
          break;

        case 34:
          /* topn */
          opt_topn = args_getlong(optarg);
          break;

        case 35:
          /* maxseqlength */
          opt_maxseqlength = args_getlong(optarg);
          break;

        case 36:
          /* sizein */
          opt_sizein = 1;
          break;

        case 37:
          /* sortbylength */
          opt_sortbylength = optarg;
          break;

        case 38:
          /* matched */
          opt_matched = optarg;
          break;

        case 39:
          /* notmatched */
          opt_notmatched = optarg;
          break;

        case 40:
          /* dbmatched */
          opt_dbmatched = optarg;
          break;

        case 41:
          /* dbnotmatched */
          opt_dbnotmatched = optarg;
          break;

        case 42:
          /* fastapairs */
          opt_fastapairs = optarg;
          break;

        case 43:
          /* sizein */
          opt_output_no_hits = 1;
          break;

        case 44:
          /* maxhits */
          opt_maxhits = args_getlong(optarg);
          break;

        case 45:
          /* top_hits_only */
          opt_top_hits_only = 1;
          break;

        case 46:
          /* fasta_width */
          opt_fasta_width = args_getlong(optarg);
          break;

        case 47:
          /* query_cov */
          opt_query_cov = atof(optarg);
          break;

        case 48:
          /* target_cov */
          opt_target_cov = atof(optarg);
          break;

        case 49:
          /* idprefix */
          opt_idprefix = args_getlong(optarg);
          break;

        case 50:
          /* idsuffix */
          opt_idsuffix = args_getlong(optarg);
          break;

        case 51:
          /* minqt */
          opt_minqt = atof(optarg);
          break;

        case 52:
          /* maxqt */
          opt_maxqt = atof(optarg);
          break;

        case 53:
          /* minsl */
          opt_minsl = atof(optarg);
          break;

        case 54:
          /* maxsl */
          opt_maxsl = atof(optarg);
          break;

        case 55:
          /* leftjust */
          opt_leftjust = 1;
          break;

        case 56:
          /* rightjust */
          opt_rightjust = 1;
          break;

        case 57:
          /* selfid */
          opt_selfid = 1;
          break;

        case 58:
          /* maxid */
          opt_maxid = atof(optarg);
          break;

        case 59:
          /* minsizeratio */
          opt_minsizeratio = atof(optarg);
          break;

        case 60:
          /* maxsizeratio */
          opt_maxsizeratio = atof(optarg);
          break;

        case 61:
          /* maxdiffs */
          opt_maxdiffs = args_getlong(optarg);
          break;

        case 62:
          /* maxsubs */
          opt_maxsubs = args_getlong(optarg);
          break;

        case 63:
          /* maxgaps */
          opt_maxgaps = args_getlong(optarg);
          break;

        case 64:
          /* mincols */
          opt_mincols = args_getlong(optarg);
          break;

        case 65:
          /* maxqsize */
          opt_maxqsize = args_getlong(optarg);
          break;

        case 66:
          /* mintsize */
          opt_mintsize = args_getlong(optarg);
          break;

        case 67:
          /* mid */
          opt_mid = atof(optarg);
          break;

        case 68:
          /* shuffle */
          opt_shuffle = optarg;
          break;

        case 69:
          /* randseed */
          opt_randseed = args_getlong(optarg);
          break;

        case 70:
          /* mask */
          opt_maskfasta = optarg;
          break;

        case 71:
          /* hardmask */
          opt_hardmask = 1;
          break;

        case 72:
          /* qmask */
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
          /* dbmask */
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
          /* cluster_smallmem */
          opt_cluster_smallmem = optarg;
          break;

        case 75:
          /* cluster_fast */
          opt_cluster_fast = optarg;
          break;

        case 76:
          /* centroids */
          opt_centroids = optarg;
          break;

        case 77:
          /* clusters */
          opt_clusters = optarg;
          break;

        case 78:
          /* consout */
          opt_consout = optarg;
          break;

        case 79:
          /* cons_truncate */
          fprintf(stderr, "WARNING: Option --cons_truncate is ignored\n");
          opt_cons_truncate = 1;
          break;

        case 80:
          /* msaout */
          opt_msaout = optarg;
          break;

        case 81:
          /* usersort */
          opt_usersort = 1;
          break;

        case 82:
          /* xn */
          opt_xn = atof(optarg);
          break;

        case 83:
          /* iddef */
          opt_iddef = args_getlong(optarg);
          break;

        case 84:
          /* slots */
          fprintf(stderr, "WARNING: Option --slots is ignored\n");
          opt_slots = args_getlong(optarg);
          break;

        case 85:
          /* pattern */
          fprintf(stderr, "WARNING: Option --pattern is ignored\n");
          opt_pattern = optarg;
          break;

        case 86:
          /* maxuniquesize */
          opt_maxuniquesize = args_getlong(optarg);
          break;

        case 87:
          /* abskew */
          opt_abskew = atof(optarg);
          break;
          
        case 88:
          /* chimeras */
          opt_chimeras = optarg;
          break;
          
        case 89:
          /* dn */
          opt_dn = atof(optarg);
          break;
          
        case 90:
          /* mindiffs */
          opt_mindiffs = args_getlong(optarg);
          break;
          
        case 91:
          /* mindiv */
          opt_mindiv = atof(optarg);
          break;
          
        case 92:
          /* minh */
          opt_minh = atof(optarg);
          break;
          
        case 93:
          /* nonchimeras */
          opt_nonchimeras = optarg;
          break;
          
        case 94:
          /* uchime_denovo */
          opt_uchime_denovo = optarg;
          break;
          
        case 95:
          /* uchime_ref */
          opt_uchime_ref = optarg;
          break;
          
        case 96:
          /* uchimealns */
          opt_uchimealns = optarg;
          break;
          
        case 97:
          /* uchimeout */
          opt_uchimeout = optarg;
          break;
          
        case 98:
          /* uchimeout5 */
          opt_uchimeout5 = 1;
          break;
          
        case 99:
          /* alignwidth */
          opt_alignwidth = args_getlong(optarg);
          break;
          
        case 100:
          /* allpairs_global */
          opt_allpairs_global = optarg;
          break;
          
        case 101:
          /* acceptall */
          opt_acceptall = 1;
          break;
          
        case 102:
          /* cluster_size */
          opt_cluster_size = optarg;
          break;

        case 103:
          /* samout */
          opt_samout = optarg;
          break;

        case 104:
          /* log */
          opt_log = optarg;
          break;

        case 105:
          /* quiet */
          opt_quiet = true;
          break;

        case 106:
          /* fastx_subsample */
          opt_fastx_subsample = optarg;
          break;

        case 107:
          /* sample_pct */
          opt_sample_pct = atof(optarg);
          break;

        case 108:
          /* fastq_chars */
          opt_fastq_chars = optarg;
          break;

        case 109:
          /* profile */
          opt_profile = optarg;
          break;

        case 110:
          /* sample_size */
          opt_sample_size = args_getlong(optarg);
          break;
          
        case 111:
          /* fastaout */
          opt_fastaout = optarg;
          break;

        case 112:
          /* xsize */
          opt_xsize = 1;
          break;

        case 113:
          /* cluster_id */
          opt_cluster_id = 1;
          break;

        case 114:
          /* cluster_sort */
          opt_cluster_sort = 1;
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
  if (opt_usearch_global)
    commands++;
  if (opt_sortbysize)
    commands++;
  if (opt_sortbylength)
    commands++;
  if (opt_derep_fulllength)
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

  if ((opt_wordlength < 3) || (opt_wordlength > 15))
    fatal("The argument to --wordlength must be in the range 3 to 15");

  if ((opt_iddef < 0) || (opt_iddef > 4))
    fatal("The argument to --iddef must in the range 0 to 4");

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

  if (opt_threads == 0)
    opt_threads = sysconf(_SC_NPROCESSORS_ONLN);

  /* set default opt_minseqlength depending on command */

  if (opt_minseqlength == 0)
    {
      if (opt_cluster_smallmem || opt_cluster_fast || opt_cluster_size ||
          opt_usearch_global || opt_derep_fulllength )
        opt_minseqlength = 32;
      else
        opt_minseqlength = 1;
    }
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
              "  --help                      display help information\n"
              "  --log FILENAME              write messages, timing and memory info to file\n"
              "  --maxseqlength INT          maximum sequence length (50000)\n"
              "  --minseqlength INT          min seq length (clust/derep/search: 32, other:1)\n"
              "  --notrunclabels             do not truncate labels at first space\n"
              "  --quiet                     output just warnings and fatal errors to stderr\n"
              "  --threads INT               number of threads to use, zero for all cores (0)\n"
              "  --version                   display version information\n"
              "\n"
              "Alignment options (most searching options also apply)\n"
              "  --allpairs_global FILENAME  perform global alignment of all sequence pairs\n"
              "  --alnout FILENAME           filename for human-readable alignment output\n"
              "  --acceptall                 output all pairwise alignments\n"
              "\n"
              "Chimera detection options\n"
              "  --abskew REAL               min abundance ratio of parent vs chimera (2.0)\n"
              "  --alignwidth INT            width of alignment in uchimealn output (80)\n"
              "  --chimeras FILENAME         output chimeric sequences to file\n"
              "  --db FILENAME               reference database for --uchime_ref\n"
              "  --dn REAL                   'no' vote pseudo-count (1.4)\n"
              "  --mindiffs INT              minimum number of differences in segment (3)\n"
              "  --mindiv REAL               minimum divergence from closest parent (0.8)\n"
              "  --minh REAL                 minimum score (0.28)\n"
              "  --nonchimeras FILENAME      output non-chimeric sequences to file\n"
              "  --self                      exclude identical labels for --uchime_ref\n"
              "  --selfid                    exclude identical sequences for --uchime_ref\n"
              "  --uchime_denovo FILENAME    detect chimeras de novo\n"
              "  --uchime_ref FILENAME       detect chimeras using a reference database\n"
              "  --uchimealns FILENAME       output chimera alignments to file\n"
              "  --uchimeout FILENAME        output to chimera info to tab-separated file\n"
              "  --uchimeout5                make output compatible with uchime version 5\n"
              "  --xn REAL                   'no' vote weight (8.0)\n"
              "\n"
              "Clustering options (most searching options also apply)\n"
              "  --centroids FILENAME        output centroid sequences to FASTA file\n"
              "  --cluster_fast FILENAME     cluster sequences after sorting by length\n"
              "  --cluster_id                add cluster id info to consout and profile files\n"
              "  --cluster_size FILENAME     cluster sequences after sorting by abundance\n"
              "  --cluster_smallmem FILENAME cluster already sorted sequences (see -usersort)\n"
              "  --cluster_sort              order msaout, consout, profile by decr abundance\n"
              "  --clusters STRING           output each cluster to a separate FASTA file\n"
              "  --consout FILENAME          output cluster consensus sequences to FASTA file\n"
              "  --cons_truncate             do not ignore terminal gaps in MSA for consensus\n"
              "  --id REAL                   reject if identity lower\n"
              "  --iddef INT                 id definition, 0-4=CD-HIT,all,int,MBL,BLAST (2)\n"
              "  --msaout FILENAME           output multiple seq. alignments to FASTA file\n"
              "  --profile FILENAME          output sequence profile of each cluster to file\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --sizeout                   write cluster abundances to centroid file\n"
              "  --strand plus|both          cluster using plus or both strands (plus)\n"
              "  --uc FILENAME               filename for UCLUST-like output\n"
              "  --usersort                  indicate sequences not presorted by length\n"
              "\n"
              "Dereplication options\n"
              "  --derep_fulllength FILENAME dereplicate sequences in the given FASTA file\n"
              "  --maxuniquesize INT         maximum abundance for output from dereplication\n"
              "  --minuniquesize INT         minimum abundance for output from dereplication\n"
              "  --output FILENAME           output FASTA file\n"
              "  --sizein                    propagate abundance annotation from input\n"
              "  --sizeout                   write abundance annotation to output\n"
              "  --strand plus|both          dereplicate plus or both strands (plus)\n"
              "  --topn INT                  output just the n most abundant sequences\n"
              "  --uc FILENAME               filename for UCLUST-like output\n"
              "\n"
              "FASTQ file processing\n"
              "  --fastq_chars FILENAME      Analyse FASTQ file for version and quality range\n"
              "\n"
              "Masking options\n"
              "  --hardmask                  mask by replacing with N instead of lower case\n"
              "  --maskfasta FILENAME        mask sequences in the given FASTA file\n"
              "  --output FILENAME           output to specified FASTA file\n"
              "  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)\n"
              "\n"
              "Searching options\n"
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
              "  --mismatch INT              score for mismatch (-4)\n"
              "  --notmatched FILENAME       FASTA file for non-matching query sequences\n"
              "  --output_no_hits            output non-matching queries to output files\n"
              "  --qmask none|dust|soft      mask query with dust, soft or no method (dust)\n"
              "  --query_cov REAL            reject if fraction of query seq. aligned lower\n"
              "  --rightjust                 reject if terminal gaps at alignment right end\n"
              "  --rowlen INT                width of alignment lines in alnout output (64)\n"
              "  --samout FILENAME           filename for SAM format output\n"
              "  --self                      reject if labels identical\n"
              "  --selfid                    reject if sequences identical\n"
              "  --sizeout                   write abundance annotation to dbmatched file\n"
              "  --strand plus|both          search plus or both strands (plus)\n"
              "  --target_cov REAL           reject if fraction of target seq. aligned lower\n"
              "  --top_hits_only             output only hits with identity equal to the best\n"
              "  --uc FILENAME               filename for UCLUST-like output\n"
              "  --uc_allhits                show all, not just top hit with uc output\n"
              "  --usearch_global FILENAME   filename of queries for global alignment search\n"
              "  --userfields STRING         fields to output in userout file\n"
              "  --userout FILENAME          filename for user-defined tab-separated output\n"
              "  --weak_id REAL              include aligned hits with >= id; continue search\n"
              "  --wordlength INT            length of words for database index 3-15 (8)\n"
              "\n"
              "Shuffling options\n"
              "  --output FILENAME           output to specified FASTA file\n"
              "  --randseed INT              seed for PRNG, zero to use random data source (0)\n"
              "  --shuffle FILENAME          shuffle order of sequences pseudo-randomly\n"
              "  --topn INT                  output just first n sequences\n"
              "\n"
              "Sorting options\n"
              "  --maxsize INT               maximum abundance for sortbysize\n"
              "  --minsize INT               minimum abundance for sortbysize\n"
              "  --output FILENAME           output FASTA file\n"
              "  --relabel STRING            relabel with this prefix string after sorting\n"
              "  --sortbylength FILENAME     sort sequences by length in given FASTA file\n"
              "  --sortbysize FILENAME       abundance sort sequences in given FASTA file\n"
              "  --topn INT                  output just top n seqs after sorting\n"
              "\n"
              "Subsampling options\n"
              "  --fastaout FILENAME         output FASTA file for subsamples\n"
              "  --fastx_subsample FILENAME  subsample sequences from given file\n"
              "  --randseed INT              seed for PRNG, zero to use random data source (0)\n"
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

void cmd_derep_fulllength()
{
  if ((!opt_output) && (!opt_uc))
    fatal("Output file for derepl_fulllength must be specified with --output or --uc");
  
  derep_fulllength();
}

void cmd_shuffle()
{
  if (!opt_output)
    fatal("Output file for shuffling must be specified with --output");
  
  shuffle();
}

void cmd_subsample()
{
  if (!opt_fastaout)
    fatal("Output file for subsampling must be specified with --fastaout");

  if ((opt_sample_pct > 0) && (opt_sample_size > 0))
    fatal("Specify either --sample_pct or --sample_size, not both");

  subsample();
}

void cmd_maskfasta()
{
  if (!opt_output)
    fatal("Output file for masking must be specified with --output");
  
  maskfasta();
}

void cmd_none()
{
  if (! opt_quiet)
    fprintf(stderr,
            "For help, please enter: %s --help\n"
            "\n"
            "For further details, please see the manual by entering: man vsearch\n"
            "\n"
            "Some basic command examples:\n"
            "\n"
            "vsearch --allpairs_global FILENAME --id 0.5 --alnout FILENAME\n"
            "vsearch --cluster_fast FILENAME --id 0.97 --centroids FILENAME\n"
            "vsearch --cluster_size FILENAME --id 0.97 --centroids FILENAME\n"
            "vsearch --cluster_smallmem FILENAME --usersort --id 0.97 --centroids FILENAME\n"
            "vsearch --derep_fulllength FILENAME --output FILENAME\n"
            "vsearch --maskfasta FILENAME --output FILENAME\n"
            "vsearch --shuffle FILENAME --output FILENAME\n"
            "vsearch --sortbylength FILENAME --output FILENAME\n"
            "vsearch --sortbysize FILENAME --output FILENAME\n"
            "vsearch --uchime_denovo FILENAME --nonchimeras FILENAME\n"
            "vsearch --uchime_ref FILENAME --db FILENAME --nonchimeras FILENAME\n"
            "vsearch --usearch_global FILENAME --db FILENAME --id 0.97 --alnout FILENAME\n"
            "\n",
            progname);
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

void fillheader()
{
  snprintf(progheader, 80, 
           "%s %s_%s, %.1fGB RAM, %ld cores",
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

  if (opt_log)
    {
      fp_log = fopen(opt_log, "w");
      if (!fp_log)
        fatal("Unable to open log file for writing");
      fprintf(fp_log, "%s\n", progheader);
      fprintf(fp_log, "%s\n", cmdline);
      char time_string[26];
      time_start = time(0);
      fprintf(fp_log, "Started %s", ctime_r(&time_start, time_string));
    }

  show_header();

  if (!sse2_present)
    fatal("Sorry, this program requires a cpu with SSE2.");

  if (opt_help)
    {
      cmd_help();
    }
  else if (opt_allpairs_global)
    {
      cmd_allpairs_global();
    }
  else if (opt_usearch_global)
    {
      cmd_usearch_global();
    }
  else if (opt_sortbysize)
    {
      cmd_sortbysize();
    }
  else if (opt_sortbylength)
    {
      cmd_sortbylength();
    }
  else if (opt_derep_fulllength)
    {
      cmd_derep_fulllength();
    }
  else if (opt_shuffle)
    {
      cmd_shuffle();
    }
  else if (opt_fastx_subsample)
    {
      cmd_subsample();
    }
  else if (opt_maskfasta)
    {
      cmd_maskfasta();
    }
  else if (opt_cluster_smallmem || opt_cluster_fast || opt_cluster_size)
    {
      cmd_cluster();
    }
  else if (opt_uchime_denovo || opt_uchime_ref)
    {
      cmd_uchime();
    }
  else if (opt_fastq_chars)
    {
      fastq_chars();
    }
  else if (opt_version)
    {
    }
  else
    {
      cmd_none();
    }
  
  if (opt_log)
    {
      time_finish = time(0);
      char time_string[26];
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Finished %s", ctime_r(&time_finish, time_string));
      time_t time_diff = time_finish - time_start;
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
}
