/*
    Copyright (C) 2014 Torbjorn Rognes

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

/* ARGUMENTS AND THEIR DEFAULTS */

static char * progname;
char * databasefilename;
char * opt_usearch_global;

char * alnoutfilename;
char * useroutfilename;
char * blast6outfilename;
char * ucfilename;
char * opt_fastapairs;
char * opt_matched;
char * opt_notmatched;
char * opt_dbmatched;
char * opt_dbnotmatched;

long opt_threads;

long wordlength;
long maxrejects;
long maxaccepts;
long match_score;
long mismatch_score;
long rowlen;

double opt_id;
double opt_query_cov;
double opt_target_cov;
long opt_idprefix;
long opt_idsuffix;
double opt_minqt;
double opt_maxqt;
double opt_minsl;
double opt_maxsl;
long opt_leftjust;
long opt_rightjust;
long opt_self;
long opt_selfid;
double opt_maxid;
double opt_minsizeratio;
double opt_maxsizeratio;
long opt_maxdiffs;
long opt_maxsubs;
long opt_maxgaps;
long opt_mincols;
long opt_maxqsize;
long opt_mintsize;
double opt_mid;

double opt_weak_id;
long opt_strand;
long opt_uc_allhits;
long opt_notrunclabels;
char * opt_sortbysize;
char * opt_sortbylength;
char * opt_output;
long opt_minsize;
long opt_maxsize;
char * opt_relabel;
long opt_sizein;
long opt_sizeout;
char * opt_derep_fulllength;
long opt_minseqlength;
long opt_maxseqlength;
long opt_minuniquesize;
long opt_topn;
long opt_help;
long opt_version;
long opt_output_no_hits;
long opt_maxhits;
long opt_top_hits_only;
long opt_fasta_width;
char * opt_shuffle;
long opt_seed;

int opt_gap_open_query_left;
int opt_gap_open_target_left;
int opt_gap_open_query_interior;
int opt_gap_open_target_interior;
int opt_gap_open_query_right;
int opt_gap_open_target_right;
int opt_gap_extension_query_left;
int opt_gap_extension_target_left;
int opt_gap_extension_query_interior;
int opt_gap_extension_target_interior;
int opt_gap_extension_query_right;
int opt_gap_extension_target_right;

/* Other variables */

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

static char progheader[80];
static char * cmdline;

#define cpuid(f,a,b,c,d)                                                \
  __asm__ __volatile__ ("cpuid":                                        \
                        "=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (f));

void cpu_features_detect()
{
  unsigned int a, b, c, d;

  cpuid(0,a,b,c,d);
  unsigned int maxlevel = a & 0xff;

  if (maxlevel >= 1)
  {
    cpuid(1,a,b,c,d);
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
      cpuid(7,a,b,c,d);
      avx2_present   = (b >>  5) & 1;
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

  databasefilename = 0;
  alnoutfilename = 0;
  useroutfilename = 0;
  opt_fastapairs = 0;
  opt_matched = 0;
  opt_notmatched = 0;
  opt_dbmatched = 0;
  opt_dbnotmatched = 0;
  
  wordlength = 8;

  maxrejects = 32;
  maxaccepts = 1;

  opt_weak_id = 10.0;
  opt_strand = 1;
  opt_threads = 0;
  rowlen = 64;
  opt_uc_allhits = 0;
  opt_notrunclabels = 0;
  opt_maxhits = LONG_MAX;

  opt_help = 0;
  opt_version = 0;
  opt_usearch_global = 0;
  opt_sortbysize = 0;
  opt_sortbylength = 0;
  opt_derep_fulllength = 0;
  opt_shuffle = 0;

  opt_seed = 0;
  opt_output = 0;
  opt_minsize = 0;
  opt_maxsize = LONG_MAX;
  opt_relabel = 0;
  opt_sizein = 0;
  opt_sizeout = 0;
  opt_minseqlength = 0;
  opt_maxseqlength = 50000;
  opt_minuniquesize = 0;
  opt_topn = LONG_MAX;
  opt_output_no_hits = 0;
  opt_top_hits_only = 0;
  opt_fasta_width = 80;

  opt_id = -1.0;
  opt_query_cov = 0.0;
  opt_target_cov = 0.0;
  opt_idprefix = 0;
  opt_idsuffix = 0;
  opt_minqt = 0.0;
  opt_maxqt = 1.0e37;
  opt_minsl = 0.0;
  opt_maxsl = 1.0e37;
  opt_leftjust = 0;
  opt_rightjust = 0;
  opt_self = 0;
  opt_selfid = 0;
  opt_maxid = 1.0;
  opt_minsizeratio = 0.0;
  opt_maxsizeratio = 1.0e37;
  opt_maxdiffs = INT_MAX;
  opt_maxsubs = INT_MAX;
  opt_maxgaps = INT_MAX;
  opt_mincols = 0;
  opt_maxqsize = INT_MAX;
  opt_mintsize = 0;
  opt_mid = 0.0;

  match_score = 2;
  mismatch_score = -4;

  opt_gap_open_query_left=2;
  opt_gap_open_target_left=2;
  opt_gap_open_query_interior=20;
  opt_gap_open_target_interior=20;
  opt_gap_open_query_right=2;
  opt_gap_open_target_right=2;

  opt_gap_extension_query_left=1;
  opt_gap_extension_target_left=1;
  opt_gap_extension_query_interior=2;
  opt_gap_extension_target_interior=2;
  opt_gap_extension_query_right=1;
  opt_gap_extension_target_right=1;

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
    {"seed",                  required_argument, 0, 0 },
    { 0, 0, 0, 0 }
  };
  
  int option_index = 0;
  int c;
  
  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
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
      alnoutfilename = optarg;
      break;
          
    case 3:
      /* usearch_global */
      opt_usearch_global = optarg;
      break;

    case 4:
      /* db */
      databasefilename = optarg;
      break;

    case 5:
      /* id */
      opt_id = atof(optarg);
      break;

    case 6:
      /* maxaccepts */
      maxaccepts = args_getlong(optarg);
      break;

    case 7:
      /* maxrejects */
      maxrejects = args_getlong(optarg);
      break;

    case 8:
      /* wordlength */
      wordlength = args_getlong(optarg);
      break;

    case 9:
      /* match */
      match_score = args_getlong(optarg);
      break;

    case 10:
      /* mismatch */
      mismatch_score = args_getlong(optarg);
      break;

    case 11:
      /* fulldp */
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
      rowlen = args_getlong(optarg);
      break;

    case 17:
      /* userfields */
      if (!parse_userfields_arg(optarg))
	fatal("Unrecognized userfield argument");
      break;

    case 18:
      /* userout */
      useroutfilename = optarg;
      break;
      
    case 19:
      /* self */
      opt_self = 1;
      break;
      
    case 20:
      /* blast6out */
      blast6outfilename = optarg;
      break;
      
    case 21:
      /* uc */
      ucfilename = optarg;
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
      /* seed */
      opt_seed = args_getlong(optarg);
      break;

    default:
      fatal("Internal error in option parsing");
    }
  }
  
  if (c != -1)
    exit(1);
  
  int commands = 0;
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

  if (commands == 0)
    opt_version = 1;

  if (commands > 1)
    fatal("More than one command specified");

  if (opt_weak_id > opt_id)
    opt_weak_id = opt_id;

  if (opt_minseqlength < 0)
    fatal("The argument to --minseqlength must be positive");

  if (maxaccepts < 0)
    fatal("The argument to --maxaccepts must not be negative");

  if (maxrejects < 0)
    fatal("The argument to --maxrejects must not be negative");

  if ((opt_threads < 0) || (opt_threads > 32))
    fatal("The argument to --threads must be in the range 0 (default) to 32");

  if ((wordlength < 3) || (wordlength > 15))
    fatal("The argument to --wordlength must be in the range 3 to 15");

  if (match_score <= 0)
    fatal("The argument to --match must be positive");

  if (mismatch_score >= 0)
    fatal("The argument to --mismatch must be negative");

  if(rowlen < 0)
    fatal("The argument to --rowlen must not be negative");

  if (opt_strand < 1)
    fatal("The argument to --strand must be plus or both");
  
  /* TODO: check valid range of gap penalties */

  /* adapt/adjust parameters */

#if 1

  /* adjust gap open penalty according to convention */

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
      if (opt_sortbylength || opt_sortbysize || opt_shuffle)
	opt_minseqlength = 1;
      else
	opt_minseqlength = 32;
    }

  if (rowlen == 0)
    rowlen = 60;
}


void cmd_help()
{
  /*       0         1         2         3         4         5         6         7          */
  /*       01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr, 
	  "Usage: %s [OPTIONS] [filename]\n", progname);
  fprintf(stderr, 
	  "\n"
	  "General options:\n"
	  "  --help                      display help information\n"
	  "  --version                   display version information\n"
	  "  --fasta_width INT           width of FASTA seq lines, 0 for no wrap (80)\n"
	  "  --maxseqlength INT          maximum sequence length (50000)\n"
	  "  --minseqlength INT          min seq length (sort/shuffle:1, search/derep: 32)\n"
	  "  --notrunclabels             do not truncate labels at first space\n"
	  "  --strand plus|both          search / dereplicate plus strand or both strands\n"
	  "  --threads INT               number of threads to use, zero for all cores (0)\n"
	  "  --uc FILENAME               UCLUST-like output filename for search / derepl.\n"
	  "  --uc_allhits                show all, not just top hit with uc output\n"
	  "\n"
	  "Search options:\n"
	  "  --usearch_global FILENAME   filename of queries for global alignment search\n"
	  "  --alnout FILENAME           filename for human-readable alignment output\n"
	  "  --blast6out FILENAME        filename for blast-like tab-separated output\n"
	  "  --db FILENAME               filename for FASTA formatted database for search\n"
	  "  --dbmatched FILENAME        FASTA file for matching database sequences\n"
	  "  --dbnotmatched FILENAME     FASTA file for non-matching database sequences\n"
	  "  --fastapairs FILENAME       FASTA file with pairs of query and target\n"
	  "  --fulldp                    full dynamic programming alignment for all hits\n"
	  "  --gapext STRING             penalties for gap extension (2I/1E)\n"
	  "  --gapopen STRING            penalties for gap opening (20I/2E)\n"
	  "  --id REAL                   reject if identity lower\n"
	  "  --idprefix INT              reject if first n nucleotides do not match\n"
	  "  --idsuffix INT              reject if last n nucleotides do not match\n"
	  "  --leftjust                  reject if terminal gaps at alignment left end\n"
	  "  --match INT                 score for match (2)\n"
	  "  --matched FILENAME          FASTA file for matching query sequences\n"
	  "  --maxaccepts INT            number of hits to accept and show (1)\n"
	  "  --maxdiffs INT              reject if more substitutions or indels\n"
	  "  --maxgaps INT               reject if more indels\n"
	  "  --maxhits INT               maximum number of hits to show\n"
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
	  "  --query_cov REAL            reject if fraction of query aligned lower\n"
	  "  --rightjust                 reject if terminal gaps at alignment right end\n"
	  "  --rowlen INT                width of alignment lines in alnout output (64)\n"
	  "  --self                      reject if labels identical\n"
	  "  --selfid                    reject if sequences identical\n"
	  "  --target_cov REAL           reject if fraction of target aligned lower\n"
	  "  --top_hits_only             output only hits with identity equal to the best\n"
	  "  --userfields STRING         fields to output in userout file\n"
	  "  --userout FILENAME          filename for user-defined tab-separated output\n"
	  "  --weak_id REAL              show hits with at least this id; continue search\n"
	  "  --wordlength INT            length of words (kmers) for database index (8)\n"
	  "\n"
	  "Dereplication, sorting and shuffling options\n"
	  "  --derep_fulllength FILENAME dereplicate sequences in the given FASTA file\n"
	  "  --shuffle FILENAME          shuffle order of sequences pseudo-randomly\n"
	  "  --sortbylength FILENAME     sort sequences by length in given FASTA file\n"
	  "  --sortbysize FILENAME       abundance sort sequences in given FASTA file\n"
	  "  --maxsize INT               maximum abundance for sortbysize\n"
	  "  --minsize INT               minimum abundance for sortbysize\n"
	  "  --minuniquesize INT         minimum abundance for output from dereplication\n"
	  "  --output FILENAME           output FASTA file for derepl./sort/shuffle\n"
	  "  --relabel STRING            relabel with this prefix string after sorting\n"
	  "  --seed INT                  seed for shuffle; zero for random seed (0)\n"
	  "  --sizein                    read abundance annotation from input\n"
	  "  --sizeout                   add abundance annotation to output\n"
	  "  --topn INT                  output just top n seqs from derepl./shuffle/sort\n"
	  );
}

void cmd_usearch_global()
{
  /* check options */

  if ((!alnoutfilename) && (!useroutfilename) &&
      (!ucfilename) && (!blast6outfilename) &&
      (!opt_matched) && (!opt_notmatched) &&
      (!opt_dbmatched) && (!opt_dbnotmatched))
    fatal("No output files specified");
  
  if (!databasefilename)
    fatal("Database filename not specified with --db");
  
  if ((opt_id < 0.0) || (opt_id > 1.0))
    fatal("Identity between 0.0 and 1.0 must be specified with --id");

  search(cmdline, progheader);
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
  if ((!opt_output) && (!ucfilename))
    fatal("Output file for derepl_fulllength must be specified with --output or --uc");
  
  derep_fulllength();
}

void cmd_shuffle()
{
  if (!opt_output)
    fatal("Output file for randomization must be specified with --output");
  
  shuffle();
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

  for (int i=0; i<argc; i++)
    {
      strcat(cmdline, argv[i]);
      strcat(cmdline, " ");
    }
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "Copyright (C) 2014 Torbjorn Rognes. License: AGPL 3.0\n");
  fprintf(stdout, "https://github.com/torognes/vsearch\n");
  fprintf(stdout, "\n");
}

int main(int argc, char** argv)
{
  fillheader();
  getentirecommandline(argc, argv);
  cpu_features_detect();

  args_init(argc, argv);

  show_header();

  if (opt_help)
    {
      cmd_help();
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
}
