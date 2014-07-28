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
char * opt_matched;
char * opt_notmatched;
char * opt_dbmatched;
char * opt_dbnotmatched;

static long threads;

long wordlength;
long maxrejects;
long maxaccepts;
long match_score;
long mismatch_score;
double identity;
long rowlen;

double opt_weak_id;
int opt_strand;
int opt_self;
int opt_uc_allhits;
int opt_notrunclabels;
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


void show_header()
{
  char title[] = PROG_NAME " " PROG_VERSION "_" PROG_ARCH;
  char ref[] = "Copyright (C) 2014 Torbjorn Rognes, all rights reserved. License: AGPL 3.0";
  //  fprintf(stderr, "%s [%s %s]\n%s\n", title, __DATE__, __TIME__, ref);
  fprintf(stderr, "%s, %ldMB RAM, %ld cores\n%s\n", 
	  title,
	  arch_get_memtotal() / 1024 / 1024,
	  sysconf(_SC_NPROCESSORS_ONLN),
	  ref);
  fprintf(stderr, "https://github.com/torognes/vsearch\n");
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
  opt_matched = 0;
  opt_notmatched = 0;
  opt_dbmatched = 0;
  opt_dbnotmatched = 0;
  
  wordlength = 8;

  maxrejects = 32;
  maxaccepts = 1;

  identity = -1.0;
  opt_weak_id = 10.0;
  opt_strand = 1;
  threads = 1;
  rowlen = 64;
  opt_self = 0;
  opt_uc_allhits = 0;
  opt_notrunclabels = 0;

  opt_help = 0;
  opt_version = 0;
  opt_usearch_global = 0;
  opt_sortbysize = 0;
  opt_sortbylength = 0;
  opt_derep_fulllength = 0;

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
      identity = atof(optarg);
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
      threads = args_getlong(optarg);
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

  if (commands == 0)
    opt_version = 1;

  if (commands > 1)
    fatal("More than one command specified");

  if (opt_weak_id > identity)
    opt_weak_id = identity;

  if (opt_minseqlength < 0)
    fatal("The argument to --minseqlength must be positive");

  if (maxaccepts < 0)
    fatal("The argument to --maxaccepts must not be negative");

  if (maxrejects < 0)
    fatal("The argument to --maxrejects must not be negative");

  if ((threads < 1) || (threads > 256))
    fatal("The argument to --threads must be in the range 1 to 256");

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

  /* currently unsupported option arguments */
  
  if (threads != 1)
    fatal("Only one thread is currently not supported");


  /* adapt/adjust parameters */

#if 0

  /* adjust gap open penalty according to convention */

  opt_gap_open_query_left -= opt_gap_extension_query_left;
  opt_gap_open_target_left -= opt_gap_extension_target_left;
  opt_gap_open_query_interior -= opt_gap_extension_query_interior;
  opt_gap_open_target_interior -= opt_gap_extension_target_interior;
  opt_gap_open_query_right -= opt_gap_extension_query_right;
  opt_gap_open_target_right -= opt_gap_extension_target_right;

#endif
  
  /* set default opt_minseqlength depending on command */

  if (opt_minseqlength == 0)
    {
      if (opt_sortbysize || opt_derep_fulllength)
	opt_minseqlength = 32;
      else
	opt_minseqlength = 1;
    }

  if (rowlen == 0)
    rowlen = 60;
}


void cmd_help()
{
  /*               0         1         2         3         4         5         6         7          */
  /*               01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr, "Usage: %s [OPTIONS] [filename]\n", progname);
  fprintf(stderr, "  --help                      display this help and exit\n");
  fprintf(stderr, "  --version                   display version information and exit\n");
  fprintf(stderr, "  --usearch_global FILENAME   query filename (FASTA) - global alignment search\n");
  fprintf(stderr, "  --db FILENAME               database filename (FASTA)\n");
  fprintf(stderr, "  --id REAL                   minimum sequence identity accepted\n");
  fprintf(stderr, "  --weak_id REAL              show hits with this id, but do not terminate\n");
  fprintf(stderr, "  --alnout FILENAME           alignment output filename\n");
  fprintf(stderr, "  --userout FILENAME          user-defined tab-separated output filename\n");
  fprintf(stderr, "  --userfields STRINGS        fields to output with --userout file\n");
  fprintf(stderr, "  --blast6out FILENAME        filename for output similar to blast -m6 option\n");
  fprintf(stderr, "  --uc FILENAME               UCLUST-like output filename for search / derepl.\n");
  fprintf(stderr, "  --uc_allhits                show all, not just top hit with uc output\n");
  fprintf(stderr, "  --matched FILENAME          FASTA file for matching query sequences\n");
  fprintf(stderr, "  --notmatched FILENAME       FASTA file for non-matching query sequences\n");
  fprintf(stderr, "  --dbmatched FILENAME        FASTA file for matching database sequences\n");
  fprintf(stderr, "  --dbnotmatched FILENAME     FASTA file for non-matching database sequences\n");
  fprintf(stderr, "  --strand plus|both          search plus strand or both\n");
  fprintf(stderr, "  --maxaccepts INT            maximum number of hits to show (1)\n");
  fprintf(stderr, "  --maxrejects INT            number of non-matching hits to consider (32)\n");
  fprintf(stderr, "  --match INT                 score for match (2)\n");
  fprintf(stderr, "  --mismatch INT              score for mismatch (-4)\n");
  fprintf(stderr, "  --gapopen STRING            penalties for gap opening (20I/2E)\n");
  fprintf(stderr, "  --gapext STRING             penalties for gap extension (2I/1E)\n");
  fprintf(stderr, "  --wordlength INT            length of words (kmers) for database index (8)\n");
  fprintf(stderr, "  --fulldp                    full dynamic programming for all alignments\n");
  fprintf(stderr, "  --threads INT               number of threads to use (1)\n");
  fprintf(stderr, "  --rowlen INT                width of alignment lines in alnout output(64)\n");
  fprintf(stderr, "  --self                      ignore hits with identical label as the query\n");
  fprintf(stderr, "  --notrunclabels             do not truncate labels at first space\n");
  fprintf(stderr, "  --sortbysize FILENAME       abundance sort sequences in FASTA file specified\n");
  fprintf(stderr, "  --output FILENAME           output FASTA file for abundance sort / derepl. \n");
  fprintf(stderr, "  --minuniquesize INT         minimum abundance for output from dereplication\n");
  fprintf(stderr, "  --minsize INT               minimum abundance for sortbysize\n");
  fprintf(stderr, "  --maxsize INT               maximum abundance for sortbysize\n");
  fprintf(stderr, "  --relabel STRING            relabel sequences with this prefix string\n");
  fprintf(stderr, "  --sizein                    read abundance annotation from input\n");
  fprintf(stderr, "  --sizeout                   add abundance annotation to output\n");
  fprintf(stderr, "  --derep_fulllength FILENAME dereplicate sequences in the given FASTA file\n");
  fprintf(stderr, "  --minseqlength INT          minimum sequence length for sort and derep\n");
  fprintf(stderr, "  --maxseqlength INT          maximum sequence length (50000)\n");
  fprintf(stderr, "  --topn INT                  output only most abundant sequences from derepl.\n");

}

void cmd_usearch_global()
{
  /* check options */

  if ((!alnoutfilename) && (!useroutfilename) &&
      (!ucfilename) && (!blast6outfilename))
    fatal("No output files specified");
  
  if (!databasefilename)
    fatal("Database filename not specified with --db");
  
  if ((identity < 0.0) || (identity > 1.0))
    fatal("Identity between 0.0 and 1.0 must be specified with --id");


  search();
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

int main(int argc, char** argv)
{
  //  setlocale(LC_ALL, "en_US");

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
}
