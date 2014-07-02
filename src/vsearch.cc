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

char * progname;
char * databasefilename;
char * queryfilename;

char * alnoutfilename;
char * useroutfilename;
char * blast6outfilename;
char * ucfilename;

int wordlength;
int maxrejects;
int maxaccepts;
int match_score;
int mismatch_score;
int gapopen_penalty;
int gapextend_penalty;
double identity;
int strand;
int threads;
int mismatch_cost;
int gapopen_cost;
int gapextend_cost;
int rowlen;
int opt_self;

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

unsigned long dbsequencecount = 0;

FILE * alnoutfile;
FILE * useroutfile;
FILE * blast6outfile;
FILE * ucfile;

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
  char title[] = PROG_NAME " " PROG_VERSION;
  char ref[] = "Copyright (C) 2014 Torbjorn Rognes";
  fprintf(stderr, "%s [%s %s]\n%s\n\n", title, __DATE__, __TIME__, ref);
}


void args_usage()
{
  /*               0         1         2         3         4         5         6         7          */
  /*               01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr, "Usage: %s [OPTIONS] [filename]\n", progname);
  fprintf(stderr, "  --help                     display this help and exit\n");
  fprintf(stderr, "  --version                  display version information and exit\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --usearch_global FILENAME  query filename (FASTA) for global alignment search\n");
  fprintf(stderr, "  --db FILENAME              database filename (FASTA)\n");
  fprintf(stderr, "  --id REAL                  minimum sequence identity accepted\n");
  fprintf(stderr, "  --alnout FILENAME          alignment output filename\n");
  fprintf(stderr, "  --userout FILENAME         user-defined tab-separated output filename\n");
  fprintf(stderr, "  --userfields STRINGS       fields to output to file specified with --userout\n");
  fprintf(stderr, "  --blast6out FILENAME       blast6-like output filename\n");
  fprintf(stderr, "  --uc FILENAME              UCLUST-like output filename\n");
  fprintf(stderr, "  --strand plus|both         search plus strand or both\n");
  fprintf(stderr, "  --maxaccepts INT           maximum number of hits to show (1)\n");
  fprintf(stderr, "  --maxrejects INT           number of non-matching hits to consider (32)\n");
  fprintf(stderr, "  --match INT                score for match (1)\n");
  fprintf(stderr, "  --mismatch INT             score for mismatch (-2)\n");
  fprintf(stderr, "  --gapopen INT              penalty for gap opening (10)\n");
  fprintf(stderr, "  --gapext INT               penalty for gap extension (1)\n");
  fprintf(stderr, "  --wordlength INT           length of words (kmers) used for database index (8)\n");
  fprintf(stderr, "  --fulldp                   full dynamic programming for all alignments\n");
  fprintf(stderr, "  --threads INT              number of threads to use (1)\n");
  fprintf(stderr, "  --rowlen INT               width of sequence alignment lines (64)\n");
  fprintf(stderr, "  --self                     ignore hits with identical label as the query\n");
}

void args_get_gap_penalty_string(char * optarg, int is_open)
{
  /* See http://www.drive5.com/usearch/manual/aln_params.html
     
     --gapopen *E/10I/1E/2L/3RQ/4RT/1IQ
     
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
  
  char *p = optarg;

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
	  pen = 1000;
	  p++;
	}
      else
	fatal("Invalid gap open penalty argument");

      char * q = p;

      int set_E = 0;
      int set_I = 0;
      int set_L = 0;
      int set_R = 0;
      int set_Q = 0;
      int set_T = 0;

      while(*p)
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
	    case '\0':
	    case '/':
	      break;
	    default:
	      fatal("Invalid char '%.1s' in gap penalty string", p);
	      break;
	    }
	  p++;
	}
      
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


long args_getint(char * optarg)
{
  int len = 0;
  long temp = 0;
  int ret = sscanf(optarg, "%ld%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(optarg)))
    fatal("Illegal option argument");
  return temp;
}

void args_init(int argc, char **argv)
{
  /* Set defaults */

  progname = argv[0];

  databasefilename = NULL;
  alnoutfilename = NULL;
  useroutfilename = NULL;
  queryfilename = NULL;
  
  wordlength = 8;
  maxrejects = 32;
  maxaccepts = 1;
  match_score = 1;
  mismatch_score = -2;
  identity = -1.0;
  strand = -1;
  threads = 1;
  gapopen_penalty = 10;
  gapextend_penalty = 1;
  rowlen = 64;
  opt_self = 0;

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
      args_usage();
      exit(1);
      break;

    case 1:
      /* version */
      show_header();
      exit(1);
      break;

    case 2:
      /* alnout */
      alnoutfilename = optarg;
      break;
          
    case 3:
      /* usearch_global */
      queryfilename = optarg;
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
      maxaccepts = args_getint(optarg);
      break;

    case 7:
      /* maxrejects */
      maxrejects = args_getint(optarg);
      break;

    case 8:
      /* wordlength */
      wordlength = args_getint(optarg);
      break;

    case 9:
      /* match */
      match_score = args_getint(optarg);
      break;

    case 10:
      /* mismatch */
      mismatch_score = args_getint(optarg);
      break;

    case 11:
      /* fulldp */
      break;

    case 12:
      /* strand */
      if (strcasecmp(optarg, "plus") == 0)
        strand = 1;
      else if (strcasecmp(optarg, "both") == 0)
        strand = 2;
      break;

    case 13:
      /* threads */
      threads = args_getint(optarg);
      break;

    case 14:
      /* gapopen */
      gapopen_penalty = args_getint(optarg);
      break;

    case 15:
      /* gapext */
      gapextend_penalty = args_getint(optarg);
      break;

    case 16:
      /* rowlen */
      rowlen = args_getint(optarg);
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
      
    default:
      args_usage();
      exit(1);
      break;
    }
  }
  
  if (c != -1)
    exit(1);
  
  if (!queryfilename)
    fatal("Query filename not specified with --usearch_global");

  if (!databasefilename)
    fatal("Database filename not specified with --db");

  if ((!alnoutfilename) && (!useroutfilename) &&
      (!ucfilename) && (!blast6outfilename))
    fatal("No output files specified");
  
  if ((identity < 0.0) || (identity > 100.0))
    fatal("Identity between 0.0 and 1.0 must be specified with --id");
  
  if (maxaccepts < 0)
    fatal("The argument to --maxaccepts must not be negative");

  if (maxrejects < 0)
    fatal("The argument to --maxrejects must not be negative");

  if ((threads < 1) || (threads > 256))
    fatal("The argument to --threads must be in the range 1 to 256");

  if ((wordlength < 3) || (wordlength > 15))
    fatal("The argument to --wordlength must be in the range 3 to 15");

  if (strand < 0)
    fatal("The argument to required option --strand must be plus or both");

  if (gapopen_penalty < 0)
    fatal("The argument to --gapopen must not be negative");

  if (gapextend_penalty < 0)
    fatal("The argument to --gapext must not be negative");

  if (match_score <= 0)
    fatal("The argument to --match must be positive");

  if (mismatch_score >= 0)
    fatal("The argument to --mismatch must be negative");

  if(rowlen < 0)
    fatal("The argument to --rowlen must not be negative");

  /* currently unsupported option arguments */
  
  if (threads != 1)
    fatal("Only one thread is currently not supported");

  if (strand != 1)
    fatal("Searching on both strands currently not supported");
}


int main(int argc, char** argv)
{
  setlocale(LC_ALL, "en_US");

  cpu_features_detect();

  /*
  if (!(sse2_present && popcnt_present))
    fatal("This program requires a processor with SSE2 and POPCNT instructions.\n");
  */
  

  args_init(argc, argv);


  /* open output files */

  if (alnoutfilename)
    {
      alnoutfile = fopen(alnoutfilename, "w");
      if (! alnoutfile)
	fatal("Unable to open alignment output file for writing");
    }

  if (useroutfilename)
    {
      useroutfile = fopen(useroutfilename, "w");
      if (! useroutfile)
	fatal("Unable to open user-defined output file for writing");
    }

  if (blast6outfilename)
    {
      blast6outfile = fopen(blast6outfilename, "w");
      if (! blast6outfile)
	fatal("Unable to open blast6-like output file for writing");
    }

  if (ucfilename)
    {
      ucfile = fopen(ucfilename, "w");
      if (! ucfile)
	fatal("Unable to open uclust-like output file for writing");
    }

  show_header();
  

  /* convert score/penalty similary system to cost distance system (Sellers) */

  mismatch_cost = match_score - mismatch_score;
  gapopen_cost = gapopen_penalty;
  gapextend_cost = match_score + gapextend_penalty;

  int cost_factor = gcd(gcd(mismatch_cost, gapopen_cost), gapextend_cost);
  
  mismatch_cost /= cost_factor;
  gapopen_cost /= cost_factor;
  gapextend_cost /= cost_factor;

  if (rowlen == 0)
    rowlen = 60;


  //  cpu_features_show();

  db_read(databasefilename);

  dbsequencecount = db_getsequencecount();

  query_open(queryfilename);


  long start_time, stop_time, difference;

  start_time = getusec();

  dbindex_build();

  search();

  dbindex_free();

  stop_time = getusec();

  difference = stop_time - start_time;


  query_close();

  db_free();


  /* close output files */

  if (ucfilename)
    fclose(ucfile);

  if (blast6outfilename)
    fclose(blast6outfile);

  if (useroutfilename)
    fclose(useroutfile);

  if (alnoutfilename)
    fclose(alnoutfile);
}
