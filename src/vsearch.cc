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
char * alnoutfilename;
char * databasefilename;
char * queryfilename;

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


void args_usage()
{
  /*               0         1         2         3         4         5         6         7          */
  /*               01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr, "Usage: %s [OPTIONS] [filename]\n", progname);
  fprintf(stderr, "  --help                     display this help and exit\n");
  fprintf(stderr, "  --version                  display version information and exit\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --usearch_global FILENAME  search for global alignments with given query\n");
  fprintf(stderr, "  --id REAL                  minimum sequence identity accepted\n");
  fprintf(stderr, "  --alnout FILENAME          alignment output filename\n");
  fprintf(stderr, "  --strand plus|both         search plus strand or both\n");
  fprintf(stderr, "\n");
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
}

void show_header()
{
  char title[] = PROG_NAME " " PROG_VERSION;
  char ref[] = "Copyright (C) 2014 Torbjorn Rognes";
  fprintf(stderr, "%s [%s %s]\n%s\n\n", title, __DATE__, __TIME__, ref);
}

long args_getint(char * optarg)
{
  int len = 0;
  long temp = 0;
  int ret = sscanf(optarg, "%ld%n", &temp, &len);
  if ((ret == 0) || (len < strlen(optarg)))
    fatal("Illegal option argument");
  return temp;
}

void args_init(int argc, char **argv)
{
  /* Set defaults */

  progname = argv[0];

  databasefilename = NULL;
  alnoutfilename = NULL;
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

  if (!alnoutfilename)
    fatal("Alignment output filename not specified with --alnout");

  alnoutfile = fopen(alnoutfilename, "w");
  if (! alnoutfile)
    fatal("Unable to open alignment output file for writing");

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

  if (alnoutfilename)
    fclose(alnoutfile);
}
