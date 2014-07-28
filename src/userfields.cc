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

static const char * userfields_names[] =
  {
    "query",  // 0
    "target", // 1
    "evalue", // 2
    "id",     // 3
    "pctpv",
    "pctgaps",
    "pairs",
    "gaps",
    "qlo",
    "qhi",
    "tlo",
    "thi",
    "pv",
    "ql",
    "tl",
    "qs",
    "ts",
    "alnlen",
    "opens",
    "exts",
    "raw",
    "bits",
    "aln",
    "caln",
    "qstrand",
    "tstrand",
    "qrow",
    "trow",
    "qframe",
    "tframe",
    "mism",
    "ids",
    "qcov",
    "tcov",  // 33
    0
  };

int * userfields_requested = 0;
int userfields_requested_count = 0;

int parse_userfields_arg(char * arg)
{
  // Parses the userfields option argument, e.g. query+target+id+alnlen+mism
  // and returns 1 if it is ok or 0 if not.

  char * p = arg;
  char * e = p + strlen(p); // pointer to end of string

  userfields_requested_count = 1;
  while(p<e)
    if (*p++ == '+')
      userfields_requested_count++;	

  userfields_requested = (int*) xmalloc(sizeof(int) * (unsigned long)userfields_requested_count);
  
  p = arg;

  char * q;
  unsigned long n;

  int fields = 0;

  while(1)
    {
      q = strchr(p, '+');
      if (!q)
	q = e;
      
      n = (unsigned long)(q - p);

      char ** u = (char**) userfields_names;

      while (*u)
	{
	  if ((strncmp(p, *u, n) == 0) && (strlen(*u) == n))
	    break;
	  u++;
	}

      if (!*u)    // reached end of list -> unrecognized field
	return 0; // bad argument

      int i = (int)(((const char**)u) - userfields_names);
      userfields_requested[fields++] = i;

      p = q;
      
      if (p == e)  // reached end of argument
	{
#if 0
	  printf("Userfields requested (%d):\n", userfields_requested_count);
	  for(int j=0; j<userfields_requested_count; j++)
	    {
	      printf("Field: %d\n", userfields_requested[j]);
	    }
#endif
	  return 1;  // ok
	}

      p++;
    }
}
