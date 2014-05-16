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


long line_pos;

char * q_seq;
char * d_seq;

long q_start;
long d_start;

long q_pos;
long d_pos;

long alignlen;

char * q_line;
char * a_line;
char * d_line;

char ntsymbols[] = "ACGT";

FILE * out;

int poswidth = 3;
int headwidth = 5;

const char * q_name;
const char * d_name;

long q_len;
long d_len;


void putop(char c, long len)
{
  long count = len;
  while(count)
    {
      if (line_pos == 0)
	{
	  q_start = q_pos;
	  d_start = d_pos;
	}

      char qs;
      char ds;

      switch(c)
	{
	case 'M':
	  qs = q_seq[q_pos++];
	  ds = d_seq[d_pos++];
	  q_line[line_pos] = ntsymbols[(int)(qs)];
	  a_line[line_pos] = (qs == ds) ? '|' : ' ';
	  d_line[line_pos] = ntsymbols[(int)(ds)];
	  line_pos++;
	  break;

	case 'D':
	  qs = q_seq[q_pos++];
	  q_line[line_pos] = ntsymbols[(int)(qs)];
	  a_line[line_pos] = ' ';
	  d_line[line_pos] = '-';
	  line_pos++;
	  break;

	case 'I':
	  ds = d_seq[d_pos++];
	  q_line[line_pos] = '-';
	  a_line[line_pos] = ' ';
	  d_line[line_pos] = ntsymbols[(int)(ds)];
	  line_pos++;
	  break;
	}

      if ((line_pos == alignlen) || ((c == 0) && (line_pos > 0)))
	{
	  q_line[line_pos] = 0;
	  a_line[line_pos] = 0;
	  d_line[line_pos] = 0;

	  long q1 = q_start + 1;
	  if (q1 > q_len)
	    q1 = q_len;

	  long q2 = q_pos;

	  long d1 = d_start + 1;
	  if (d1 > d_len)
	    d1 = d_len;

	  long d2 = d_pos;

	  fprintf(out, "\n");
	  fprintf(out, "%*s %*ld + %s %ld\n", headwidth, q_name, poswidth,
		  q1, q_line, q2);
	  fprintf(out, "%*s %*s   %s\n",      headwidth, "",     poswidth,
		  "", a_line);
	  fprintf(out, "%*s %*ld + %s %ld\n", headwidth, d_name, poswidth,
		  d1, d_line, d2);

	  line_pos = 0;
	}
      count--;
    }
}

void showalign(FILE * f,
	       char * seq1,
	       long seq1len,
	       const char * seq1name,
	       char * seq2,
	       long seq2len,
	       const char * seq2name,
	       char * cigar,
	       int numwidth,
	       int namewidth,
	       int alignwidth)
{
  out = f;

  q_seq = seq1;
  q_len = seq1len;
  q_name = seq1name;
  
  d_seq = seq2;
  d_len = seq2len;
  d_name = seq2name;
  
  char * p = cigar;
  char * e = p + strlen(p);
  
  poswidth = numwidth;
  headwidth = namewidth;
  alignlen = alignwidth;

  q_line = (char*) xmalloc(alignwidth+1);
  a_line = (char*) xmalloc(alignwidth+1);
  d_line = (char*) xmalloc(alignwidth+1);

  q_pos = 0;
  d_pos = 0;
  
  line_pos = 0;

  while(p < e)
    {
      long len;
      int n;
      if (!sscanf(p, "%ld%n", & len, & n))
	{
	  n = 0;
	  len = 1;
	}
      p += n;
      char op = *p++;
      putop(op, len);
    }
  
  putop(0, 1);

  free(q_line);
  free(a_line);
  free(d_line);
}
	       
