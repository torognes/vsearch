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

void results_show_uncompressed_alignment(FILE * f, char * cigar)
{
  char * p = cigar;
  while(*p)
    {
      if (*p > '9')
	fprintf(f, "%c", *p++);
      else
	{
	  int n = 0;
	  char c = 0;
	  int x = 0;
	  if (sscanf(p, "%d%c%n", &n, &c, &x) == 2)
	    {
	      for(int i = 0; i<n; i++)
		fprintf(f, "%c", c);
	      p += x;
	    }
	  else
	    fatal("bad cigar");
	}
    }
}
	    
void results_show_blast6out(struct hit * hits, int accepts,
			   char * query_head,
			   char * qsequence, long qseqlen)
{

  /* http://www.drive5.com/usearch/manual/blast6out.html */

  for(int t = 0; t < accepts; t++)
    {
      struct hit * hp = hits + t;
	  
      /*
	query label
	target label
	percent identity
	alignment length
	number of mismatches
	number of gap opens
	1-based position of start in query
	1-based position of end in query
	1-based position of start in target
	1-based position of end in target
	E-value
	bit score
      */
	     
      fprintf(blast6outfile,
	      "%s\t%s\t%.1f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t*\t*\n",
	      query_head,
	      db_getheader(hp->target),
	      hp->internal_id,
	      hp->internal_alignmentlength,
	      hp->mismatches,
	      hp->internal_gaps,
	      hp->trim_q_left + 1,
	      qseqlen - hp->trim_q_right,
	      hp->trim_t_left + 1,
	      db_getsequencelen(hp->target) - hp->trim_t_right);
    }
}

void results_show_uc(struct hit * hits, int accepts,
		    char * query_head,
		    char * qsequence, long qseqlen,
		    int allhits)
{

  /* http://www.drive5.com/usearch/manual/ucout.html */

  if (accepts == 0)
    {
      fprintf(ucfile, "N\t*\t*\t*\t.\t*\t*\t*\t%s\t*\n", query_head);
    }
  else
    {
      int limit = 1;
      if (allhits)
	limit = accepts;

      for(int t = 0; t < limit; t++)
	{
	  struct hit * hp = hits + t;
	  
	  /*
	    H
	    cluster no (0-based) (target sequence no)
	    sequence length (query)
	    percent identity
	    strand: + or -
	    0
	    0
	    compressed alignment, e.g. 9I92M14D
	    query label
	    target label
	  */
	     
	  fprintf(ucfile,
		  "H\t%ld\t%ld\t%.1f\t%c\t0\t0\t%s\t%s\t%s\n",
		  hp->target,
		  qseqlen,
		  hp->internal_id,
		  hp->strand,
		  hp->nwalignment,
		  query_head,
		  db_getheader(hp->target));
	}
    }
}


void results_show_userout(struct hit * hits, int accepts,
			 char * query_head,
			 char * qsequence, long qseqlen)
{

  /* http://drive5.com/usearch/manual/userout.html */

  for(int t = 0; t < accepts; t++)
    {
      struct hit * hp = hits + t;

      for (int c = 0; c < userfields_requested_count; c++)
	{
	  if (c>0)
	    fprintf(useroutfile, "\t");

	  int field = userfields_requested[c];
	  
	  char * tsequence;
	  long tseqlen;
          db_getsequenceandlength(hp->target, & tsequence, & tseqlen);

	  char * t_head = db_getheader(hp->target);

	  char * qrow;
	  char * trow;

	  /*
	    some modifications necessary:
	    - gaps at both ends of the alignment must be removed
	    - counting of gaps must be updated accordingly
	    - testing
	  */

	  switch (field)
	    {
	    case 0: /* query */
	      fprintf(useroutfile, "%s", query_head);
	      break;
	    case 1: /* target */
	      fprintf(useroutfile, "%s", t_head);
	      break;
	    case 2: /* evalue */
	      fprintf(useroutfile, "-1");
	      break;
	    case 3: /* id */
	      fprintf(useroutfile, "%.1f", hp->internal_id);
	      break;
	    case 4: /* pctpv */
	      fprintf(useroutfile, "%.1f", 100.0 * hp->matches / hp->internal_alignmentlength);
	      break;
	    case 5: /* pctgaps */
	      fprintf(useroutfile, "%.1f", 100.0 * hp->internal_indels / hp->internal_alignmentlength);
	      break;
	    case 6: /* pairs */
	      fprintf(useroutfile, "%ld", hp->matches + hp->mismatches);
	      break;
	    case 7: /* gaps */
	      fprintf(useroutfile, "%ld", hp->internal_indels);
	      break;
	    case 8: /* qlo */
	      fprintf(useroutfile, "%d", 1);
	      break;
	    case 9: /* qhi */
	      fprintf(useroutfile, "%ld", qseqlen);
	      break;
	    case 10: /* tlo */
	      fprintf(useroutfile, "%d", 1);
	      break;
	    case 11: /* thi */
	      fprintf(useroutfile, "%ld", tseqlen);
	      break;
	    case 12: /* pv */
	      fprintf(useroutfile, "%ld", hp->matches);
	      break;
	    case 13: /* ql */
	      fprintf(useroutfile, "%ld", qseqlen);
	      break;
	    case 14: /* tl */
	      fprintf(useroutfile, "%ld", tseqlen);
	      break;
	    case 15: /* qs */
	      fprintf(useroutfile, "%ld", qseqlen);
	      break;
	    case 16: /* ts */
	      fprintf(useroutfile, "%ld", tseqlen);
	      break;
	    case 17: /* alnlen */
	      fprintf(useroutfile, "%ld", hp->internal_alignmentlength);
	      break;
	    case 18: /* opens */
	      fprintf(useroutfile, "%ld", hp->internal_gaps);
	      break;
	    case 19: /* exts */
	      fprintf(useroutfile, "%ld", hp->internal_indels - hp->internal_gaps);
	      break;
	    case 20: /* raw */
	      fprintf(useroutfile, "%d", 0);
	      break;
	    case 21: /* bits */
	      fprintf(useroutfile, "%d", 0);
	      break;
	    case 22: /* aln */ 
	      results_show_uncompressed_alignment(useroutfile, hp->nwalignment);
	      break;
	    case 23: /* caln */
	      fprintf(useroutfile, "%s", hp->nwalignment);
	      break;
	    case 24: /* qstrand */
	      fprintf(useroutfile, "%c", hp->strand);
	      break;
	    case 25: /* tstrand */
	      fprintf(useroutfile, "%c", '+');
	      break;
	    case 26: /* qrow */
	      qrow = align_getrow(qsequence,
				  hp->nwalignment,
				  hp->nwalignmentlength,
				  0);
	      fprintf(useroutfile, "%.*s",
		      (int)(hp->internal_alignmentlength),
		      qrow + hp->trim_q_left + hp->trim_t_left);
	      free(qrow);
	      break;
	    case 27: /* trow */
	      trow = align_getrow(tsequence,
				  hp->nwalignment,
				  hp->nwalignmentlength,
				  1);
	      fprintf(useroutfile, "%.*s",
		      (int)(hp->internal_alignmentlength),
		      trow + hp->trim_q_left + hp->trim_t_left);
	      free(trow);
	      break;
	    case 28: /* qframe */
	      fprintf(useroutfile, "+0");
	      break;
	    case 29: /* tframe */
	      fprintf(useroutfile, "+0");
	      break;
	    case 30: /* mism */
	      fprintf(useroutfile, "%ld", hp->mismatches);
	      break;
	    case 31: /* ids */
	      fprintf(useroutfile, "%ld", hp->matches);
	      break;
	    case 32: /* qcov */
	      fprintf(useroutfile, "%.0f",
		      100.0 * (hp->matches + hp->mismatches) / qseqlen);
	      break;
	    case 33: /* tcov */
	      fprintf(useroutfile, "%.0f",
		      100.0 * (hp->matches + hp->mismatches) / tseqlen);
	      break;
	    }
	}
      fprintf(useroutfile, "\n");
    }
}


void results_show_alnout(struct hit * hits, int accepts,
                        char * query_head,
                        char * qsequence, long qseqlen)
{

  /* http://drive5.com/usearch/manual/alnout.html */

  if (accepts > 0)
    {
      fprintf(alnoutfile,"Query >%s\n", query_head);
      fprintf(alnoutfile," %%Id   TLen  Target\n");
      
      for(int t = 0; t < accepts; t++)
        fprintf(alnoutfile,"%3.0f%% %6lu  %s\n",
                hits[t].internal_id,
                db_getsequencelen(hits[t].target),
                db_getheader(hits[t].target));
      
      fprintf(alnoutfile,"\n");
      
      for(int t = 0; t < accepts; t++)
        {
          unsigned int target = hits[t].target;
          unsigned int count = hits[t].count;
          char * dseq;
          long dseqlen;
          
          db_getsequenceandlength(target, & dseq, & dseqlen);
          
          unsigned long nwscore = hits[t].nwscore;
          unsigned long nwdiff = hits[t].nwdiff;
          unsigned long nwgaps = hits[t].nwgaps;
          unsigned long nwindels = hits[t].nwindels;
          unsigned long nwalignmentlength = hits[t].nwalignmentlength;
          char * nwalignment = hits[t].nwalignment;
          double nwid = hits[t].nwid;
          char * thead = db_getheader(target);
              
          fprintf(alnoutfile," Query %ldnt >%s\n", qseqlen, query_head);
          fprintf(alnoutfile,"Target %ldnt >%s\n", dseqlen, thead);

          align_show(alnoutfile,
		     qsequence,
		     qseqlen,
		     "Qry",
		     dseq,
		     dseqlen,
		     "Tgt",
		     nwalignment,
		     3,
		     3,
		     rowlen);
              
          fprintf(alnoutfile,"\n%ld cols, %ld ids (%3.1f%%), %ld gaps (%3.1f%%)",
                  nwalignmentlength,
                  nwalignmentlength - nwdiff,
                  nwid,
                  nwindels,
                  100.0 * nwindels / nwalignmentlength);

#if 1
          fprintf(alnoutfile," [%u kmers, %lu costs, %lu gap opens]\n",
                  count, nwscore, nwgaps);
#endif
          
          fprintf(alnoutfile, "\n");
	}
    }
}
