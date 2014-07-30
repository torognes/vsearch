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

struct topscore
{
  unsigned int count;
  unsigned int seqno;
#if 0
  unsigned int length;
#endif
};
      
int hit_compare(const void * a, const void * b)
{
  struct hit * x = (struct hit *) a;
  struct hit * y = (struct hit *) b;

  // high id, then low id
  // early target, then late target

  if (x->internal_id > y->internal_id)
    return -1;
  else if (x->internal_id < y->internal_id)
    return +1;
  else
    if (x->target < y->target)
      return -1;
    else if (x->target > y->target)
      return +1;
    else
      return 0;
}

#define MAXSAMPLES 255

static long scorematrix[16][16];
static long tophits;

static unsigned char * hitcount = 0;
static struct topscore * topscores = 0;
static struct hit * hits = 0;

static unsigned int * targetlist = 0;
static unsigned int targetcount = 0;

static unsigned int kmersample[256];
static unsigned int kmersamplecount;

static unsigned long nwalignments = 0;

static FILE * fp_alnout = 0;
static FILE * fp_userout = 0;
static FILE * fp_blast6out = 0;
static FILE * fp_uc = 0;
static FILE * fp_fastapairs = 0;
static FILE * fp_matched = 0;
static FILE * fp_notmatched = 0;
static FILE * fp_dbmatched = 0;
static FILE * fp_dbnotmatched = 0;

void search_get_query_samples(char * qsequence, unsigned int qseqlen, 
                              unsigned int wl, unsigned int samples)
{
  count_kmers(wl, qsequence, qseqlen);
  unsigned int unique = count_kmers_unique();
  kmersamplecount = 0;

  unsigned int pos = 0;
  unsigned int u = 0;

  unsigned int kmer = 0;
  unsigned int mask = (1<<(2*wl)) - 1;
  char * s = qsequence;
  char * e1 = s + wl - 1;
  char * e2 = s + qseqlen;
  if (e2 < e1)
    e1 = e2;

  while (s < e1)
    {
      kmer <<= 2;
      kmer |= chrmap_2bit[(int)(*s++)];
    }

#if 0
  fprintf(stderr, 
          "Sequence length: %d  Unique kmers: %d  Samples: %d\n", 
          qseqlen, unique, samples);
#endif

  int z = 0;

  while (s < e2)
    {
      kmer <<= 2;
      kmer |= chrmap_2bit[(int)(*s++)];
      kmer &= mask;

      if (count_kmers_getcount(wl, kmer) == 1)
        {
          if (z>=0)
            {
              kmersample[kmersamplecount++] = kmer;
#if 0
              fprintf(stderr, "Query kmer sample %d at pos %d u %d: ", kmersamplecount, pos, u);
              fprint_kmer(stderr, wl, kmer);
              fprintf(stderr, "\n");
#endif
              z -= unique+1;
            }
          z += samples;
          u++;
        }
      pos++;
    }
}


unsigned int search_topscores(unsigned int seqcount, 
                              unsigned int th)
{
  /*
    Count the kmer hits in each database sequence and
    make a sorted list of a given number (th)
    of the database sequences with the highest number of matching kmers.

    These are stored in the topscores array.
    
    The number to keep: th = min(seqcount, maxaccepts + maxrejects)

    This is a partial sort.
  */


  

  /* count kmer hits in the database sequences */
  /* compute total hits */

  memset(hitcount, 0, seqcount);
  
  unsigned int totalhits = 0;

  for(unsigned int i=0; i<kmersamplecount; i++)
    {
      unsigned int kmer = kmersample[i];
      totalhits += kmerhash[kmer+1] - kmerhash[kmer];
    }

  double hitspertarget = totalhits * 1.0 / seqcount;

#if 0
  printf("Total hits: %u SeqCount: %u Hits/Target: %.2f\n", 
         totalhits, seqcount, hitspertarget);
#endif

  
  targetcount = 0;

  unsigned int sparse;

  if (hitspertarget < 0.1)
    sparse = 1;
  else
    sparse = 0;

  //  sparse = 1;

  if (! sparse)
    {
      /* dense hit distribution - check all targets - no need for a list*/

      for(unsigned int i=0; i<kmersamplecount; i++)
        {
          unsigned int kmer = kmersample[i];
          unsigned int * a = kmerindex + kmerhash[kmer];
          unsigned int * b = kmerindex + kmerhash[kmer+1];
          for(unsigned int * j=a; j<b; j++)
            hitcount[*j]++;
        }
    }
  else
    {
      /* sparse hits - check only a list of targets */
      
      /* create list with targets (no duplications) */

      for(unsigned int i=0; i<kmersamplecount; i++)
        {
          unsigned int kmer = kmersample[i];
          unsigned int * a = kmerindex + kmerhash[kmer];
          unsigned int * b = kmerindex + kmerhash[kmer+1];
          
          for(unsigned int * j=a; j<b; j++)
            {
              /* append to target list */
              if (hitcount[*j] == 0)
                targetlist[targetcount++] = *j;
              
              hitcount[*j]++;
            }
        }
      
#if 0
      printf("Unique targets: %u\n", targetcount);
#endif

    }

  unsigned int topcount = 0;
  
  if (sparse)
    {
      for(unsigned int z=0; z < targetcount; z++)
        {
          unsigned int i = targetlist[z];

          unsigned int count = hitcount[i];

#if 0
	  unsigned long seqlen = db_getsequencelen(i);
#endif
          
          /* find insertion point */
          
          unsigned int p = topcount;
          
#if 0
          while ((p > 0) &&
		 ((count > topscores[p-1].count) ||
		  ((count == topscores[p-1].count) &&
		   (seqlen < topscores[p-1].length))))
            p--;
#else
          while ((p > 0) &&
		 (count >= topscores[p-1].count))
            p--;
#endif
          
          
          /* p = index in array where new data should be placed */
          
          if (p < th)
            {
              
              /* find new bottom of list */
              
              int bottom = topcount;
              if (topcount == th)
                bottom--;
              
              
              /* shift lower counts down */
              
              for(unsigned int j = bottom; j > p; j--)
                topscores[j] = topscores[j-1];
              
              
              /* insert or overwrite */
              
              topscores[p].count  = count;
              topscores[p].seqno  = i;
#if 0
              topscores[p].length = seqlen;
#endif
              
              if (topcount < th)
                topcount++;
            }
        }
    }
  else
    {
      for(unsigned int i=0; i < seqcount; i++)
        {
          unsigned int count = hitcount[i];
#if 0
	  unsigned long seqlen = db_getsequencelen(i);
#endif
          
          /* find insertion point */
          
          unsigned int p = topcount;
          
#if 0
          while ((p > 0) &&
		 ((count > topscores[p-1].count) ||
		  ((count == topscores[p-1].count) &&
		   (seqlen < topscores[p-1].length))))
            p--;
#else
          while ((p > 0) &&
		 (count >= topscores[p-1].count))
            p--;
#endif
          
          /* p = index in array where new data should be placed */
          
          if (p < th)
            {
              
              /* find new bottom of list */
              
              int bottom = topcount;
              if (topcount == th)
                bottom--;
              
              
              /* shift lower counts down */
              
              for(unsigned int j = bottom; j > p; j--)
                topscores[j] = topscores[j-1];
              
              
              /* insert or overwrite */
              
              topscores[p].count = count;
              topscores[p].seqno = i;
#if 0
              topscores[p].length = seqlen;
#endif
              
              if (topcount < th)
                topcount++;
            }
        }
    }
  
  return topcount;
}


int search_onequery(char * query_head, char * qsequence, 
		    long qseqlen, long query_no, char * rc)
{
  int seqcount = db_getsequencecount();


  /* compute necessary number of samples */

  unsigned int minsamples = 1;
  unsigned int maxsamples = MIN(qseqlen-wordlength+1, MAXSAMPLES);
  //  double default_samples = 32.0;
  double default_samples = 255.0;
  unsigned int samples;

  //samples = ceil(default_samples * exp(wordlength * (1.0 - identity)));
  samples = round(default_samples);

  if (samples < minsamples)
    samples = minsamples;
  if (samples > maxsamples)
    samples = maxsamples;

  //  printf("seqlen: %ld samples: %d\n", qseqlen, samples);

  int hit_count = 0;
  
  for(int s = 0; s < opt_strand; s++)
    {
      /* check plus or both strands*/
      /* s=0: plus; s=1: minus */

      char * qseq;

      if (s)
	qseq = rc;
      else
	qseq = qsequence;

      
      /* extract unique kmer samples from query*/
      
      search_get_query_samples(qseq, qseqlen, wordlength, samples);
  
      
      /* find database sequences with the most kmer hits */
      
      unsigned int topcount = search_topscores(seqcount, tophits);
      
      //  printf("topcount: %u\n", topcount);
  
      /* analyse targets with the highest number of kmer hits */

      int accepts = 0;
      int rejects = 0;

      for(unsigned int t = 0; (accepts < maxaccepts) && (rejects <= maxrejects) && (t<topcount); t++)
	{
	  unsigned int target = topscores[t].seqno;
	  char * dlabel = db_getheader(target);
      
	  /* Test these accept/reject criteria without alignment:
	     --idprefix
	     --idsuffix
	     --minqt
	     --maxqt
	     --minsl
	     --maxsl
	     --self (implemented)
	     --selfid
	     --minsizeratio
	     --maxsizeratio
	     --maxqsize
	     --mintsize
	  */

	  if ((opt_self) && (strcmp(query_head, dlabel) == 0))
	    {
	      rejects++;
	    }
	  else
	    {

	      int count = topscores[t].count;

	      char * dseq;
	      long dseqlen;
      
	      db_getsequenceandlength(target, & dseq, & dseqlen);
      
	      /* compute global alignment */
      
	      long nwscore;
	      long nwdiff;
	      long nwgaps;
	      long nwindels;
	      long nwalignmentlength;
	      char * nwalignment;
	  
	      nw_align(dseq,
		       dseq + dseqlen,
		       qseq,
		       qseq + qseqlen,
		       (long*) scorematrix,
		       opt_gap_open_query_left,
		       opt_gap_open_query_interior,
		       opt_gap_open_query_right,
		       opt_gap_open_target_left,
		       opt_gap_open_target_interior,
		       opt_gap_open_target_right,
		       opt_gap_extension_query_left,
		       opt_gap_extension_query_interior,
		       opt_gap_extension_query_right,
		       opt_gap_extension_target_left,
		       opt_gap_extension_target_interior,
		       opt_gap_extension_target_right,
		       & nwscore,
		       & nwdiff,
		       & nwgaps,
		       & nwindels,
		       & nwalignmentlength,
		       & nwalignment,
		       query_no,
		       target);

	      nwalignments++;
	      
	      double nwid = (nwalignmentlength - nwdiff) * 100.0 / nwalignmentlength;

	      /* info for semi-global alignment (without gaps at ends) */
	  
	      long trim_aln_left = 0;
	      long trim_q_left = 0;
	      long trim_t_left = 0;
	      long trim_aln_right = 0;
	      long trim_q_right = 0;
	      long trim_t_right = 0;


	      /* left trim alignment */
	  
	      char * p = nwalignment;
	      long run = 1;
	      int scanlength = 0;
	      sscanf(p, "%ld%n", &run, &scanlength);
	      char op = *(p+scanlength);
	      if (op != 'M')
		{
		  trim_aln_left = 1 + scanlength;
		  if (op == 'D')
		    trim_q_left = run;
		  else
		    trim_t_left = run;
		}

	      /* right trim alignment */
	  
	      char * e = nwalignment + strlen(nwalignment);
	      p = e - 1;
	      op = *p;
	      if (op != 'M')
		{
		  while (*(p-1) <= '9')
		    p--;
		  run = 1;
		  sscanf(p, "%ld", &run);
		  trim_aln_right = e - p;
		  if (op == 'D')
		    trim_q_right = run;
		  else
		    trim_t_right = run;
		}

#if 0
	      printf("Alignment string: %s\n", nwalignment);
	      printf("Trim aln: %ld,%ld q: %ld,%ld t: %ld,%ld\n",
		     trim_aln_left, trim_aln_right,
		     trim_q_left, trim_q_right,
		     trim_t_left, trim_t_right);
#endif

	      long mismatches = nwdiff - nwindels;
	      long matches = nwalignmentlength - nwdiff;
	      long internal_alignmentlength = nwalignmentlength
		- trim_q_left - trim_t_left - trim_q_right - trim_t_right;
	      long internal_gaps = nwgaps
		- (trim_q_left  + trim_t_left  > 0 ? 1 : 0)
		- (trim_q_right + trim_t_right > 0 ? 1 : 0);
	      long internal_indels = nwindels
		- trim_q_left - trim_t_left - trim_q_right - trim_t_right;
	      double internal_id = 100.0 * matches / internal_alignmentlength;

	      /*

		More accept criteria:

		 --query_cov
		 --target_cov
		 --leftjust
		 --rightjust
		 --maxid
		 --maxdiffs
		 --maxsubs
		 --maxgaps
		 --mincols
		 --mid
		 ...

	      */

	      if (internal_id >= 100.0 * opt_weak_id)
		{
		  hits[hit_count].target = target;
		  hits[hit_count].count = count;
		  hits[hit_count].nwscore = nwscore;
		  hits[hit_count].nwdiff = nwdiff;
		  hits[hit_count].nwgaps = nwgaps;
		  hits[hit_count].nwindels = nwindels;
		  hits[hit_count].nwalignmentlength = nwalignmentlength;
		  hits[hit_count].nwalignment = nwalignment;
		  hits[hit_count].nwid = nwid;
		  hits[hit_count].strand = s;
		  hits[hit_count].matches = matches;
		  hits[hit_count].mismatches = mismatches;
		  hits[hit_count].trim_q_left = trim_q_left;
		  hits[hit_count].trim_q_right = trim_q_right;
		  hits[hit_count].trim_t_left = trim_t_left;
		  hits[hit_count].trim_t_right = trim_t_right;
		  hits[hit_count].trim_aln_left = trim_aln_left;
		  hits[hit_count].trim_aln_right = trim_aln_right;
		  hits[hit_count].internal_alignmentlength = internal_alignmentlength;
		  hits[hit_count].internal_gaps = internal_gaps;
		  hits[hit_count].internal_indels = internal_indels;
		  hits[hit_count].internal_id = internal_id;
		  hit_count++;
		  
		  if (internal_id >= 100.0 * identity)
		    accepts++;
		}
	      else
		{
		  free(nwalignment);
		  rejects++;
		}
	    }
	}  
    }
    
  /* sort and return hits */
  
  qsort(hits, hit_count, sizeof(struct hit), hit_compare);

  return hit_count;
}

void search(char * cmdline, char * progheader)
{
  /* check options */

  if ((!alnoutfilename) && (!useroutfilename) &&
      (!ucfilename) && (!blast6outfilename) &&
      (!opt_matched) && (!opt_notmatched) &&
      (!opt_dbmatched) && (!opt_dbnotmatched))
    fatal("No output files specified");

  if (!databasefilename)
    fatal("Database filename not specified with --db");

  if ((identity < 0.0) || (identity > 1.0))
    fatal("Identity between 0.0 and 1.0 must be specified with --id");


  /* open output files */

  if (alnoutfilename)
    {
      fp_alnout = fopen(alnoutfilename, "w");
      if (! fp_alnout)
	fatal("Unable to open alignment output file for writing");
    }

  if (useroutfilename)
    {
      fp_userout = fopen(useroutfilename, "w");
      if (! fp_userout)
        fatal("Unable to open user-defined output file for writing");
    }

  if (blast6outfilename)
    {
      fp_blast6out = fopen(blast6outfilename, "w");
      if (! fp_blast6out)
        fatal("Unable to open blast6-like output file for writing");
    }

  if (ucfilename)
    {
      fp_uc = fopen(ucfilename, "w");
      if (! fp_uc)
        fatal("Unable to open uclust-like output file for writing");
    }

  if (opt_fastapairs)
    {
      fp_fastapairs = fopen(opt_fastapairs, "w");
      if (! fp_fastapairs)
        fatal("Unable to open fastapairs output file for writing");
    }

  if (opt_matched)
    {
      fp_matched = fopen(opt_matched, "w");
      if (! fp_matched)
        fatal("Unable to open matched output file for writing");
    }

  if (opt_notmatched)
    {
      fp_notmatched = fopen(opt_notmatched, "w");
      if (! fp_notmatched)
        fatal("Unable to open notmatched output file for writing");
    }

  if (opt_dbmatched)
    {
      fp_dbmatched = fopen(opt_dbmatched, "w");
      if (! fp_dbmatched)
        fatal("Unable to open dbmatched output file for writing");
    }

  if (opt_dbnotmatched)
    {
      fp_dbnotmatched = fopen(opt_dbnotmatched, "w");
      if (! fp_dbnotmatched)
        fatal("Unable to open dbnotmatched output file for writing");
    }

  if (alnoutfilename)
    {
      fprintf(fp_alnout, "%s\n", cmdline);
      fprintf(fp_alnout, "%s\n", progheader);
    }

  db_read(databasefilename);

  dbindex_build();

  for(int i=0; i<16; i++)
    for(int j=0; j<16; j++)
      if (i==j)
        scorematrix[i][j] = match_score;
      else
        scorematrix[i][j] = mismatch_score;
  
  int seqcount = db_getsequencecount();
  tophits = maxrejects + maxaccepts;
  if (tophits > seqcount)
    tophits = seqcount;
  if (maxaccepts > seqcount)
    maxaccepts = seqcount;

  nw_init();

  count_kmers_init();

  hitcount = (unsigned char *) xmalloc(seqcount);
  topscores = (struct topscore *) xmalloc(sizeof(struct topscore) * tophits);
  hits = (struct hit *) xmalloc(sizeof(struct hit) * (maxaccepts+maxrejects) * opt_strand);
  targetlist = (unsigned int*) xmalloc(sizeof(unsigned int)*seqcount);

  char * query_head;
  long query_head_len;
  char * qsequence;
  long qseqlen;
  long query_no;
  long qmatches = 0;

  query_open(opt_usearch_global);

  progress_init("Searching", query_getfilesize());

  long * dbmatched = (long *) xmalloc(seqcount * sizeof(long));
  memset(dbmatched, 0, seqcount * sizeof(long));

  long rc_seq_alloc = 50001;
  char * rc = 0;

  while(query_getnext(& query_head, & query_head_len,
		      & qsequence, & qseqlen,
		      & query_no))
    {

      if (qseqlen + 1 > rc_seq_alloc)
	{
	  rc_seq_alloc += 10000;
	  rc = (char *) xrealloc(rc, (size_t)rc_seq_alloc);
	}

      if (opt_strand > 1)
	reverse_complement(rc, qsequence, qseqlen);

      /* perform search */
      
      int hit_count = search_onequery(query_head,
				      qsequence, qseqlen, 
				      query_no, rc);
      
      if (hit_count > 0)
	qmatches++;

      /* show results */
      
      long toreport = MIN(opt_maxhits, hit_count);

      if (toreport)
	{
	  if (fp_alnout)
	    results_show_alnout(fp_alnout,
				hits, toreport, query_head,
				qsequence, qseqlen, rc);
      
	  double top_hit_id = hits[0].internal_id;
	  
	  for(int t = 0; t < toreport; t++)
	    {
	      struct hit * hp = hits + t;
	      
	      if (opt_top_hits_only && (hp->internal_id < top_hit_id))
		break;
	      
	      if (fp_fastapairs)
		results_show_fastapairs_one(fp_fastapairs, hp, query_head,
					    qsequence, qseqlen, rc);

	      if (fp_uc)
		results_show_uc_one(fp_uc, hp, query_head,
				    qsequence, qseqlen, rc);
	      
	      if (fp_userout)
		results_show_userout_one(fp_userout, hp, query_head, 
					 qsequence, qseqlen, rc);
	      
	      if (fp_blast6out)
		results_show_blast6out_one(fp_blast6out, hp, query_head,
					   qsequence, qseqlen, rc);
	    }
	}
      else
	{
	  if (fp_uc)
	    results_show_uc_one(fp_uc, 0, query_head,
				qsequence, qseqlen, rc);
	  
	  if (opt_output_no_hits)
	    {
	      if (fp_userout)
		results_show_userout_one(fp_userout, 0, query_head, 
					 qsequence, qseqlen, rc);
	      
	      if (fp_blast6out)
		results_show_blast6out_one(fp_blast6out, 0, query_head,
					   qsequence, qseqlen, rc);
	    }
	}


      if (hit_count)
	{
	  if (opt_matched)
	    {
	      fprintf(fp_matched, ">%s\n", query_head);
	      fprint_fasta_seq_only(fp_matched, qsequence, qseqlen, 80);
	    }
	}
      else
	{
	  if (opt_notmatched)
	    {
	      fprintf(fp_notmatched, ">%s\n", query_head);
	      fprint_fasta_seq_only(fp_notmatched, qsequence, qseqlen, 80);
	    }
	}


      /* update matching db sequences */

      for (long i=0; i<hit_count; i++)
	if (hits[i].internal_id >= 100.0 * identity)
	  dbmatched[hits[i].target]++;
      
      
      /* free memory for alignment strings */
      
      for(int i=0; i<hit_count; i++)
	free(hits[i].nwalignment);
      
      if(rc)
	free(rc);
      
      progress_update(query_getfilepos());
    }
  progress_done();

  if (rc)
    free(rc);

  free(targetlist);
  free(hits);
  free(topscores);
  free(hitcount);

  count_kmers_exit();

  nw_exit();

  if (opt_dbmatched || opt_dbnotmatched)
    {
      for(long i=0; i<seqcount; i++)
	if (dbmatched[i])
	  {
	    if (opt_dbmatched)
	      db_fprint_fasta_with_size(fp_dbmatched, i, dbmatched[i]);
	  }
	else
	  {
	    if (opt_dbnotmatched)
	      db_fprint_fasta(fp_dbnotmatched, i);
	  }
    }



#if 1
  fprintf(stderr, "Matching query sequences: %ld (%.1f%%)\n", qmatches, 100.0 * qmatches / (query_no+1));
  fprintf(stderr, "Alignments computed: %ld\n", nwalignments);
#endif

  query_close();
  dbindex_free();
  db_free();

  /* close output files */

  if (opt_matched)
    fclose(fp_matched);

  if (opt_notmatched)
    fclose(fp_notmatched);

  if (opt_dbmatched)
    fclose(fp_dbmatched);

  if (opt_dbnotmatched)
    fclose(fp_dbnotmatched);

  if (opt_fastapairs)
    fclose(fp_fastapairs);

  if (fp_uc)
    fclose(fp_uc);

  if (fp_blast6out)
    fclose(fp_blast6out);

  if (fp_userout)
    fclose(fp_userout);

  if (fp_alnout)
    fclose(fp_alnout);

  show_rusage();
}


