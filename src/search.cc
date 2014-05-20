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
};
      
struct hit
{
  long target;
  long count;
  long nwscore;
  long nwdiff;
  long nwgaps;
  long nwindels;
  long nwalignmentlength;
  double nwid;
  char * nwalignment;
};

int hit_compare(const void * a, const void * b)
{
  struct hit * x = (struct hit *) a;
  struct hit * y = (struct hit *) b;

  // high id, then low id
  // early target, then late target

  if (x->nwid > y->nwid)
    return -1;
  else if (x->nwid < y->nwid)
    return +1;
  else
    if (x->target < y->target)
      return -1;
    else if (x->target > y->target)
      return +1;
    else
      return 0;
}

void search_showresults(struct hit * hits, int accepts,
                        char * query_head,
                        char * qsequence, long qseqlen)
{
  if (accepts > 0)
    {
      fprintf(alnoutfile,"Query >%s\n", query_head);
      fprintf(alnoutfile," %%Id   TLen  Target\n");
      
      for(int t = 0; t < accepts; t++)
        fprintf(alnoutfile,"%3.0f%% %6lu  %s\n",
                hits[t].nwid,
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
          
#if 0
          fprintf(alnoutfile,"\nCIGAR: %s\n", nwalignment);
#endif

          showalign(alnoutfile,
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
          
          free(hits[t].nwalignment);
        }
    }
}

long scorematrix[32][32];
long tophits;

unsigned char * hitcount = 0;
struct topscore * topscores = 0;
struct hit * hits = 0;

#define MAXSAMPLES 255

unsigned int hit_distr[MAXSAMPLES+1];

unsigned int * targetlist = 0;
unsigned int targetcount = 0;

unsigned int kmersample[256];
unsigned int kmersamplecount;

void search_get_query_samples(char * qsequence, unsigned int qseqlen, 
                              unsigned int wordlength, unsigned int samples)
{
  count_kmers(wordlength, qsequence, qseqlen);
  unsigned int unique = count_kmers_unique();
  kmersamplecount = 0;

  unsigned int pos = 0;
  unsigned int u = 0;

  unsigned int kmer = 0;
  unsigned int mask = (1<<(2*wordlength)) - 1;
  char * s = qsequence;
  char * e1 = s + wordlength - 1;
  char * e2 = s + qseqlen;
  if (e2 < e1)
    e1 = e2;

  while (s < e1)
    {
      kmer <<= 2;
      kmer |= *s++;
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
      kmer |= *s++;
      kmer &= mask;

      if (count_kmers_getcount(wordlength, kmer) == 1)
        {
          if (z>=0)
            {
              kmersample[kmersamplecount++] = kmer;
#if 0
              fprintf(stderr, "Query kmer sample %d at pos %d u %d: ", kmersamplecount, pos, u);
              fprint_kmer(stderr, wordlength, kmer);
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


unsigned int search_topscores(unsigned int samples,
                              unsigned int seqcount, 
                              unsigned int tophits)
{
  /*
    Count the kmer hits in each database sequence and
    make a sorted list of a given number (tophits)
    of the database sequences with the highest number of matching kmers.

    These are stored in the topscores array.
    
    The number to keep: tophits = min(seqcount, maxaccepts + maxrejects)

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

  if (hitspertarget < 0.5)
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

  
#if 0
  /* Start Experimental code */

  /* start by finding the distribution of kmer hits */

  for(unsigned int i=0; i <= samples; i++)
    hit_distr[i] = 0;

  if (sparse)
    {
      hit_distr[0] = seqcount - targetcount;
      for(unsigned int i=0; i<targetcount; i++)
        hit_distr[hitcount[targetlist[i]]]++;
    }
  else
    {
      for(unsigned int i=0; i < seqcount; i++)
        hit_distr[hitcount[i]]++;
    }

#if 0
  printf("Kmer count distribution:\n");
  for (unsigned int i=0; i <= samples; i++)
    printf("%u: %u\n", i, hit_distr[i]);
#endif

  /* End Experimental alternative code */
#endif
  

  if (sparse)
    {
      for(unsigned int z=0; z < targetcount; z++)
        {
          unsigned int i = targetlist[z];

          unsigned int count = hitcount[i];
          
          
          /* find insertion point */
          
          unsigned int p = topcount;
          
          while ((p > 0) && (count > topscores[p-1].count))
            p--;
          
          
          /* p = index in array where new data should be placed */
          
          if (p < tophits)
            {
              
              /* find new bottom of list */
              
              int bottom = topcount;
              if (topcount == tophits)
                bottom--;
              
              
              /* shift lower counts down */
              
              for(unsigned int j = bottom; j > p; j--)
                topscores[j] = topscores[j-1];
              
              
              /* insert or overwrite */
              
              topscores[p].count = count;
              topscores[p].seqno = i;
              
              if (topcount < tophits)
                topcount++;
            }
        }
    }
  else
    {
      for(unsigned int i=0; i < seqcount; i++)
        {
          unsigned int count = hitcount[i];
          
          
          /* find insertion point */
          
          unsigned int p = topcount;
          
          while ((p > 0) && (count > topscores[p-1].count))
            p--;
          
          
          /* p = index in array where new data should be placed */
          
          if (p < tophits)
            {
              
              /* find new bottom of list */
              
              int bottom = topcount;
              if (topcount == tophits)
                bottom--;
              
              
              /* shift lower counts down */
              
              for(unsigned int j = bottom; j > p; j--)
                topscores[j] = topscores[j-1];
              
              
              /* insert or overwrite */
              
              topscores[p].count = count;
              topscores[p].seqno = i;
              
              if (topcount < tophits)
                topcount++;
            }
        }
    }
  
  return topcount;
}


int search_onequery(char * query_head, long query_head_len,
                    char * qsequence, long qseqlen, long query_no)
{
  int seqcount = db_getsequencecount();

  /* compute necessary number of samples */

  unsigned int minsamples = 1;
  unsigned int maxsamples = MIN(qseqlen / 2, MAXSAMPLES);
  double default_samples = 24.0; // 18.0;
  unsigned int samples;

  if (identity > 0.0)
    samples = ceil(default_samples * exp(- wordlength * log(identity)));
  else
    samples = maxsamples;

  if (samples < minsamples)
    samples = minsamples;
  if (samples > maxsamples)
    samples = maxsamples;

  // printf("seqlen: %ld samples: %d\n", qseqlen, samples);


  /* extract unique kmer samples from query*/

  search_get_query_samples(qsequence, qseqlen, wordlength, samples);
  
  
  /* find database sequences with the most kmer hits */
  
  unsigned int topcount = search_topscores(samples, seqcount, tophits);

  
  /* analyze targets with the highest number of kmer hits */
  
  unsigned int t = 0;
  int accepts = 0;
  int rejects = 0;
  

#define FILTER 1
#ifdef FILTER
  /* compute query kmer vector */
  unsigned char query_kmervector[KMERVECTORBYTES];
  findkmers((unsigned char *)qsequence, qseqlen, query_kmervector);
#endif


  while((accepts < maxaccepts) && (rejects <= maxrejects) && (t<topcount))
    {
      unsigned int target = topscores[t].seqno;
      unsigned int count = topscores[t].count;
      char * dseq;
      long dseqlen;
      
      db_getsequenceandlength(target, & dseq, & dseqlen);
      
#ifdef FILTER
      /* perform kmer vector similarity test */

      unsigned char target_kmervector[KMERVECTORBYTES];
      findkmers((unsigned char*) dseq, qseqlen, target_kmervector);
      unsigned long diffkmers = comparekmervectors(query_kmervector,
                                                   target_kmervector);
      unsigned long mindiff = (diffkmers + 2*KMERLENGTH - 1)/(2*KMERLENGTH);
      double maxid = double(qseqlen - mindiff) / double(qseqlen);
      
      if (maxid < identity)
        rejects++;
      else
#endif
        {
          /* compute global alignment */
          
          unsigned long nwscore;
          unsigned long nwdiff;
          unsigned long nwgaps;
          unsigned long nwindels;
          unsigned long nwalignmentlength;
          char * nwalignment;
          
          nw_align(dseq,
                   dseq + dseqlen,
                   qsequence,
                   qsequence + qseqlen,
                   (long*) scorematrix,
                   gapopen_cost,
                   gapextend_cost,
                   & nwscore,
                   & nwdiff,
                   & nwgaps,
                   & nwindels,
                   & nwalignmentlength,
                   & nwalignment,
                   query_no,
                   target);
          
          double nwid = (nwalignmentlength - nwdiff) * 100.0 / nwalignmentlength;
          
          if (nwid >= 100.0 * identity)
            {
              hits[accepts].target = target;
              hits[accepts].count = count;
              hits[accepts].nwscore = nwscore;
              hits[accepts].nwdiff = nwdiff;
              hits[accepts].nwgaps = nwgaps;
              hits[accepts].nwindels = nwindels;
              hits[accepts].nwalignmentlength = nwalignmentlength;
              hits[accepts].nwalignment = nwalignment;
              hits[accepts].nwid = nwid;
              accepts++;
            }
          else
            {
              free(nwalignment);
              rejects++;
            }
        }
      t++;
    }

  
  /* sort accepted targets */
  
  qsort(hits, accepts, sizeof(struct hit), hit_compare);
  
  
  return accepts;

}

void search()
{
  fprintf(stderr,"Searching\n");

  for(int i=0; i<32; i++)
    for(int j=0; j<32; j++)
      if (i==j)
        scorematrix[i][j] = 0;
      else
        scorematrix[i][j] = mismatch_cost;
  
  int seqcount = db_getsequencecount();
  tophits = maxrejects + maxaccepts;
  if (tophits > seqcount)
    tophits = seqcount;
  if (maxaccepts > seqcount)
    maxaccepts = seqcount;

  if (identity == 1.0)
    maxrejects = 0;
  
  nw_init();
  count_kmers_init();

  hitcount = (unsigned char *) xmalloc(seqcount);
  topscores = (struct topscore *) xmalloc(sizeof(struct topscore) * tophits);
  hits = (struct hit *) xmalloc(sizeof(struct hit) * maxaccepts);
  targetlist = (unsigned int*) xmalloc(sizeof(unsigned int)*seqcount);

  while(1)
    {
      char * query_head;
      long query_head_len;
      char * qsequence;
      long qseqlen;
      long query_no;
      
      if (query_getnext(& query_head, & query_head_len,
                        & qsequence, & qseqlen,
                        & query_no))
        {

          /* perform search */

          int accepts = search_onequery(query_head, query_head_len, 
                                        qsequence, qseqlen, 
                                        query_no);
          
          /* show results */
          
          search_showresults(hits, accepts, query_head, qsequence, qseqlen);

        }
      else
        break;
    }
  
  free(targetlist);
  free(hits);
  free(topscores);
  free(hitcount);

  count_kmers_exit();
  nw_exit();

  show_rusage();
}


