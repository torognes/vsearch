/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.


  GNU General Public License version 3

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The BSD 2-Clause License

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "vsearch.h"

static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */

typedef struct clusterinfo_s
{
  int seqno;
  int clusterno;
  char * cigar;
  int strand;
} clusterinfo_t;

static clusterinfo_t * clusterinfo = 0;
static int clusters = 0;

static int count_matched = 0;
static int count_notmatched = 0;

static int64_t * cluster_abundance;

static FILE * fp_centroids = 0;
static FILE * fp_uc = 0;
static FILE * fp_alnout = 0;
static FILE * fp_samout = 0;
static FILE * fp_userout = 0;
static FILE * fp_blast6out = 0;
static FILE * fp_fastapairs = 0;
static FILE * fp_matched = 0;
static FILE * fp_notmatched = 0;
static FILE * fp_otutabout = 0;
static FILE * fp_mothur_shared_out = 0;
static FILE * fp_biomout = 0;

static pthread_attr_t attr;

static struct searchinfo_s * si_plus;
static struct searchinfo_s * si_minus;

typedef struct thread_info_s
{
  pthread_t thread;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  int work;
  int query_first;
  int query_count;
} thread_info_t;

static thread_info_t * ti;

inline int compare_byclusterno(const void * a, const void * b)
{
  clusterinfo_t * x = (clusterinfo_t *) a;
  clusterinfo_t * y = (clusterinfo_t *) b;
  if (x->clusterno < y->clusterno)
    return -1;
  else if (x->clusterno > y->clusterno)
    return +1;
  else if (x->seqno < y->seqno)
    return -1;
  else if (x->seqno > y->seqno)
    return +1;
  else
    return 0;
}

inline int compare_byclusterabundance(const void * a, const void * b)
{
  clusterinfo_t * x = (clusterinfo_t *) a;
  clusterinfo_t * y = (clusterinfo_t *) b;
  if (cluster_abundance[x->clusterno] > cluster_abundance[y->clusterno])
    return -1;
  else if (cluster_abundance[x->clusterno] < cluster_abundance[y->clusterno])
    return +1;
  else if (x->clusterno < y->clusterno)
    return -1;
  else if (x->clusterno > y->clusterno)
    return +1;
  else if (x->seqno < y->seqno)
    return -1;
  else if (x->seqno > y->seqno)
    return +1;
  else
    return 0;
}


inline void cluster_query_core(struct searchinfo_s * si)
{
  /* the main core function for clustering */
  
  /* get sequence etc */
  int seqno = si->query_no;
  si->query_head_len = db_getheaderlen(seqno);
  si->query_head = db_getheader(seqno);
  si->qsize = db_getabundance(seqno);
  si->qseqlen = db_getsequencelen(seqno);
  if (si->strand)
    reverse_complement(si->qsequence, db_getsequence(seqno), si->qseqlen);
  else
    strcpy(si->qsequence, db_getsequence(seqno));
  
  /* perform search */
  search_onequery(si, opt_qmask);
}

inline void cluster_worker(int64_t t)
{
  /* wrapper for the main threaded core function for clustering */
  for (int q = 0; q < ti[t].query_count; q++)
    {
      cluster_query_core(si_plus + ti[t].query_first + q);
      if (opt_strand>1)
        cluster_query_core(si_minus + ti[t].query_first + q);
    }
}

void * threads_worker(void * vp)
{
  int64_t t = (int64_t) vp;
  thread_info_s * tip = ti + t;
  pthread_mutex_lock(&tip->mutex);
  /* loop until signalled to quit */
  while (tip->work >= 0)
    {
      /* wait for work available */
      if (tip->work == 0)
        pthread_cond_wait(&tip->cond, &tip->mutex);
      if (tip->work > 0)
        {
          cluster_worker(t);
          tip->work = 0;
          pthread_cond_signal(&tip->cond);
        }
    }
  pthread_mutex_unlock(&tip->mutex);
  return 0;
}

void threads_wakeup(int queries)
{
  int threads = queries > opt_threads ? opt_threads : queries;
  int queries_rest = queries;
  int threads_rest = threads;
  int query_next = 0;

  /* tell the threads that there is work to do */
  for(int t=0; t < threads; t++)
    {
      thread_info_t * tip = ti + t;

      tip->query_first = query_next;
      tip->query_count = (queries_rest + threads_rest - 1) / threads_rest;
      queries_rest -= tip->query_count;
      query_next += tip->query_count;
      threads_rest--;

      pthread_mutex_lock(&tip->mutex);
      tip->work = 1;
      pthread_cond_signal(&tip->cond);
      pthread_mutex_unlock(&tip->mutex);
    }
  
  /* wait for theads to finish their work */
  for(int t=0; t < threads; t++)
    {
      thread_info_t * tip = ti + t;
      pthread_mutex_lock(&tip->mutex);
      while (tip->work > 0)
        pthread_cond_wait(&tip->cond, &tip->mutex);
      pthread_mutex_unlock(&tip->mutex);
    }
}

void threads_init()
{
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  /* allocate memory for thread info */
  ti = (thread_info_t *) xmalloc(opt_threads * sizeof(thread_info_t));
  
  /* init and create worker threads */
  for(int t=0; t < opt_threads; t++)
    {
      thread_info_t * tip = ti + t;
      tip->work = 0;
      pthread_mutex_init(&tip->mutex, 0);
      pthread_cond_init(&tip->cond, 0);
      if (pthread_create(&tip->thread, &attr, threads_worker, (void*)(int64_t)t))
        fatal("Cannot create thread");
    }
}

void threads_exit()
{
  /* finish and clean up worker threads */
  for(int t=0; t<opt_threads; t++)
    {
      struct thread_info_s * tip = ti + t;
      
      /* tell worker to quit */
      pthread_mutex_lock(&tip->mutex);
      tip->work = -1;
      pthread_cond_signal(&tip->cond);
      pthread_mutex_unlock(&tip->mutex);

      /* wait for worker to quit */
      if (pthread_join(tip->thread, 0))
        fatal("Cannot join thread");

      pthread_cond_destroy(&tip->cond);
      pthread_mutex_destroy(&tip->mutex);
    }
  xfree(ti);
  pthread_attr_destroy(&attr);
}

void cluster_query_init(struct searchinfo_s * si)
{
  /* initialisation of data for one thread; run once for each thread */
  /* thread specific initialiation */

  si->qsize = 1;
  si->nw = 0;
  si->hit_count = 0;

  /* allocate memory for sequence */

  si->seq_alloc = db_getlongestsequence() + 1;
  si->qsequence = (char *) xmalloc(si->seq_alloc);

  si->kmers = (count_t *) xmalloc(seqcount * sizeof(count_t) + 32);
  si->hits = (struct hit *) xmalloc(sizeof(struct hit) * tophits);

  si->uh = unique_init();
  si->m = minheap_init(tophits);
  si->s = search16_init(opt_match,
                        opt_mismatch,
                        opt_gap_open_query_left,
                        opt_gap_open_target_left,
                        opt_gap_open_query_interior,
                        opt_gap_open_target_interior,
                        opt_gap_open_query_right,
                        opt_gap_open_target_right,
                        opt_gap_extension_query_left,
                        opt_gap_extension_target_left,
                        opt_gap_extension_query_interior,
                        opt_gap_extension_target_interior,
                        opt_gap_extension_query_right,
                        opt_gap_extension_target_right);
  si->nw = nw_init();
}

void cluster_query_exit(struct searchinfo_s * si)
{
  /* clean up after thread execution; called once per thread */

  search16_exit(si->s);
  unique_exit(si->uh);
  minheap_exit(si->m);
  nw_exit(si->nw);
  
  if (si->qsequence)
    xfree(si->qsequence);
  if (si->hits)
    xfree(si->hits);
  if (si->kmers)
    xfree(si->kmers);
}

char * relabel_otu(int clusterno, char * sequence, int seqlen)
{
  char * label = 0;
  if (opt_relabel)
    {
      label = (char*) xmalloc(strlen(opt_relabel) + 21);
      sprintf(label, "%s%d", opt_relabel, clusterno+1);
    }
  else if (opt_relabel_sha1)
    {
      label = (char*) xmalloc(LEN_HEX_DIG_SHA1);
      get_hex_seq_digest_sha1(label, sequence, seqlen);
    }
  else if (opt_relabel_md5)
    {
      label = (char*) xmalloc(LEN_HEX_DIG_MD5);
      get_hex_seq_digest_md5(label, sequence, seqlen);
    }
  return label;
}

void cluster_core_results_hit(struct hit * best,
                              int clusterno,
                              char * query_head,
                              int qseqlen,
                              char * qsequence,
                              char * qsequence_rc,
                              int qsize)
{
  count_matched++;

  if (opt_otutabout || opt_mothur_shared_out || opt_biomout)
    {
      if (opt_relabel || opt_relabel_sha1 || opt_relabel_md5)
        {
          char * label = relabel_otu(clusterno,
                                     db_getsequence(best->target),
                                     db_getsequencelen(best->target));
          otutable_add(query_head, label, qsize);
          xfree(label);
        }
      else
        otutable_add(query_head,
                     db_getheader(best->target),
                     qsize);
    }

  if (fp_uc)
    results_show_uc_one(fp_uc,
                        best, query_head,
                        qsequence, qseqlen, qsequence_rc,
                        clusterno);
  
  if (fp_alnout)
    results_show_alnout(fp_alnout,
                        best, 1, query_head,
                        qsequence, qseqlen, qsequence_rc);
  
  if (fp_samout)
    results_show_samout(fp_samout,
                        best, 1, query_head,
                        qsequence, qseqlen, qsequence_rc);
  
  if (fp_fastapairs)
    results_show_fastapairs_one(fp_fastapairs,
                                best,
                                query_head,
                                qsequence,
                                qseqlen,
                                qsequence_rc);
  
  if (fp_userout)
    results_show_userout_one(fp_userout, best, query_head, 
                             qsequence, qseqlen, qsequence_rc);
  
  if (fp_blast6out)
    results_show_blast6out_one(fp_blast6out, best, query_head,
                               qsequence, qseqlen, qsequence_rc);
  
  if (opt_matched)
    fasta_print_general(fp_matched,
                        0,
                        qsequence,
                        qseqlen,
                        query_head,
                        strlen(query_head),
                        qsize,
                        count_matched,
                        -1, -1, 0, 0.0);
}

void cluster_core_results_nohit(int clusterno,
                                char * query_head,
                                int qseqlen,
                                char * qsequence,
                                char * qsequence_rc,
                                int qsize)
{
  count_notmatched++;

  if (opt_otutabout || opt_mothur_shared_out || opt_biomout)
    {
      if (opt_relabel || opt_relabel_sha1 || opt_relabel_md5)
        {
          char * label = relabel_otu(clusterno, qsequence, qseqlen);
          otutable_add(query_head, label, qsize);
          xfree(label);
        }
      else
        otutable_add(query_head, query_head, qsize);
    }

  if (opt_uc)
    {
      fprintf(fp_uc, "S\t%d\t%d\t*\t*\t*\t*\t*\t%s\t*\n",
              clusters, qseqlen, query_head);
    }
  
  if (opt_output_no_hits)
    {
      if (fp_userout)
        results_show_userout_one(fp_userout, 0, query_head, 
                                 qsequence, qseqlen, qsequence_rc);
      
      if (fp_blast6out)
        results_show_blast6out_one(fp_blast6out, 0, query_head,
                                   qsequence, qseqlen, qsequence_rc);
    }
  
  if (opt_notmatched)
    fasta_print_general(fp_notmatched,
                        0,
                        qsequence,
                        qseqlen,
                        query_head,
                        strlen(query_head),
                        qsize,
                        count_notmatched,
                        -1, -1, 0, 0.0);
}

int compare_kmersample(const void * a, const void * b)
{
  unsigned int x = * (unsigned int *) a;
  unsigned int y = * (unsigned int *) b;
  
  if (x < y)
    return -1;
  else if (x > y)
    return +1;
  else
    return 0;
}

void cluster_core_parallel()
{
  /* create threads and set them in stand-by mode */
  threads_init();

  const int queries_per_thread = 1;
  int max_queries = queries_per_thread * opt_threads;

  /* allocate memory for the search information for each query;
     and initialize it */
  si_plus  = (struct searchinfo_s *) xmalloc(max_queries *
                                             sizeof(struct searchinfo_s));
  if (opt_strand>1)
    si_minus = (struct searchinfo_s *) xmalloc(max_queries *
                                               sizeof(struct searchinfo_s));
  for(int i = 0; i < max_queries; i++)
    {
      cluster_query_init(si_plus+i);
      si_plus[i].strand = 0;
      if (opt_strand > 1)
        {
          cluster_query_init(si_minus+i);
          si_minus[i].strand = 1;
        }
    }

  int * extra_list = (int*) xmalloc(max_queries*sizeof(int));

  LinearMemoryAligner lma;
  int64_t * scorematrix = lma.scorematrix_create(opt_match, opt_mismatch);
  lma.set_parameters(scorematrix,
                     opt_gap_open_query_left,
                     opt_gap_open_target_left,
                     opt_gap_open_query_interior,
                     opt_gap_open_target_interior,
                     opt_gap_open_query_right,
                     opt_gap_open_target_right,
                     opt_gap_extension_query_left,
                     opt_gap_extension_target_left,
                     opt_gap_extension_query_interior,
                     opt_gap_extension_target_interior,
                     opt_gap_extension_query_right,
                     opt_gap_extension_target_right);
  
  int aligncount = 0;

  int lastlength = INT_MAX;

  int seqno = 0;

  int64_t sum_nucleotides = 0;

  progress_init("Clustering", db_getnucleotidecount());

  while(seqno < seqcount)
    {
      /* prepare work for the threads in sia[i] */
      /* read query sequences into the search info (si) for each thread */

      int queries = 0;
      
      for(int i = 0; i < max_queries; i++)
        {
          if (seqno < seqcount)
            {
              int length = db_getsequencelen(seqno);

#if 1
              if (opt_cluster_smallmem && (!opt_usersort) && (length > lastlength))
                fatal("Sequences not sorted by length and --usersort not specified.");
#endif

              lastlength = length;

              si_plus[i].query_no = seqno;
              si_plus[i].strand = 0;

              if (opt_strand > 1)
                {
                  si_minus[i].query_no = seqno;
                  si_minus[i].strand = 1;
                }

              queries++;
              seqno++;
            }
        }

      /* perform work in threads */
      threads_wakeup(queries);
      
      /* analyse results */
      int extra_count = 0;

      for(int i=0; i < queries; i++)
        {
          struct searchinfo_s * si_p = si_plus + i;
          struct searchinfo_s * si_m = opt_strand > 1 ? si_minus + i : 0;
          
          for(int s = 0; s < opt_strand; s++)
            {
              struct searchinfo_s * si = s ? si_m : si_p;

              int added = 0;

              if (extra_count)
                {
                  /* Check if there is a hit with one of the non-matching
                     extra sequences just analysed in this round */
                  
                  for (int j=0; j<extra_count; j++)
                    {
                      struct searchinfo_s * sic = si_plus + extra_list[j];
                      
                      /* find the number of shared unique kmers */
                      unsigned int shared
                        = unique_count_shared(si->uh,
                                              opt_wordlength,
                                              sic->kmersamplecount,
                                              sic->kmersample);

                      /* check if min number of shared kmers is satisfied */
                      if (search_enough_kmers(si, shared))
                        {
                          unsigned int length = sic->qseqlen;
                          
                          /* Go through the list of hits and see if the current
                             match is better than any on the list in terms of
                             more shared kmers (or shorter length if equal
                             no of kmers). Determine insertion point (x). */
                          
                          int x = si->hit_count;
                          while ((x > 0) &&
                                 ((si->hits[x-1].count < shared) ||
                                  ((si->hits[x-1].count == shared) &&
                                   (db_getsequencelen(si->hits[x-1].target)
                                    > length))))
                            x--;
                          
                          if (x < opt_maxaccepts + opt_maxrejects - 1)
                            {
                              /* insert into list at position x */
                              
                              /* trash bottom element if no more space */
                              if (si->hit_count >= opt_maxaccepts + opt_maxrejects - 1)
                                {
                                  if (si->hits[si->hit_count-1].aligned)
                                    xfree(si->hits[si->hit_count-1].nwalignment);
                                  si->hit_count--;
                                }
                              
                              /* move the rest down */
                              for(int z = si->hit_count; z > x; z--)
                                si->hits[z] = si->hits[z-1];
                              
                              /* init new hit */
                              struct hit * hit = si->hits + x;
                              si->hit_count++;
                              
                              hit->target = sic->query_no;
                              hit->strand = si->strand;
                              hit->count = shared;
                              hit->accepted = 0;
                              hit->rejected = 0;
                              hit->aligned = 0;
                              hit->weak = 0;
                              hit->nwalignment = 0;
                              
                              added++;
                            }
                        }
                    }
                }
              
              /* now go through the hits and determine final status of each */

              if (added)
                {
                  si->rejects = 0;
                  si->accepts = 0;
                  
                  /* set all statuses to undetermined */
                  
                  for(int t=0; t< si->hit_count; t++)
                    {
                      si->hits[t].accepted = 0;
                      si->hits[t].rejected = 0;
                    }
                  
                  for(int t = 0;
                      (si->accepts < opt_maxaccepts) && 
                        (si->rejects < opt_maxrejects) &&
                        (t < si->hit_count);
                      t++)
                    {
                      struct hit * hit = si->hits + t;

                      if (! hit->aligned)
                        {
                          /* Test accept/reject criteria before alignment */
                          unsigned int target = hit->target;
                          if (search_acceptable_unaligned(si, target))
                            {
                              aligncount++;
                              
                              /* perform vectorized alignment */
                              /* but only using 1 sequence ! */

                              unsigned int nwtarget = target;
                              
                              int64_t nwscore;
                              int64_t nwalignmentlength;
                              int64_t nwmatches;
                              int64_t nwmismatches;
                              int64_t nwgaps;
                              char * nwcigar = 0;

                              /* short variants for simd aligner */
                              CELL snwscore;
                              unsigned short snwalignmentlength;
                              unsigned short snwmatches;
                              unsigned short snwmismatches;
                              unsigned short snwgaps;
                              
                              search16(si->s,
                                       1,
                                       & nwtarget,
                                       & snwscore,
                                       & snwalignmentlength,
                                       & snwmatches,
                                       & snwmismatches,
                                       & snwgaps,
                                       & nwcigar);
                              
                              int64_t tseqlen = db_getsequencelen(target);

                              if (snwscore == SHRT_MAX)
                                {
                                  /* In case the SIMD aligner cannot align,
                                     perform a new alignment with the
                                     linear memory aligner */
                                  
                                  char * tseq = db_getsequence(target);
                                  
                                  if (nwcigar)
                                    xfree(nwcigar);
                                  
                                  nwcigar = xstrdup(lma.align(si->qsequence,
                                                             tseq,
                                                             si->qseqlen,
                                                             tseqlen));

                                  lma.alignstats(nwcigar,
                                                 si->qsequence,
                                                 tseq,
                                                 & nwscore,
                                                 & nwalignmentlength,
                                                 & nwmatches,
                                                 & nwmismatches,
                                                 & nwgaps);
                                }
                              else
                                {
                                  nwscore = snwscore;
                                  nwalignmentlength = snwalignmentlength;
                                  nwmatches = snwmatches;
                                  nwmismatches = snwmismatches;
                                  nwgaps = snwgaps;
                                }
                              

                              int64_t nwdiff = nwalignmentlength - nwmatches;
                              int64_t nwindels = nwdiff - nwmismatches;
                              
                              hit->aligned = 1;
                              hit->nwalignment = nwcigar;
                              hit->nwscore = nwscore;
                              hit->nwdiff = nwdiff;
                              hit->nwgaps = nwgaps;
                              hit->nwindels = nwindels;
                              hit->nwalignmentlength = nwalignmentlength;
                              hit->matches = nwmatches;
                              hit->mismatches = nwmismatches;

                              hit->nwid = 100.0 *
                                (nwalignmentlength - hit->nwdiff) /
                                nwalignmentlength;
                              
                              hit->shortest = MIN(si->qseqlen, tseqlen);
                              hit->longest = MAX(si->qseqlen, tseqlen);
                              
                              /* trim alignment and compute numbers
                                 excluding terminal gaps */
                              align_trim(hit);
                            }
                          else
                            {
                              /* rejection without alignment */
                              hit->rejected = 1;
                              si->rejects++;
                            }
                        }
                          
                      if (! hit->rejected)
                        {
                          /* test accept/reject criteria after alignment */
                          if (search_acceptable_aligned(si, hit))
                            si->accepts++;
                          else
                            si->rejects++;
                        }
                    }

                  /* delete all undetermined hits */
                  
                  int new_hit_count = si->hit_count;
                  for(int t=si->hit_count-1; t>=0; t--)
                    {
                      struct hit * hit = si->hits + t;
                      if (!hit->accepted && !hit->rejected)
                        {
                          new_hit_count = t;
                          if (hit->aligned)
                            xfree(hit->nwalignment);
                        }
                    }
                  si->hit_count = new_hit_count;
                }
            }

          /* find best hit */
          struct hit * best = 0;
          if (opt_sizeorder)
            best = search_findbest2_bysize(si_p, si_m);
          else
            best = search_findbest2_byid(si_p, si_m);
            
          int myseqno = si_p->query_no;
          
          if (best)
            {
              /* a hit was found, cluster current sequence with hit */
              int target = best->target;
              
              /* output intermediate results to uc etc */
              cluster_core_results_hit(best,
                                       clusterinfo[target].clusterno,
                                       si_p->query_head,
                                       si_p->qseqlen,
                                       si_p->qsequence,
                                       best->strand ? si_m->qsequence : 0,
                                       si_p->qsize);

              /* update cluster info about this sequence */
              clusterinfo[myseqno].seqno = myseqno;
              clusterinfo[myseqno].clusterno = clusterinfo[target].clusterno;
              clusterinfo[myseqno].cigar = best->nwalignment;
              clusterinfo[myseqno].strand = best->strand;
              best->nwalignment = 0;
            }
          else
            {
              /* no hit found; add it to the list of extra sequences
                 that must be considered by the coming queries in this 
                 round */
              extra_list[extra_count++] = i;
              
              /* update cluster info about this sequence */
              clusterinfo[myseqno].seqno = myseqno;
              clusterinfo[myseqno].clusterno = clusters;
              clusterinfo[myseqno].cigar = 0;
              clusterinfo[myseqno].strand = 0;
              
              /* add current sequence to database */
              dbindex_addsequence(myseqno, opt_qmask);
              
              /* output intermediate results to uc etc */
              cluster_core_results_nohit(clusters,
                                         si_p->query_head,
                                         si_p->qseqlen,
                                         si_p->qsequence,
                                         0,
                                         si_p->qsize);
              clusters++;
            }
          
          /* free alignments */
          for (int s = 0; s < opt_strand; s++)
            {
              struct searchinfo_s * si = s ? si_m : si_p;
              for(int j=0; j<si->hit_count; j++)
                if (si->hits[j].aligned)
                  if (si->hits[j].nwalignment)
                    xfree(si->hits[j].nwalignment);
            }

          sum_nucleotides += si_p->qseqlen;
        }
      
      progress_update(sum_nucleotides);
    }
  progress_done();

#if 0
  if (!opt_quiet)
    fprintf(stderr, "Extra alignments computed: %d\n", aligncount);
#endif

  /* clean up search info */
  for(int i = 0; i < max_queries; i++)
    {
      cluster_query_exit(si_plus+i);
      if (opt_strand > 1)
        cluster_query_exit(si_minus+i);
    }

  xfree(extra_list);

  xfree(si_plus);
  if (opt_strand>1)
    xfree(si_minus);

  /* terminate threads and clean up */
  threads_exit();

  xfree(scorematrix);
}

void cluster_core_serial()
{
  struct searchinfo_s si_p[1];
  struct searchinfo_s si_m[1];

  cluster_query_init(si_p);
  if (opt_strand > 1)
    cluster_query_init(si_m);

  int lastlength = INT_MAX;

  progress_init("Clustering", seqcount);
  for (int seqno=0; seqno<seqcount; seqno++)
    {
      int length = db_getsequencelen(seqno);

#if 1
      if (opt_cluster_smallmem && (!opt_usersort) && (length > lastlength))
        fatal("Sequences not sorted by length and --usersort not specified.");
#endif

      lastlength = length;

      si_p->query_no = seqno;
      si_p->strand = 0;
      cluster_query_core(si_p);

      if (opt_strand > 1)
        {
          si_m->query_no = seqno;
          si_m->strand = 1;
          cluster_query_core(si_m);
        }

      struct hit * best = 0;
      if (opt_sizeorder)
        best = search_findbest2_bysize(si_p, si_m);
      else
        best = search_findbest2_byid(si_p, si_m);
      
      if (best)
        {
          int target = best->target;
          cluster_core_results_hit(best,
                                   clusterinfo[target].clusterno,
                                   si_p->query_head,
                                   si_p->qseqlen,
                                   si_p->qsequence,
                                   best->strand ? si_m->qsequence : 0,
                                   si_p->qsize);
          clusterinfo[seqno].seqno = seqno;
          clusterinfo[seqno].clusterno = clusterinfo[target].clusterno;
          clusterinfo[seqno].cigar = best->nwalignment;
          clusterinfo[seqno].strand = best->strand;
          best->nwalignment = 0;
        }
      else
        {
          clusterinfo[seqno].seqno = seqno;
          clusterinfo[seqno].clusterno = clusters;
          clusterinfo[seqno].cigar = 0;
          clusterinfo[seqno].strand = 0;
          dbindex_addsequence(seqno, opt_qmask);
          cluster_core_results_nohit(clusters,
                                     si_p->query_head,
                                     si_p->qseqlen,
                                     si_p->qsequence,
                                     0,
                                     si_p->qsize);
          clusters++;
        }
      
      /* free alignments */
      for (int s = 0; s < opt_strand; s++)
        {
          struct searchinfo_s * si = s ? si_m : si_p;
          for(int i=0; i<si->hit_count; i++)
            if (si->hits[i].aligned)
              if (si->hits[i].nwalignment)
                xfree(si->hits[i].nwalignment);
        }

      progress_update(seqno);
    }
  progress_done();

  cluster_query_exit(si_p);
  if (opt_strand>1)
    cluster_query_exit(si_m);
}


void cluster(char * dbname, 
             char * cmdline,
             char * progheader)
{
  if (opt_centroids)
    {
      fp_centroids = fopen_output(opt_centroids);
      if (!fp_centroids)
        fatal("Unable to open centroids file for writing");
    }

  if (opt_uc)
    {
      fp_uc = fopen_output(opt_uc);
      if (!fp_uc)
        fatal("Unable to open uc file for writing");
    }

  if (opt_alnout)
    {
      fp_alnout = fopen_output(opt_alnout);
      if (! fp_alnout)
        fatal("Unable to open alignment output file for writing");

      fprintf(fp_alnout, "%s\n", cmdline);
      fprintf(fp_alnout, "%s\n", progheader);
    }

  if (opt_samout)
    {
      fp_samout = fopen_output(opt_samout);
      if (! fp_samout)
        fatal("Unable to open SAM output file for writing");
    }

  if (opt_userout)
    {
      fp_userout = fopen_output(opt_userout);
      if (! fp_userout)
        fatal("Unable to open user-defined output file for writing");
    }

  if (opt_blast6out)
    {
      fp_blast6out = fopen_output(opt_blast6out);
      if (! fp_blast6out)
        fatal("Unable to open blast6-like output file for writing");
    }

  if (opt_fastapairs)
    {
      fp_fastapairs = fopen_output(opt_fastapairs);
      if (! fp_fastapairs)
        fatal("Unable to open fastapairs output file for writing");
    }

  if (opt_matched)
    {
      fp_matched = fopen_output(opt_matched);
      if (! fp_matched)
        fatal("Unable to open matched output file for writing");
    }

  if (opt_notmatched)
    {
      fp_notmatched = fopen_output(opt_notmatched);
      if (! fp_notmatched)
        fatal("Unable to open notmatched output file for writing");
    }

  if (opt_otutabout)
    {
      fp_otutabout = fopen_output(opt_otutabout);
      if (! fp_otutabout)
        fatal("Unable to open OTU table (text format) output file for writing");
    }

  if (opt_mothur_shared_out)
    {
      fp_mothur_shared_out = fopen_output(opt_mothur_shared_out);
      if (! fp_mothur_shared_out)
        fatal("Unable to open OTU table (mothur format) output file for writing");
    }

  if (opt_biomout)
    {
      fp_biomout = fopen_output(opt_biomout);
      if (! fp_biomout)
        fatal("Unable to open OTU table (biom 1.0 format) output file for writing");
    }

  db_read(dbname, 0);

  otutable_init();

  results_show_samheader(fp_samout, cmdline, dbname);

  if (opt_qmask == MASK_DUST)
    dust_all();
  else if ((opt_qmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();
  
  show_rusage();
  
  seqcount = db_getsequencecount();
  
  if (opt_cluster_fast)
    db_sortbylength();
  else if (opt_cluster_size || opt_cluster_unoise)
    db_sortbyabundance();

  dbindex_prepare(1, opt_qmask);
  
  /* tophits = the maximum number of hits we need to store */

  if ((opt_maxrejects == 0) || (opt_maxrejects > seqcount))
    opt_maxrejects = seqcount;

  if ((opt_maxaccepts == 0) || (opt_maxaccepts > seqcount))
    opt_maxaccepts = seqcount;

  tophits = opt_maxrejects + opt_maxaccepts + MAXDELAYED;

  if (tophits > seqcount)
    tophits = seqcount;

  clusterinfo = (clusterinfo_t *) xmalloc(seqcount * sizeof(clusterinfo_t));

  if (opt_log)
    {
      uint64_t slots = 1ULL << (opt_wordlength << 1ULL);
      fprintf(fp_log, "\n");
      fprintf(fp_log, "      Alphabet  nt\n");
      fprintf(fp_log, "    Word width  %" PRId64 "\n", opt_wordlength);
      fprintf(fp_log, "     Word ones  %" PRId64 "\n", opt_wordlength);
      fprintf(fp_log, "        Spaced  No\n");
      fprintf(fp_log, "        Hashed  No\n");
      fprintf(fp_log, "         Coded  No\n");
      fprintf(fp_log, "       Stepped  No\n");
      fprintf(fp_log, "         Slots  %" PRIu64 " (%.1fk)\n", slots, slots/1000.0);
      fprintf(fp_log, "       DBAccel  100%%\n");
      fprintf(fp_log, "\n");
    }      

  if (opt_threads == 1)
    cluster_core_serial();
  else
    cluster_core_parallel();


  /* find size and abundance of each cluster and save stats */

  cluster_abundance = (int64_t *) xmalloc(clusters * sizeof(int64_t));
  int * cluster_size = (int *) xmalloc(clusters * sizeof(int));

  memset(cluster_abundance, 0, clusters * sizeof(int64_t));
  memset(cluster_size, 0, clusters * sizeof(int));

  for(int i=0; i<seqcount; i++)
    {
      int seqno = clusterinfo[i].seqno;
      int clusterno = clusterinfo[i].clusterno;
      cluster_abundance[clusterno] += opt_sizein ? db_getabundance(seqno) : 1;
      cluster_size[clusterno]++;
    }
      
  int64_t abundance_min = LONG_MAX;
  int64_t abundance_max = 0;
  int size_max = 0;
  int singletons = 0;

  for(int z=0; z<clusters; z++)
    {
      int64_t abundance = cluster_abundance[z];
      if (abundance < abundance_min)
        abundance_min = abundance;
      if (abundance > abundance_max)
        abundance_max = abundance;

      if (abundance == 1)
        singletons++;

      int size = cluster_size[z];
      if (size > size_max)
        size_max = size;
    }


  /* Sort clusters */
  /* Sequences in same cluster must always come right after each other. */
  /* The centroid sequence must be the first in each cluster. */

  progress_init("Sorting clusters", clusters);
  qsort(clusterinfo, seqcount, sizeof(clusterinfo_t), compare_byclusterno);
  progress_done();

  progress_init("Writing clusters", seqcount);

  /* allocate memory for full file name of the clusters files */
  FILE * fp_clusters = 0;
  char * fn_clusters = 0;
  if (opt_clusters)
    fn_clusters = (char *) xmalloc(strlen(opt_clusters) + 25);

  int lastcluster = -1;
  int ordinal = 0;

  for(int i=0; i<seqcount; i++)
    {
      int seqno = clusterinfo[i].seqno;
      int clusterno = clusterinfo[i].clusterno;

      if (clusterno != lastcluster)
        {
          /* prepare for new cluster */
          /* performed with first sequence only in each cluster */
          /* the first sequence is always the centroid */

          if (opt_centroids)
            fasta_print_general(fp_centroids,
                                0,
                                db_getsequence(seqno),
                                db_getsequencelen(seqno),
                                db_getheader(seqno),
                                db_getheaderlen(seqno),
                                cluster_abundance[clusterno],
                                clusterno+1,
                                -1, -1, 0, 0.0);

          if (opt_uc)
            fprintf(fp_uc, "C\t%d\t%" PRId64 "\t*\t*\t*\t*\t*\t%s\t*\n",
                    clusterno,
                    cluster_abundance[clusterno],
                    db_getheader(seqno));
          
          if (opt_clusters)
            {
              /* close previous (except for first time) and open new file */
              if (lastcluster != -1)
                fclose(fp_clusters);
              
              ordinal = 0;
              sprintf(fn_clusters, "%s%d", opt_clusters, clusterno);
              fp_clusters = fopen_output(fn_clusters);
              if (!fp_clusters)
                fatal("Unable to open clusters file for writing");
            }
          
          lastcluster = clusterno;
        }
      
      /* performed for all sequences */

      if (opt_clusters)
        {
          ordinal++;
          fasta_print_db_relabel(fp_clusters, seqno, ordinal);
        }
      
      progress_update(i);
    }

  if (lastcluster != -1)
    {
      /* performed with the last sequence */
      if (opt_clusters)
        {
          fclose(fp_clusters);
          if (fn_clusters)
            xfree(fn_clusters);
        }
    }

  progress_done();

  if (clusters < 1)
    {
      if (!opt_quiet)
        {
          fprintf(stderr, "Clusters: 0\n");
          fprintf(stderr, "Singletons: 0\n");
        }
      if (opt_log)
        {
          fprintf(fp_log, "Clusters: 0\n");
          fprintf(fp_log, "Singletons: 0\n");
        }
    }
  else
    {
      if (!opt_quiet)
        {
          fprintf(stderr,
                  "Clusters: %d Size min %" PRId64 ", max %" PRId64 ", avg %.1f\n",
                  clusters,
                  abundance_min,
                  abundance_max,
                  1.0 * seqcount / clusters);
          fprintf(stderr,
                  "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
                  singletons,
                  100.0 * singletons / seqcount,
                  100.0 * singletons / clusters);
        }

      if (opt_log)
        {
          fprintf(fp_log,
                  "Clusters: %d Size min %" PRId64 ", max %" PRId64 ", avg %.1f\n",
                  clusters,
                  abundance_min,
                  abundance_max,
                  1.0 * seqcount / clusters);
          fprintf(fp_log,
                  "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
                  singletons,
                  100.0 * singletons / seqcount,
                  100.0 * singletons / clusters);
          fprintf(fp_log, "\n");
        }
    }

  if (opt_clusterout_sort)
    {
      /* Optionally sort clusters by abundance */
      progress_init("Sorting clusters by abundance", clusters);
      qsort(clusterinfo, seqcount, sizeof(clusterinfo_t),
            compare_byclusterabundance);
      progress_done();
    }

  if (opt_msaout || opt_consout || opt_profile)
    {
      int msa_target_count = 0;
      struct msa_target_s * msa_target_list =
        (struct msa_target_s *) xmalloc(sizeof(struct msa_target_s) * size_max);
      progress_init("Multiple alignments", seqcount);

      FILE * fp_msaout = 0;
      FILE * fp_consout = 0;
      FILE * fp_profile = 0;

      if (opt_msaout)
        if (!(fp_msaout = fopen_output(opt_msaout)))
          fatal("Unable to open msaout file");

      if (opt_consout)
        if (!(fp_consout = fopen_output(opt_consout)))
          fatal("Unable to open consout file");

      if (opt_profile)
        if (!(fp_profile = fopen_output(opt_profile)))
          fatal("Unable to open profile file");

      lastcluster = -1;

      for(int i=0; i<seqcount; i++)
        {
          int clusterno = clusterinfo[i].clusterno;
          int seqno = clusterinfo[i].seqno;
          char * cigar = clusterinfo[i].cigar;
          int strand = clusterinfo[i].strand;

          if (clusterno != lastcluster)
            {
              if (lastcluster != -1)
                {
                  /* compute msa & consensus */
                  msa(fp_msaout, fp_consout, fp_profile,
                      lastcluster,
                      msa_target_count, msa_target_list,
                      cluster_abundance[lastcluster]);
                }

              /* start new cluster */
              msa_target_count = 0;
              lastcluster = clusterno;
            }

          /* add current sequence to the cluster */
          msa_target_list[msa_target_count].seqno = seqno;
          msa_target_list[msa_target_count].cigar = cigar;
          msa_target_list[msa_target_count].strand = strand;
          msa_target_count++;

          progress_update(i);
        }

      if (lastcluster != -1)
        {
          /* compute msa & consensus */
          msa(fp_msaout, fp_consout, fp_profile,
              lastcluster,
              msa_target_count, msa_target_list,
              cluster_abundance[lastcluster]);
        }

      progress_done();

      if (fp_profile)
        fclose(fp_profile);

      if (fp_msaout)
        fclose(fp_msaout);

      if (fp_consout)
        fclose(fp_consout);

      xfree(msa_target_list);
    }

  xfree(cluster_abundance);
  xfree(cluster_size);

  /* free cigar strings for all aligned sequences */

  for(int i=0; i<seqcount; i++)
    if (clusterinfo[i].cigar)
      xfree(clusterinfo[i].cigar);

  xfree(clusterinfo);

  if (fp_biomout)
    {
      otutable_print_biomout(fp_biomout);
      fclose(fp_biomout);
    }

  if (fp_otutabout)
    {
      otutable_print_otutabout(fp_otutabout);
      fclose(fp_otutabout);
    }

  if (fp_mothur_shared_out)
    {
      otutable_print_mothur_shared_out(fp_mothur_shared_out);
      fclose(fp_mothur_shared_out);
    }

  otutable_done();

  if (opt_matched)
    fclose(fp_matched);
  if (opt_notmatched)
    fclose(fp_notmatched);
  if (opt_fastapairs)
    fclose(fp_fastapairs);
  if (fp_blast6out)
    fclose(fp_blast6out);
  if (fp_userout)
    fclose(fp_userout);
  if (fp_alnout)
    fclose(fp_alnout);
  if (fp_samout)
    fclose(fp_samout);
  if (fp_uc)
    fclose(fp_uc);
  if (fp_centroids)
    fclose(fp_centroids);
  
  dbindex_free();
  db_free();
  show_rusage();
}

void cluster_fast(char * cmdline, char * progheader)
{
  cluster(opt_cluster_fast, cmdline, progheader);
}

void cluster_smallmem(char * cmdline, char * progheader)
{
  cluster(opt_cluster_smallmem, cmdline, progheader);
}

void cluster_size(char * cmdline, char * progheader)
{
  cluster(opt_cluster_size, cmdline, progheader);
}

void cluster_unoise(char * cmdline, char * progheader)
{
  cluster(opt_cluster_unoise, cmdline, progheader);
}
