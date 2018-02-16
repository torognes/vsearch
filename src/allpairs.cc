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

static pthread_t * pthread;

/* global constants/data, no need for synchronization */
static int seqcount; /* number of database sequences */
static pthread_attr_t attr;

/* global data protected by mutex */
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static int qmatches;
static int queries;
static int64_t progress = 0;
static FILE * fp_alnout = 0;
static FILE * fp_samout = 0;
static FILE * fp_userout = 0;
static FILE * fp_blast6out = 0;
static FILE * fp_uc = 0;
static FILE * fp_fastapairs = 0;
static FILE * fp_matched = 0;
static FILE * fp_notmatched = 0;

static int count_matched = 0;
static int count_notmatched = 0;

inline int allpairs_hit_compare_typed(struct hit * x, struct hit * y)
{
  // high id, then low id
  // early target, then late target

  if (x->id > y->id)
    return -1;
  else
    if (x->id < y->id)
      return +1;
    else
      if (x->target < y->target)
        return -1;
      else
        if (x->target > y->target)
          return +1;
        else
          return 0;
}

int allpairs_hit_compare(const void * a, const void * b)
{
  return allpairs_hit_compare_typed((struct hit *) a, (struct hit *) b);
}

void allpairs_output_results(int hit_count,
                             struct hit * hits,
                             char * query_head,
                             int qseqlen,
                             char * qsequence,
                             char * qsequence_rc)
{
  /* show results */
  int64_t toreport = MIN(opt_maxhits, hit_count);

  if (fp_alnout)
    results_show_alnout(fp_alnout,
                        hits,
                        toreport,
                        query_head,
                        qsequence,
                        qseqlen, 
                        qsequence_rc);

  if (fp_samout)
    results_show_samout(fp_samout,
                        hits,
                        toreport,
                        query_head,
                        qsequence,
                        qseqlen, 
                        qsequence_rc);

  if (toreport)
    {
      double top_hit_id = hits[0].id;
      
      for(int t = 0; t < toreport; t++)
        {
          struct hit * hp = hits + t;

          if (opt_top_hits_only && (hp->id < top_hit_id))
            break;
              
          if (fp_fastapairs)
            results_show_fastapairs_one(fp_fastapairs,
                                        hp, 
                                        query_head,
                                        qsequence,
                                        qseqlen,
                                        qsequence_rc);

          if (fp_uc)
            if ((t==0) || opt_uc_allhits)
              results_show_uc_one(fp_uc,
                                  hp,
                                  query_head,
                                  qsequence,
                                  qseqlen,
                                  qsequence_rc,
                                  hp->target);
              
          if (fp_userout)
            results_show_userout_one(fp_userout,
                                     hp,
                                     query_head, 
                                     qsequence,
                                     qseqlen,
                                     qsequence_rc);
              
          if (fp_blast6out)
            results_show_blast6out_one(fp_blast6out,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
        }
    }
  else
    {
      if (fp_uc)
        results_show_uc_one(fp_uc,
                            0,
                            query_head,
                            qsequence,
                            qseqlen,
                            qsequence_rc,
                            0);
      
      if (opt_output_no_hits)
        {
          if (fp_userout)
            results_show_userout_one(fp_userout,
                                     0,
                                     query_head,
                                     qsequence,
                                     qseqlen,
                                     qsequence_rc);

          if (fp_blast6out)
            results_show_blast6out_one(fp_blast6out,
                                       0,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
        }
    }

  if (hit_count)
    {
      count_matched++;
      if (opt_matched)
        fasta_print_general(fp_matched,
                            0,
                            qsequence,
                            qseqlen,
                            query_head,
                            strlen(query_head),
                            0,
                            count_matched,
                            -1, -1, 0, 0.0);
    }
  else
    {
      count_notmatched++;
      if (opt_notmatched)
        fasta_print_general(fp_notmatched,
                            0,
                            qsequence,
                            qseqlen,
                            query_head,
                            strlen(query_head),
                            0,
                            count_notmatched,
                            -1, -1, 0, 0.0);
    }
}

void allpairs_thread_run(int64_t t)
{
  struct searchinfo_s sia;

  struct searchinfo_s * si = & sia;
  
  si->strand = 0;
  si->query_head_alloc = 0;
  si->seq_alloc = 0;
  si->kmersamplecount = 0;
  si->kmers = 0;
  si->m = 0;
  si->finalized = 0;

  si->hits = (struct hit *) xmalloc(sizeof(struct hit) * seqcount);

  struct nwinfo_s * nw = nw_init();

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

  /* allocate memory for alignment results */
  unsigned int maxhits = seqcount;
  unsigned int * pseqnos = 
    (unsigned int *) xmalloc(sizeof(unsigned int) * maxhits);
  CELL * pscores = 
    (CELL*) xmalloc(sizeof(CELL) * maxhits);
  unsigned short * paligned =
    (unsigned short*) xmalloc(sizeof(unsigned short) * maxhits);
  unsigned short * pmatches =
    (unsigned short*) xmalloc(sizeof(unsigned short) * maxhits);
  unsigned short * pmismatches =
    (unsigned short*) xmalloc(sizeof(unsigned short) * maxhits);
  unsigned short * pgaps =
    (unsigned short*) xmalloc(sizeof(unsigned short) * maxhits);
  char** pcigar = (char**) xmalloc(sizeof(char*) * maxhits);

  struct hit * finalhits
    = (struct hit *) xmalloc(sizeof(struct hit) * seqcount);
  
  bool cont = 1;

  while (cont)
    {
      pthread_mutex_lock(&mutex_input);
      
      int query_no = queries;

      if (query_no < seqcount)
        {
          queries++;

          /* let other threads read input */
          pthread_mutex_unlock(&mutex_input);

          /* init search info */
          si->query_no = query_no;
          si->qsize = db_getabundance(query_no);
          si->query_head_len = db_getheaderlen(query_no);
          si->query_head = db_getheader(query_no);
          si->qseqlen = db_getsequencelen(query_no);
          si->qsequence = db_getsequence(query_no);
          si->rejects = 0;
          si->accepts = 0;
          si->hit_count = 0;

          for(int target = si->query_no + 1; 
              target < seqcount; target++)
            {
              if (opt_acceptall || search_acceptable_unaligned(si, target))
                pseqnos[si->hit_count++] = target;
            }
          
          if (si->hit_count)
            {
              /* perform alignments */

              search16_qprep(si->s, si->qsequence, si->qseqlen);
              
              search16(si->s,
                       si->hit_count,
                       pseqnos,
                       pscores,
                       paligned,
                       pmatches,
                       pmismatches,
                       pgaps,
                       pcigar);
              
              /* convert to hit structure */
              for (int h = 0; h < si->hit_count; h++)
                {
                  struct hit * hit = si->hits + h;
                  
                  unsigned int target = pseqnos[h];
                  int64_t nwscore = pscores[h];
                  
                  char * nwcigar;
                  int64_t nwalignmentlength;
                  int64_t nwmatches;
                  int64_t nwmismatches;
                  int64_t nwgaps;
                  
                  if (nwscore == SHRT_MAX)
                    {
                      /* In case the SIMD aligner cannot align,
                         perform a new alignment with the
                         linear memory aligner */
                      
                      char * tseq = db_getsequence(target);
                      int64_t tseqlen = db_getsequencelen(target);
                      
                      if (pcigar[h])
                        xfree(pcigar[h]);
                      
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
                      nwcigar = pcigar[h];
                      nwalignmentlength = paligned[h];
                      nwmatches = pmatches[h];
                      nwmismatches = pmismatches[h];
                      nwgaps = pgaps[h];
                    }
                  
                  hit->target = target;
                  hit->strand = 0;
                  hit->count = 0;
                  
                  hit->accepted = 0;
                  hit->rejected = 0;
                  hit->aligned = 1;
                  hit->weak = 0;
                  
                  hit->nwscore = nwscore;
                  hit->nwdiff = nwalignmentlength - nwmatches;
                  hit->nwgaps = nwgaps;
                  hit->nwindels = nwalignmentlength - nwmatches - nwmismatches;
                  hit->nwalignmentlength = nwalignmentlength;
                  hit->nwid = 100.0 * (nwalignmentlength - hit->nwdiff) /
                    nwalignmentlength;
                  hit->nwalignment = nwcigar;
                  hit->matches = nwalignmentlength - hit->nwdiff;
                  hit->mismatches = hit->nwdiff - hit->nwindels;
                  
                  int64_t dseqlen = db_getsequencelen(target);
                  hit->shortest = MIN(si->qseqlen, dseqlen);
                  hit->longest = MAX(si->qseqlen, dseqlen);
                  
                  /* trim alignment, compute numbers excluding terminal gaps */
                  align_trim(hit);
                  
                  /* test accept/reject criteria after alignment */
                  if (opt_acceptall || search_acceptable_aligned(si, hit))
                    finalhits[si->accepts++] = *hit;
                }
              
              /* sort hits */
              qsort(finalhits, si->accepts,
                    sizeof(struct hit), allpairs_hit_compare);
            }
          
          /* lock mutex for update of global data and output */
          pthread_mutex_lock(&mutex_output);
          
          /* output results */
          allpairs_output_results(si->accepts,
                                  finalhits,
                                  si->query_head,
                                  si->qseqlen,
                                  si->qsequence,
                                  0);
          
          /* update stats */
          if (si->accepts)
            qmatches++;
          
          /* show progress */
          progress += seqcount - query_no - 1;
          progress_update(progress);
          
          pthread_mutex_unlock(&mutex_output);
          
          /* free memory for alignment strings */
          for(int i=0; i < si->hit_count; i++)
            if (si->hits[i].aligned)
              xfree(si->hits[i].nwalignment);
        }
      else
        {
          /* let other threads read input */
          pthread_mutex_unlock(&mutex_input);

          cont = 0;
        }
    }

  xfree(finalhits);

  xfree(pcigar);
  xfree(pgaps);
  xfree(pmismatches);
  xfree(pmatches);
  xfree(paligned);
  xfree(pscores);
  xfree(pseqnos);

  search16_exit(si->s);

  nw_exit(nw);

  xfree(scorematrix);

  xfree(si->hits);
}

void * allpairs_thread_worker(void * vp)
{
  int64_t t = (int64_t) vp;
  allpairs_thread_run(t);
  return 0;
}

void allpairs_thread_worker_run()
{
  /* initialize threads, start them, join them and return */

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  /* init and create worker threads, put them into stand-by mode */
  for(int t=0; t<opt_threads; t++)
    {
      if (pthread_create(pthread+t, &attr,
                         allpairs_thread_worker, (void*)(int64_t)t))
        fatal("Cannot create thread");
    }

  /* finish and clean up worker threads */
  for(int t=0; t<opt_threads; t++)
    {
      if (pthread_join(pthread[t], NULL))
        fatal("Cannot join thread");
    }

  pthread_attr_destroy(&attr);
}


void allpairs_global(char * cmdline, char * progheader)
{
  opt_strand = 1;
  opt_uc_allhits = 1;

  /* open output files */

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

  if (opt_uc)
    {
      fp_uc = fopen_output(opt_uc);
      if (! fp_uc)
        fatal("Unable to open uc output file for writing");
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

  db_read(opt_allpairs_global, 0);

  results_show_samheader(fp_samout, cmdline, opt_allpairs_global);

  if (opt_qmask == MASK_DUST)
    dust_all();
  else if ((opt_qmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();

  show_rusage();

  seqcount = db_getsequencecount();

  /* prepare reading of queries */
  qmatches = 0;
  queries = 0;

  pthread = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));

  /* init mutexes for input and output */
  pthread_mutex_init(&mutex_input, NULL);
  pthread_mutex_init(&mutex_output, NULL);

  progress = 0;
  progress_init("Aligning", MAX(0,((int64_t)seqcount)*((int64_t)seqcount-1))/2);
  allpairs_thread_worker_run();
  progress_done();
  
  if (!opt_quiet)
    fprintf(stderr, "Matching query sequences: %d of %d (%.2f%%)\n", 
            qmatches, queries, 100.0 * qmatches / queries);

  if (opt_log)
    {
      fprintf(fp_log, "Matching query sequences: %d of %d (%.2f%%)\n\n", 
              qmatches, queries, 100.0 * qmatches / queries);
    }

  pthread_mutex_destroy(&mutex_output);
  pthread_mutex_destroy(&mutex_input);

  xfree(pthread);

  /* clean up, global */
  db_free();
  if (opt_matched)
    fclose(fp_matched);
  if (opt_notmatched)
    fclose(fp_notmatched);
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
  if (fp_samout)
    fclose(fp_samout);
  show_rusage();
}
