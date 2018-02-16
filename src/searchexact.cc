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

static struct searchinfo_s * si_plus;
static struct searchinfo_s * si_minus;
static pthread_t * pthread;

/* global constants/data, no need for synchronization */
static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */
static pthread_attr_t attr;
static fastx_handle query_fasta_h;

/* global data protected by mutex */
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static int qmatches;
static int queries;
static int * dbmatched;
static FILE * fp_samout = 0;
static FILE * fp_alnout = 0;
static FILE * fp_userout = 0;
static FILE * fp_blast6out = 0;
static FILE * fp_uc = 0;
static FILE * fp_fastapairs = 0;
static FILE * fp_matched = 0;
static FILE * fp_notmatched = 0;
static FILE * fp_dbmatched = 0;
static FILE * fp_dbnotmatched = 0;
static FILE * fp_otutabout = 0;
static FILE * fp_mothur_shared_out = 0;
static FILE * fp_biomout = 0;

static int count_matched = 0;
static int count_notmatched = 0;

void add_hit(struct searchinfo_s * si, uint64_t seqno)
{
  if (search_acceptable_unaligned(si, seqno))
    {
      struct hit * hp = si->hits + si->hit_count;
      si->hit_count++;
      
      hp->target = seqno;
      hp->strand = si->strand;
      
      hp->count = 0;
      
      hp->nwscore = si->qseqlen * opt_match;
      hp->nwdiff = 0;
      hp->nwgaps = 0;
      hp->nwindels = 0;
      hp->nwalignmentlength = si->qseqlen;
      hp->nwid = 100.0;
      hp->matches = si->qseqlen;
      hp->mismatches = 0;
      
      int ret = xsprintf(&hp->nwalignment, "%dM", si->qseqlen);
      if ((ret == -1) || (!hp->nwalignment))
        fatal("Out of memory");
      
      hp->internal_alignmentlength = si->qseqlen;
      hp->internal_gaps = 0;
      hp->internal_indels = 0;
      hp->trim_q_left = 0;
      hp->trim_q_right = 0;
      hp->trim_t_left = 0;
      hp->trim_t_right = 0;
      hp->trim_aln_left = 0;
      hp->trim_aln_right = 0;
      
      hp->id = 100.0;
      hp->id0 = 100.0;
      hp->id1 = 100.0;
      hp->id2 = 100.0;
      hp->id3 = 100.0;
      hp->id4 = 100.0;
      
      hp->shortest = si->qseqlen;
      hp->longest = si->qseqlen;
      
      hp->aligned = 1;

      hp->accepted = 0;
      hp->rejected = 0;
      hp->weak = 0;
      (void) search_acceptable_aligned(si, hp);
    }
}

void search_exact_onequery(struct searchinfo_s * si)
{
  dbhash_search_info_s info;

  char * seq = si->qsequence;
  uint64_t seqlen = si->qseqlen;
  char * normalized = (char*) xmalloc(seqlen+1);
  string_normalize(normalized, seq, seqlen);

  si->hit_count = 0;

  int64_t ret = dbhash_search_first(normalized, seqlen, & info);
  while (ret >= 0)
    {
      add_hit(si, ret);
      ret = dbhash_search_next(&info);
    }
  xfree(normalized);
}

void search_exact_output_results(int hit_count,
                                 struct hit * hits,
                                 char * query_head,
                                 int qseqlen,
                                 char * qsequence,
                                 char * qsequence_rc,
                                 int qsize)
{
  pthread_mutex_lock(&mutex_output);

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

      if (opt_otutabout || opt_mothur_shared_out || opt_biomout)
        otutable_add(query_head,
                     db_getheader(hits[0].target),
                     qsize);

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
                            qsize,
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
                            qsize,
                            count_notmatched,
                            -1, -1, 0, 0.0);
    }

  /* update matching db sequences */
  for (int i=0; i < hit_count; i++)
    if (hits[i].accepted)
      dbmatched[hits[i].target]++;
  
  pthread_mutex_unlock(&mutex_output);
}

int search_exact_query(int64_t t)
{
  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = s ? si_minus+t : si_plus+t;

      /* mask query */
      if (opt_qmask == MASK_DUST)
        {
          dust(si->qsequence, si->qseqlen);
        }
      else if ((opt_qmask == MASK_SOFT) && (opt_hardmask))
        {
          hardmask(si->qsequence, si->qseqlen);
        }

      /* perform search */
      search_exact_onequery(si);
    }

  struct hit * hits;
  int hit_count;

  search_joinhits(si_plus + t,
                  opt_strand > 1 ? si_minus + t : 0,
                  & hits,
                  & hit_count);

  search_exact_output_results(hit_count,
                              hits,
                              si_plus[t].query_head,
                              si_plus[t].qseqlen,
                              si_plus[t].qsequence,
                              opt_strand > 1 ? si_minus[t].qsequence : 0,
                              si_plus[t].qsize);

  /* free memory for alignment strings */
  for(int i=0; i<hit_count; i++)
    if (hits[i].aligned)
      xfree(hits[i].nwalignment);

  xfree(hits);

  return hit_count;
}

void search_exact_thread_run(int64_t t)
{
  while (1)
    {
      pthread_mutex_lock(&mutex_input);

      if (fasta_next(query_fasta_h, ! opt_notrunclabels, chrmap_no_change))
        {
          char * qhead = fasta_get_header(query_fasta_h);
          int query_head_len = fasta_get_header_length(query_fasta_h);
          char * qseq = fasta_get_sequence(query_fasta_h);
          int qseqlen = fasta_get_sequence_length(query_fasta_h);
          int query_no = fasta_get_seqno(query_fasta_h);
          int qsize = fasta_get_abundance(query_fasta_h);
          
          for (int s = 0; s < opt_strand; s++)
            {
              struct searchinfo_s * si = s ? si_minus+t : si_plus+t;
              
              si->query_head_len = query_head_len;
              si->qseqlen = qseqlen;
              si->query_no = query_no;
              si->qsize = qsize;
              si->strand = s;
              
              /* allocate more memory for header and sequence, if necessary */
              
              if (si->query_head_len + 1 > si->query_head_alloc)
                {
                  si->query_head_alloc = si->query_head_len + 2001;
                  si->query_head = (char*)
                    xrealloc(si->query_head, (size_t)(si->query_head_alloc));
                }
              
              if (si->qseqlen + 1 > si->seq_alloc)
                {
                  si->seq_alloc = si->qseqlen + 2001;
                  si->qsequence = (char*)
                    xrealloc(si->qsequence, (size_t)(si->seq_alloc));
                }
            }
              
          /* plus strand: copy header and sequence */
          strcpy(si_plus[t].query_head, qhead);
          strcpy(si_plus[t].qsequence, qseq);
          
          /* get progress as amount of input file read */
          uint64_t progress = fasta_get_position(query_fasta_h);

          /* let other threads read input */
          pthread_mutex_unlock(&mutex_input);
          
          /* minus strand: copy header and reverse complementary sequence */
          if (opt_strand > 1)
            {
              strcpy(si_minus[t].query_head, si_plus[t].query_head);
              reverse_complement(si_minus[t].qsequence,
                                 si_plus[t].qsequence,
                                 si_plus[t].qseqlen);
            }
          
          int match = search_exact_query(t);
          
          /* lock mutex for update of global data and output */
          pthread_mutex_lock(&mutex_output);

          /* update stats */
          queries++;

          if (match)
            qmatches++;

          /* show progress */
          progress_update(progress);

          pthread_mutex_unlock(&mutex_output);
        }
      else
        {
          pthread_mutex_unlock(&mutex_input);
          break;
        }
    }
}

void search_exact_thread_init(struct searchinfo_s * si)
{
  /* thread specific initialiation */
  si->uh = 0;
  si->kmers = 0;
  si->m = 0;
  si->hits = (struct hit *) xmalloc
    (sizeof(struct hit) * (tophits) * opt_strand);
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = 0;
  si->seq_alloc = 0;
  si->qsequence = 0;
  si->nw = 0;
  si->s = 0;
}

void search_exact_thread_exit(struct searchinfo_s * si)
{
  /* thread specific clean up */
  xfree(si->hits);
  if (si->query_head)
    xfree(si->query_head);
  if (si->qsequence)
    xfree(si->qsequence);
}

void * search_exact_thread_worker(void * vp)
{
  int64_t t = (int64_t) vp;
  search_exact_thread_run(t);
  return 0;
}

void search_exact_thread_worker_run()
{
  /* initialize threads, start them, join them and return */

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  /* init and create worker threads, put them into stand-by mode */
  for(int t=0; t<opt_threads; t++)
    {
      search_exact_thread_init(si_plus+t);
      if (si_minus)
        search_exact_thread_init(si_minus+t);
      if (pthread_create(pthread+t, &attr,
                         search_exact_thread_worker, (void*)(int64_t)t))
        fatal("Cannot create thread");
    }

  /* finish and clean up worker threads */
  for(int t=0; t<opt_threads; t++)
    {
      if (pthread_join(pthread[t], NULL))
        fatal("Cannot join thread");
      search_exact_thread_exit(si_plus+t);
      if (si_minus)
        search_exact_thread_exit(si_minus+t);
    }

  pthread_attr_destroy(&attr);
}

void search_exact_prep(char * cmdline, char * progheader)
{
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

  if (opt_dbmatched)
    {
      fp_dbmatched = fopen_output(opt_dbmatched);
      if (! fp_dbmatched)
        fatal("Unable to open dbmatched output file for writing");
    }

  if (opt_dbnotmatched)
    {
      fp_dbnotmatched = fopen_output(opt_dbnotmatched);
      if (! fp_dbnotmatched)
        fatal("Unable to open dbnotmatched output file for writing");
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

  db_read(opt_db, 0);

  results_show_samheader(fp_samout, cmdline, opt_db);

  if (opt_dbmask == MASK_DUST)
    dust_all();
  else if ((opt_dbmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();

  show_rusage();

  seqcount = db_getsequencecount();

  /* tophits = the maximum number of hits we need to store */
  tophits = seqcount;

  dbmatched = (int*) xmalloc(seqcount * sizeof(int*));
  memset(dbmatched, 0, seqcount * sizeof(int*));

  dbhash_open(seqcount);
  dbhash_add_all();
}

void search_exact_done()
{
  /* clean up, global */
  dbhash_close();

  db_free();
  xfree(dbmatched);

  if (opt_dbmatched)
    fclose(fp_dbmatched);
  if (opt_dbnotmatched)
    fclose(fp_dbnotmatched);
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

void search_exact(char * cmdline, char * progheader)
{
  opt_id = 1.0;

  search_exact_prep(cmdline, progheader);

  otutable_init();

  /* prepare reading of queries */
  qmatches = 0;
  queries = 0;
  query_fasta_h = fasta_open(opt_search_exact);

  /* allocate memory for thread info */
  si_plus = (struct searchinfo_s *) xmalloc(opt_threads * 
                                            sizeof(struct searchinfo_s));
  if (opt_strand > 1)
    si_minus = (struct searchinfo_s *) xmalloc(opt_threads * 
                                               sizeof(struct searchinfo_s));
  else
    si_minus = 0;
  
  pthread = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));

  /* init mutexes for input and output */
  pthread_mutex_init(&mutex_input, NULL);
  pthread_mutex_init(&mutex_output, NULL);

  progress_init("Searching", fasta_get_size(query_fasta_h));
  search_exact_thread_worker_run();
  progress_done();
  
  pthread_mutex_destroy(&mutex_output);
  pthread_mutex_destroy(&mutex_input);

  xfree(pthread);
  xfree(si_plus);
  if (si_minus)
    xfree(si_minus);

  fasta_close(query_fasta_h);

  if (!opt_quiet)
    fprintf(stderr, "Matching query sequences: %d of %d (%.2f%%)\n", 
            qmatches, queries, 100.0 * qmatches / queries);

  if (opt_log)
    fprintf(fp_log, "Matching query sequences: %d of %d (%.2f%%)\n", 
            qmatches, queries, 100.0 * qmatches / queries);

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

  int count_dbmatched = 0;
  int count_dbnotmatched = 0;

  if (opt_dbmatched || opt_dbnotmatched)
    {
      for(int64_t i=0; i<seqcount; i++)
        if (dbmatched[i])
          {
            count_dbmatched++;
            if (opt_dbmatched)
              fasta_print_general(fp_dbmatched,
                                  0,
                                  db_getsequence(i),
                                  db_getsequencelen(i),
                                  db_getheader(i),
                                  db_getheaderlen(i),
                                  dbmatched[i],
                                  count_dbmatched,
                                  -1, -1, 0, 0.0);
          }
        else
          {
            count_dbnotmatched++;
            if (opt_dbnotmatched)
              fasta_print_general(fp_dbnotmatched,
                                  0,
                                  db_getsequence(i),
                                  db_getsequencelen(i),
                                  db_getheader(i),
                                  db_getheaderlen(i),
                                  0,
                                  count_dbnotmatched,
                                  -1, -1, 0, 0.0);
          }
    }

  search_exact_done();
}
