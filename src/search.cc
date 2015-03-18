/*
  Copyright (C) 2014-2015 Torbjorn Rognes

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

static struct searchinfo_s * si_plus;
static struct searchinfo_s * si_minus;
static pthread_t * pthread;

/* global constants/data, no need for synchronization */
static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */
static pthread_attr_t attr;

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

void search_output_results(int hit_count,
                           struct hit * hits,
                           char * query_head,
                           int qseqlen,
                           char * qsequence,
                           char * qsequence_rc)
{
  pthread_mutex_lock(&mutex_output);

  /* show results */
  long toreport = MIN(opt_maxhits, hit_count);

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
                                  qsequence_rc);
              
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
  else if (opt_output_no_hits)
    {
      if (fp_uc)
        results_show_uc_one(fp_uc,
                            0,
                            query_head,
                            qsequence,
                            qseqlen,
                            qsequence_rc);
      
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

  if (hit_count)
    {
      if (opt_matched)
        {
          fprintf(fp_matched,
                  ">%s\n",
                  query_head);
          fprint_fasta_seq_only(fp_matched,
                                qsequence,
                                qseqlen,
                                opt_fasta_width);
        }
    }
  else
    {
      if (opt_notmatched)
        {
          fprintf(fp_notmatched,
                  ">%s\n",
                  query_head);
          fprint_fasta_seq_only(fp_notmatched,
                                qsequence,
                                qseqlen,
                                opt_fasta_width);
        }
    }

  /* update matching db sequences */
  for (int i=0; i < hit_count; i++)
    if (hits[i].accepted)
      dbmatched[hits[i].target]++;
  
  pthread_mutex_unlock(&mutex_output);
}

int search_query(long t)
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
      search_onequery(si);
    }

  struct hit * hits;
  int hit_count;

  search_joinhits(si_plus + t,
                  opt_strand > 1 ? si_minus + t : 0,
                  & hits,
                  & hit_count);

  search_output_results(hit_count,
                        hits,
                        si_plus[t].query_head,
                        si_plus[t].qseqlen,
                        si_plus[t].qsequence,
                        opt_strand > 1 ? si_minus[t].qsequence : 0);

  /* free memory for alignment strings */
  for(int i=0; i<hit_count; i++)
    if (hits[i].aligned)
      free(hits[i].nwalignment);

  free(hits);

  return hit_count;
}

void search_thread_run(long t)
{
  while (1)
    {
      pthread_mutex_lock(&mutex_input);
      
      char * qhead;
      int query_head_len;
      char * qseq;
      int qseqlen;
      int query_no;
      int qsize;

      if (query_getnext(& qhead, & query_head_len,
                        & qseq, & qseqlen,
                        & query_no, & qsize,
                        opt_qmask != MASK_SOFT))
        {

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
          unsigned long progress = query_getfilepos();

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
          
          int match = search_query(t);
          
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

void search_thread_init(struct searchinfo_s * si)
{
  /* thread specific initialiation */
  si->uh = unique_init();
  si->kmers = (count_t *) xmalloc(seqcount * sizeof(count_t) + 32);
  si->m = minheap_init(tophits);
  si->hits = (struct hit *) xmalloc
    (sizeof(struct hit) * (tophits) * opt_strand);
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = 0;
  si->seq_alloc = 0;
  si->qsequence = 0;
#ifdef COMPARENONVECTORIZED
  si->nw = nw_init();
#else
  si->nw = 0;
#endif
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
}

void search_thread_exit(struct searchinfo_s * si)
{
  /* thread specific clean up */
  search16_exit(si->s);
#ifdef COMPARENONVECTORIZED
  nw_exit(si->nw);
#endif
  unique_exit(si->uh);
  free(si->hits);
  minheap_exit(si->m);
  free(si->kmers);
  if (si->query_head)
    free(si->query_head);
  if (si->qsequence)
    free(si->qsequence);
}



void * search_thread_worker(void * vp)
{
  long t = (long) vp;
  search_thread_run(t);
  return 0;
}

void search_thread_worker_run()
{
  /* initialize threads, start them, join them and return */

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  /* init and create worker threads, put them into stand-by mode */
  for(int t=0; t<opt_threads; t++)
    {
      search_thread_init(si_plus+t);
      if (si_minus)
        search_thread_init(si_minus+t);
      if (pthread_create(pthread+t, &attr,
                         search_thread_worker, (void*)(long)t))
        fatal("Cannot create thread");
    }

  /* finish and clean up worker threads */
  for(int t=0; t<opt_threads; t++)
    {
      if (pthread_join(pthread[t], NULL))
        fatal("Cannot join thread");
      search_thread_exit(si_plus+t);
      if (si_minus)
        search_thread_exit(si_minus+t);
    }

  pthread_attr_destroy(&attr);
}



void search_prep(char * cmdline, char * progheader)
{
  /* open output files */

  if (opt_alnout)
    {
      fp_alnout = fopen(opt_alnout, "w");
      if (! fp_alnout)
        fatal("Unable to open alignment output file for writing");

      fprintf(fp_alnout, "%s\n", cmdline);
      fprintf(fp_alnout, "%s\n", progheader);
    }

  if (opt_samout)
    {
      fp_samout = fopen(opt_samout, "w");
      if (! fp_samout)
        fatal("Unable to open SAM output file for writing");
    }

  if (opt_userout)
    {
      fp_userout = fopen(opt_userout, "w");
      if (! fp_userout)
        fatal("Unable to open user-defined output file for writing");
    }

  if (opt_blast6out)
    {
      fp_blast6out = fopen(opt_blast6out, "w");
      if (! fp_blast6out)
        fatal("Unable to open blast6-like output file for writing");
    }

  if (opt_uc)
    {
      fp_uc = fopen(opt_uc, "w");
      if (! fp_uc)
        fatal("Unable to open uc output file for writing");
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

  db_read(opt_db, opt_dbmask != MASK_SOFT);

  if (opt_dbmask == MASK_DUST)
    dust_all();
  else if ((opt_dbmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();

  show_rusage();

  seqcount = db_getsequencecount();

  dbindex_prepare(1);
  dbindex_addallsequences();

  /* tophits = the maximum number of hits we need to store */

  if ((opt_maxrejects == 0) || (opt_maxrejects > seqcount))
    opt_maxrejects = seqcount;

  if ((opt_maxaccepts == 0) || (opt_maxaccepts > seqcount))
    opt_maxaccepts = seqcount;

  tophits = opt_maxrejects + opt_maxaccepts + MAXDELAYED;

  if (tophits > seqcount)
    tophits = seqcount;
}

void search_done()
{
  /* clean up, global */
  dbindex_free();
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



void usearch_global(char * cmdline, char * progheader)
{
  search_prep(cmdline, progheader);

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

  dbmatched = (int*) xmalloc(seqcount * sizeof(int*));
  memset(dbmatched, 0, seqcount * sizeof(int*));

  /* prepare reading of queries */
  qmatches = 0;
  queries = 0;
  query_open(opt_usearch_global);

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

  progress_init("Searching", query_getfilesize());
  search_thread_worker_run();
  progress_done();
  
  pthread_mutex_destroy(&mutex_output);
  pthread_mutex_destroy(&mutex_input);

  free(pthread);
  free(si_plus);
  if (si_minus)
    free(si_minus);

  query_close();

  if (!opt_quiet)
    fprintf(stderr, "Matching query sequences: %d of %d (%.2f%%)\n", 
            qmatches, queries, 100.0 * qmatches / queries);

  if (opt_log)
    fprintf(fp_log, "Matching query sequences: %d of %d (%.2f%%)\n", 
            qmatches, queries, 100.0 * qmatches / queries);

  if (opt_dbmatched || opt_dbnotmatched)
    {
      for(long i=0; i<seqcount; i++)
        if (dbmatched[i])
          {
            if (opt_dbmatched)
              {
                if (opt_sizeout)
                  db_fprint_fasta_with_size(fp_dbmatched, i, dbmatched[i]);
                else
                  db_fprint_fasta(fp_dbmatched, i);
              }
          }
        else
          {
            if (opt_dbnotmatched)
              db_fprint_fasta(fp_dbnotmatched, i);
          }
    }

  free(dbmatched);

  if (opt_dbmatched)
    fclose(fp_dbmatched);
  if (opt_dbnotmatched)
    fclose(fp_dbnotmatched);

  search_done();
}
