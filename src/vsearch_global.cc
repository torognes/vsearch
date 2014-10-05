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

static struct searchinfo_s * sia;

/* global constants/data, no need for synchronization */
static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */
static pthread_attr_t attr;

/* global data protected by mutex */
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static long qmatches;
static long queries;
static int * dbmatched;
static FILE * fp_alnout = 0;
static FILE * fp_userout = 0;
static FILE * fp_blast6out = 0;
static FILE * fp_uc = 0;
static FILE * fp_fastapairs = 0;
static FILE * fp_matched = 0;
static FILE * fp_notmatched = 0;
static FILE * fp_dbmatched = 0;
static FILE * fp_dbnotmatched = 0;

void search_output_results(struct searchinfo_s * si)
{
  /* show results */
  long toreport = MIN(opt_maxhits, si->hit_count);
  if (opt_weak_id == opt_id)
    toreport = MIN(toreport, maxaccepts);

  if (toreport)
    {
      if (fp_alnout)
        results_show_alnout(fp_alnout,
                            si->hits, toreport, si->query_head,
                            si->qsequence, si->qseqlen, si->rc);
      
      double top_hit_id = si->hits[0].internal_id;
          
      for(int t = 0; t < toreport; t++)
        {
          struct hit * hp = si->hits + t;
              
          if (opt_top_hits_only && (hp->internal_id < top_hit_id))
            break;
              
          if (fp_fastapairs)
            results_show_fastapairs_one(fp_fastapairs, hp, 
                                        si->query_head,
                                        si->qsequence,
                                        si->qseqlen,
                                        si->rc);

          if (fp_uc)
            if (opt_uc_allhits || (t==0))
              results_show_uc_one(fp_uc, hp, si->query_head,
                                  si->qsequence, si->qseqlen, si->rc);
              
          if (fp_userout)
            results_show_userout_one(fp_userout, hp, si->query_head, 
                                     si->qsequence, si->qseqlen, si->rc);
              
          if (fp_blast6out)
            results_show_blast6out_one(fp_blast6out, hp, si->query_head,
                                       si->qsequence, si->qseqlen, si->rc);
        }
    }
  else
    {
      if (fp_uc)
        results_show_uc_one(fp_uc, 0, si->query_head,
                            si->qsequence, si->qseqlen, si->rc);
          
      if (opt_output_no_hits)
        {
          if (fp_userout)
            results_show_userout_one(fp_userout, 0, si->query_head, 
                                     si->qsequence, si->qseqlen, si->rc);
              
          if (fp_blast6out)
            results_show_blast6out_one(fp_blast6out, 0, si->query_head,
                                       si->qsequence, si->qseqlen, si->rc);
        }
    }

  if (si->hit_count)
    {
      if (opt_matched)
        {
          fprintf(fp_matched, ">%s\n", si->query_head);
          fprint_fasta_seq_only(fp_matched, si->qsequence, si->qseqlen,
                                opt_fasta_width);
        }
    }
  else
    {
      if (opt_notmatched)
        {
          fprintf(fp_notmatched, ">%s\n", si->query_head);
          fprint_fasta_seq_only(fp_notmatched, si->qsequence, si->qseqlen,
                                opt_fasta_width);
        }
    }
}

void search_query(struct searchinfo_s * si,
                  char * qhead, int query_head_len,
                  char * qseq, int qseqlen,
                  int query_no, int qsize,
                  pthread_mutex_t * mutex_output)
{
  /* mask query */
  if (opt_qmask == MASK_DUST)
    {
      dust(si->qsequence, si->qseqlen);
    }
  else if ((opt_qmask == MASK_SOFT) && (opt_hardmask))
    {
      hardmask(si->qsequence, si->qseqlen);
    }
  
  /* compute reverse complement query sequence */
  if (opt_strand > 1)
    {
      if (si->qseqlen + 1 > si->rc_seq_alloc)
        {
          si->rc_seq_alloc = si->qseqlen + 2001;
          si->rc = (char *) xrealloc(si->rc, (size_t)(si->rc_seq_alloc));
        }
      reverse_complement(si->rc, si->qsequence, si->qseqlen);
    }

  /* perform search */
  search_onequery(si);
  
  pthread_mutex_lock(mutex_output);
  search_output_results(si);
  pthread_mutex_unlock(mutex_output);

  /* free memory for alignment strings */
  for(int i=0; i<si->hit_count; i++)
    free(si->hits[i].nwalignment);
}

void search_thread_run(struct searchinfo_s * si)
{
  while (1)
    {
      pthread_mutex_lock(&mutex_input);
      
      char * qhead;
      long query_head_len;
      char * qseq;
      long qseqlen;
      long query_no;
      long qsize;

      if (query_getnext(& qhead, & query_head_len,
                        & qseq, & qseqlen,
                        & query_no, & qsize,
                        opt_qmask != MASK_SOFT))
        {

          si->query_head_len = query_head_len;
          si->qseqlen = qseqlen;
          si->query_no = query_no;
          si->qsize = qsize;
          
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
          
          /* copy header and sequence */
          strcpy(si->query_head, qhead);
          strcpy(si->qsequence, qseq);
          
          /* let other threads read input */
          pthread_mutex_unlock(&mutex_input);

          search_query(si,
                       qhead, query_head_len,
                       qseq, qseqlen,
                       query_no, qsize,
                       & mutex_output);
          
          queries++;
          
          if (si->hit_count > 0)
            qmatches++;
          
          /* update matching db sequences */
          for (long i=0; i<si->hit_count; i++)
            if (si->hits[i].internal_id >= 100.0 * opt_id)
              dbmatched[si->hits[i].target]++;
          
          /* show progress */
          pthread_mutex_lock(&mutex_input);
          pthread_mutex_lock(&mutex_output);
          int progress = query_getfilepos();
          progress_update(progress);
          pthread_mutex_unlock(&mutex_input);
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
  si->kmers = (count_t *) xmalloc(seqcount * sizeof(count_t));
  si->m = minheap_init(tophits);
  si->targetlist = (unsigned int*) xmalloc
    (sizeof(unsigned int)*seqcount);
  si->hits = (struct hit *) xmalloc
    (sizeof(struct hit) * (tophits) * opt_strand);
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = 0;
  si->seq_alloc = 0;
  si->qsequence = 0;
  si->rc_seq_alloc = 0;
  si->rc = 0;
#ifdef COMPARENONVECTORIZED
  si->nw = nw_init();
#else
  si->nw = 0;
#endif
  si->s = search16_init(match_score,
                        mismatch_score,
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
  free(si->targetlist);
  free(si->hits);
  minheap_exit(si->m);
  free(si->kmers);
  if (si->query_head)
    free(si->query_head);
  if (si->qsequence)
    free(si->qsequence);
  if (si->rc)
    free(si->rc);
}



void * search_thread_worker(void * vp)
{
  long t = (long) vp;
  struct searchinfo_s * si = sia + t;
  search_thread_run(si);
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
      struct searchinfo_s * si = sia + t;
      search_thread_init(si);
      if (pthread_create(&si->pthread, &attr, search_thread_worker,
                         (void*)(long)t))
        fatal("Cannot create thread");
    }

  /* finish and clean up worker threads */
  for(int t=0; t<opt_threads; t++)
    {
      struct searchinfo_s * si = sia + t;
      
      /* wait for worker to quit */
      if (pthread_join(si->pthread, NULL))
        fatal("Cannot join thread");

      search_thread_exit(si);
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

  if (useroutfilename)
    {
      fp_userout = fopen(useroutfilename, "w");
      if (! fp_userout)
        fatal("Unable to open user-defined output file for writing");
    }

  if (opt_blast6out)
    {
      fp_blast6out = fopen(opt_blast6out, "w");
      if (! fp_blast6out)
        fatal("Unable to open blast6-like output file for writing");
    }

  if (ucfilename)
    {
      fp_uc = fopen(ucfilename, "w");
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

  dbindex_prepare(1);
  dbindex_addallsequences();

  seqcount = db_getsequencecount();

  if ((maxrejects == 0) || (maxrejects > seqcount))
    maxrejects = seqcount;

  if (maxaccepts > seqcount)
    maxaccepts = seqcount;

  tophits = maxrejects + maxaccepts + 8;
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
  show_rusage();
}



void search(char * cmdline, char * progheader)
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
  query_open(opt_vsearch_global);

  /* allocate memory for thread info */
  sia = (struct searchinfo_s *) xmalloc(opt_threads * 
                                       sizeof(struct searchinfo_s));

  /* init mutexes for input and output */
  pthread_mutex_init(&mutex_input, NULL);
  pthread_mutex_init(&mutex_output, NULL);

  progress_init("Searching", query_getfilesize());
  search_thread_worker_run();
  progress_done();

  pthread_mutex_destroy(&mutex_output);
  pthread_mutex_destroy(&mutex_input);

  free(sia);

  query_close();

  fprintf(stderr, "Matching query sequences: %ld of %ld (%.2f%%)\n", 
          qmatches, queries, 100.0 * qmatches / queries);

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

  free(dbmatched);

  if (opt_dbmatched)
    fclose(fp_dbmatched);
  if (opt_dbnotmatched)
    fclose(fp_dbnotmatched);

  search_done();
}
