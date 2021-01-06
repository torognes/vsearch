/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2021, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

/*

  Implements the Sintax algorithm as desribed in Robert Edgar's preprint:

  Robert Edgar (2016)
  SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences
  BioRxiv, 074161
  doi: https://doi.org/10.1101/074161

  Further details:

  https://www.drive5.com/usearch/manual/cmd_sintax.html

*/

#include "vsearch.h"

static struct searchinfo_s * si_plus;
static struct searchinfo_s * si_minus;
static pthread_t * pthread;

/* global constants/data, no need for synchronization */
static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */
static pthread_attr_t attr;
static fastx_handle query_fastx_h;

const int tax_levels = 8;
const char * tax_letters = "dkpcofgs";
const int subset_size = 32;
const int bootstrap_count = 100;

/* global data protected by mutex */
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static FILE * fp_tabbedout;
static int queries = 0;
static int classified = 0;

bool sintax_parse_tax(const char * header,
                      int header_length,
                      int * tax_start,
                      int * tax_end)
{
  /*
    Identify the first occurence of the pattern (^|;)tax=([^;]*)(;|$)
  */

  if (! header)
    return false;

  const char * attribute = "tax=";

  int hlen = header_length;
  int alen = strlen(attribute);

  int i = 0;

  while (i < hlen - alen)
    {
      char * r = (char *) strstr(header + i, attribute);

      /* no match */
      if (r == NULL)
        break;

      i = r - header;

      /* check for ';' in front */
      if ((i > 0) && (header[i-1] != ';'))
        {
          i += alen + 1;
          continue;
        }

      * tax_start = i;

      /* find end (semicolon or end of header) */
      const char * s = strchr(header+i+alen, ';');
      if (s == 0)
        * tax_end = hlen;
      else
        * tax_end = s - header;

      return true;
    }
  return false;
}

void sintax_split(int seqno, int * level_start, int * level_len)
{
  /* Parse taxonomy string into the following parts
     d domain
     k kingdom
     p phylum
     c class
     o order
     f family
     g genus
     s species
  */

  for (int i = 0; i < tax_levels; i++)
    {
      level_start[i] = 0;
      level_len[i] = 0;
    }

  int tax_start, tax_end;
  char * h = db_getheader(seqno);
  int hlen = db_getheaderlen(seqno);
  if (sintax_parse_tax(h, hlen, & tax_start, & tax_end))
    {
      int t = tax_start + 4;

      while (t < tax_end)
        {
          /* Is the next char a recogized tax level letter? */
          const char * r = strchr(tax_letters, tolower(h[t]));
          if (r)
            {
              int level = r - tax_letters;

              /* Is there a colon after it? */
              if (h[t + 1] == ':')
                {
                  level_start[level] = t + 2;

                  char * z = strchr(h + t + 2, ',');
                  if (z)
                    level_len[level] = z - h - t - 2;
                  else
                    level_len[level] = tax_end - t - 2;
                }
            }

          /* skip past next comma */
          char * x = strchr(h + t, ',');
          if (x)
            t = x - h + 1;
          else
            t = tax_end;
        }
    }
}

void sintax_analyse(char * query_head,
                    int strand,
                    int best_seqno,
                    int best_count,
                    int * all_seqno,
                    int count)
{
  int best_level_start[tax_levels];
  int best_level_len[tax_levels];
  int level_match[tax_levels];

  /* check number of successful bootstraps */
  if (count >= (bootstrap_count+1) / 2)
    {
      char * best_h = db_getheader(best_seqno);

      sintax_split(best_seqno, best_level_start, best_level_len);

      for (int j = 0; j < tax_levels; j++)
        level_match[j] = 0;

      for (int i = 0; i < count; i++)
        {
          /* For each bootstrap experiment */

          int level_start[tax_levels];
          int level_len[tax_levels];
          sintax_split(all_seqno[i], level_start, level_len);

          char * h = db_getheader(all_seqno[i]);

          for (int j = 0; j < tax_levels; j++)
            {
              /* For each taxonomic level */

              if ((level_len[j] == best_level_len[j]) &&
                  (strncmp(best_h + best_level_start[j],
                           h + level_start[j],
                           level_len[j]) == 0))
                {
                  level_match[j]++;
                }
            }
        }
    }

  /* write to tabbedout file */
  xpthread_mutex_lock(&mutex_output);
  fprintf(fp_tabbedout, "%s\t", query_head);

  queries++;

  if (count >= bootstrap_count / 2)
    {
      char * best_h = db_getheader(best_seqno);

      classified++;

      bool comma = false;
      for (int j = 0; j < tax_levels; j++)
        {
          if (best_level_len[j] > 0)
            {
              fprintf(fp_tabbedout,
                      "%s%c:%.*s(%.2f)",
                      (comma ? "," : ""),
                      tax_letters[j],
                      best_level_len[j],
                      best_h + best_level_start[j],
                      1.0 * level_match[j] / count);
              comma = true;
            }
        }

      fprintf(fp_tabbedout, "\t%c", strand ? '-' : '+');

      if (opt_sintax_cutoff > 0.0)
        {
          fprintf(fp_tabbedout, "\t");
          bool comma = false;
          for (int j = 0; j < tax_levels; j++)
            {
              if ((best_level_len[j] > 0) &&
                  (1.0 * level_match[j] / count >= opt_sintax_cutoff))
                {
                  fprintf(fp_tabbedout,
                          "%s%c:%.*s",
                          (comma ? "," : ""),
                          tax_letters[j],
                          best_level_len[j],
                          best_h + best_level_start[j]);
                  comma = true;
                }
            }
        }
    }
  else
    {
      if (opt_sintax_cutoff > 0.0)
        fprintf(fp_tabbedout, "\t\t\t");
      else
        fprintf(fp_tabbedout, "\t\t");
    }

#if 0
  fprintf(fp_tabbedout, "\t%d\t%d", best_count, count);
#endif

  fprintf(fp_tabbedout, "\n");
  xpthread_mutex_unlock(&mutex_output);
}

void sintax_query(int64_t t)
{
  int all_seqno[2][bootstrap_count];
  int best_seqno[2];
  int boot_count[2];
  unsigned int best_count[2];

  best_count[0] = 0;
  best_count[1] = 0;
  best_seqno[0] = 0;
  best_seqno[1] = 0;
  boot_count[0] = 0;
  boot_count[1] = 0;

  int qseqlen = si_plus[t].qseqlen;
  char * query_head = si_plus[t].query_head;

  bitmap_t * b = bitmap_init(qseqlen);

  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = s ? si_minus+t : si_plus+t;

      /* perform search */

      unsigned int kmersamplecount;
      unsigned int * kmersample;

      /* find unique kmers */
      unique_count(si->uh, opt_wordlength,
                   si->qseqlen, si->qsequence,
                   & kmersamplecount, & kmersample, MASK_NONE);

      /* perform 100 bootstraps */

      if (kmersamplecount >= subset_size)
        {
          for (int i = 0; i < bootstrap_count ; i++)
            {
              /* subsample 32 kmers */
              unsigned int kmersample_subset[subset_size];
              int subsamples = 0;
              bitmap_reset_all(b);
              for(int j = 0; j < subset_size ; j++)
                {
                  int64_t x = random_int(kmersamplecount);
                  if (! bitmap_get(b, x))
                    {
                      kmersample_subset[subsamples++] = kmersample[x];
                      bitmap_set(b, x);
                    }
                }

              si->kmersamplecount = subsamples;
              si->kmersample = kmersample_subset;

              search_topscores(si);

              while(!minheap_isempty(si->m))
                {
                  elem_t e = minheap_poplast(si->m);

                  all_seqno[s][boot_count[s]++] = e.seqno;

                  if (e.count > best_count[s])
                    {
                      best_count[s] = e.count;
                      best_seqno[s] = e.seqno;
                    }
                }
            }
        }
    }

  int best_strand;

  if (opt_strand == 1)
    best_strand = 0;
  else
    {
      if (best_count[0] > best_count[1])
        best_strand = 0;
      else if (best_count[1] > best_count[0])
        best_strand = 1;
      else
        {
          if (boot_count[0] >= boot_count[1])
            best_strand = 0;
          else
            best_strand = 1;
        }
    }

  sintax_analyse(query_head,
                 best_strand,
                 best_seqno[best_strand],
                 best_count[best_strand],
                 all_seqno[best_strand],
                 boot_count[best_strand]);

  bitmap_free(b);
}

void sintax_thread_run(int64_t t)
{
  while (1)
    {
      xpthread_mutex_lock(&mutex_input);

      if (fastx_next(query_fastx_h,
                     ! opt_notrunclabels,
                     chrmap_no_change))
        {
          char * qhead = fastx_get_header(query_fastx_h);
          int query_head_len = fastx_get_header_length(query_fastx_h);
          char * qseq = fastx_get_sequence(query_fastx_h);
          int qseqlen = fastx_get_sequence_length(query_fastx_h);
          int query_no = fastx_get_seqno(query_fastx_h);
          int qsize = fastx_get_abundance(query_fastx_h);

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
          uint64_t progress = fastx_get_position(query_fastx_h);

          /* let other threads read input */
          xpthread_mutex_unlock(&mutex_input);

          /* minus strand: copy header and reverse complementary sequence */
          if (opt_strand > 1)
            {
              strcpy(si_minus[t].query_head, si_plus[t].query_head);
              reverse_complement(si_minus[t].qsequence,
                                 si_plus[t].qsequence,
                                 si_plus[t].qseqlen);
            }

          sintax_query(t);

          /* lock mutex for update of global data and output */
          xpthread_mutex_lock(&mutex_output);

          /* show progress */
          progress_update(progress);

          xpthread_mutex_unlock(&mutex_output);
        }
      else
        {
          xpthread_mutex_unlock(&mutex_input);
          break;
        }
    }
}

void sintax_thread_init(struct searchinfo_s * si)
{
  /* thread specific initialiation */
  si->uh = unique_init();
  si->kmers = (count_t *) xmalloc(seqcount * sizeof(count_t) + 32);
  si->m = minheap_init(tophits);
  si->hits = 0;
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = 0;
  si->seq_alloc = 0;
  si->qsequence = 0;
  si->nw = 0;
  si->s = 0;
}

void sintax_thread_exit(struct searchinfo_s * si)
{
  /* thread specific clean up */
  unique_exit(si->uh);
  minheap_exit(si->m);
  xfree(si->kmers);
  if (si->query_head)
    xfree(si->query_head);
  if (si->qsequence)
    xfree(si->qsequence);
}

void * sintax_thread_worker(void * vp)
{
  int64_t t = (int64_t) vp;
  sintax_thread_run(t);
  return 0;
}

void sintax_thread_worker_run()
{
  /* initialize threads, start them, join them and return */

  xpthread_attr_init(&attr);
  xpthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* init and create worker threads, put them into stand-by mode */
  for(int t=0; t<opt_threads; t++)
    {
      sintax_thread_init(si_plus+t);
      if (si_minus)
        sintax_thread_init(si_minus+t);
      xpthread_create(pthread+t, &attr,
                      sintax_thread_worker, (void*)(int64_t)t);
    }

  /* finish and clean up worker threads */
  for(int t=0; t<opt_threads; t++)
    {
      xpthread_join(pthread[t], NULL);
      sintax_thread_exit(si_plus+t);
      if (si_minus)
        sintax_thread_exit(si_minus+t);
    }

  xpthread_attr_destroy(&attr);
}

void sintax()
{
  /* tophits = the maximum number of hits we need to store */

  tophits = 1;

  /* open output files */

  if (! opt_db)
    fatal("No database file specified with --db");

  if (opt_tabbedout)
    {
      fp_tabbedout = fopen_output(opt_tabbedout);
      if (! fp_tabbedout)
        fatal("Unable to open tabbedout output file for writing");
    }
  else
    fatal("No output file specified with --tabbedout");

  /* check if db may be an UDB file */

  bool is_udb = udb_detect_isudb(opt_db);

  if (is_udb)
    udb_read(opt_db, 1, 1);
  else
    db_read(opt_db, 0);

  seqcount = db_getsequencecount();

  if (!is_udb)
    {
      dbindex_prepare(1, opt_dbmask);
      dbindex_addallsequences(opt_dbmask);
    }

  /* prepare reading of queries */

  query_fastx_h = fastx_open(opt_sintax);

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
  xpthread_mutex_init(&mutex_input, NULL);
  xpthread_mutex_init(&mutex_output, NULL);

  /* run */

  progress_init("Classifying sequences", fastx_get_size(query_fastx_h));
  sintax_thread_worker_run();
  progress_done();

  if (! opt_quiet)
    fprintf(stderr, "Classified %d of %d sequences (%.2f%%)\n",
            classified, queries, 100.0 * classified / queries);

  if (opt_log)
    fprintf(fp_log, "Classified %d of %d sequences (%.2f%%)\n",
            classified, queries, 100.0 * classified / queries);

  /* clean up */

  xpthread_mutex_destroy(&mutex_output);
  xpthread_mutex_destroy(&mutex_input);

  xfree(pthread);
  xfree(si_plus);
  if (si_minus)
    xfree(si_minus);

  fastx_close(query_fastx_h);
  fclose(fp_tabbedout);

  dbindex_free();
  db_free();
}
