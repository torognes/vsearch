/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

  Implements the Sintax algorithm as described in Robert Edgar's preprint:

  Robert Edgar (2016)
  SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences
  BioRxiv, 074161
  doi: https://doi.org/10.1101/074161

  Further details:

  https://www.drive5.com/usearch/manual/cmd_sintax.html


  Note that due to the lack of details in the description, this implementation
  in vsearch is surely somewhat different from the one in usearch.

*/

#include "vsearch.h"
#include "minheap.h"
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t
#include <cstring>  // std::memset, std::strncmp, std::strcpy
#include <pthread.h>


static struct searchinfo_s * si_plus;
static struct searchinfo_s * si_minus;
static pthread_t * pthread;

/* global constants/data, no need for synchronization */
static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */
static pthread_attr_t attr;
static fastx_handle query_fastx_h;

const int subset_size = 32;
const int bootstrap_count = 100;

/* global data protected by mutex */
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static FILE * fp_tabbedout;
static int queries = 0;
static int classified = 0;


auto sintax_analyse(char * query_head,
                    int strand,
                    int * all_seqno,
                    int count) -> void
{
  int level_match[tax_levels];
  int level_best[tax_levels];

  char * cand_level_start[bootstrap_count][tax_levels];
  int cand_level_len[bootstrap_count][tax_levels];
  int cand_level_matches[bootstrap_count][tax_levels];

  /* Check number of successful bootstraps, must be at least half */

  bool const enough = count >= (bootstrap_count + 1) / 2;

  if (enough)
    {
      /* Find the most common name at each taxonomic rank,
         but with the same names at higher ranks. */

      for (int i = 0; i < count ; i++)
        {
          /* Split headers of all candidates by taxonomy ranks */

          int const seqno = all_seqno[i];
          int new_level_start[tax_levels];
          int new_level_len[tax_levels];
          tax_split(seqno, new_level_start, new_level_len);
          for (int k = 0; k < tax_levels; k++)
            {
              cand_level_start[i][k] = db_getheader(seqno)+new_level_start[k];
              cand_level_len[i][k] = new_level_len[k];
            }
        }

      /* Count matching names among candidates */

      for (int i = 0; i < count ; i++) {
        for (int k = 0; k < tax_levels; k++) {
          cand_level_matches[i][k] = 0;
        }
      }

      for (int k = 0; k < tax_levels; k++) {
        for (int i = 0; i < count ; i++) {
          for (int c = 0; c <= i ; c++)
            {
              /* check match at current and all higher levels */
              bool match = true;
              for (int j = 0; j <= k; j++)
                {
                  if ((cand_level_len[i][j] != cand_level_len[c][j]) ||
                      (strncmp(cand_level_start[i][j],
                               cand_level_start[c][j],
                               cand_level_len[i][j]) != 0))
                    {
                      match = false;
                      break;
                    }
                }
              if (match)
                {
                  cand_level_matches[c][k]++;
                  break; /* stop at first match */
                }
            }
        }
      }

      /* Find most common name at each taxonomic level */

      for (int k = 0; k < tax_levels; k++)
        {
          level_best[k] = -1;
          level_match[k] = 0;
          int m = 0;
          for (int i = 0; i < count ; i++)
            {
              m += cand_level_matches[i][k];
              if (cand_level_matches[i][k] > level_match[k])
                {
                  level_best[k] = i;
                  level_match[k] = cand_level_matches[i][k];
                }
              if (m >= count) {
                break;
              }
            }
        }
    }

  /* write to tabbedout file */
  xpthread_mutex_lock(&mutex_output);
  fprintf(fp_tabbedout, "%s\t", query_head);

  queries++;

  if (enough)
    {
      classified++;

      bool comma = false;
      for (int j = 0; j < tax_levels; j++)
        {
          int const best = level_best[j];
          if (cand_level_len[best][j] > 0)
            {
              fprintf(fp_tabbedout,
                      "%s%c:%.*s(%.2f)",
                      (comma ? "," : ""),
                      tax_letters[j],
                      cand_level_len[best][j],
                      cand_level_start[best][j],
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
              int const best = level_best[j];
              if ((cand_level_len[best][j] > 0) &&
                  (1.0 * level_match[j] / count >= opt_sintax_cutoff))
                {
                  fprintf(fp_tabbedout,
                          "%s%c:%.*s",
                          (comma ? "," : ""),
                          tax_letters[j],
                          cand_level_len[best][j],
                          cand_level_start[best][j]);
                  comma = true;
                }
            }
        }
    }
  else
    {
      if (opt_sintax_cutoff > 0.0)
        {
          fprintf(fp_tabbedout, "\t\t");
        }
      else
        {
          fprintf(fp_tabbedout, "\t");
        }
    }

  fprintf(fp_tabbedout, "\n");
  xpthread_mutex_unlock(&mutex_output);
}

auto sintax_search_topscores(struct searchinfo_s * si) -> void
{
  /*
    Count the number of kmer hits in each database sequence and select
    the database sequence with the highest number of matching kmers.
    If several sequences have equally many kmer matches, choose one of
    them according to the following rules: By default, choose the
    shortest. If two are equally short, choose the one that comes
    first in the database.  If the sintax_random option is in effect,
    ties will instead be chosen randomly.
  */

  /* count kmer hits in the database sequences */
  const int indexed_count = dbindex_getcount();

  /* zero counts */
  memset(si->kmers, 0, indexed_count * sizeof(count_t));

  for(unsigned int i = 0; i < si->kmersamplecount; i++)
    {
      unsigned int const kmer = si->kmersample[i];
      unsigned char * bitmap = dbindex_getbitmap(kmer);

      if (bitmap)
        {
#ifdef __x86_64__
          if (ssse3_present)
            {
              increment_counters_from_bitmap_ssse3(si->kmers,
                                                   bitmap, indexed_count);
            }
          else
            {
              increment_counters_from_bitmap_sse2(si->kmers,
                                                  bitmap, indexed_count);
            }
#else
          increment_counters_from_bitmap(si->kmers, bitmap, indexed_count);
#endif
        }
      else
        {
          unsigned int * list = dbindex_getmatchlist(kmer);
          unsigned int const count = dbindex_getmatchcount(kmer);
          for(unsigned int j = 0; j < count; j++)
            {
              si->kmers[list[j]]++;
            }
        }
    }

  unsigned int tophits = 0;

  elem_t best;
  best.count = 0;
  best.seqno = 0;
  best.length = 0;

  for (int i = 0; i < indexed_count; i++)
    {
      count_t const count = si->kmers[i];
      unsigned int const seqno = dbindex_getmapping(i);
      unsigned int const length = db_getsequencelen(best.seqno);

      if (count > best.count)
        {
          best.count = count;
          best.seqno = seqno;
          best.length = length;
          tophits = 1;
        }
      else if (count == best.count)
        {
          if (opt_sintax_random)
            {
              tophits++;
              if (random_int(tophits) == 0)
                {
                  best.seqno = seqno;
                  best.length = length;
                }
            }
          else
            {
              if (length < best.length)
                {
                  best.seqno = seqno;
                  best.length = length;
                }
              else if (length == best.length)
                {
                  if (seqno < best.seqno)
                    {
                      best.seqno = seqno;
                    }
                }
            }
        }
    }

  minheap_empty(si->m);
  if (best.count > 1) {
    minheap_add(si->m, &best);
  }
}

auto sintax_query(int64_t t) -> void
{
  int all_seqno[2][bootstrap_count];
  int boot_count[2] = {0, 0};
  unsigned int best_count[2] = {0, 0};
  int const qseqlen = si_plus[t].qseqlen;
  char * query_head = si_plus[t].query_head;

  bitmap_t * b = bitmap_init(qseqlen);

  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = s ? si_minus+t : si_plus+t;

      /* perform search */

      unsigned int kmersamplecount = 0;
      unsigned int * kmersample = nullptr;

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
              for (int j = 0; j < subset_size ; j++)
                {
                  int64_t const x = random_int(kmersamplecount);
                  if (! bitmap_get(b, x))
                    {
                      kmersample_subset[subsamples++] = kmersample[x];
                      bitmap_set(b, x);
                    }
                }

              si->kmersamplecount = subsamples;
              si->kmersample = kmersample_subset;

              sintax_search_topscores(si);

              if (! minheap_isempty(si->m))
                {
                  elem_t const e = minheap_poplast(si->m);

                  all_seqno[s][boot_count[s]++] = e.seqno;

                  if (e.count > best_count[s])
                    {
                      best_count[s] = e.count;
                    }
                }
            }
        }
    }

  int best_strand = 0;

  if (opt_strand == 1)
    {
      best_strand = 0;
    }
  else
    {
      if (best_count[0] > best_count[1])
        {
          best_strand = 0;
        }
      else if (best_count[1] > best_count[0])
        {
          best_strand = 1;
        }
      else
        {
          if (boot_count[0] >= boot_count[1])
            {
              best_strand = 0;
            }
          else
            {
              best_strand = 1;
            }
        }
    }

  sintax_analyse(query_head,
                 best_strand,
                 all_seqno[best_strand],
                 boot_count[best_strand]);

  bitmap_free(b);
}

auto sintax_thread_run(int64_t t) -> void
{
  while (true)
    {
      xpthread_mutex_lock(&mutex_input);

      if (fastx_next(query_fastx_h,
                     ! opt_notrunclabels,
                     chrmap_no_change))
        {
          char * qhead = fastx_get_header(query_fastx_h);
          int const query_head_len = fastx_get_header_length(query_fastx_h);
          char * qseq = fastx_get_sequence(query_fastx_h);
          int const qseqlen = fastx_get_sequence_length(query_fastx_h);
          int const query_no = fastx_get_seqno(query_fastx_h);
          int const qsize = fastx_get_abundance(query_fastx_h);

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
                  si->query_head = (char *)
                    xrealloc(si->query_head, (size_t)(si->query_head_alloc));
                }

              if (si->qseqlen + 1 > si->seq_alloc)
                {
                  si->seq_alloc = si->qseqlen + 2001;
                  si->qsequence = (char *)
                    xrealloc(si->qsequence, (size_t)(si->seq_alloc));
                }
            }

          /* plus strand: copy header and sequence */
          strcpy(si_plus[t].query_head, qhead);
          strcpy(si_plus[t].qsequence, qseq);

          /* get progress as amount of input file read */
          uint64_t const progress = fastx_get_position(query_fastx_h);

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

auto sintax_thread_init(struct searchinfo_s * si) -> void
{
  /* thread specific initialiation */
  si->uh = unique_init();
  si->kmers = (count_t *) xmalloc((seqcount * sizeof(count_t)) + 32);
  si->m = minheap_init(tophits);
  si->hits = nullptr;
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = nullptr;
  si->seq_alloc = 0;
  si->qsequence = nullptr;
  si->nw = nullptr;
  si->s = nullptr;
}

auto sintax_thread_exit(struct searchinfo_s * si) -> void
{
  /* thread specific clean up */
  unique_exit(si->uh);
  minheap_exit(si->m);
  xfree(si->kmers);
  if (si->query_head)
    {
      xfree(si->query_head);
    }
  if (si->qsequence)
    {
      xfree(si->qsequence);
    }
}

auto sintax_thread_worker(void * vp) -> void *
{
  auto t = (int64_t) vp;
  sintax_thread_run(t);
  return nullptr;
}

auto sintax_thread_worker_run() -> void
{
  /* initialize threads, start them, join them and return */

  xpthread_attr_init(&attr);
  xpthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* init and create worker threads, put them into stand-by mode */
  for (int t = 0; t < opt_threads; t++)
    {
      sintax_thread_init(si_plus + t);
      if (si_minus)
        {
          sintax_thread_init(si_minus + t);
        }
      xpthread_create(pthread + t, &attr,
                      sintax_thread_worker, (void *) (int64_t)t);
    }

  /* finish and clean up worker threads */
  for (int t = 0; t < opt_threads; t++)
    {
      xpthread_join(pthread[t], nullptr);
      sintax_thread_exit(si_plus + t);
      if (si_minus)
        {
          sintax_thread_exit(si_minus + t);
        }
    }

  xpthread_attr_destroy(&attr);
}

auto sintax() -> void
{
  /* tophits = the maximum number of hits we need to store */

  tophits = 1;

  /* open output files */

  if (! opt_db)
    {
      fatal("No database file specified with --db");
    }

  if (opt_tabbedout)
    {
      fp_tabbedout = fopen_output(opt_tabbedout);
      if (! fp_tabbedout)
        {
          fatal("Unable to open tabbedout output file for writing");
        }
    }
  else
    {
      fatal("No output file specified with --tabbedout");
    }

  /* check if db may be an UDB file */

  bool const is_udb = udb_detect_isudb(opt_db);

  if (is_udb)
    {
      udb_read(opt_db, true, true);
    }
  else
    {
      db_read(opt_db, 0);
    }

  seqcount = db_getsequencecount();

  if (! is_udb)
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
    {
      si_minus = (struct searchinfo_s *) xmalloc(opt_threads *
                                                 sizeof(struct searchinfo_s));
    }
  else
    {
      si_minus = nullptr;
    }

  pthread = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));

  /* init mutexes for input and output */
  xpthread_mutex_init(&mutex_input, nullptr);
  xpthread_mutex_init(&mutex_output, nullptr);

  /* run */

  progress_init("Classifying sequences", fastx_get_size(query_fastx_h));
  sintax_thread_worker_run();
  progress_done();

  if (! opt_quiet)
    {
      fprintf(stderr, "Classified %d of %d sequences", classified, queries);
      if (queries > 0)
        {
          fprintf(stderr, " (%.2f%%)", 100.0 * classified / queries);
        }
      fprintf(stderr, "\n");
    }

  if (opt_log)
    {
      fprintf(fp_log, "Classified %d of %d sequences", classified, queries);
      if (queries > 0)
        {
          fprintf(fp_log, " (%.2f%%)", 100.0 * classified / queries);
        }
      fprintf(fp_log, "\n");
    }

  /* clean up */

  xpthread_mutex_destroy(&mutex_output);
  xpthread_mutex_destroy(&mutex_input);

  xfree(pthread);
  xfree(si_plus);
  if (si_minus)
    {
      xfree(si_minus);
    }

  fastx_close(query_fastx_h);
  fclose(fp_tabbedout);

  dbindex_free();
  db_free();
}
