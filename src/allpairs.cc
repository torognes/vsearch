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

#include "vsearch.h"
#include "align_simd.h"
#include "mask.h"
#include <algorithm>  // std::min, std::max
#include <cstdint>  // int64_t
#include <cstdio>  // std::fprintf, std::FILE, std:fclose, std::size_t
#include <cstdlib>  // std::qsort
#include <cstring>  // std::strlen
#include <limits>
#include <pthread.h>
#include <vector>


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
static FILE * fp_alnout = nullptr;
static FILE * fp_samout = nullptr;
static FILE * fp_userout = nullptr;
static FILE * fp_blast6out = nullptr;
static FILE * fp_uc = nullptr;
static FILE * fp_fastapairs = nullptr;
static FILE * fp_matched = nullptr;
static FILE * fp_notmatched = nullptr;
static FILE * fp_qsegout = nullptr;
static FILE * fp_tsegout = nullptr;

static int count_matched = 0;
static int count_notmatched = 0;

inline auto allpairs_hit_compare_typed(struct hit * x, struct hit * y) -> int
{
  // high id, then low id
  // early target, then late target

  if (x->id > y->id)
    {
      return -1;
    }
  else if (x->id < y->id)
    {
      return +1;
    }
  else if (x->target < y->target)
    {
      return -1;
    }
  else if (x->target > y->target)
    {
      return +1;
    }
  else
    {
      return 0;
    }
}

auto allpairs_hit_compare(const void * a, const void * b) -> int
{
  return allpairs_hit_compare_typed((struct hit *) a, (struct hit *) b);
}

auto allpairs_output_results(int hit_count,
                             struct hit * hits,
                             char * query_head,
                             int qseqlen,
                             char * qsequence,
                             char * qsequence_rc) -> void
{
  /* show results */
  auto const toreport = std::min(opt_maxhits, static_cast<int64_t>(hit_count));

  if (fp_alnout)
    {
      results_show_alnout(fp_alnout,
                          hits,
                          toreport,
                          query_head,
                          qsequence,
                          qseqlen);
    }

  if (fp_samout)
    {
      results_show_samout(fp_samout,
                          hits,
                          toreport,
                          query_head,
                          qsequence,
                          qsequence_rc);
    }

  if (toreport)
    {
      double const top_hit_id = hits[0].id;

      for (int t = 0; t < toreport; t++)
        {
          struct hit * hp = hits + t;

          if (opt_top_hits_only and (hp->id < top_hit_id))
            {
              break;
            }

          if (fp_fastapairs)
            {
              results_show_fastapairs_one(fp_fastapairs,
                                          hp,
                                          query_head,
                                          qsequence,
                                          qsequence_rc);
            }

          if (fp_qsegout)
            {
              results_show_qsegout_one(fp_qsegout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (fp_tsegout)
            {
              results_show_tsegout_one(fp_tsegout,
                                       hp);
            }

          if (fp_uc)
            {
              if ((t == 0) or opt_uc_allhits)
                {
                  results_show_uc_one(fp_uc,
                                      hp,
                                      query_head,
                                      qseqlen,
                                      hp->target);
                }
            }

          if (fp_userout)
            {
              results_show_userout_one(fp_userout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (fp_blast6out)
            {
              results_show_blast6out_one(fp_blast6out,
                                         hp,
                                         query_head,
                                         qseqlen);
            }
        }
    }
  else
    {
      if (fp_uc)
        {
          results_show_uc_one(fp_uc,
                              nullptr,
                              query_head,
                              qseqlen,
                              0);
        }

      if (opt_output_no_hits)
        {
          if (fp_userout)
            {
              results_show_userout_one(fp_userout,
                                       nullptr,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (fp_blast6out)
            {
              results_show_blast6out_one(fp_blast6out,
                                         nullptr,
                                         query_head,
                                         qseqlen);
            }
        }
    }

  if (hit_count)
    {
      ++count_matched;
      if (opt_matched)
        {
          fasta_print_general(fp_matched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              strlen(query_head),
                              0,
                              count_matched,
                              -1.0,
                              -1, -1, nullptr, 0.0);
        }
    }
  else
    {
      ++count_notmatched;
      if (opt_notmatched)
        {
          fasta_print_general(fp_notmatched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              strlen(query_head),
                              0,
                              count_notmatched,
                              -1.0,
                              -1, -1, nullptr, 0.0);
        }
    }
}

auto allpairs_thread_run(int64_t t) -> void
{
  (void) t;

  struct searchinfo_s searchinfo;

  struct searchinfo_s * si = & searchinfo;
  searchinfo.hits_v.resize(seqcount);
  searchinfo.hits = searchinfo.hits_v.data();

  searchinfo.s = search16_init(opt_match,
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
  auto const maxhits = static_cast<std::size_t>(seqcount);
  std::vector<unsigned int> pseqnos(maxhits);
  std::vector<CELL> pscores(maxhits);
  std::vector<unsigned short> paligned(maxhits);
  std::vector<unsigned short> pmatches(maxhits);
  std::vector<unsigned short> pmismatches(maxhits);
  std::vector<unsigned short> pgaps(maxhits);
  std::vector<char *> pcigar(maxhits);
  std::vector<struct hit> finalhits(maxhits);

  auto cont = true;

  while (cont)
    {
      xpthread_mutex_lock(&mutex_input);

      int const query_no = queries;

      if (query_no < seqcount)
        {
          ++queries;

          /* let other threads read input */
          xpthread_mutex_unlock(&mutex_input);

          /* init search info */
          searchinfo.query_no = query_no;
          searchinfo.qsize = db_getabundance(query_no);
          searchinfo.query_head_len = db_getheaderlen(query_no);
          searchinfo.query_head = db_getheader(query_no);
          searchinfo.qseqlen = db_getsequencelen(query_no);
          searchinfo.qsequence = db_getsequence(query_no);
          searchinfo.rejects = 0;
          searchinfo.accepts = 0;
          searchinfo.hit_count = 0;

          for (int target = searchinfo.query_no + 1; target < seqcount; target++)
            {
              if (opt_acceptall or search_acceptable_unaligned(si, target))
                {
                  pseqnos[searchinfo.hit_count] = target;
                  ++searchinfo.hit_count;
                }
            }

          if (searchinfo.hit_count)
            {
              /* perform alignments */

              search16_qprep(searchinfo.s, searchinfo.qsequence, searchinfo.qseqlen);

              search16(searchinfo.s,
                       searchinfo.hit_count,
                       pseqnos.data(),
                       pscores.data(),
                       paligned.data(),
                       pmatches.data(),
                       pmismatches.data(),
                       pgaps.data(),
                       pcigar.data());

              /* convert to hit structure */
              for (int h = 0; h < searchinfo.hit_count; h++)
                {
                  struct hit * hit = &searchinfo.hits_v[h];

                  unsigned int const target = pseqnos[h];
                  int64_t nwscore = pscores[h];

                  char * nwcigar {nullptr};
                  int64_t nwalignmentlength {0};
                  int64_t nwmatches {0};
                  int64_t nwmismatches {0};
                  int64_t nwgaps {0};

                  if (nwscore == std::numeric_limits<short>::max())
                    {
                      /* In case the SIMD aligner cannot align,
                         perform a new alignment with the
                         linear memory aligner */

                      char * tseq = db_getsequence(target);
                      int64_t const tseqlen = db_getsequencelen(target);

                      if (pcigar[h] != nullptr)
                        {
                          xfree(pcigar[h]);
                        }

                      nwcigar = xstrdup(lma.align(searchinfo.qsequence,
                                                  tseq,
                                                  searchinfo.qseqlen,
                                                  tseqlen));
                      lma.alignstats(nwcigar,
                                     searchinfo.qsequence,
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

                  hit->accepted = false;
                  hit->rejected = false;
                  hit->aligned = true;
                  hit->weak = false;

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

                  auto const dseqlen = static_cast<int>(db_getsequencelen(target));
                  hit->shortest = std::min(searchinfo.qseqlen, dseqlen);
                  hit->longest = std::max(searchinfo.qseqlen, dseqlen);

                  /* trim alignment, compute numbers excluding terminal gaps */
                  align_trim(hit);

                  /* test accept/reject criteria after alignment */
                  if (opt_acceptall or search_acceptable_aligned(si, hit))
                    {
                      finalhits[searchinfo.accepts] = *hit;
                      ++searchinfo.accepts;
                    }
                }

              /* sort hits */
              qsort(finalhits.data(), searchinfo.accepts,
                    sizeof(struct hit), allpairs_hit_compare);
            }

          /* lock mutex for update of global data and output */
          xpthread_mutex_lock(&mutex_output);

          /* output results */
          allpairs_output_results(searchinfo.accepts,
                                  finalhits.data(),
                                  searchinfo.query_head,
                                  searchinfo.qseqlen,
                                  searchinfo.qsequence,
                                  nullptr);

          /* update stats */
          if (searchinfo.accepts)
            {
              ++qmatches;
            }

          /* show progress */
          progress += seqcount - query_no - 1;
          progress_update(progress);

          xpthread_mutex_unlock(&mutex_output);

          /* free memory for alignment strings */
          for (int i = 0; i < searchinfo.hit_count; i++)
            {
              if (searchinfo.hits_v[i].aligned)
                {
                  xfree(searchinfo.hits_v[i].nwalignment);
                }
            }
        }
      else
        {
          /* let other threads read input */
          xpthread_mutex_unlock(&mutex_input);

          cont = false;
        }
    }

  search16_exit(searchinfo.s);

  xfree(scorematrix);
}

auto allpairs_thread_worker(void * void_ptr) -> void *
{
  auto const nth_thread = reinterpret_cast<int64_t>(void_ptr);
  allpairs_thread_run(nth_thread);
  return nullptr;
}

auto allpairs_thread_worker_run() -> void
{
  /* initialize threads, start them, join them and return */

  xpthread_attr_init(&attr);
  xpthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* init and create worker threads, put them into stand-by mode */
  for (int t = 0; t < opt_threads; t++)
    {
      xpthread_create(pthread + t, &attr,
                      allpairs_thread_worker, (void *) (int64_t) t);
    }

  /* finish and clean up worker threads */
  for (int t = 0; t < opt_threads; t++)
    {
      xpthread_join(pthread[t], nullptr);
    }

  xpthread_attr_destroy(&attr);
}


auto allpairs_global(char * cmdline, char * progheader) -> void
{
  opt_strand = 1;
  opt_uc_allhits = 1;

  /* open output files */

  if (opt_alnout)
    {
      fp_alnout = fopen_output(opt_alnout);
      if (not fp_alnout)
        {
          fatal("Unable to open alignment output file for writing");
        }

      fprintf(fp_alnout, "%s\n", cmdline);
      fprintf(fp_alnout, "%s\n", progheader);
    }

  if (opt_samout)
    {
      fp_samout = fopen_output(opt_samout);
      if (not fp_samout)
        {
          fatal("Unable to open SAM output file for writing");
        }
    }

  if (opt_userout)
    {
      fp_userout = fopen_output(opt_userout);
      if (not fp_userout)
        {
          fatal("Unable to open user-defined output file for writing");
        }
    }

  if (opt_blast6out)
    {
      fp_blast6out = fopen_output(opt_blast6out);
      if (not fp_blast6out)
        {
          fatal("Unable to open blast6-like output file for writing");
        }
    }

  if (opt_uc)
    {
      fp_uc = fopen_output(opt_uc);
      if (not fp_uc)
        {
          fatal("Unable to open uc output file for writing");
        }
    }

  if (opt_fastapairs)
    {
      fp_fastapairs = fopen_output(opt_fastapairs);
      if (not fp_fastapairs)
        {
          fatal("Unable to open fastapairs output file for writing");
        }
    }

  if (opt_qsegout)
    {
      fp_qsegout = fopen_output(opt_qsegout);
      if (not fp_qsegout)
        {
          fatal("Unable to open qsegout output file for writing");
        }
    }

  if (opt_tsegout)
    {
      fp_tsegout = fopen_output(opt_tsegout);
      if (not fp_tsegout)
        {
          fatal("Unable to open tsegout output file for writing");
        }
    }

  if (opt_matched)
    {
      fp_matched = fopen_output(opt_matched);
      if (not fp_matched)
        {
          fatal("Unable to open matched output file for writing");
        }
    }

  if (opt_notmatched)
    {
      fp_notmatched = fopen_output(opt_notmatched);
      if (not fp_notmatched)
        {
          fatal("Unable to open notmatched output file for writing");
        }
    }

  db_read(opt_allpairs_global, 0);

  results_show_samheader(fp_samout, cmdline, opt_allpairs_global);

  if (opt_qmask == MASK_DUST)
    {
      dust_all();
    }
  else if ((opt_qmask == MASK_SOFT) and (opt_hardmask))
    {
      hardmask_all();
    }

  show_rusage();

  seqcount = db_getsequencecount();

  /* prepare reading of queries */
  qmatches = 0;
  queries = 0;

  std::vector<pthread_t> pthread_v(opt_threads);
  pthread = pthread_v.data();

  /* init mutexes for input and output */
  xpthread_mutex_init(&mutex_input, nullptr);
  xpthread_mutex_init(&mutex_output, nullptr);

  progress = 0;
  progress_init("Aligning", MAX(0, ((int64_t) seqcount) * ((int64_t) seqcount - 1)) / 2);  // refactoring: issue with parenthesis?
  allpairs_thread_worker_run();
  progress_done();

  if (not opt_quiet)
    {
      fprintf(stderr, "Matching query sequences: %d of %d",
              qmatches, queries);
      if (queries > 0)
        {
          fprintf(stderr, " (%.2f%%)", 100.0 * qmatches / queries);
        }
      fprintf(stderr, "\n");
    }

  if (opt_log)
    {
      fprintf(fp_log, "Matching query sequences: %d of %d",
              qmatches, queries);
      if (queries > 0)
        {
          fprintf(fp_log, " (%.2f%%)", 100.0 * qmatches / queries);
        }
      fprintf(fp_log, "\n\n");
    }

  xpthread_mutex_destroy(&mutex_output);
  xpthread_mutex_destroy(&mutex_input);

  // pthread_v not used after this point

  /* clean up, global */
  db_free();
  if (opt_matched)
    {
      fclose(fp_matched);
    }
  if (opt_notmatched)
    {
      fclose(fp_notmatched);
    }
  if (opt_fastapairs)
    {
      fclose(fp_fastapairs);
    }
  if (opt_qsegout)
    {
      fclose(fp_qsegout);
    }
  if (opt_tsegout)
    {
      fclose(fp_tsegout);
    }
  if (fp_uc)
    {
      fclose(fp_uc);
    }
  if (fp_blast6out)
    {
      fclose(fp_blast6out);
    }
  if (fp_userout)
    {
      fclose(fp_userout);
    }
  if (fp_alnout)
    {
      fclose(fp_alnout);
    }
  if (fp_samout)
    {
      fclose(fp_samout);
    }
  show_rusage();
}
