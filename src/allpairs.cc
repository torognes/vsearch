/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "linmemalign.h"
#include "mask.h"
#include "utils/fatal.hpp"
#include "utils/threads.hpp"
#include <algorithm>  // std::min, std::max
#include <cstdint>  // int64_t
#include <cstdio>  // std::fprintf, std::FILE, std:fclose, std::size_t
#include <cstdlib>  // std::qsort
#include <cstring>  // std::strlen
#include <limits>
#include <mutex>  // std::mutex, std::lock_guard, std::unique_lock
#include <vector>


/* global constants/data, no need for synchronization */
static int seqcount; /* number of database sequences */

/* global data protected by mutex */
static std::mutex mutex_input;
static std::mutex mutex_output;
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


inline auto allpairs_hit_compare_typed(struct hit const * lhs, struct hit const * rhs) -> int
{
  // high id, then low id
  // early target, then late target

  if (lhs->id > rhs->id)
    {
      return -1;
    }
  if (lhs->id < rhs->id)
    {
      return +1;
    }
  if (lhs->target < rhs->target)
    {
      return -1;
    }
  if (lhs->target > rhs->target)
    {
      return +1;
    }
  return 0;
}


auto allpairs_hit_compare(const void * lhs, const void * rhs) -> int
{
  return allpairs_hit_compare_typed(static_cast<struct hit const *>(lhs), static_cast<struct hit const *>(rhs));
}


auto allpairs_output_results(int hit_count,
                             struct hit * hits,
                             char const * query_head,
                             int qseqlen,
                             char const * qsequence,
                             char const * qsequence_rc) -> void
{
  /* show results */
  auto const toreport = std::min(opt_maxhits, static_cast<int64_t>(hit_count));

  if (fp_alnout != nullptr)
    {
      results_show_alnout(fp_alnout,
                          hits,
                          static_cast<int>(toreport),
                          query_head,
                          qsequence,
                          qseqlen);
    }

  if (fp_samout != nullptr)
    {
      results_show_samout(fp_samout,
                          hits,
                          static_cast<int>(toreport),
                          query_head,
                          qsequence,
                          qsequence_rc);
    }

  if (toreport != 0)
    {
      double const top_hit_id = hits[0].id;

      for (int t = 0; t < toreport; t++)
        {
          struct hit const * hp = hits + t;

          if ((opt_top_hits_only != 0) and (hp->id < top_hit_id))
            {
              break;
            }

          if (fp_fastapairs != nullptr)
            {
              results_show_fastapairs_one(fp_fastapairs,
                                          hp,
                                          query_head,
                                          qsequence,
                                          qsequence_rc);
            }

          if (fp_qsegout != nullptr)
            {
              results_show_qsegout_one(fp_qsegout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (fp_tsegout != nullptr)
            {
              results_show_tsegout_one(fp_tsegout,
                                       hp);
            }

          if (fp_uc != nullptr)
            {
              if ((t == 0) or (opt_uc_allhits != 0))
                {
                  results_show_uc_one(fp_uc,
                                      hp,
                                      query_head,
                                      qseqlen,
                                      hp->target);
                }
            }

          if (fp_userout != nullptr)
            {
              results_show_userout_one(fp_userout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (fp_blast6out != nullptr)
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
      if (fp_uc != nullptr)
        {
          results_show_uc_one(fp_uc,
                              nullptr,
                              query_head,
                              qseqlen,
                              0);
        }

      if (opt_output_no_hits != 0)
        {
          if (fp_userout != nullptr)
            {
              results_show_userout_one(fp_userout,
                                       nullptr,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (fp_blast6out != nullptr)
            {
              results_show_blast6out_one(fp_blast6out,
                                         nullptr,
                                         query_head,
                                         qseqlen);
            }
        }
    }

  if (hit_count != 0)
    {
      ++count_matched;
      if (opt_matched != nullptr)
        {
          fasta_print_general(fp_matched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              static_cast<int>(std::strlen(query_head)),
                              0,
                              count_matched,
                              -1.0,
                              -1, -1, nullptr, 0.0,
                              0);
        }
    }
  else
    {
      ++count_notmatched;
      if (opt_notmatched != nullptr)
        {
          fasta_print_general(fp_notmatched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              static_cast<int>(std::strlen(query_head)),
                              0,
                              count_notmatched,
                              -1.0,
                              -1, -1, nullptr, 0.0,
                              0);
        }
    }
}


auto allpairs_thread_run(uint64_t t) -> void
{
  (void) t;

  struct searchinfo_s searchinfo;

  searchinfo.hits_v.resize(static_cast<std::size_t>(seqcount));
  searchinfo.hits = searchinfo.hits_v.data();

  searchinfo.s = search16_init(static_cast<CELL>(opt_match),
                        static_cast<CELL>(opt_mismatch),
                        static_cast<CELL>(opt_gap_open_query_left),
                        static_cast<CELL>(opt_gap_open_target_left),
                        static_cast<CELL>(opt_gap_open_query_interior),
                        static_cast<CELL>(opt_gap_open_target_interior),
                        static_cast<CELL>(opt_gap_open_query_right),
                        static_cast<CELL>(opt_gap_open_target_right),
                        static_cast<CELL>(opt_gap_extension_query_left),
                        static_cast<CELL>(opt_gap_extension_target_left),
                        static_cast<CELL>(opt_gap_extension_query_interior),
                        static_cast<CELL>(opt_gap_extension_target_interior),
                        static_cast<CELL>(opt_gap_extension_query_right),
                        static_cast<CELL>(opt_gap_extension_target_right));


  struct Scoring scoring;
  scoring.match = opt_match;
  scoring.mismatch = opt_mismatch;
  scoring.gap_open_query_left = opt_gap_open_query_left;
  scoring.gap_open_target_left = opt_gap_open_target_left;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_open_target_interior = opt_gap_open_target_interior;
  scoring.gap_open_query_right = opt_gap_open_query_right;
  scoring.gap_open_target_right = opt_gap_open_target_right;
  scoring.gap_extension_query_left = opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opt_gap_extension_target_left;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_extension_target_interior = opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opt_gap_extension_target_right;


  LinearMemoryAligner lma(scoring);


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
      std::unique_lock<std::mutex> input_lock(mutex_input);

      int const query_no = queries;

      if (query_no < seqcount)
        {
          ++queries;

          /* let other threads read input */
          input_lock.unlock();

          /* init search info */
          auto const query_no_u = static_cast<uint64_t>(query_no);
          searchinfo.query_no = query_no;
          searchinfo.qsize = static_cast<int64_t>(db_getabundance(query_no_u));
          searchinfo.query_head_len = static_cast<int>(db_getheaderlen(query_no_u));
          searchinfo.query_head = db_getheader(query_no_u);
          searchinfo.qseqlen = static_cast<int>(db_getsequencelen(query_no_u));
          searchinfo.qsequence = db_getsequence(query_no_u);
          searchinfo.rejects = 0;
          searchinfo.accepts = 0;
          searchinfo.hit_count = 0;

          for (int target = searchinfo.query_no + 1; target < seqcount; target++)
            {
              if ((opt_acceptall != 0) or search_acceptable_unaligned(searchinfo, target))
                {
                  pseqnos[static_cast<std::size_t>(searchinfo.hit_count)] = static_cast<unsigned int>(target);
                  ++searchinfo.hit_count;
                }
            }

          if (searchinfo.hit_count != 0)
            {
              /* perform alignments */

              search16_qprep(searchinfo.s, searchinfo.qsequence, searchinfo.qseqlen);

              search16(searchinfo.s,
                       static_cast<unsigned int>(searchinfo.hit_count),
                       pseqnos.data(),
                       pscores.data(),
                       paligned.data(),
                       pmatches.data(),
                       pmismatches.data(),
                       pgaps.data(),
                       pcigar.data());

              /* convert to hit structure */
              for (std::size_t h = 0; h < static_cast<std::size_t>(searchinfo.hit_count); h++)
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
                      int64_t const tseqlen = static_cast<int64_t>(db_getsequencelen(target));

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

                  hit->target = static_cast<int>(target);
                  hit->strand = 0;
                  hit->count = 0;

                  hit->accepted = false;
                  hit->rejected = false;
                  hit->aligned = true;
                  hit->weak = false;

                  hit->nwscore = static_cast<int>(nwscore);
                  hit->nwdiff = static_cast<int>(nwalignmentlength - nwmatches);
                  hit->nwgaps = static_cast<int>(nwgaps);
                  hit->nwindels = static_cast<int>(nwalignmentlength - nwmatches - nwmismatches);
                  hit->nwalignmentlength = static_cast<int>(nwalignmentlength);
                  hit->nwid = 100.0 * static_cast<double>(nwalignmentlength - hit->nwdiff) /
                    static_cast<double>(nwalignmentlength);
                  hit->nwalignment = nwcigar;
                  hit->matches = static_cast<int>(nwalignmentlength - hit->nwdiff);
                  hit->mismatches = hit->nwdiff - hit->nwindels;

                  auto const dseqlen = static_cast<int>(db_getsequencelen(target));
                  hit->shortest = std::min(searchinfo.qseqlen, dseqlen);
                  hit->longest = std::max(searchinfo.qseqlen, dseqlen);

                  /* trim alignment, compute numbers excluding terminal gaps */
                  align_trim(hit);

                  /* test accept/reject criteria after alignment */
                  if ((opt_acceptall != 0) or search_acceptable_aligned(searchinfo, hit))
                    {
                      finalhits[static_cast<std::size_t>(searchinfo.accepts)] = *hit;
                      ++searchinfo.accepts;
                    }
                }

              /* sort hits (skip when empty: qsort requires a non-null
                 pointer even for zero elements) */
              if (searchinfo.accepts > 0)
                {
                  std::qsort(finalhits.data(), static_cast<std::size_t>(searchinfo.accepts),
                        sizeof(struct hit), allpairs_hit_compare);
                }
            }

          /* lock mutex for update of global data and output */
          std::unique_lock<std::mutex> output_lock(mutex_output);

          /* output results */
          allpairs_output_results(searchinfo.accepts,
                                  finalhits.data(),
                                  searchinfo.query_head,
                                  searchinfo.qseqlen,
                                  searchinfo.qsequence,
                                  nullptr);

          /* update stats */
          if (searchinfo.accepts != 0)
            {
              ++qmatches;
            }

          /* show progress */
          progress += seqcount - query_no - 1;
          progress_update(static_cast<uint64_t>(progress));

          output_lock.unlock();

          /* free memory for alignment strings */
          for (std::size_t i = 0; i < static_cast<std::size_t>(searchinfo.hit_count); i++)
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
          input_lock.unlock();

          cont = false;
        }
    }

  search16_exit(searchinfo.s);
}


auto allpairs_thread_worker_run() -> void
{
  /* run the worker pool; each worker keeps its own search state and
     processes queries until the shared counter is exhausted */
  ThreadRunner threadrunner(static_cast<std::size_t>(opt_threads),
                            allpairs_thread_run);
  threadrunner.run();
}


auto allpairs_global(struct Parameters const & parameters, char const * cmdline, char const * progheader) -> void
{
  /* open output files */

  if (opt_alnout != nullptr)
    {
      fp_alnout = fopen_output(opt_alnout);
      if (fp_alnout == nullptr)
        {
          fatal("Unable to open alignment output file for writing");
        }

      std::fprintf(fp_alnout, "%s\n", parameters.command_line.c_str());
      std::fprintf(fp_alnout, "%s\n", progheader);
    }

  if (opt_samout != nullptr)
    {
      fp_samout = fopen_output(opt_samout);
      if (fp_samout == nullptr)
        {
          fatal("Unable to open SAM output file for writing");
        }
    }

  if (opt_userout != nullptr)
    {
      fp_userout = fopen_output(opt_userout);
      if (fp_userout == nullptr)
        {
          fatal("Unable to open user-defined output file for writing");
        }
    }

  if (opt_blast6out != nullptr)
    {
      fp_blast6out = fopen_output(opt_blast6out);
      if (fp_blast6out == nullptr)
        {
          fatal("Unable to open blast6-like output file for writing");
        }
    }

  if (opt_uc != nullptr)
    {
      fp_uc = fopen_output(opt_uc);
      if (fp_uc == nullptr)
        {
          fatal("Unable to open uc output file for writing");
        }
    }

  if (opt_fastapairs != nullptr)
    {
      fp_fastapairs = fopen_output(opt_fastapairs);
      if (fp_fastapairs == nullptr)
        {
          fatal("Unable to open fastapairs output file for writing");
        }
    }

  if (opt_qsegout != nullptr)
    {
      fp_qsegout = fopen_output(opt_qsegout);
      if (fp_qsegout == nullptr)
        {
          fatal("Unable to open qsegout output file for writing");
        }
    }

  if (opt_tsegout != nullptr)
    {
      fp_tsegout = fopen_output(opt_tsegout);
      if (fp_tsegout == nullptr)
        {
          fatal("Unable to open tsegout output file for writing");
        }
    }

  if (opt_matched != nullptr)
    {
      fp_matched = fopen_output(opt_matched);
      if (fp_matched == nullptr)
        {
          fatal("Unable to open matched output file for writing");
        }
    }

  if (opt_notmatched != nullptr)
    {
      fp_notmatched = fopen_output(opt_notmatched);
      if (fp_notmatched == nullptr)
        {
          fatal("Unable to open notmatched output file for writing");
        }
    }

  db_read(parameters.opt_allpairs_global, 0);

  results_show_samheader(fp_samout, cmdline, parameters.opt_allpairs_global);

  if (parameters.opt_qmask == MASK_DUST)
    {
      dust_all();
    }
  else if ((parameters.opt_qmask == MASK_SOFT) and parameters.opt_hardmask)
    {
      hardmask_all();
    }

  show_rusage();

  seqcount = static_cast<int>(db_getsequencecount());

  /* prepare reading of queries */
  qmatches = 0;
  queries = 0;

  progress = 0;
  progress_init("Aligning", static_cast<uint64_t>(std::max(int64_t{0}, (static_cast<int64_t>(seqcount)) * (static_cast<int64_t>(seqcount) - 1)) / 2));  // refactoring: issue with parenthesis?
  allpairs_thread_worker_run();
  progress_done();

  if (not parameters.opt_quiet)
    {
      std::fprintf(stderr, "Matching query sequences: %d of %d",
              qmatches, queries);
      if (queries > 0)
        {
          std::fprintf(stderr, " (%.2f%%)", 100.0 * qmatches / queries);
        }
      std::fprintf(stderr, "\n");
    }

  if (parameters.opt_log != nullptr)
    {
      std::fprintf(fp_log, "Matching query sequences: %d of %d",
              qmatches, queries);
      if (queries > 0)
        {
          std::fprintf(fp_log, " (%.2f%%)", 100.0 * qmatches / queries);
        }
      std::fprintf(fp_log, "\n\n");
    }

  /* clean up, global */
  db_free();
  if (opt_matched != nullptr)
    {
      std::fclose(fp_matched);
    }
  if (opt_notmatched != nullptr)
    {
      std::fclose(fp_notmatched);
    }
  if (opt_fastapairs != nullptr)
    {
      std::fclose(fp_fastapairs);
    }
  if (opt_qsegout != nullptr)
    {
      std::fclose(fp_qsegout);
    }
  if (opt_tsegout != nullptr)
    {
      std::fclose(fp_tsegout);
    }
  if (fp_uc != nullptr)
    {
      std::fclose(fp_uc);
    }
  if (fp_blast6out != nullptr)
    {
      std::fclose(fp_blast6out);
    }
  if (fp_userout != nullptr)
    {
      std::fclose(fp_userout);
    }
  if (fp_alnout != nullptr)
    {
      std::fclose(fp_alnout);
    }
  if (fp_samout != nullptr)
    {
      std::fclose(fp_samout);
    }
  show_rusage();

  if (fp_userout != nullptr) {
    clean_up(); // free userfields allocation
  }
}
