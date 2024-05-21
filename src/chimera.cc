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
#include <limits>
#include <vector>


/*
  This code implements the method described in this paper:

  Robert C. Edgar, Brian J. Haas, Jose C. Clemente, Christopher Quince
  and Rob Knight (2011)
  UCHIME improves sensitivity and speed of chimera detection
  Bioinformatics, 27, 16, 2194-2200
  https://doi.org/10.1093/bioinformatics/btr381
*/

/* global constants/data, no need for synchronization */
static int parts = 0;
const int maxparts = 100;
const int window = 64;
const int few = 4;
const int maxcandidates = few * maxparts;
const int rejects = 16;
const double chimera_id = 0.55;
static int tophits;
static pthread_attr_t attr;
static pthread_t * pthread;
static fastx_handle query_fasta_h;

/* mutexes and global data protected by mutex */
static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static unsigned int seqno = 0;
static uint64_t progress = 0;
static int chimera_count = 0;
static int nonchimera_count = 0;
static int borderline_count = 0;
static int total_count = 0;
static int64_t chimera_abundance = 0;
static int64_t nonchimera_abundance = 0;
static int64_t borderline_abundance = 0;
static int64_t total_abundance = 0;
static FILE * fp_chimeras = nullptr;
static FILE * fp_nonchimeras = nullptr;
static FILE * fp_uchimealns = nullptr;
static FILE * fp_uchimeout = nullptr;
static FILE * fp_borderline = nullptr;

/* information for each query sequence to be checked */
struct chimera_info_s
{
  int query_alloc; /* the longest query sequence allocated memory for */
  int head_alloc; /* the longest header allocated memory for */

  int query_no;
  char * query_head;
  int query_head_len;
  int query_size;
  char * query_seq;
  int query_len;

  struct searchinfo_s si[maxparts];

  unsigned int cand_list[maxcandidates];
  int cand_count;

  struct s16info_s * s;
  CELL snwscore[maxcandidates];
  unsigned short snwalignmentlength[maxcandidates];
  unsigned short snwmatches[maxcandidates];
  unsigned short snwmismatches[maxcandidates];
  unsigned short snwgaps[maxcandidates];
  int64_t nwscore[maxcandidates];
  int64_t nwalignmentlength[maxcandidates];
  int64_t nwmatches[maxcandidates];
  int64_t nwmismatches[maxcandidates];
  int64_t nwgaps[maxcandidates];
  char * nwcigar[maxcandidates];

  int match_size;
  int * match;
  int * insert;
  int * smooth;
  int * maxsmooth;

  double * scan_p;
  double * scan_q;

  int parents_found;
  int best_parents[maxparents];
  int best_start[maxparents];
  int best_len[maxparents];

  int best_target;
  char * best_cigar;

  int * maxi;
  char * paln[maxparents];
  char * qaln;
  char * diffs;
  char * votes;
  char * model;
  char * ignore;

  struct hit * all_hits;
  double best_h;
};

static struct chimera_info_s * cia;

void realloc_arrays(struct chimera_info_s * ci)
{
  if (opt_chimeras_denovo)
    {
      if (opt_chimeras_parts == 0) {
        parts = (ci->query_len + maxparts - 1) / maxparts;
      }
      else {
        parts = opt_chimeras_parts;
      }
      if (parts < 2) {
        parts = 2;
      }
      else if (parts > maxparts) {
        parts = maxparts;
      }
    }
  else
    {
      /* default for uchime, uchime2, and uchime3 */
      parts = 4;
    }

  const int maxhlen = MAX(ci->query_head_len, 1);
  if (maxhlen > ci->head_alloc)
    {
      ci->head_alloc = maxhlen;
      ci->query_head = (char *) xrealloc(ci->query_head, maxhlen + 1);
    }

  /* realloc arrays based on query length */

  const int maxqlen = MAX(ci->query_len, 1);
  const int maxpartlen = (maxqlen + parts - 1) / parts;

  if (maxqlen > ci->query_alloc)
    {
      ci->query_alloc = maxqlen;

      ci->query_seq = (char *) xrealloc(ci->query_seq, maxqlen + 1);

      for(auto & i: ci->si)
        {
          i.qsequence = (char *) xrealloc(i.qsequence, maxpartlen + 1);
        }

      ci->maxi = (int *) xrealloc(ci->maxi, (maxqlen + 1) * sizeof(int));
      ci->maxsmooth = (int *) xrealloc(ci->maxsmooth, maxqlen * sizeof(int));
      ci->match = (int *) xrealloc(ci->match,
                                  maxcandidates * maxqlen * sizeof(int));
      ci->insert = (int *) xrealloc(ci->insert,
                                   maxcandidates * maxqlen * sizeof(int));
      ci->smooth = (int *) xrealloc(ci->smooth,
                                   maxcandidates * maxqlen * sizeof(int));

      ci->scan_p = (double *) xrealloc(ci->scan_p,
                                       (maxqlen + 1) * sizeof(double));
      ci->scan_q = (double *) xrealloc(ci->scan_q,
                                       (maxqlen + 1) * sizeof(double));

      const int maxalnlen = maxqlen + 2 * db_getlongestsequence();
      for (int f = 0; f < maxparents ; f++)
        {
          ci->paln[f] = (char *) xrealloc(ci->paln[f], maxalnlen + 1);
        }
      ci->qaln = (char *) xrealloc(ci->qaln, maxalnlen + 1);
      ci->diffs = (char *) xrealloc(ci->diffs, maxalnlen + 1);
      ci->votes = (char *) xrealloc(ci->votes, maxalnlen + 1);
      ci->model = (char *) xrealloc(ci->model, maxalnlen + 1);
      ci->ignore = (char *) xrealloc(ci->ignore, maxalnlen + 1);
    }
}

void find_matches(struct chimera_info_s * ci)
{
  /* find the positions with matches for each potential parent */
  /* also note the positions with inserts in front */

  char * qseq = ci->query_seq;

  for (int i = 0; i < ci->cand_count; i++)
    for (int j = 0; j < ci->query_len; j++)
      {
        int x = i * ci->query_len + j;
        ci->match[x] = 0;
        ci->insert[x] = 0;
      }

  for(int i = 0; i < ci->cand_count; i++)
    {
      char * tseq = db_getsequence(ci->cand_list[i]);

      int qpos = 0;
      int tpos = 0;

      char * p = ci->nwcigar[i];
      char * e = p + strlen(p);

      while (p < e)
        {
          int run = 1;
          int scanlength = 0;
          sscanf(p, "%d%n", &run, &scanlength);
          p += scanlength;
          char op = *p++;
          switch (op)
            {
            case 'M':
              for(int k = 0; k < run; k++)
                {
                  if (chrmap_4bit[(int) (qseq[qpos])] &
                      chrmap_4bit[(int) (tseq[tpos])])
                    {
                      ci->match[i * ci->query_len + qpos] = 1;
                    }
                  ++qpos;
                  ++tpos;
                }
              break;

            case 'I':
              ci->insert[i * ci->query_len + qpos] = run;
              tpos += run;
              break;

            case 'D':
              qpos += run;
              break;
            }
        }
    }
}

struct parents_info_s
{
  int cand;
  int start;
  int len;
};

int compare_positions(const void * a, const void * b)
{
  const int x = ((const parents_info_s *) a)->start;
  const int y = ((const parents_info_s *) b)->start;

  if (x < y)
    return -1;
  else if (x > y)
    return +1;
  else
    return 0;
}

bool scan_matches(struct chimera_info_s * ci,
                  int * matches,
                  int len,
                  double percentage,
                  int * best_start,
                  int * best_len)
{
  /*
    Scan matches array of zeros and ones, and find the longest subsequence
    having a match fraction above or equal to the given percentage (e.g. 2%).
    Based on an idea of finding the longest positive sum substring:
    https://stackoverflow.com/questions/28356453/longest-positive-sum-substring
    If the percentage is 2%, matches are given a score of 2 and mismatches -98.
  */

  double score_match = percentage;
  double score_mismatch = percentage - 100.0;

  double * p = ci->scan_p;
  double * q = ci->scan_q;

  p[0] = 0.0;
  for (int i = 0; i < len; i++)
    p[i + 1] = p[i] + (matches[i] ? score_match : score_mismatch);

  q[len] = p[len];
  for (int i = len - 1; i >= 0; i--)
    q[i] = MAX(q[i + 1], p[i]);

  int best_i = 0;
  int best_d = -1;
  double best_c = -1.0;
  int i = 1;
  int j = 1;
  while (j <= len)
    {
      double c = q[j] - p[i - 1];
      if (c >= 0.0)
        {
          int d = j - i + 1;
          if (d > best_d)
            {
              best_i = i;
              best_d = d;
              best_c = c;
            }
          j += 1;
        }
      else
        {
          i += 1;
        }
    }

  if (best_c >= 0.0)
    {
      * best_start = best_i - 1;
      * best_len = best_d;
      return true;
    }
  else
    return false;
}

int find_best_parents_long(struct chimera_info_s * ci)
{
  /* Find parents with longest matching regions, without indels, allowing
     a given percentage of mismatches (specified with --chimeras_diff_pct),
     and excluding regions matched by previously identified parents. */

  find_matches(ci);

  struct parents_info_s best_parents[maxparents];

  for (int f = 0; f < maxparents; f++)
    {
      best_parents[f].cand = -1;
      best_parents[f].start = -1;
    }

  std::vector<bool> position_used(ci->query_len, false);

  int pos_remaining = ci->query_len;
  int parents_found = 0;

  for (int f = 0; f < opt_chimeras_parents_max; f++)
    {
      /* scan each candidate and find longest matching region */

      int best_start = 0;
      int best_len = 0;
      int best_cand = -1;

      for (int i = 0; i < ci->cand_count; i++)
        {
          int start = 0;
          int len = 0;
          int j = 0;
          while (j < ci->query_len)
            {
              start = j;
              len = 0;
              while ((j < ci->query_len) &&
                     (not position_used[j]) &&
                     ((len == 0) or (ci->insert[i * ci->query_len + j] == 0)))
                {
                  len++;
                  j++;
                }
              if (len > best_len)
                {
                  int scan_best_start = 0;
                  int scan_best_len = 0;
                  if (scan_matches(ci,
                                   ci->match + i * ci->query_len + start,
                                   len,
                                   opt_chimeras_diff_pct,
                                   & scan_best_start,
                                   & scan_best_len))
                    {
                      if (scan_best_len > best_len)
                        {
                          best_cand = i;
                          best_start = start + scan_best_start;
                          best_len = scan_best_len;
                        }
                    }
                }
              j++;
            }
        }

      if (best_len >= opt_chimeras_length_min)
        {
          best_parents[f].cand = best_cand;
          best_parents[f].start = best_start;
          best_parents[f].len = best_len;
          ++parents_found;

#if 0
          if (f == 0)
            printf("\n");
          printf("Best parents long: %d %d %d %d %s %s\n",
                 f,
                 best_cand,
                 best_start,
                 best_len,
                 ci->query_head,
                 db_getheader(ci->cand_list[best_cand]));
#endif

          /* mark positions used */
          for (int j = best_start; j < best_start + best_len; j++)
            {
              position_used[j] = true;
            }
          pos_remaining -= best_len;
        }
      else
        break;
    }

  /* sort parents by position */
  qsort(best_parents,
        parents_found,
        sizeof(struct parents_info_s),
        compare_positions);

  ci->parents_found = parents_found;

  for (int f = 0; f < parents_found; f++)
    {
      ci->best_parents[f] = best_parents[f].cand;
      ci->best_start[f] = best_parents[f].start;
      ci->best_len[f] = best_parents[f].len;
    }

#if 0
  if (pos_remaining == 0)
    printf("Fully covered!\n");
  else
    printf("Not covered completely (%d).\n", pos_remaining);
#endif

  return (parents_found > 1) and (pos_remaining == 0);
}

int find_best_parents(struct chimera_info_s * ci)
{
  find_matches(ci);

  int best_parent_cand[maxparents];

  for (int f = 0; f < 2; f++)
    {
      best_parent_cand[f] = -1;
      ci->best_parents[f] = -1;
    }

  std::vector<bool> cand_selected(ci->cand_count, false);

  for (int f = 0; f < 2; f++)
    {
      if (f > 0)
        {
          /* for all parents except the first */

          /* wipe out matches for all candidates in positions
             covered by the previous parent */

          for(int qpos = window - 1; qpos < ci->query_len; qpos++)
            {
              int z = best_parent_cand[f - 1] * ci->query_len + qpos;
              if (ci->smooth[z] == ci->maxsmooth[qpos])
                {
                  for(int i = qpos + 1 - window; i <= qpos; i++)
                    {
                      for(int j = 0; j < ci->cand_count; j++)
                        {
                          ci->match[j * ci->query_len + i] = 0;
                        }
                    }
                }
            }
        }


      /* Compute smoothed score in a 32bp window for each candidate. */
      /* Record max smoothed score for each position among candidates left. */

      for (int j = 0; j < ci->query_len; j++)
        ci->maxsmooth[j] = 0;

      for(int i = 0; i < ci->cand_count; i++)
        {
          if (not cand_selected[i])
            {
              int sum = 0;
              for(int qpos = 0; qpos < ci->query_len; qpos++)
                {
                  int z = i * ci->query_len + qpos;
                  sum += ci->match[z];
                  if (qpos >= window)
                    {
                      sum -= ci->match[z - window];
                    }
                  if (qpos >= window - 1)
                    {
                      ci->smooth[z] = sum;
                      if (ci->smooth[z] > ci->maxsmooth[qpos])
                        {
                          ci->maxsmooth[qpos] = ci->smooth[z];
                        }
                    }
                }
            }
        }


      /* find parent with the most wins */

      std::vector<int> wins(ci->cand_count, 0);

      for(int qpos = window - 1; qpos < ci->query_len; qpos++)
        {
          if (ci->maxsmooth[qpos] != 0)
            {
              for(int i = 0; i < ci->cand_count; i++)
                {
                  if (not cand_selected[i])
                    {
                      int z = i * ci->query_len + qpos;
                      if (ci->smooth[z] == ci->maxsmooth[qpos])
                        {
                          wins[i]++;
                        }
                    }
                }
            }
        }

      /* select best parent based on most wins */

      int maxwins = 0;
      for(int i = 0; i < ci->cand_count; i++)
        {
          int w = wins[i];
          if (w > maxwins)
            {
              maxwins = w;
              best_parent_cand[f] = i;
            }
        }

      /* terminate loop if no parent found */

      if (best_parent_cand[f] < 0) {
        break;
      }

#if 0
      printf("Query %d: Best parent (%d) candidate: %d. Wins: %d\n",
             ci->query_no, f, best_parent_cand[f], maxwins);
#endif

      ci->best_parents[f] = best_parent_cand[f];
      cand_selected[best_parent_cand[f]] = true;
    }

  /* Check if at least 2 candidates selected */

  return (best_parent_cand[0] >= 0) and (best_parent_cand[1] >= 0);
}


int find_max_alignment_length(struct chimera_info_s * ci)
{
  /* find max insertions in front of each position in the query sequence */

  for (int i = 0; i <= ci->query_len; i++)
    ci->maxi[i] = 0;

  for (int f = 0; f < ci->parents_found; f++)
    {
      int best_parent = ci->best_parents[f];
      char * p = ci->nwcigar[best_parent];
      char * e = p + strlen(p);
      int pos = 0;
      while (p < e)
        {
          int run = 1;
          int scanlength = 0;
          sscanf(p, "%d%n", &run, &scanlength);
          p += scanlength;
          char op = *p++;
          switch (op)
            {
            case 'M':
            case 'D':
              pos += run;
              break;

            case 'I':
              if (run > ci->maxi[pos])
                {
                  ci->maxi[pos] = run;
                }
              break;
            }
        }
    }

  /* find total alignment length */
  int alnlen = 0;
  for(int i = 0; i < ci->query_len + 1; i++)
    {
      alnlen += ci->maxi[i];
    }
  alnlen += ci->query_len;

  return alnlen;
}

void fill_alignment_parents(struct chimera_info_s * ci)
{
  /* fill in alignment strings for the parents */

  for(int j = 0; j < ci->parents_found; j++)
    {
      int cand = ci->best_parents[j];
      int target_seqno = ci->cand_list[cand];
      char * target_seq = db_getsequence(target_seqno);

      int inserted = 0;
      int qpos = 0;
      int tpos = 0;

      char * t = ci->paln[j];
      char * p = ci->nwcigar[cand];
      char * e = p + strlen(p);

      while (p < e)
        {
          int run = 1;
          int scanlength = 0;
          sscanf(p, "%d%n", &run, &scanlength);
          p += scanlength;
          char op = *p++;

          if (op == 'I')
            {
              for(int x=0; x < ci->maxi[qpos]; x++)
                {
                  if (x < run)
                    {
                      *t++ = chrmap_upcase[(int) (target_seq[tpos++])];
                    }
                  else
                    {
                      *t++ = '-';
                    }
                }
              inserted = 1;
            }
          else
            {
              for(int x = 0; x < run; x++)
                {
                  if (not inserted)
                    {
                      for(int y = 0; y < ci->maxi[qpos]; y++)
                        {
                          *t++ = '-';
                        }
                    }

                  if (op == 'M')
                    {
                      *t++ = chrmap_upcase[(int) (target_seq[tpos++])];
                    }
                  else
                    {
                      *t++ = '-';
                    }

                  ++qpos;
                  inserted = 0;
                }
            }
        }

      /* add any gaps at the end */

      if (not inserted)
        {
          for(int x=0; x < ci->maxi[qpos]; x++)
            {
              *t++ = '-';
            }
        }

      /* end of sequence string */
      *t = 0;
    }
}


int eval_parents_long(struct chimera_info_s * ci)
{
  /* always chimeric if called */
  int status = 4;

  int alnlen = find_max_alignment_length(ci);

  fill_alignment_parents(ci);

  /* fill in alignment string for query */

  char * pm = ci->model;
  int m = 0;
  char * q = ci->qaln;
  int qpos = 0;
  for (int i = 0; i < ci->query_len; i++)
    {
      if (qpos >= (ci->best_start[m] + ci->best_len[m]))
        m++;
      for (int j = 0; j < ci->maxi[i]; j++)
        {
          *q++ = '-';
          *pm++ = 'A' + m;
        }
      *q++ = chrmap_upcase[(int)(ci->query_seq[qpos++])];
      *pm++ = 'A' + m;
    }
  for (int j = 0; j < ci->maxi[ci->query_len]; j++)
    {
      *q++ = '-';
      *pm++ = 'A' + m;
    }
  *q = 0;
  *pm = 0;

  for(int i = 0; i < alnlen; i++)
    {
      unsigned int qsym = chrmap_4bit[(int) (ci->qaln[i])];
      unsigned int psym[maxparents];
      for (int f = 0; f < maxparents; f++)
        psym[f] = 0;
      for (int f = 0; f < ci->parents_found; f++)
        psym[f] = chrmap_4bit[(int) (ci->paln[f][i])];

      /* lower case parent symbols that differ from query */

      for (int f = 0; f < ci->parents_found; f++)
        if (psym[f] and (psym[f] != qsym))
          ci->paln[f][i] = tolower(ci->paln[f][i]);

      /* compute diffs */

      char diff = ' ';

      bool all_defined = qsym;
      for (int f = 0; f < ci->parents_found; f++)
        if (not psym[f])
          all_defined = false;

      if (all_defined)
        {
          int z = 0;
          for (int f = 0; f < ci->parents_found; f++)
            if (psym[f] == qsym)
              {
                diff = 'A' + f;
                z++;
              }
          if (z > 1)
            diff = ' ';
        }

      ci->diffs[i] = diff;
    }

  ci->diffs[alnlen] = 0;


  /* count matches */

  int match_QP[maxparents];
  int cols = 0;

  for(int f = 0; f < ci->parents_found; f++)
    match_QP[f] = 0;

  for(int i = 0; i < alnlen; i++)
    {
      cols++;

      char qsym = chrmap_4bit[(int) (ci->qaln[i])];

      for(int f = 0; f < ci->parents_found; f++)
        {
          char psym = chrmap_4bit[(int) (ci->paln[f][i])];
          if (qsym == psym)
            match_QP[f]++;
        }
    }


  int seqno_a = ci->cand_list[ci->best_parents[0]];
  int seqno_b = ci->cand_list[ci->best_parents[1]];
  int seqno_c = -1;
  if (ci->parents_found > 2)
    seqno_c = ci->cand_list[ci->best_parents[2]];

  double QP[maxparents];
  double QT = 0.0;

  for (int f = 0; f < maxparents; f++)
    {
      if (f < ci->parents_found)
        QP[f] = 100.0 * match_QP[f] / cols;
      else
        QP[f] = 0.0;
      if (QP[f] > QT)
        QT = QP[f];
    }

  double QA = QP[0];
  double QB = QP[1];
  double QC = ci->parents_found > 2 ? QP[2] : 0.00;
  double QM = 100.00;
  double divfrac = 100.00 * (QM - QT) / QT;

  xpthread_mutex_lock(&mutex_output);

  if (opt_alnout and (status == 4))
    {
      fprintf(fp_uchimealns, "\n");
      fprintf(fp_uchimealns, "----------------------------------------"
              "--------------------------------\n");
      fprintf(fp_uchimealns, "Query   (%5d nt) ",
              ci->query_len);
      header_fprint_strip(fp_uchimealns,
                          ci->query_head,
                          ci->query_head_len,
                          opt_xsize,
                          opt_xee,
                          opt_xlength);

      for (int f = 0; f < ci->parents_found; f++)
        {
          int seqno = ci->cand_list[ci->best_parents[f]];
          fprintf(fp_uchimealns, "\nParent%c (%5" PRIu64 " nt) ",
                  'A' + f,
                  db_getsequencelen(seqno));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(seqno),
                              db_getheaderlen(seqno),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
        }

      fprintf(fp_uchimealns, "\n\n");


      int width = opt_alignwidth > 0 ? opt_alignwidth : alnlen;
      qpos = 0;
      int ppos[maxparents];
      for (int f = 0; f < ci->parents_found; f++)
        ppos[f] = 0;
      int rest = alnlen;

      for(int i = 0; i < alnlen; i += width)
        {
          /* count non-gap symbols on current line */

          int qnt = 0;
          int pnt[maxparents];
          for (int f = 0; f < ci->parents_found; f++)
            pnt[f] = 0;

          int w = MIN(rest, width);

          for(int j = 0; j < w; j++)
            {
              if (ci->qaln[i + j] != '-')
                {
                  qnt++;
                }

              for (int f = 0; f < ci->parents_found; f++)
                if (ci->paln[f][i + j] != '-')
                  {
                    pnt[f]++;
                  }
            }

          fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                  qpos + 1, w, ci->qaln + i, qpos + qnt);

          for (int f = 0; f < ci->parents_found; f++)
            {
              fprintf(fp_uchimealns, "%c %5d %.*s %d\n",
                      'A' + f,
                      ppos[f] + 1, w, ci->paln[f] + i, ppos[f] + pnt[f]);
            }

          fprintf(fp_uchimealns, "Diffs   %.*s\n", w, ci->diffs+i);
          fprintf(fp_uchimealns, "Model   %.*s\n", w, ci->model+i);
          fprintf(fp_uchimealns, "\n");

          rest -= width;
          qpos += qnt;
          for (int f = 0; f < ci->parents_found; f++)
            ppos[f] += pnt[f];
        }

      fprintf(fp_uchimealns, "Ids.  QA %.2f%%, QB %.2f%%, QC %.2f%%, "
              "QT %.2f%%, QModel %.2f%%, Div. %+.2f%%\n",
              QA, QB, QC, QT, QM, divfrac);
    }

  if (opt_tabbedout)
    {
      fprintf(fp_uchimeout, "%.4f\t", 99.9999);

      header_fprint_strip(fp_uchimeout,
                          ci->query_head,
                          ci->query_head_len,
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(fp_uchimeout, "\t");
      header_fprint_strip(fp_uchimeout,
                          db_getheader(seqno_a),
                          db_getheaderlen(seqno_a),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(fp_uchimeout, "\t");
      header_fprint_strip(fp_uchimeout,
                          db_getheader(seqno_b),
                          db_getheaderlen(seqno_b),
                          opt_xsize,
                          opt_xee,
                          opt_xlength);
      fprintf(fp_uchimeout, "\t");
      if (seqno_c >= 0)
        {
          header_fprint_strip(fp_uchimeout,
                              db_getheader(seqno_c),
                              db_getheaderlen(seqno_c),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
        }
      else
        {
          fprintf(fp_uchimeout, "*");
        }
      fprintf(fp_uchimeout, "\t");

      fprintf(fp_uchimeout,
              "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t"
              "%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%c\n",
              QM,
              QA,
              QB,
              QC,
              QT,
              0, /* ignore, left yes */
              0, /* ignore, left no */
              0, /* ignore, left abstain */
              0, /* ignore, right yes */
              0, /* ignore, right no */
              0, /* ignore, right abstain */
              0.00,
              status == 4 ? 'Y' : (status == 2 ? 'N' : '?'));
    }

  xpthread_mutex_unlock(&mutex_output);

  return status;
}

int eval_parents(struct chimera_info_s * ci)
{
  int status = 1;
  ci->parents_found = 2;

  int alnlen = find_max_alignment_length(ci);

  fill_alignment_parents(ci);

  /* fill in alignment string for query */

  char * q = ci->qaln;
  int qpos = 0;
  for (int i = 0; i < ci->query_len; i++)
    {
      for (int j=0; j < ci->maxi[i]; j++)
        {
          *q++ = '-';
        }
      *q++ = chrmap_upcase[(int) (ci->query_seq[qpos++])];
    }
  for (int j = 0; j < ci->maxi[ci->query_len]; j++)
    {
      *q++ = '-';
    }
  *q = 0;

  /* mark positions to ignore in voting */

  for (int i = 0; i < alnlen; i++)
    ci->ignore[i] = 0;

  for(int i = 0; i < alnlen; i++)
    {
      unsigned int qsym  = chrmap_4bit[(int) (ci->qaln   [i])];
      unsigned int p1sym = chrmap_4bit[(int) (ci->paln[0][i])];
      unsigned int p2sym = chrmap_4bit[(int) (ci->paln[1][i])];

      /* ignore gap positions and those next to the gap */
      if ((not qsym) or (not p1sym) or (not p2sym))
        {
          ci->ignore[i] = 1;
          if (i > 0)
            {
              ci->ignore[i - 1] = 1;
            }
          if (i < alnlen - 1)
            {
              ci->ignore[i + 1] = 1;
            }
        }

      /* ignore ambiguous symbols */
      if ((ambiguous_4bit[qsym]) or
          (ambiguous_4bit[p1sym]) or
          (ambiguous_4bit[p2sym]))
        {
          ci->ignore[i] = 1;
        }

      /* lower case parent symbols that differ from query */

      if (p1sym and (p1sym != qsym))
        {
          ci->paln[0][i] = tolower(ci->paln[0][i]);
        }

      if (p2sym and (p2sym != qsym))
        {
          ci->paln[1][i] = tolower(ci->paln[1][i]);
        }

      /* compute diffs */

      char diff;

      if (qsym and p1sym and p2sym)
        {
          if (p1sym == p2sym)
            {
              if (qsym == p1sym)
                {
                  diff = ' ';
                }
              else
                {
                  diff = 'N';
                }
            }
          else
            {
              if (qsym == p1sym)
                {
                  diff = 'A';
                }
              else if (qsym == p2sym)
                {
                  diff = 'B';
                }
              else
                {
                  diff = '?';
                }
            }
        }
      else
        {
          diff = ' ';
        }

      ci->diffs[i] = diff;
    }

  ci->diffs[alnlen] = 0;

  /* compute score */

  int sumA = 0;
  int sumB = 0;
  int sumN = 0;

  for (int i = 0; i < alnlen; i++)
    {
      if (not ci->ignore[i])
        {
          char diff = ci->diffs[i];

          if (diff == 'A')
            {
              ++sumA;
            }
          else if (diff == 'B')
            {
              ++sumB;
            }
          else if (diff != ' ')
            {
              ++sumN;
            }
        }

    }

  int left_n = 0;
  int left_a = 0;
  int left_y = 0;
  int right_n = sumA;
  int right_a = sumN;
  int right_y = sumB;

  double best_h = -1;
  int best_i = -1;
  int best_reverse = 0;

  int best_left_y = 0;
  int best_right_y = 0;
  int best_left_n = 0;
  int best_right_n = 0;
  int best_left_a = 0;
  int best_right_a = 0;

  for (int i = 0; i < alnlen; i++)
    {
      if(not ci->ignore[i])
        {
          char diff = ci->diffs[i];
          if (diff != ' ')
            {
              if (diff == 'A')
                {
                  ++left_y;
                  --right_n;
                }
              else if (diff == 'B')
                {
                  ++left_n;
                  --right_y;
                }
              else
                {
                  ++left_a;
                  --right_a;
                }

              double left_h = 0;
              double right_h = 0;
              double h = 0;

              if ((left_y > left_n) and (right_y > right_n))
                {
                  left_h = left_y / (opt_xn * (left_n + opt_dn) + left_a);
                  right_h = right_y / (opt_xn * (right_n + opt_dn) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_reverse = 0;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_n;
                      best_left_y = left_y;
                      best_left_a = left_a;
                      best_right_n = right_n;
                      best_right_y = right_y;
                      best_right_a = right_a;
                    }
                }
              else if ((left_n > left_y) and (right_n > right_y))
                {
                  /* swap left/right and yes/no */

                  left_h = left_n / (opt_xn * (left_y + opt_dn) + left_a);
                  right_h = right_n / (opt_xn * (right_y + opt_dn) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_reverse = 1;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_y;
                      best_left_y = left_n;
                      best_left_a = left_a;
                      best_right_n = right_y;
                      best_right_y = right_n;
                      best_right_a = right_a;
                    }
                }
            }
        }
    }

  ci->best_h = best_h > 0 ? best_h : 0.0;

  if (best_h >= 0.0)
    {
      status = 2;

      /* flip A and B if necessary */

      if (best_reverse)
        {
          for(int i = 0; i < alnlen; i++)
            {
              char diff = ci->diffs[i];
              if (diff == 'A')
                {
                  ci->diffs[i] = 'B';
                }
              else if (diff == 'B')
                {
                  ci->diffs[i] = 'A';
                }
            }
        }

      /* fill in votes and model */

      for(int i = 0; i < alnlen; i++)
        {
          char m = i <= best_i ? 'A' : 'B';
          ci->model[i] = m;

          char v = ' ';
          if (not ci->ignore[i])
            {
              char d = ci->diffs[i];

              if ((d == 'A') or (d == 'B'))
                {
                  if (d == m)
                    {
                      v = '+';
                    }
                  else
                    {
                      v = '!';
                    }
                }
              else if ((d == 'N') or (d == '?'))
                {
                  v = '0';
                }
            }
          ci->votes[i] = v;

          /* lower case diffs for no votes */
          if (v == '!')
            {
              ci->diffs[i] = tolower(ci->diffs[i]);
            }
        }

      /* fill in crossover region */

      for(int i = best_i + 1; i < alnlen; i++)
        {
          if ((ci->diffs[i] == ' ') or (ci->diffs[i] == 'A'))
            {
              ci->model[i] = 'x';
            }
          else
            {
              break;
            }
        }

      ci->votes[alnlen] = 0;
      ci->model[alnlen] = 0;

      /* count matches */

      int index_a = best_reverse ? 1 : 0;
      int index_b = best_reverse ? 0 : 1;

      int match_QA = 0;
      int match_QB = 0;
      int match_AB = 0;
      int match_QM = 0;
      int cols = 0;

      for(int i = 0; i < alnlen; i++)
        {
          if (not ci->ignore[i])
            {
              cols++;

              char qsym = chrmap_4bit[(int) (ci->qaln[i])];
              char asym = chrmap_4bit[(int) (ci->paln[index_a][i])];
              char bsym = chrmap_4bit[(int) (ci->paln[index_b][i])];
              char msym = (i <= best_i) ? asym : bsym;

              if (qsym == asym)
                {
                  match_QA++;
                }

              if (qsym == bsym)
                {
                  match_QB++;
                }

              if (asym == bsym)
                {
                  match_AB++;
                }

              if (qsym == msym)
                {
                  match_QM++;
                }
            }
        }

      int seqno_a = ci->cand_list[ci->best_parents[index_a]];
      int seqno_b = ci->cand_list[ci->best_parents[index_b]];

      double QA = 100.0 * match_QA / cols;
      double QB = 100.0 * match_QB / cols;
      double AB = 100.0 * match_AB / cols;
      double QT = MAX(QA, QB);
      double QM = 100.0 * match_QM / cols;
      double divdiff = QM - QT;
      double divfrac = 100.0 * divdiff / QT;

      int sumL = best_left_n + best_left_a + best_left_y;
      int sumR = best_right_n + best_right_a + best_right_y;

      if (opt_uchime2_denovo or opt_uchime3_denovo)
        {
          // fix -Wfloat-equal: if match_QM == cols, then QM == 100.0
          if ((match_QM == cols) and (QT < 100.0))
            {
              status = 4;
            }
        }
      else
        if (best_h >= opt_minh)
          {
            status = 3;
            if ((divdiff >= opt_mindiv) and
                (sumL >= opt_mindiffs) and
                (sumR >= opt_mindiffs))
              {
                status = 4;
              }
          }

      /* print alignment */

      xpthread_mutex_lock(&mutex_output);

      if (opt_uchimealns and (status == 4))
        {
          fprintf(fp_uchimealns, "\n");
          fprintf(fp_uchimealns, "----------------------------------------"
                  "--------------------------------\n");
          fprintf(fp_uchimealns, "Query   (%5d nt) ",
                  ci->query_len);

          header_fprint_strip(fp_uchimealns,
                              ci->query_head,
                              ci->query_head_len,
                              opt_xsize,
                              opt_xee,
                              opt_xlength);

          fprintf(fp_uchimealns, "\nParentA (%5" PRIu64 " nt) ",
                  db_getsequencelen(seqno_a));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(seqno_a),
                              db_getheaderlen(seqno_a),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);

          fprintf(fp_uchimealns, "\nParentB (%5" PRIu64 " nt) ",
                  db_getsequencelen(seqno_b));
          header_fprint_strip(fp_uchimealns,
                              db_getheader(seqno_b),
                              db_getheaderlen(seqno_b),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimealns, "\n\n");

          int width = opt_alignwidth > 0 ? opt_alignwidth : alnlen;
          qpos = 0;
          int p1pos = 0;
          int p2pos = 0;
          int rest = alnlen;

          for(int i = 0; i < alnlen; i += width)
            {
              /* count non-gap symbols on current line */

              int qnt = 0;
              int p1nt = 0;
              int p2nt = 0;

              int w = MIN(rest,width);

              for(int j = 0; j < w; j++)
                {
                  if (ci->qaln[i + j] != '-')
                    {
                      qnt++;
                    }
                  if (ci->paln[0][i + j] != '-')
                    {
                      p1nt++;
                    }
                  if (ci->paln[1][i + j] != '-')
                    {
                      p2nt++;
                    }
                }

              if (not best_reverse)
                {
                  fprintf(fp_uchimealns, "A %5d %.*s %d\n",
                          p1pos + 1, w, ci->paln[0] + i, p1pos + p1nt);
                  fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos + 1, w, ci->qaln + i, qpos + qnt);
                  fprintf(fp_uchimealns, "B %5d %.*s %d\n",
                          p2pos + 1, w, ci->paln[1] + i, p2pos + p2nt);
                }
              else
                {
                  fprintf(fp_uchimealns, "A %5d %.*s %d\n",
                          p2pos + 1, w, ci->paln[1] + i, p2pos + p2nt);
                  fprintf(fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos + 1, w, ci->qaln + i, qpos + qnt);
                  fprintf(fp_uchimealns, "B %5d %.*s %d\n",
                          p1pos + 1, w, ci->paln[0] + i, p1pos + p1nt);
                }

              fprintf(fp_uchimealns, "Diffs   %.*s\n", w, ci->diffs + i);
              fprintf(fp_uchimealns, "Votes   %.*s\n", w, ci->votes + i);
              fprintf(fp_uchimealns, "Model   %.*s\n", w, ci->model + i);
              fprintf(fp_uchimealns, "\n");

              qpos += qnt;
              p1pos += p1nt;
              p2pos += p2nt;
              rest -= width;
            }

          fprintf(fp_uchimealns, "Ids.  QA %.1f%%, QB %.1f%%, AB %.1f%%, "
                  "QModel %.1f%%, Div. %+.1f%%\n",
                  QA, QB, AB, QM, divfrac);

          fprintf(fp_uchimealns, "Diffs Left %d: N %d, A %d, Y %d (%.1f%%); "
                  "Right %d: N %d, A %d, Y %d (%.1f%%), Score %.4f\n",
                  sumL, best_left_n, best_left_a, best_left_y,
                  100.0 * best_left_y / sumL,
                  sumR, best_right_n, best_right_a, best_right_y,
                  100.0 * best_right_y / sumR,
                  best_h);
        }

      if (opt_uchimeout)
        {
          fprintf(fp_uchimeout, "%.4f\t", best_h);

          header_fprint_strip(fp_uchimeout,
                              ci->query_head,
                              ci->query_head_len,
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimeout, "\t");
          header_fprint_strip(fp_uchimeout,
                              db_getheader(seqno_a),
                              db_getheaderlen(seqno_a),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimeout, "\t");
          header_fprint_strip(fp_uchimeout,
                              db_getheader(seqno_b),
                              db_getheaderlen(seqno_b),
                              opt_xsize,
                              opt_xee,
                              opt_xlength);
          fprintf(fp_uchimeout, "\t");

          if(not opt_uchimeout5)
            {
              if (QA >= QB)
                {
                  header_fprint_strip(fp_uchimeout,
                                      db_getheader(seqno_a),
                                      db_getheaderlen(seqno_a),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                }
              else
                {
                  header_fprint_strip(fp_uchimeout,
                                      db_getheader(seqno_b),
                                      db_getheaderlen(seqno_b),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                }
              fprintf(fp_uchimeout, "\t");
            }

          fprintf(fp_uchimeout,
                  "%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t"
                  "%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%c\n",
                  QM,
                  QA,
                  QB,
                  AB,
                  QT,
                  best_left_y,
                  best_left_n,
                  best_left_a,
                  best_right_y,
                  best_right_n,
                  best_right_a,
                  divdiff,
                  status == 4 ? 'Y' : (status == 2 ? 'N' : '?'));
        }
      xpthread_mutex_unlock(&mutex_output);
    }

  return status;
}

// refactoring: enum struct status {};
/*
  new chimeric status:
  0: no parents, non-chimeric
  1: score < 0 (no alignment), non-chimeric
  2: score < minh, non-chimeric
  3: score >= minh, suspicious -> not available with uchime2_denovo and uchime3_denovo
  4: score >= minh && (divdiff >= opt_mindiv) && ..., chimeric
*/

void query_init(struct searchinfo_s * si)
{
  si->qsequence = nullptr;
  si->kmers = nullptr;
  si->hits = (struct hit *) xmalloc(sizeof(struct hit) * tophits);
  si->kmers = (count_t *) xmalloc(db_getsequencecount() * sizeof(count_t) + 32);
  si->hit_count = 0;
  si->uh = unique_init();
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
  si->m = minheap_init(tophits);
}

void query_exit(struct searchinfo_s * si)
{
  search16_exit(si->s);
  unique_exit(si->uh);
  minheap_exit(si->m);
  nw_exit(si->nw);

  if (si->qsequence)
    {
      xfree(si->qsequence);
      si->qsequence = nullptr;
    }
  if (si->hits)
    {
      xfree(si->hits);
      si->hits = nullptr;
    }
  if (si->kmers)
    {
      xfree(si->kmers);
      si->kmers = nullptr;
    }
}

void partition_query(struct chimera_info_s * ci)
{
  int rest = ci->query_len;
  char * p = ci->query_seq;
  for (int i = 0; i < parts; i++)
    {
      int len = (rest + (parts - i - 1)) / (parts - i);

      struct searchinfo_s * si = ci->si + i;

      si->query_no = ci->query_no;
      si->strand = 0;
      si->qsize = ci->query_size;
      si->query_head_len = ci->query_head_len;
      si->query_head = ci->query_head;
      si->qseqlen = len;
      strncpy(si->qsequence, p, len);
      si->qsequence[len] = 0;

      rest -= len;
      p += len;
    }
}

void chimera_thread_init(struct chimera_info_s * ci)
{
  ci->query_alloc = 0;
  ci->head_alloc = 0;
  ci->query_head = nullptr;
  ci->query_seq = nullptr;
  ci->maxi = nullptr;
  ci->maxsmooth = nullptr;
  ci->match = nullptr;
  ci->insert = nullptr;
  ci->smooth = nullptr;
  ci->qaln = nullptr;
  ci->diffs = nullptr;
  ci->votes = nullptr;
  ci->model = nullptr;
  ci->ignore = nullptr;
  ci->scan_p = nullptr;
  ci->scan_q = nullptr;

  for (int f = 0; f < maxparents; f++)
    {
      ci->paln[f] = nullptr;
    }

  for(int i = 0; i < maxparts; i++)
    {
      query_init(ci->si + i);
    }

  ci->s = search16_init(opt_match,
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

void chimera_thread_exit(struct chimera_info_s * ci)
{
  search16_exit(ci->s);

  for(int i = 0; i < maxparts; i++)
    {
      query_exit(ci->si + i);
    }

  if (ci->maxsmooth)
    {
      xfree(ci->maxsmooth);
    }
  if (ci->match)
    {
      xfree(ci->match);
    }
  if (ci->insert)
    {
      xfree(ci->insert);
    }
  if (ci->smooth)
    {
      xfree(ci->smooth);
    }
  if (ci->diffs)
    {
      xfree(ci->diffs);
    }
  if (ci->votes)
    {
      xfree(ci->votes);
    }
  if (ci->model)
    {
      xfree(ci->model);
    }
  if (ci->ignore)
    {
      xfree(ci->ignore);
    }
  if (ci->maxi)
    {
      xfree(ci->maxi);
    }
  if (ci->qaln)
    {
      xfree(ci->qaln);
    }
  if (ci->query_seq)
    {
      xfree(ci->query_seq);
    }
  if (ci->query_head)
    {
      xfree(ci->query_head);
    }
  if (ci->scan_p)
    {
      xfree(ci->scan_p);
    }
  if (ci->scan_q)
    {
      xfree(ci->scan_q);
    }

  for (int f = 0; f < maxparents; f++)
    if (ci->paln[f])
      {
        xfree(ci->paln[f]);
      }
}

uint64_t chimera_thread_core(struct chimera_info_s * ci)
{
  chimera_thread_init(ci);

  auto * allhits_list = (struct hit *) xmalloc(maxcandidates *
                                               sizeof(struct hit));

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

  while(true)
    {
      /* get next sequence */

      xpthread_mutex_lock(&mutex_input);

      if (opt_uchime_ref)
        {
          if (fasta_next(query_fasta_h, not opt_notrunclabels,
                         chrmap_no_change))
            {
              ci->query_head_len = fasta_get_header_length(query_fasta_h);
              ci->query_len = fasta_get_sequence_length(query_fasta_h);
              ci->query_no = fasta_get_seqno(query_fasta_h);
              ci->query_size = fasta_get_abundance(query_fasta_h);

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              /* copy the data locally (query seq, head) */
              strcpy(ci->query_head, fasta_get_header(query_fasta_h));
              strcpy(ci->query_seq, fasta_get_sequence(query_fasta_h));
            }
          else
            {
              xpthread_mutex_unlock(&mutex_input);
              break; /* end while loop */
            }
        }
      else
        {
          if (seqno < db_getsequencecount())
            {
              ci->query_no = seqno;
              ci->query_head_len = db_getheaderlen(seqno);
              ci->query_len = db_getsequencelen(seqno);
              ci->query_size = db_getabundance(seqno);

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              strcpy(ci->query_head, db_getheader(seqno));
              strcpy(ci->query_seq, db_getsequence(seqno));
            }
          else
            {
              xpthread_mutex_unlock(&mutex_input);
              break; /* end while loop */
            }
        }

      xpthread_mutex_unlock(&mutex_input);



      int status = 0;

      /* partition query */
      partition_query(ci);

      /* perform searches and collect candidate parents */
      ci->cand_count = 0;
      int allhits_count = 0;

      if (ci->query_len >= parts)
        {
          for (int i = 0; i < parts; i++)
            {
              struct hit * hits;
              int hit_count;
              search_onequery(ci->si + i, opt_qmask);
              search_joinhits(ci->si + i, nullptr, & hits, & hit_count);
              for(int j = 0; j < hit_count; j++)
                {
                  if (hits[j].accepted)
                    {
                      allhits_list[allhits_count++] = hits[j];
                    }
                }
              xfree(hits);
            }
        }

      for(int i = 0; i < allhits_count; i++)
        {
          unsigned int target = allhits_list[i].target;

          /* skip duplicates */
          int k {0};
          for(k = 0; k < ci->cand_count; k++)
            {
              if (ci->cand_list[k] == target)
                {
                  break;
                }
            }

          if (k == ci->cand_count)
            {
              ci->cand_list[ci->cand_count++] = target;
            }

          /* deallocate cigar */
          if (allhits_list[i].nwalignment)
            {
              xfree(allhits_list[i].nwalignment);
              allhits_list[i].nwalignment = nullptr;
            }
        }


      /* align full query to each candidate */

      search16_qprep(ci->s, ci->query_seq, ci->query_len);

      search16(ci->s,
               ci->cand_count,
               ci->cand_list,
               ci->snwscore,
               ci->snwalignmentlength,
               ci->snwmatches,
               ci->snwmismatches,
               ci->snwgaps,
               ci->nwcigar);

      for(int i = 0; i < ci->cand_count; i++)
        {
          int64_t target = ci->cand_list[i];
          int64_t nwscore = ci->snwscore[i];
          char * nwcigar;
          int64_t nwalignmentlength;
          int64_t nwmatches;
          int64_t nwmismatches;
          int64_t nwgaps;

          if (nwscore == std::numeric_limits<short>::max())
            {
              /* In case the SIMD aligner cannot align,
                 perform a new alignment with the
                 linear memory aligner */

              char * tseq = db_getsequence(target);
              int64_t tseqlen = db_getsequencelen(target);

              if (ci->nwcigar[i])
                {
                  xfree(ci->nwcigar[i]);
                }

              nwcigar = xstrdup(lma.align(ci->query_seq,
                                          tseq,
                                          ci->query_len,
                                          tseqlen));
              lma.alignstats(nwcigar,
                             ci->query_seq,
                             tseq,
                             & nwscore,
                             & nwalignmentlength,
                             & nwmatches,
                             & nwmismatches,
                             & nwgaps);

              ci->nwcigar[i] = nwcigar;
              ci->nwscore[i] = nwscore;
              ci->nwalignmentlength[i] = nwalignmentlength;
              ci->nwmatches[i] = nwmatches;
              ci->nwmismatches[i] = nwmismatches;
              ci->nwgaps[i] = nwgaps;
            }
          else
            {
              ci->nwscore[i] = ci->snwscore[i];
              ci->nwalignmentlength[i] = ci->snwalignmentlength[i];
              ci->nwmatches[i] = ci->snwmatches[i];
              ci->nwmismatches[i] = ci->snwmismatches[i];
              ci->nwgaps[i] = ci->snwgaps[i];
            }
        }


      /* find the best pair of parents, then compute score for them */

      if (opt_chimeras_denovo)
        {
          /* long high-quality reads */
          if (find_best_parents_long(ci))
            {
              status = eval_parents_long(ci);
            }
          else
            {
              status = 0;
            }
        }
      else
        {
          if (find_best_parents(ci))
            {
              status = eval_parents(ci);
            }
          else
            {
              status = 0;
            }
        }

      /* output results */

      xpthread_mutex_lock(&mutex_output);

      total_count++;
      total_abundance += ci->query_size;

      if (status == 4)
        {
          chimera_count++;
          chimera_abundance += ci->query_size;

          if (opt_chimeras)
            {
              fasta_print_general(fp_chimeras,
                                  nullptr,
                                  ci->query_seq,
                                  ci->query_len,
                                  ci->query_head,
                                  ci->query_head_len,
                                  ci->query_size,
                                  chimera_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( opt_uchime_ref ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h);

            }
        }

      if (status == 3)
        {
          borderline_count++;
          borderline_abundance += ci->query_size;

          if (opt_borderline)
            {
              fasta_print_general(fp_borderline,
                                  nullptr,
                                  ci->query_seq,
                                  ci->query_len,
                                  ci->query_head,
                                  ci->query_head_len,
                                  ci->query_size,
                                  borderline_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( opt_uchime_ref ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h);

            }
        }

      if (status < 3)
        {
          nonchimera_count++;
          nonchimera_abundance += ci->query_size;

          /* output no parents, no chimeras */
          if ((status < 2) and opt_uchimeout)
            {
              fprintf(fp_uchimeout, "0.0000\t");

              header_fprint_strip(fp_uchimeout,
                                  ci->query_head,
                                  ci->query_head_len,
                                  opt_xsize,
                                  opt_xee,
                                  opt_xlength);

              if (opt_uchimeout5)
                {
                  fprintf(fp_uchimeout,
                          "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
                }
              else
                {
                  fprintf(fp_uchimeout,
                          "\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
                }
            }

          if (opt_nonchimeras)
            {
              fasta_print_general(fp_nonchimeras,
                                  nullptr,
                                  ci->query_seq,
                                  ci->query_len,
                                  ci->query_head,
                                  ci->query_head_len,
                                  ci->query_size,
                                  nonchimera_count,
                                  -1.0,
                                  -1,
                                  -1,
                                  opt_fasta_score ?
                                  ( opt_uchime_ref ?
                                    "uchime_ref" : "uchime_denovo" ) : nullptr,
                                  ci->best_h);
            }
        }

      if (status < 3)
        {
          /* uchime_denovo: add non-chimeras to db */
          if (opt_uchime_denovo or opt_uchime2_denovo or opt_uchime3_denovo or opt_chimeras_denovo)
            {
              dbindex_addsequence(seqno, opt_qmask);
            }
        }

      for (int i = 0; i < ci->cand_count; i++)
        {
          if (ci->nwcigar[i])
            {
              xfree(ci->nwcigar[i]);
            }
        }

      if (opt_uchime_ref)
        {
          progress = fasta_get_position(query_fasta_h);
        }
      else
        {
          progress += db_getsequencelen(seqno);
        }

      progress_update(progress);

      seqno++;

      xpthread_mutex_unlock(&mutex_output);
    }

  if (allhits_list)
    {
      xfree(allhits_list);
    }

  chimera_thread_exit(ci);

  xfree(scorematrix);

  return 0;
}

void * chimera_thread_worker(void * vp)
{
  return (void *) chimera_thread_core(cia + (int64_t) vp);
}

void chimera_threads_run()
{
  xpthread_attr_init(&attr);
  xpthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* create worker threads */
  for(int64_t t = 0; t < opt_threads; t++)
    {
      xpthread_create(pthread + t, & attr,
                      chimera_thread_worker, (void*)t);
    }

  /* finish worker threads */
  for(int t = 0; t < opt_threads; t++)
    {
      xpthread_join(pthread[t], nullptr);
    }

  xpthread_attr_destroy(&attr);
}

void open_chimera_file(FILE ** f, char * name)
{
  if (name)
    {
      *f = fopen_output(name);
      if (not *f)
        {
          fatal("Unable to open file %s for writing", name);
        }
    }
  else
    {
      *f = nullptr;
    }
}

void close_chimera_file(FILE * f)
{
  if (f)
    {
      fclose(f);
    }
}

void chimera()
{
  open_chimera_file(&fp_chimeras, opt_chimeras);
  open_chimera_file(&fp_nonchimeras, opt_nonchimeras);
  open_chimera_file(&fp_borderline, opt_borderline);

  if (opt_chimeras_denovo)
    {
      open_chimera_file(&fp_uchimealns, opt_alnout);
      open_chimera_file(&fp_uchimeout, opt_tabbedout);
    }
  else
    {
      open_chimera_file(&fp_uchimealns, opt_uchimealns);
      open_chimera_file(&fp_uchimeout, opt_uchimeout);
    }


  /* override any options the user might have set */
  opt_maxaccepts = few;
  opt_maxrejects = rejects;
  opt_id = chimera_id;

  if (opt_strand != 1)
    {
      fatal("Only --strand plus is allowed with uchime_ref.");
    }

  if (not opt_uchime_ref)
    {
      opt_self = 1;
      opt_selfid = 1;
      opt_threads = 1;
      opt_maxsizeratio = 1.0 / opt_abskew;
    }

  tophits = opt_maxaccepts + opt_maxrejects;

  uint64_t progress_total;
  chimera_count = 0;
  nonchimera_count = 0;
  progress = 0;
  seqno = 0;

  /* prepare threads */
  pthread = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));
  cia = (struct chimera_info_s *) xmalloc(opt_threads *
                                          sizeof(struct chimera_info_s));

  /* init mutexes for input and output */
  xpthread_mutex_init(&mutex_input, nullptr);
  xpthread_mutex_init(&mutex_output, nullptr);

  char * denovo_dbname = nullptr;

  /* prepare queries / database */
  if (opt_uchime_ref)
    {
      /* check if the reference database may be an UDB file */

      bool is_udb = udb_detect_isudb(opt_db);

      if (is_udb)
        {
          udb_read(opt_db, true, true);
        }
      else
        {
          db_read(opt_db, 0);
          if (opt_dbmask == MASK_DUST)
            {
              dust_all();
            }
          else if ((opt_dbmask == MASK_SOFT) and (opt_hardmask))
            {
              hardmask_all();
            }
          dbindex_prepare(1, opt_dbmask);
          dbindex_addallsequences(opt_dbmask);
        }

      query_fasta_h = fasta_open(opt_uchime_ref);
      progress_total = fasta_get_size(query_fasta_h);
    }
  else
    {

      if (opt_uchime_denovo)
        {
          denovo_dbname = opt_uchime_denovo;
        }
      else if (opt_uchime2_denovo)
        {
          denovo_dbname = opt_uchime2_denovo;
        }
      else if (opt_uchime3_denovo)
        {
          denovo_dbname = opt_uchime3_denovo;
        }
      else if (opt_chimeras_denovo)
        {
          denovo_dbname = opt_chimeras_denovo;
        }
      else
        fatal("Internal error");

      db_read(denovo_dbname, 0);

      if (opt_qmask == MASK_DUST)
        {
          dust_all();
        }
      else if ((opt_qmask == MASK_SOFT) and (opt_hardmask))
        {
          hardmask_all();
        }

      db_sortbyabundance();
      dbindex_prepare(1, opt_qmask);
      progress_total = db_getnucleotidecount();
    }

  if (opt_log)
    {
      if (opt_uchime_ref or opt_uchime_denovo)
        {
          fprintf(fp_log, "%8.2f  minh\n", opt_minh);
        }

      if (opt_uchime_ref or
          opt_uchime_denovo or
          opt_uchime2_denovo or
          opt_uchime3_denovo)
        {
          fprintf(fp_log, "%8.2f  xn\n", opt_xn);
          fprintf(fp_log, "%8.2f  dn\n", opt_dn);
          fprintf(fp_log, "%8.2f  xa\n", 1.0);
        }

      if (opt_uchime_ref or opt_uchime_denovo)
        {
          fprintf(fp_log, "%8.2f  mindiv\n", opt_mindiv);
        }

      fprintf(fp_log, "%8.2f  id\n", opt_id);

      if (opt_uchime_ref or
          opt_uchime_denovo or
          opt_uchime2_denovo or
          opt_uchime3_denovo)
        {
          fprintf(fp_log, "%8d  maxp\n", 2);
        }

      fprintf(fp_log, "\n");
    }


  progress_init("Detecting chimeras", progress_total);

  chimera_threads_run();

  progress_done();

  if (not opt_quiet)
    {
      if (total_count > 0)
        {
          if (opt_chimeras_denovo)
            {
              fprintf(stderr,
                      "Found %d (%.1f%%) chimeras and "
                      "%d (%.1f%%) non-chimeras "
                      "in %u unique sequences.\n",
                      chimera_count,
                      100.0 * chimera_count / total_count,
                      nonchimera_count,
                      100.0 * nonchimera_count / total_count,
                      total_count);
            }
          else
            {
              fprintf(stderr,
                      "Found %d (%.1f%%) chimeras, "
                      "%d (%.1f%%) non-chimeras,\n"
                      "and %d (%.1f%%) borderline sequences "
                      "in %u unique sequences.\n",
                      chimera_count,
                      100.0 * chimera_count / total_count,
                      nonchimera_count,
                      100.0 * nonchimera_count / total_count,
                      borderline_count,
                      100.0 * borderline_count / total_count,
                      total_count);
            }
        }
      else
        {
          if (opt_chimeras_denovo)
            {
              fprintf(stderr,
                      "Found %d chimeras and "
                      "%d non-chimeras "
                      "in %u unique sequences.\n",
                      chimera_count,
                      nonchimera_count,
                      total_count);
            }
          else
            {
              fprintf(stderr,
                      "Found %d chimeras, "
                      "%d non-chimeras,\n"
                      "and %d borderline sequences "
                      "in %u unique sequences.\n",
                      chimera_count,
                      nonchimera_count,
                      borderline_count,
                      total_count);
            }
        }

      if (total_abundance > 0)
        {
          if (opt_chimeras_denovo)
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " (%.1f%%) chimeras and "
                      "%" PRId64 " (%.1f%%) non-chimeras "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      100.0 * chimera_abundance / total_abundance,
                      nonchimera_abundance,
                      100.0 * nonchimera_abundance / total_abundance,
                      total_abundance);
            }
          else
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " (%.1f%%) chimeras, "
                      "%" PRId64 " (%.1f%%) non-chimeras,\n"
                      "and %" PRId64 " (%.1f%%) borderline sequences "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      100.0 * chimera_abundance / total_abundance,
                      nonchimera_abundance,
                      100.0 * nonchimera_abundance / total_abundance,
                      borderline_abundance,
                      100.0 * borderline_abundance / total_abundance,
                      total_abundance);
            }
        }
      else
        {
          if (opt_chimeras_denovo)
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " chimeras, "
                      "%" PRId64 " non-chimeras "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      nonchimera_abundance,
                      total_abundance);
            }
          else
            {
              fprintf(stderr,
                      "Taking abundance information into account, "
                      "this corresponds to\n"
                      "%" PRId64 " chimeras, "
                      "%" PRId64 " non-chimeras,\n"
                      "and %" PRId64 " borderline sequences "
                      "in %" PRId64 " total sequences.\n",
                      chimera_abundance,
                      nonchimera_abundance,
                      borderline_abundance,
                      total_abundance);
            }
        }
    }

  if (opt_log)
    {
      if (opt_uchime_ref)
        {
          fprintf(fp_log, "%s", opt_uchime_ref);
        }
      else
        {
          fprintf(fp_log, "%s", denovo_dbname);
        }

      if (seqno > 0)
        {
          fprintf(fp_log, ": %d/%u chimeras (%.1f%%)\n",
                  chimera_count,
                  seqno,
                  100.0 * chimera_count / seqno);
        }
      else
        {
          fprintf(fp_log, ": %d/%u chimeras\n",
                  chimera_count,
                  seqno);
        }
    }


  if (opt_uchime_ref)
    {
      fasta_close(query_fasta_h);
    }

  dbindex_free();
  db_free();

  xpthread_mutex_destroy(&mutex_output);
  xpthread_mutex_destroy(&mutex_input);

  xfree(cia);
  xfree(pthread);

  close_chimera_file(fp_borderline);
  close_chimera_file(fp_uchimeout);
  close_chimera_file(fp_uchimealns);
  close_chimera_file(fp_nonchimeras);
  close_chimera_file(fp_chimeras);

  show_rusage();
}
