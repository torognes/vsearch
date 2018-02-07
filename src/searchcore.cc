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

/* per thread data */

inline int hit_compare_byid_typed(struct hit * x, struct hit * y)
{
  // high id, then low id
  // early target, then late target

  if (x->rejected < y->rejected)
    return -1;
  else
    if (x->rejected > y->rejected)
      return +1;
    else
      if (x->rejected == 1)
        return 0;
      else
        if (x->aligned > y->aligned)
          return -1;
        else
          if (x->aligned < y->aligned)
            return +1;
          else
            if (x->aligned == 0)
              return 0;
            else
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

inline int hit_compare_bysize_typed(struct hit * x, struct hit * y)
{
  // high abundance, then low abundance
  // high id, then low id
  // early target, then late target

  if (x->rejected < y->rejected)
    return -1;
  else
    if (x->rejected > y->rejected)
      return +1;
    else
      if (x->rejected == 1)
        return 0;
      else
        if (x->aligned > y->aligned)
          return -1;
        else
          if (x->aligned < y->aligned)
            return +1;
          else
            if (x->aligned == 0)
              return 0;
            else
              if (db_getabundance(x->target) > db_getabundance(y->target))
                return -1;
              else
                if (db_getabundance(x->target) < db_getabundance(y->target))
                  return +1;
                else
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

int hit_compare_byid(const void * a, const void * b)
{
  return hit_compare_byid_typed((struct hit *) a, (struct hit *) b);
}

int hit_compare_bysize(const void * a, const void * b)
{
  return hit_compare_bysize_typed((struct hit *) a, (struct hit *) b);
}

bool search_enough_kmers(struct searchinfo_s * si,
                         unsigned int count)
{
  return (count >= opt_minwordmatches) || (count >= si->kmersamplecount);
}

void search_topscores(struct searchinfo_s * si)
{
  /*
    Count the kmer hits in each database sequence and
    make a sorted list of a given number (th)
    of the database sequences with the highest number of matching kmers.
    These are stored in the min heap array.
  */

  /* count kmer hits in the database sequences */
  int indexed_count = dbindex_getcount();
  
  /* zero counts */
  memset(si->kmers, 0, indexed_count * sizeof(count_t));
  
  minheap_empty(si->m);

  for(unsigned int i=0; i<si->kmersamplecount; i++)
    {
      unsigned int kmer = si->kmersample[i];
      unsigned char * bitmap = dbindex_getbitmap(kmer);
      
      if (bitmap)
        {
#ifdef __PPC__
	  increment_counters_from_bitmap(si->kmers, bitmap, indexed_count);
#else
          if (ssse3_present)
            increment_counters_from_bitmap_ssse3(si->kmers, 
                                                 bitmap, indexed_count);
          else
            increment_counters_from_bitmap_sse2(si->kmers,
                                                bitmap, indexed_count);
#endif
        }
      else
        {
          unsigned int * list = dbindex_getmatchlist(kmer);
          unsigned int count = dbindex_getmatchcount(kmer);
          for(unsigned int j=0; j < count; j++)
            si->kmers[list[j]]++;
        }
    }

  int minmatches = MIN(opt_minwordmatches, si->kmersamplecount);

  for(int i=0; i < indexed_count; i++)
    {
      count_t count = si->kmers[i];
      if (count >= minmatches)
        {
          unsigned int seqno = dbindex_getmapping(i);
          unsigned int length = db_getsequencelen(seqno);

          elem_t novel;
          novel.count = count;
          novel.seqno = seqno;
          novel.length = length;

          minheap_add(si->m, & novel);
        }
    }

  minheap_sort(si->m);
}

int seqncmp(char * a, char * b, uint64_t n)
{
  for(unsigned int i = 0; i<n; i++)
    {
      int x = chrmap_4bit[(int)(a[i])];
      int y = chrmap_4bit[(int)(b[i])];
      if (x < y)
        return -1;
      else if (x > y)
        return +1;
    }
  return 0;
}


void align_trim(struct hit * hit)
{
  /* trim alignment and fill in info */
  /* assumes that the hit has been aligned */

  /* info for semi-global alignment (without gaps at ends) */
  
  hit->trim_aln_left = 0;
  hit->trim_q_left = 0;
  hit->trim_t_left = 0;
  hit->trim_aln_right = 0;
  hit->trim_q_right = 0;
  hit->trim_t_right = 0;
  
  /* left trim alignment */
  
  char * p = hit->nwalignment;
  char op;
  int64_t run;
  if (*p)
    {
      run = 1;
      int scanlength = 0;
      sscanf(p, "%" PRId64 "%n", &run, &scanlength);
      op = *(p+scanlength);
      if (op != 'M')
        {
          hit->trim_aln_left = 1 + scanlength;
          if (op == 'D')
            hit->trim_q_left = run;
          else
            hit->trim_t_left = run;
        }
    }
  
  /* right trim alignment */
  
  char * e = hit->nwalignment + strlen(hit->nwalignment);
  if (e > hit->nwalignment)
    {
      p = e - 1;
      op = *p;
      if (op != 'M')
        {
          while ((p > hit->nwalignment) && (*(p-1) <= '9'))
            p--;
          run = 1;
          sscanf(p, "%" PRId64, &run);
          hit->trim_aln_right = e - p;
          if (op == 'D')
            hit->trim_q_right = run;
          else
            hit->trim_t_right = run;
        }
    }

  if (hit->trim_q_left >= hit->nwalignmentlength)
    hit->trim_q_right = 0;

  if (hit->trim_t_left >= hit->nwalignmentlength)
    hit->trim_t_right = 0;

  hit->internal_alignmentlength = hit->nwalignmentlength
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;
  
  hit->internal_indels = hit->nwindels
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;

  hit->internal_gaps = hit->nwgaps
    - ((hit->trim_q_left  + hit->trim_t_left)  > 0 ? 1 : 0)
    - ((hit->trim_q_right + hit->trim_t_right) > 0 ? 1 : 0);

  /* CD-HIT */
  hit->id0 = hit->shortest > 0 ? 100.0 * hit->matches / hit->shortest : 0.0;
  /* all diffs */
  hit->id1 = hit->nwalignmentlength > 0 ?
    100.0 * hit->matches / hit->nwalignmentlength : 0.0;
  /* internal diffs */
  hit->id2 = hit->internal_alignmentlength > 0 ?
    100.0 * hit->matches / hit->internal_alignmentlength : 0.0;
  /* Marine Biology Lab */
  hit->id3 = MAX(0.0, 100.0 * (1.0 - (1.0 * (hit->mismatches + hit->nwgaps) /
                                      hit->longest)));
  /* BLAST */
  hit->id4 = hit->nwalignmentlength > 0 ?
    100.0 * hit->matches / hit->nwalignmentlength : 0.0;

  switch (opt_iddef)
    {
    case 0:
      hit->id = hit->id0;
      break;
    case 1:
      hit->id = hit->id1;
      break;
    case 2:
      hit->id = hit->id2;
      break;
    case 3:
      hit->id = hit->id3;
      break;
    case 4:
      hit->id = hit->id4;
      break;
    }
}

int search_acceptable_unaligned(struct searchinfo_s * si,
                                int target)
{
  /* consider whether a hit satisfy accept criteria before alignment */

  char * qseq = si->qsequence;
  char * dlabel = db_getheader(target);
  char * dseq = db_getsequence(target);
  int64_t dseqlen = db_getsequencelen(target);
  int64_t tsize = db_getabundance(target);

  if (
      /* maxqsize */
      (si->qsize <= opt_maxqsize)
      &&
      /* mintsize */
      (tsize >= opt_mintsize)
      &&
      /* minsizeratio */
      (si->qsize >= opt_minsizeratio * tsize)
      &&
      /* maxsizeratio */
      (si->qsize <= opt_maxsizeratio * tsize)
      &&
      /* minqt */
      (si->qseqlen >= opt_minqt * dseqlen)
      &&
      /* maxqt */
      (si->qseqlen <= opt_maxqt * dseqlen)
      &&
      /* minsl */
      (si->qseqlen < dseqlen ? 
       si->qseqlen >= opt_minsl * dseqlen : 
       dseqlen >= opt_minsl * si->qseqlen)
      &&
      /* maxsl */
      (si->qseqlen < dseqlen ? 
       si->qseqlen <= opt_maxsl * dseqlen :
       dseqlen <= opt_maxsl * si->qseqlen)
      &&
      /* idprefix */
      ((si->qseqlen >= opt_idprefix) &&
       (dseqlen >= opt_idprefix) &&
       (!seqncmp(qseq, dseq, opt_idprefix)))
      &&
      /* idsuffix */
      ((si->qseqlen >= opt_idsuffix) &&
       (dseqlen >= opt_idsuffix) &&
       (!seqncmp(qseq+si->qseqlen-opt_idsuffix,
                dseq+dseqlen-opt_idsuffix,
                opt_idsuffix)))
      &&
      /* self */
      ((!opt_self) || (strcmp(si->query_head, dlabel)))
      &&
      /* selfid */
      ((!opt_selfid) ||
       (si->qseqlen != dseqlen) ||
       (seqncmp(qseq, dseq, si->qseqlen)))
      )
    {
      /* needs further consideration */
      return 1;
    }
  else
    {
      /* reject */
      return 0;
    }
}

int search_acceptable_aligned(struct searchinfo_s * si,
                              struct hit * hit)
{
  if (/* weak_id */
      (hit->id >= 100.0 * opt_weak_id) &&
      /* maxsubs */
      (hit->mismatches <= opt_maxsubs) &&
      /* maxgaps */
      (hit->internal_gaps <= opt_maxgaps) &&
      /* mincols */
      (hit->internal_alignmentlength >= opt_mincols) &&
      /* leftjust */
      ((!opt_leftjust) || (hit->trim_q_left + 
                           hit->trim_t_left == 0)) &&
      /* rightjust */
      ((!opt_rightjust) || (hit->trim_q_right +
                            hit->trim_t_right == 0)) &&
      /* query_cov */
      (hit->matches + hit->mismatches >= opt_query_cov * si->qseqlen) &&
      /* target_cov */
      (hit->matches + hit->mismatches >=
       opt_target_cov * db_getsequencelen(hit->target)) &&
      /* maxid */
      (hit->id <= 100.0 * opt_maxid) &&
      /* mid */
      (100.0 * hit->matches / (hit->matches + hit->mismatches) >= opt_mid) &&
      /* maxdiffs */
      (hit->mismatches + hit->internal_indels <= opt_maxdiffs))
    {
      if (opt_cluster_unoise)
        {
          int d = hit->mismatches;
          double skew = 1.0 * si->qsize / db_getabundance(hit->target);
          double beta = 1.0 / pow(2, 1.0 * opt_unoise_alpha * d + 1);

          if (skew <= beta || d == 0)
            {
              /* accepted */
              hit->accepted = 1;
              hit->weak = 0;
              return 1;
            }
          else
            {
              /* rejected, but weak hit */
              hit->rejected = 1;
              hit->weak = 1;
              return 0;
            }
        }
      else
        {
          if (hit->id >= 100.0 * opt_id)
            {
              /* accepted */
              hit->accepted = 1;
              hit->weak = 0;
              return 1;
            }
          else
            {
              /* rejected, but weak hit */
              hit->rejected = 1;
              hit->weak = 1;
              return 0;
            }
        }
    }
  else
    {
      /* rejected */
      hit->rejected = 1;
      hit->weak = 0;
      return 0;
    }
}

void align_delayed(struct searchinfo_s * si)
{
  /* compute global alignment */
  
  unsigned int target_list[MAXDELAYED];
  CELL  nwscore_list[MAXDELAYED];
  unsigned short nwalignmentlength_list[MAXDELAYED];
  unsigned short nwmatches_list[MAXDELAYED];
  unsigned short nwmismatches_list[MAXDELAYED];
  unsigned short nwgaps_list[MAXDELAYED];
  char * nwcigar_list[MAXDELAYED];
  
  int target_count = 0;

  for(int x = si->finalized; x < si->hit_count; x++)
    {
      struct hit * hit = si->hits + x;
      if (! hit->rejected)
        target_list[target_count++] = hit->target;
    }

  if (target_count)
    search16(si->s,
             target_count,
             target_list,
             nwscore_list,
             nwalignmentlength_list,
             nwmatches_list,
             nwmismatches_list,
             nwgaps_list,
             nwcigar_list);
  
  int i = 0;

  for(int x = si->finalized; x < si->hit_count; x++)
    {
      /* maxrejects or maxaccepts reached - ignore remaining hits */
      if ((si->rejects < opt_maxrejects) && (si->accepts < opt_maxaccepts))
        {
          struct hit * hit = si->hits + x;
      
          if (hit->rejected)
            {
              si->rejects++;
            }
          else
            {
              int64_t target = hit->target;
              int64_t nwscore = nwscore_list[i];

              char * nwcigar;
              int64_t nwalignmentlength;
              int64_t nwmatches;
              int64_t nwmismatches;
              int64_t nwgaps;
              
              int64_t dseqlen = db_getsequencelen(target);

              if (nwscore == SHRT_MAX)
                {
                  /* In case the SIMD aligner cannot align,
                     perform a new alignment with the
                     linear memory aligner */
                  
                  char * dseq = db_getsequence(target);
                  
                  if (nwcigar_list[i])
                    xfree(nwcigar_list[i]);
                  
                  nwcigar = xstrdup(si->lma->align(si->qsequence,
                                                  dseq,
                                                  si->qseqlen,
                                                  dseqlen));

                  si->lma->alignstats(nwcigar,
                                      si->qsequence,
                                      dseq,
                                      & nwscore,
                                      & nwalignmentlength,
                                      & nwmatches,
                                      & nwmismatches,
                                      & nwgaps);
                }
              else
                {
                  nwalignmentlength = nwalignmentlength_list[i];
                  nwmatches = nwmatches_list[i];
                  nwmismatches = nwmismatches_list[i];
                  nwgaps = nwgaps_list[i];
                  nwcigar = nwcigar_list[i];
                }
              
              hit->aligned = 1;
              hit->shortest = MIN(si->qseqlen, dseqlen);
              hit->longest = MAX(si->qseqlen, dseqlen);
              hit->nwalignment = nwcigar;
              hit->nwscore = nwscore;
              hit->nwdiff = nwalignmentlength - nwmatches;
              hit->nwgaps = nwgaps;
              hit->nwindels = nwalignmentlength - nwmatches - nwmismatches;
              hit->nwalignmentlength = nwalignmentlength;
              hit->nwid = 100.0 * (nwalignmentlength - hit->nwdiff) /
                nwalignmentlength;
              hit->matches = nwalignmentlength - hit->nwdiff;
              hit->mismatches = hit->nwdiff - hit->nwindels;

              /* trim alignment and compute numbers excluding terminal gaps */
              align_trim(hit);
              
              /* test accept/reject criteria after alignment */
              if (search_acceptable_aligned(si, hit))
                si->accepts++;
              else
                si->rejects++;

              i++;
            }
        }
    }

  /* free ignored alignments */
  while (i < target_count)
    xfree(nwcigar_list[i++]);

  si->finalized = si->hit_count;
}

void search_onequery(struct searchinfo_s * si, int seqmask)
{
  si->hit_count = 0;

  search16_qprep(si->s, si->qsequence, si->qseqlen);

  si->lma = new LinearMemoryAligner;

  int64_t * scorematrix = si->lma->scorematrix_create(opt_match, opt_mismatch);

  si->lma->set_parameters(scorematrix,
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
  
  /* extract unique kmer samples from query*/
  unique_count(si->uh, opt_wordlength, 
               si->qseqlen, si->qsequence,
               & si->kmersamplecount, & si->kmersample, seqmask);
  
  /* find database sequences with the most kmer hits */
  search_topscores(si);
  
  /* analyse targets with the highest number of kmer hits */
  si->accepts = 0;
  si->rejects = 0;
  si->finalized = 0;

  int delayed = 0;

  int t = 0;
  while ((si->finalized + delayed < opt_maxaccepts + opt_maxrejects - 1) &&
         (si->rejects < opt_maxrejects) &&
         (si->accepts < opt_maxaccepts) && 
         (!minheap_isempty(si->m)))
    {
      elem_t e = minheap_poplast(si->m);
      
      struct hit * hit = si->hits + si->hit_count;

      hit->target = e.seqno;
      hit->count = e.count;
      hit->strand = si->strand;
      hit->rejected = 0;
      hit->accepted = 0;
      hit->aligned = 0;
      hit->weak = 0;
      hit->nwalignment = 0;

      /* Test some accept/reject criteria before alignment */
      if (search_acceptable_unaligned(si, e.seqno))
        {
          delayed++;
        }
      else
        {
          hit->rejected = 1;
        }

      si->hit_count++;

      if (delayed == MAXDELAYED)
        {
          align_delayed(si);
          delayed = 0;
        }
      t++;
    }  
  if (delayed > 0)
    align_delayed(si);
  
  delete si->lma;
  xfree(scorematrix);
}

struct hit * search_findbest2_byid(struct searchinfo_s * si_p,
                                   struct searchinfo_s * si_m)
{
  struct hit * best = 0;

  for(int i=0; i < si_p->hit_count; i++)
    if ((!best) || (hit_compare_byid_typed(si_p->hits + i, best) < 0))
      best = si_p->hits + i;
  
  if (opt_strand>1)
    for(int i=0; i < si_m->hit_count; i++)
      if ((!best) || (hit_compare_byid_typed(si_m->hits + i, best) < 0))
        best = si_m->hits + i;
  
  if (best && ! best->accepted)
    best = 0;

  return best;
}

struct hit * search_findbest2_bysize(struct searchinfo_s * si_p,
                                     struct searchinfo_s * si_m)
{
  struct hit * best = 0;

  for(int i=0; i < si_p->hit_count; i++)
    if ((!best) || (hit_compare_bysize_typed(si_p->hits + i, best) < 0))
      best = si_p->hits + i;
  
  if (opt_strand>1)
    for(int i=0; i < si_m->hit_count; i++)
      if ((!best) || (hit_compare_bysize_typed(si_m->hits + i, best) < 0))
        best = si_m->hits + i;
  
  if (best && ! best->accepted)
    best = 0;

  return best;
}

void search_joinhits(struct searchinfo_s * si_p,
                     struct searchinfo_s * si_m,
                     struct hit * * hitsp,
                     int * hit_count)
{
  /* join and sort accepted hits from both strands */
  /* remove and unallocate unaccepted hits */

  int a = 0;
  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = s ? si_m : si_p;
      for(int i=0; i<si->hit_count; i++)
        if (si->hits[i].accepted)
          a++;
    }
  
  struct hit * hits = (struct hit *) xmalloc(a * sizeof(struct hit));
  
  a = 0;

  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = s ? si_m : si_p;
      for(int i=0; i<si->hit_count; i++)
        {
          struct hit * h = si->hits + i;
          if (h->accepted)
            hits[a++] = *h;
          else if (h->aligned)
            xfree(h->nwalignment);
        }
    }
  
  qsort(hits, a, sizeof(struct hit), hit_compare_byid);

  *hitsp = hits;
  *hit_count = a;
}
