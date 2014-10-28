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


/* per thread data */

inline int hit_compare_typed(struct hit * x, struct hit * y)
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

int hit_compare(const void * a, const void * b)
{
  return hit_compare_typed((struct hit *) a, (struct hit *) b);
}

inline void topscore_insert(int i, struct searchinfo_s * si)
{
  count_t count = si->kmers[i];
  
  /* ignore sequences with very few kmer matches */

  if (count < MINMATCHSAMPLECOUNT)
    return;

  if (count < MINMATCHSAMPLEFREQ * si->kmersamplecount)
    return;
  
  unsigned int seqno = dbindex_getmapping(i);
  unsigned int length = db_getsequencelen(seqno);

  elem_t novel;
  novel.count = count;
  novel.seqno = seqno;
  novel.length = length;
  
#if 0
  printf("i, count, seqno, length: %d %d %d %d\n", i, count, seqno, length);
#endif

  minheap_add(si->m, & novel);
}

void _mm_print_epi8(__m128i x)
{
  unsigned char * y = (unsigned char*)&x;
  for (int i=0; i<16; i++)
    printf("%s%02x", (i>0?" ":""), y[15-i]);
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
  
#if 0
  printf("Indexed sequences: %d\n", indexed_count);
#endif

  /* zero counts */
  memset(si->kmers, 0, indexed_count * sizeof(count_t));
  
  minheap_empty(si->m);

  for(unsigned int i=0; i<si->kmersamplecount; i++)
    {
      unsigned int kmer = si->kmersample[i];

      unsigned char * bitmap = dbindex_getbitmap(kmer);
      
      if (bitmap)
        {
          /* http://stackoverflow.com/questions/21622212/how-to-perform-the-inverse-of-mm256-movemask-epi8-vpmovmskb */
          
          const __m128i c1 =
            _mm_set_epi32(0x01010101, 0x01010101, 0x00000000, 0x00000000);
          const __m128i c2 = 
            _mm_set_epi32(0x7fbfdfef, 0xf7fbfdfe, 0x7fbfdfef, 0xf7fbfdfe);
          const __m128i c3 =
            _mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);
          
          int r = (indexed_count + 15) / 16;
          unsigned short * p = (unsigned short *)(bitmap);
          __m128i * q = (__m128i *)(si->kmers);
          
          for(int j=0; j<r; j++)
            {
              __m128i xmm0, xmm1, xmm2, xmm3, xmm4, xmm5;
              xmm0 = _mm_loadu_si128((__m128i*)p++);
              xmm1 = _mm_shuffle_epi8(xmm0, c1);
              xmm2 = _mm_or_si128(xmm1, c2);
              xmm3 = _mm_cmpeq_epi8(xmm2, c3);
              xmm4 = _mm_unpacklo_epi8(xmm3, xmm3);
              xmm5 = _mm_unpackhi_epi8(xmm3, xmm3);
              *q = _mm_subs_epi16(*q, xmm4);
              q++;
              *q = _mm_subs_epi16(*q, xmm5);
              q++;
            }
        }
      else
        {
          unsigned int * list = dbindex_getmatchlist(kmer);
          unsigned int count = dbindex_getmatchcount(kmer);
          for(unsigned int j=0; j < count; j++)
            si->kmers[list[j]]++;
        }
    }

  for(int i=0; i < indexed_count; i++)
    topscore_insert(i, si);
  
  minheap_sort(si->m);
}

int seqncmp(char * a, char * b, unsigned long n)
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
  long run = 1;
  int scanlength = 0;
  sscanf(p, "%ld%n", &run, &scanlength);
  char op = *(p+scanlength);
  if (op != 'M')
    {
      hit->trim_aln_left = 1 + scanlength;
      if (op == 'D')
        hit->trim_q_left = run;
      else
        hit->trim_t_left = run;
    }
  
  /* right trim alignment */
  
  char * e = hit->nwalignment + strlen(hit->nwalignment);
  p = e - 1;
  op = *p;
  if (op != 'M')
    {
      while (*(p-1) <= '9')
        p--;
      run = 1;
      sscanf(p, "%ld", &run);
      hit->trim_aln_right = e - p;
      if (op == 'D')
        hit->trim_q_right = run;
      else
        hit->trim_t_right = run;
    }
  
  hit->internal_alignmentlength = hit->nwalignmentlength
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;
  
  hit->internal_indels = hit->nwindels
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;

  hit->internal_gaps = hit->nwgaps
    - (hit->trim_q_left  + hit->trim_t_left  > 0 ? 1 : 0)
    - (hit->trim_q_right + hit->trim_t_right > 0 ? 1 : 0);
  
  /* CD-HIT */
  hit->id0 = 100.0 * hit->matches / hit->shortest;
  /* all diffs */
  hit->id1 = 100.0 * hit->matches / hit->nwalignmentlength;
  /* internal diffs */
  hit->id2 = 100.0 * hit->matches / hit->internal_alignmentlength; 
  /* Marine Biology Lab */
  hit->id3 = MAX(0.0, 100.0 * (1.0 - (1.0 * (hit->mismatches + hit->nwgaps) /
                                      hit->shortest)));
  /* BLAST */
  hit->id4 = 100.0 * hit->matches / hit->nwalignmentlength;

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
  long dseqlen = db_getsequencelen(target);
  long tsize = db_getabundance(target);

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
      (hit->internal_alignmentlength >= opt_query_cov * si->qseqlen) &&
      /* target_cov */
      (hit->internal_alignmentlength >=
       opt_target_cov * db_getsequencelen(hit->target)) &&
      /* maxid */
      (hit->id <= 100.0 * opt_maxid) &&
      /* mid */
      (100.0 * hit->matches / (hit->matches + hit->mismatches) >= opt_mid) &&
      /* maxdiffs */
      (hit->mismatches + hit->internal_indels <= opt_maxdiffs))
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
      if ((si->rejects < maxrejects) && (si->accepts < maxaccepts))
        {
          struct hit * hit = si->hits + x;
      
          if (hit->rejected)
            {
              si->rejects++;
            }
          else
            {
              CELL  nwscore = nwscore_list[i];
              unsigned short nwalignmentlength = nwalignmentlength_list[i];
              unsigned short nwmatches = nwmatches_list[i];
              unsigned short nwmismatches = nwmismatches_list[i];
              unsigned short nwgaps = nwgaps_list[i];
              char * nwcigar = nwcigar_list[i];

              long dseqlen = db_getsequencelen(hit->target);

#ifdef COMPARENONVECTORIZED

              /* compare results with non-vectorized nw */

              char * qseq = si->qsequence;
              char * dseq = db_getsequence(hit->target);

              long xnwscore;
              long xnwdiff;
              long xnwindels;
              long xnwgaps;
              long xnwalignmentlength;
              char * xnwcigar = 0;
      
              long scorematrix[16][16];

              for(int x=0; x<16; x++)
                for(int y=0; y<16; y++)
                  if ((x==0) || (y==0) || (x>4) || (y>4))
                    scorematrix[x][y] = 0;
                  else if (x==y)
                    scorematrix[x][y] = match_score;
                  else
                    scorematrix[x][y] = mismatch_score;

              nw_align(dseq,
                       dseq + dseqlen,
                       qseq,
                       qseq + si->qseqlen,
                       (long*) scorematrix,
                       opt_gap_open_query_left,
                       opt_gap_open_query_interior,
                       opt_gap_open_query_right,
                       opt_gap_open_target_left,
                       opt_gap_open_target_interior,
                       opt_gap_open_target_right,
                       opt_gap_extension_query_left,
                       opt_gap_extension_query_interior,
                       opt_gap_extension_query_right,
                       opt_gap_extension_target_left,
                       opt_gap_extension_target_interior,
                       opt_gap_extension_target_right,
                       & xnwscore,
                       & xnwdiff,
                       & xnwgaps,
                       & xnwindels,
                       & xnwalignmentlength,
                       & xnwcigar,
                       si->query_no,
                       hit->target,
                       si->nw);

              if ((xnwscore != nwscore) || (strcmp(xnwcigar, nwcigar) != 0))
                {
                  printf("Alignment error in channel %d:\n", i);
                  printf("qlen, dlen: %d, %ld\n",si->qseqlen, dseqlen);
                  printf("qseq: [%s]\n", qseq);
                  printf("dseq: [%s]\n", dseq);
                  printf("Non-vectorized: %ld %s\n", xnwscore, xnwcigar);
                  printf("Vectorized:     %d %s\n",   nwscore,  nwcigar);
                }
      
              free(xnwcigar);

#endif

              hit->aligned = 1;
              hit->shortest = MIN(si->qseqlen, dseqlen);
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
    free(nwcigar_list[i++]);

  si->finalized = si->hit_count;
}

void search_onequery(struct searchinfo_s * si)
{
  si->hit_count = 0;

  search16_qprep(si->s, si->qsequence, si->qseqlen);
  
  /* extract unique kmer samples from query*/
  unique_count(si->uh, opt_wordlength, 
               si->qseqlen, si->qsequence,
               & si->kmersamplecount, & si->kmersample);
  
  /* find database sequences with the most kmer hits */
  search_topscores(si);
  
  /* analyse targets with the highest number of kmer hits */
  si->accepts = 0;
  si->rejects = 0;
  si->finalized = 0;

  int delayed = 0;

  int t = 0;
  while ((si->finalized + delayed < maxaccepts + maxrejects - 1) &&
         (si->rejects < maxrejects) &&
         (si->accepts < maxaccepts) && 
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
}

struct hit * search_findbest2(struct searchinfo_s * si_p,
                              struct searchinfo_s * si_m)
{
  struct hit * best = 0;

  for(int i=0; i < si_p->hit_count; i++)
    if ((!best) || (hit_compare_typed(si_p->hits + i, best) < 0))
      best = si_p->hits + i;
  
  if (opt_strand>1)
    for(int i=0; i < si_m->hit_count; i++)
      if ((!best) || (hit_compare_typed(si_m->hits + i, best) < 0))
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
            free(h->nwalignment);
        }
    }
  
  qsort(hits, a, sizeof(struct hit), hit_compare);

  *hitsp = hits;
  *hit_count = a;
}
