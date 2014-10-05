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

#define MINMATCHSAMPLECOUNT 6
#define MINMATCHSAMPLEFREQ (1/16)

//#define COMPARENONVECTORIZED

/* per thread data */

int hit_compare(const void * a, const void * b)
{
  struct hit * x = (struct hit *) a;
  struct hit * y = (struct hit *) b;

  // high id, then low id
  // early target, then late target

  if (x->internal_id > y->internal_id)
    return -1;
  else if (x->internal_id < y->internal_id)
    return +1;
  else
    if (x->target < y->target)
      return -1;
    else if (x->target > y->target)
      return +1;
    else
      return 0;
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
          /*
            http://stackoverflow.com/questions/21622212/
            how-to-perform-the-inverse-of-mm256-movemask-epi8-vpmovmskb
          */
     
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
              __m128i xmm, xmm0, xmm1;
              xmm = _mm_loadu_si128((__m128i*)p++);
              xmm = _mm_shuffle_epi8(xmm, c1);
              xmm = _mm_or_si128(xmm, c2);
              xmm = _mm_cmpeq_epi8(xmm, c3);
              xmm0 = _mm_unpacklo_epi8(xmm, xmm);
              xmm1 = _mm_unpackhi_epi8(xmm, xmm);
              *q = _mm_subs_epi16(*q, xmm0);
              q++;
              *q = _mm_subs_epi16(*q, xmm1);
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
    if (chrmap_4bit[(int)(a[i])] != chrmap_4bit[(int)(b[i])])
      return 1;
  return 0;
}


void align_candidates(struct searchinfo_s * si)
{
  /* compute global alignment */
  
  unsigned int target_list[8];
  CELL  nwscore_list[8];
  unsigned short nwalignmentlength_list[8];
  unsigned short nwmatches_list[8];
  unsigned short nwmismatches_list[8];
  unsigned short nwgaps_list[8];
  char * nwcigar_list[8];
  
  for(int i=0; i < si->candidatecount; i++)
    target_list[i] = si->candidate_target[i];
  
  search16(si->s,
           si->candidatecount,
           target_list,
           nwscore_list,
           nwalignmentlength_list,
           nwmatches_list,
           nwmismatches_list,
           nwgaps_list,
           nwcigar_list);
  
  for(int i = 0; i < si->candidatecount; i++)
    {
      unsigned int target = target_list[i];
      CELL  nwscore = nwscore_list[i];
      unsigned short nwalignmentlength = nwalignmentlength_list[i];
      unsigned short nwmatches = nwmatches_list[i];
      unsigned short nwmismatches = nwmismatches_list[i];
      unsigned short nwgaps = nwgaps_list[i];
      char * nwcigar = nwcigar_list[i];

#ifdef COMPARENONVECTORIZED
      /* compare results with non-vectorized nw */

      char * qseq;
      if (si->strand)
        qseq = si->rc;
      else
        qseq = si->qsequence;

      char * dseq = db_getsequence(target);
      long dseqlen = db_getsequencelen(target);

      long xnwscore;
      long xnwdiff;
      long xnwindels;
      long xnwgaps;
      long xnwalignmentlength;
      char * xnwcigar = 0;
      
      for(int i=0; i<16; i++)
        for(int j=0; j<16; j++)
          if (i==j)
            scorematrix[i][j] = match_score;
          else if ((i==0) || (j==0) || (i>4) || (j>4))
            scorematrix[i][j] = 0;
          else
            scorematrix[i][j] = mismatch_score;

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
               target,
               si->nw);

      if ((xnwscore != nwscore) || (strcmp(xnwcigar, nwcigar) != 0))
        {
          printf("Alignment error in channel %d:\n", i);
          printf("qlen, dlen: %ld, %ld\n",si->qseqlen, dseqlen);
          printf("qseq: [%s]\n", qseq);
          printf("dseq: [%s]\n", dseq);
          printf("Non-vectorized: %ld %s\n", xnwscore, xnwcigar);
          printf("Vectorized:     %d %s\n",   nwscore,  nwcigar);
        }
      
      free(xnwcigar);

#endif

      int kmercount = si->candidate_kmercount[i];

      int nwdiff = nwalignmentlength - nwmatches;
      int nwindels = nwalignmentlength - nwmatches - nwmismatches;
      char * nwalignment = nwcigar;
      
      double nwid = (nwalignmentlength - nwdiff) * 100.0 / nwalignmentlength;
      
      /* info for semi-global alignment (without gaps at ends) */
      
      long trim_aln_left = 0;
      long trim_q_left = 0;
      long trim_t_left = 0;
      long trim_aln_right = 0;
      long trim_q_right = 0;
      long trim_t_right = 0;
      
      /* left trim alignment */
      
      char * p = nwalignment;
      long run = 1;
      int scanlength = 0;
      sscanf(p, "%ld%n", &run, &scanlength);
      char op = *(p+scanlength);
      if (op != 'M')
        {
          trim_aln_left = 1 + scanlength;
          if (op == 'D')
            trim_q_left = run;
          else
            trim_t_left = run;
        }
      
      /* right trim alignment */
      
      char * e = nwalignment + strlen(nwalignment);
      p = e - 1;
      op = *p;
      if (op != 'M')
        {
          while (*(p-1) <= '9')
            p--;
          run = 1;
          sscanf(p, "%ld", &run);
          trim_aln_right = e - p;
          if (op == 'D')
            trim_q_right = run;
          else
            trim_t_right = run;
        }
      
      long mismatches = nwdiff - nwindels;
      long matches = nwalignmentlength - nwdiff;
      long internal_alignmentlength = nwalignmentlength
        - trim_q_left - trim_t_left - trim_q_right - trim_t_right;
      long internal_gaps = nwgaps
        - (trim_q_left  + trim_t_left  > 0 ? 1 : 0)
        - (trim_q_right + trim_t_right > 0 ? 1 : 0);
      long internal_indels = nwindels
        - trim_q_left - trim_t_left - trim_q_right - trim_t_right;
      double internal_id = 100.0 * matches / internal_alignmentlength;
      
      /* test accept/reject criteria after alignment */

      if (/* weak_id */
          (internal_id >= 100.0 * opt_weak_id) &&
          /* maxsubs */
          (mismatches <= opt_maxsubs) &&
          /* maxgaps */
          (internal_gaps <= opt_maxgaps) &&
          /* mincols */
          (internal_alignmentlength >= opt_mincols) &&
          /* leftjust */
          ((!opt_leftjust) || (trim_q_left+trim_t_left == 0)) &&
          /* rightjust */
          ((!opt_rightjust) || (trim_q_right+trim_t_right == 0)) &&
          /* query_cov */
          (internal_alignmentlength >= opt_query_cov * si->qseqlen) &&
          /* target_cov */
          (internal_alignmentlength >= opt_target_cov *  
           db_getsequencelen(target)) &&
          /* maxid */
          (internal_id <= 100.0 * opt_maxid) &&
          /* mid */
          (100.0 * matches / (matches + mismatches) >= opt_mid) &&
          /* maxdiffs */
          (mismatches + internal_indels <= opt_maxdiffs))
        {
          si->hits[si->hit_count].target = target;
          si->hits[si->hit_count].count = kmercount;
          si->hits[si->hit_count].nwscore = nwscore;
          si->hits[si->hit_count].nwdiff = nwdiff;
          si->hits[si->hit_count].nwgaps = nwgaps;
          si->hits[si->hit_count].nwindels = nwindels;
          si->hits[si->hit_count].nwalignmentlength = nwalignmentlength;
          si->hits[si->hit_count].nwalignment = nwalignment;
          si->hits[si->hit_count].nwid = nwid;
          si->hits[si->hit_count].strand = si->strand;
          si->hits[si->hit_count].matches = matches;
          si->hits[si->hit_count].mismatches = mismatches;
          si->hits[si->hit_count].trim_q_left = trim_q_left;
          si->hits[si->hit_count].trim_q_right = trim_q_right;
          si->hits[si->hit_count].trim_t_left = trim_t_left;
          si->hits[si->hit_count].trim_t_right = trim_t_right;
          si->hits[si->hit_count].trim_aln_left = trim_aln_left;
          si->hits[si->hit_count].trim_aln_right = trim_aln_right;
          si->hits[si->hit_count].internal_alignmentlength =
            internal_alignmentlength;
          si->hits[si->hit_count].internal_gaps = internal_gaps;
          si->hits[si->hit_count].internal_indels = internal_indels;
          si->hits[si->hit_count].internal_id = internal_id;

          si->hit_count++;
          
          if (internal_id >= 100.0 * opt_id)
            si->accepts++;
        }
      else
        {
          free(nwalignment);
          si->rejects++;
        }
      
#if 0
      if ((si->accepts >= maxaccepts) || (si->rejects > maxrejects))
        {
          /* TODO: free alignments */
          break;
        }
#endif
    }
  //  printf("accepts: %d rejects: %d\n", si->accepts, si->rejects);
}

void search_onequery(struct searchinfo_s * si)
{
  si->hit_count = 0;

  for(int s = 0; s < opt_strand; s++)
    {
      si->strand = s;

      /* check plus or both strands*/
      /* s=0: plus; s=1: minus */
      char * qseq;
      if (s)
        qseq = si->rc;
      else
        qseq = si->qsequence;

      search16_qprep(si->s, qseq, si->qseqlen);

      /* extract unique kmer samples from query*/
      unique_count(si->uh, opt_wordlength, 
                   si->qseqlen, qseq,
                   & si->kmersamplecount, & si->kmersample);
      
      /* find database sequences with the most kmer hits */
      search_topscores(si);
      
      /* analyse targets with the highest number of kmer hits */
      si->accepts = 0;
      si->rejects = 0;
      si->candidatecount = 0;

      for(int t = 0;
          (si->accepts < maxaccepts) && 
            (si->rejects <= maxrejects) &&
            (!minheap_isempty(si->m));
          t++)
        {
          elem_t e = minheap_poplast(si->m);

#if 0
          printf("acc=%d rej=%d popped seqno, count: %d %d\n",
                 si->accepts, si->rejects, e.seqno, e.count);
#endif

          unsigned int target = e.seqno;
          unsigned int count = e.count;
          char * dlabel = db_getheader(target);
          char * dseq = db_getsequence(target);
          long dseqlen = db_getsequencelen(target);
          long tsize = db_getabundance(target);
          
          /* Test these accept/reject criteria before alignment */
          
          if (/* maxqsize */
              (si->qsize > opt_maxqsize) ||
              /* mintsize */
              (tsize < opt_mintsize) ||
              /* minsizeratio */
              (si->qsize < opt_minsizeratio * tsize) ||
              /* maxsizeratio */
              (si->qsize > opt_maxsizeratio * tsize) ||
              /* minqt */
              (si->qseqlen < opt_minqt * dseqlen) ||
              /* maxqt */
              (si->qseqlen > opt_maxqt * dseqlen) ||
              /* minsl */
              (si->qseqlen < dseqlen ? 
               si->qseqlen < opt_minsl * dseqlen : 
               dseqlen < opt_minsl * si->qseqlen) ||
              /* maxsl */
              (si->qseqlen < dseqlen ? 
               si->qseqlen > opt_maxsl * dseqlen : 
               dseqlen > opt_maxsl * si->qseqlen) ||
              /* idprefix */
              ((si->qseqlen < opt_idprefix) || 
               (dseqlen < opt_idprefix) ||
               (seqncmp(qseq,
                        dseq,
                        opt_idprefix))) ||
              /* idsuffix */
              ((si->qseqlen < opt_idsuffix) ||
               (dseqlen < opt_idsuffix) ||
               (seqncmp(qseq+si->qseqlen-opt_idsuffix,
                        dseq+dseqlen-opt_idsuffix,
                        opt_idsuffix))) ||
              /* self */
              (opt_self && (strcmp(si->query_head, dlabel) == 0)) ||
              /* selfid */
              (opt_selfid && (si->qseqlen == dseqlen) && 
               (seqncmp(qseq, dseq, si->qseqlen) == 0)))
            {
              si->rejects++;
            }
          else
            {
              si->candidate_target[si->candidatecount] = target;
              si->candidate_kmercount[si->candidatecount] = count;
              si->candidatecount++;
              if (si->candidatecount == 8)
                {
                  align_candidates(si);
                  si->candidatecount = 0;
                }
            }
        }  

      if (si->candidatecount > 0)
        {
          align_candidates(si);
          si->candidatecount = 0;
        }
    }
    
  /* sort hits */
  
  qsort(si->hits, si->hit_count, sizeof(struct hit), hit_compare);
}
