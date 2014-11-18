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

/*
  This code implements the method described in this paper:

  Robert C. Edgar, Brian J. Haas, Jose C. Clemente, Christopher Quince
  and Rob Knight (2011)
  UCHIME improves sensitivity and speed of chimera detection
  Bioinformatics, 27, 16, 2194-2200
  http://dx.doi.org/10.1093/bioinformatics/btr381
*/

const int parts = 4;
const int few = 1;
const int maxcandidates = few * parts;

static struct chimera_info_s
{
  int qalloc; /* the longest query we allocated memory for */

  int query_no;
  char * query_head;
  int query_head_len;
  int query_size;
  char * query_seq;
  int query_len;

  struct searchinfo_s si[parts];

  unsigned int cand_list[maxcandidates];
  int cand_count;

  struct s16info_s * s;

  CELL nwscore[maxcandidates];
  unsigned short nwalignmentlength[maxcandidates];
  unsigned short nwmatches[maxcandidates];
  unsigned short nwmismatches[maxcandidates];
  unsigned short nwgaps[maxcandidates];
  char * nwcigar[maxcandidates];

  int match_size;
  int * match;
  int * smooth;
  int * maxsmooth;

  int best_parents[2];

  FILE * fp_chimeras;
  FILE * fp_nonchimeras;
  FILE * fp_uchimealns;
  FILE * fp_uchimeout;

  int chimeras;
  int tophits;
  
  int best_target;
  char * best_cigar;

  /* dynamic arrays */

  int * maxi;
  char * paln[2];
  char * qaln;
  char * diffs;
  char * votes;
  char * model;
  char * ignore;

} ci;

void realloc_arrays(int maxqlen)
{
  /* realloc arrays based on query length */

  if (maxqlen > ci.qalloc)
    {
      ci.qalloc = maxqlen;
      
      for(int i=0; i < parts; i++)
        {
          int maxpartlen = (maxqlen + parts - 1) / parts;
          ci.si[i].qsequence = (char*) xrealloc(ci.si[i].qsequence,
                                                maxpartlen + 1);
        }
      
      ci.maxi = (int *) xrealloc(ci.maxi, (maxqlen + 1) * sizeof(int));
      ci.maxsmooth = (int*) xrealloc(ci.maxsmooth, maxqlen * sizeof(int));
      ci.match = (int*) xrealloc(ci.match,
                                 maxcandidates * maxqlen * sizeof(int));
      ci.smooth = (int*) xrealloc(ci.smooth,
                                  maxcandidates * maxqlen * sizeof(int));
      
      int maxalnlen = maxqlen + db_getsequencecount();
      
      ci.paln[0] = (char*) xrealloc(ci.paln[0], maxalnlen+1);
      ci.paln[1] = (char*) xrealloc(ci.paln[1], maxalnlen+1);
      ci.qaln = (char*) xrealloc(ci.qaln, maxalnlen+1);
      ci.diffs = (char*) xrealloc(ci.diffs, maxalnlen+1);
      ci.votes = (char*) xrealloc(ci.votes, maxalnlen+1);
      ci.model = (char*) xrealloc(ci.model, maxalnlen+1);
      ci.ignore = (char*) xrealloc(ci.ignore, maxalnlen+1);
    }
}


int find_best_parents()
{
  ci.best_parents[0] = -1;
  ci.best_parents[1] = -1;

  /* find the positions with matches for each potential parent */

  char * qseq = ci.query_seq;

  memset(ci.match, 0, ci.cand_count * ci.query_len * sizeof(int));

  for(int i=0; i < ci.cand_count; i++)
    {
      char * tseq = db_getsequence(ci.cand_list[i]);

      int qpos = 0;
      int tpos = 0;
      
      char * p = ci.nwcigar[i];
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
              for(int k=0; k<run; k++)
                {
                  if (chrmap_4bit[(int)(qseq[qpos])] ==
                      chrmap_4bit[(int)(tseq[tpos])])
                    ci.match[i * ci.query_len + qpos] = 1;
                  qpos++;
                  tpos++;
                }
              break;
              
            case 'I':
              tpos += run;
              break;
              
            case 'D':
              qpos += run;
              break;
            }
        }
    }
  
  /* compute smoothed identity score in 32-wide window for each candidate */
  /* and record max smoothed score for each position among candidates     */

  memset(ci.maxsmooth, 0, ci.query_len * sizeof(int));

  const int window = 32;

  for(int i = 0; i < ci.cand_count; i++)
    {
      int sum = 0;
      for(int qpos = 0; qpos < ci.query_len; qpos++)
        {
          int z = i * ci.query_len + qpos;
          sum += ci.match[z];
          if (qpos >= window)
            sum -= ci.match[z-window];
          if (qpos >= window-1)
            {
              ci.smooth[z] = sum;
              if (ci.smooth[z] > ci.maxsmooth[qpos])
                ci.maxsmooth[qpos] = ci.smooth[z];
            }
        }
    }
  
  /* find first parent */

  int wins[ci.cand_count];

  memset(wins, 0, ci.cand_count * sizeof(int));

  for(int qpos = window-1; qpos < ci.query_len; qpos++)
    {
      if (ci.maxsmooth[qpos] != 0)
        for(int i=0; i < ci.cand_count; i++)
          {
            int z = i * ci.query_len + qpos;
            if (ci.smooth[z] == ci.maxsmooth[qpos])
              wins[i]++;
          }
    }

  int best1_w = -1;
  int best1_i = -1;
  int best2_w = -1;
  int best2_i = -1;

  for(int i=0; i < ci.cand_count; i++)
    {
      int w = wins[i];
      if (w > best1_w)
        {
          best1_w = w;
          best1_i = i;
        }
    }

  if (best1_w >= 0)
    {
      /* find second parent */

      /* wipe out matches in positions covered by first parent */
      
      for(int qpos = window - 1; qpos < ci.query_len; qpos++)
        {
          int z = best1_i * ci.query_len + qpos;
          if (ci.smooth[z] == ci.maxsmooth[qpos])
            {
              for(int i = qpos + 1 - window; i <= qpos; i++)
                for(int j = 0; j < ci.cand_count; j++)
                  {
                    int z = j * ci.query_len + i;
                    ci.match[z] = 0;
                  }
            }
        }

      /*
        recompute smoothed identity over window, and record max smoothed
        score for each position among remaining candidates
      */
      
      memset(ci.maxsmooth, 0, ci.query_len * sizeof(int));

      for(int i = 0; i < ci.cand_count; i++)
        {
          if (i != best1_i)
            {
              int sum = 0;
              for(int qpos = 0; qpos < ci.query_len; qpos++)
                {
                  int z = i * ci.query_len + qpos;
                  sum += ci.match[z];
                  if (qpos >= window)
                    sum -= ci.match[z-window];
                  if (qpos >= window-1)
                    {
                      ci.smooth[z] = sum;
                      if (ci.smooth[z] > ci.maxsmooth[qpos])
                        ci.maxsmooth[qpos] = ci.smooth[z];
                    }
                }
            }
        }
  
      /* find second parent */

      memset(wins, 0, ci.cand_count * sizeof(int));

      for(int qpos = window-1; qpos < ci.query_len; qpos++)
        {
          if (ci.maxsmooth[qpos] != 0)
            for(int i=0; i < ci.cand_count; i++)
              if (i != best1_i)
                {
                  int z = i * ci.query_len + qpos;
                  if (ci.smooth[z] == ci.maxsmooth[qpos])
                    wins[i]++;
                }
        }

      for(int i=0; i < ci.cand_count; i++)
        {
          int w = wins[i];
          if (w > best2_w)
            {
              best2_w = w;
              best2_i = i;
            }
        }
    }

  ci.best_parents[0] = best1_i;
  ci.best_parents[1] = best2_i;
  
  return (best1_w >= 0) && (best2_w >= 0);
}

int eval_parents()
{
  /* create msa */

  /* find max insertions in front of each position in the query sequence */
  memset(ci.maxi, 0, (ci.query_len + 1) * sizeof(int));

  for(int j=0; j<2; j++)
    {
      char * p = ci.nwcigar[ci.best_parents[j]];
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
              if (run > ci.maxi[pos])
                ci.maxi[pos] = run;
              break;
            }
        }
    }  

  /* find total alignment length */
  int alnlen = 0;
  for(int i=0; i < ci.query_len+1; i++)
    alnlen += ci.maxi[i];
  alnlen += ci.query_len;

  /* fill in alignment string for query */

  char * p = ci.qaln;
  int qpos = 0;
  for (int i=0; i < ci.query_len; i++)
    {
      for (int j=0; j < ci.maxi[i]; j++)
        *p++ = '-';
      *p++ = ci.query_seq[qpos++];
    }
  for (int j=0; j < ci.maxi[ci.query_len]; j++)
    *p++ = '-';
  *p = 0;
  
  /* fill in alignment strings for the 2 parents */

  for(int j=0; j<2; j++)
    {
      int cand = ci.best_parents[j];
      int target_seqno = ci.cand_list[cand];
      char * target_seq = db_getsequence(target_seqno);
      
      int inserted = 0;
      int qpos = 0;
      int tpos = 0;

      char * t = ci.paln[j];

      char * p = ci.nwcigar[cand];
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
              for(int x=0; x < ci.maxi[qpos]; x++)
                {
                  if (x < run)
                    *t++ = target_seq[tpos++];
                  else
                    *t++ = '-';
                }
              inserted = 1;
            }
          else
            {
              for(int x=0; x < run; x++)
                {
                  if (!inserted)
                    for(int x=0; x < ci.maxi[qpos]; x++)
                      *t++ = '-';
                      
                  if (op == 'M')
                    *t++ = target_seq[tpos++];
                  else
                    *t++ = '-';
                      
                  qpos++;
                  inserted = 0;
                }
            }
        }
      
      /* add any gaps at the end */

      if (!inserted)
        for(int x=0; x < ci.maxi[qpos]; x++)
          *t++ = '-';
      
      /* end of sequence string */
      *t = 0;
    }

  memset(ci.ignore, 0, alnlen);

  for(int i = 0; i < alnlen; i++)
    {
      char qsym  = chrmap_4bit[(int)(ci.qaln   [i])];
      char p1sym = chrmap_4bit[(int)(ci.paln[0][i])];
      char p2sym = chrmap_4bit[(int)(ci.paln[1][i])];
          
      /* mark positions to ignore in voting */

      /* ignore gap positions and those next to the gap */
      if ((!qsym) || (!p1sym) || (!p2sym))
        {
          ci.ignore[i] = 1;
          if (i>0)
            ci.ignore[i-1] = 1;
          if (i<alnlen-1)
            ci.ignore[i+1] = 1;
        }

      /* ignore ambigous symbols */
      if ((qsym>4) || (p1sym>4) || (p2sym>4))
        ci.ignore[i] = 1;

      /* lower case parent symbols that differ from query */

      if (p1sym && (p1sym != qsym))
        ci.paln[0][i] = tolower(ci.paln[0][i]);

      if (p2sym && (p2sym != qsym))
        ci.paln[1][i] = tolower(ci.paln[1][i]);

      /* compute diffs */

      char diff;

      if (qsym && p1sym && p2sym)
        {
          if (p1sym == p2sym)
            {
              if (qsym == p1sym)
                diff = ' ';
              else
                diff = 'N';
            }
          else
            {
              if (qsym == p1sym)
                diff = 'A';
              else if (qsym == p2sym)
                diff = 'B';
              else
                diff = '?';
            }
        }
      else
        diff = ' ';

      ci.diffs[i] = diff;
    }

  ci.diffs[alnlen] = 0;

  /* compute score */

  int sumA = 0;
  int sumB = 0;
  int sumN = 0;

  for (int i = 0; i < alnlen; i++)
    if (!ci.ignore[i])
      {
        char diff = ci.diffs[i];
        
        if (diff == 'A')
          sumA++;
        else if (diff == 'B')
          sumB++;
        else if (diff != ' ')
          sumN++;
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

  for (int i=0; i<alnlen; i++)
  {
    if(!ci.ignore[i])
      {
        char diff = ci.diffs[i];
        if (diff != ' ')
          {
            if (diff == 'A')
              {
                left_y++;
                right_n--;
              }
            else if (diff == 'B')
              {
                left_n++;
                right_y--;
              }
            else
              {
                left_a++;
                right_a--;
              }
            
            double left_h, right_h, h;
            
            if ((left_y > left_n) && (right_y > right_n))
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
            else if ((left_n > left_y) && (right_n > right_y))
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

  int ischimeric = 0;
  int isnonchimeric = 0;

  if (best_h >= 0.0)
    {
      /* flip A and B if necessary */

      if (best_reverse)
        for(int i = 0; i < alnlen; i++)
          {
            char diff = ci.diffs[i];
            if (diff == 'A')
              ci.diffs[i] = 'B';
            else if (diff == 'B')
              ci.diffs[i] = 'A';
          }

      /* fill in votes and model */

      for(int i = 0; i < alnlen; i++)
        {
          char m = i <= best_i ? 'A' : 'B';
          ci.model[i] = m;
          
          char v = ' ';
          if (!ci.ignore[i])
            {
              char d = ci.diffs[i];
          
              if ((d == 'A') || (d == 'B'))
                {
                  if (d == m)
                    v = '+';
                  else
                    v = '!';
                }
              else if ((d == 'N') || (d == '?'))
                {
                  v = '0';
                }
            }
          ci.votes[i] = v;

          /* lower case diffs for no votes */
          if (v == '!')
            ci.diffs[i] = tolower(ci.diffs[i]);
        }

      /* fill in crossover region */

      for(int i = best_i + 1; i < alnlen; i++)
        if ((ci.diffs[i] == ' ') || (ci.diffs[i] == 'A'))
          ci.model[i] = 'x';
        else
          break;

      ci.votes[alnlen] = 0;
      ci.model[alnlen] = 0;
  
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
          if (! ci.ignore[i])
            {
              cols++;

              char qsym = chrmap_4bit[(int)(ci.qaln[i])];
              char asym = chrmap_4bit[(int)(ci.paln[index_a][i])];
              char bsym = chrmap_4bit[(int)(ci.paln[index_b][i])];
              char msym = i <= best_i ? asym : bsym;

              if (qsym == asym)
                match_QA++;
              
              if (qsym == bsym)
                match_QB++;
              
              if (asym == bsym)
                match_AB++;

              if (qsym == msym)
                match_QM++;
            }
        }

      int seqno_a = ci.cand_list[ci.best_parents[index_a]];
      int seqno_b = ci.cand_list[ci.best_parents[index_b]];

      double QA = 100.0 * match_QA / cols;
      double QB = 100.0 * match_QB / cols;
      double AB = 100.0 * match_AB / cols;
      double QT = MAX(QA, QB);
      double QM = 100.0 * match_QM / cols;
      double divdiff = QM - QT;
      double divfrac = 100.0 * divdiff / QT;

      int sumL = best_left_n + best_left_a + best_left_y;
      int sumR = best_right_n + best_right_a + best_right_y;

      if (best_h >= opt_minh)
        {
          if ((divdiff >= opt_mindiv) &&
              (sumL >= opt_mindiffs) &&
              (sumR >= opt_mindiffs))
            ischimeric = 1;
        }
      else
        isnonchimeric = 1;

      /* print alignment */

      if (opt_uchimealns && ischimeric)
        {
          fprintf(ci.fp_uchimealns, "\n");
          fprintf(ci.fp_uchimealns, "----------------------------------------"
                  "--------------------------------\n");
          fprintf(ci.fp_uchimealns, "Query   (%5d nt) %s\n",
                  ci.query_len, ci.query_head);
          fprintf(ci.fp_uchimealns, "ParentA (%5lu nt) %s\n",
                  db_getsequencelen(seqno_a),
                  db_getheader(seqno_a));
          fprintf(ci.fp_uchimealns, "ParentB (%5lu nt) %s\n",
                  db_getsequencelen(seqno_b),
                  db_getheader(seqno_b));
          fprintf(ci.fp_uchimealns, "\n");

          int width = 80;
          qpos = 0;
          int p1pos = 0;
          int p2pos = 0;
          int rest = alnlen;

          for(int i = 0; i < alnlen; i += width)
            {
              /* count non-gap symbols on current line */

              int qnt, p1nt, p2nt;
              qnt = p1nt = p2nt = 0;

              int w = MIN(rest,width);

              for(int j=0; j<w; j++)
                {
                  if (ci.qaln[i+j] != '-')
                    qnt++;
                  if (ci.paln[0][i+j] != '-')
                    p1nt++;
                  if (ci.paln[1][i+j] != '-')
                    p2nt++;
                }

              if (! best_reverse)
                {
                  fprintf(ci.fp_uchimealns, "A %5d %.*s %d\n",
                          p1pos+1, w, ci.paln[0]+i, p1pos+p1nt);
                  fprintf(ci.fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos+1,  w, ci.qaln+i,    qpos+qnt);
                  fprintf(ci.fp_uchimealns, "B %5d %.*s %d\n",
                          p2pos+1, w, ci.paln[1]+i, p2pos+p2nt);
                }
              else
                {
                  fprintf(ci.fp_uchimealns, "A %5d %.*s %d\n",
                          p2pos+1, w, ci.paln[1]+i, p2pos+p2nt);
                  fprintf(ci.fp_uchimealns, "Q %5d %.*s %d\n",
                          qpos+1,  w, ci.qaln+i,    qpos+qnt);
                  fprintf(ci.fp_uchimealns, "B %5d %.*s %d\n",
                          p1pos+1, w, ci.paln[0]+i, p1pos+p1nt);
                }

              fprintf(ci.fp_uchimealns, "Diffs   %.*s\n", w, ci.diffs+i);
              fprintf(ci.fp_uchimealns, "Votes   %.*s\n", w, ci.votes+i);
              fprintf(ci.fp_uchimealns, "Model   %.*s\n", w, ci.model+i);
              fprintf(ci.fp_uchimealns, "\n");

              qpos += qnt;
              p1pos += p1nt;
              p2pos += p2nt;
              rest -= width;
            }

          fprintf(ci.fp_uchimealns, "Ids.  QA %.1f%%, QB %.1f%%, AB %.1f%%, "
                  "QModel %.1f%%, Div. %+.1f%%\n",
                  QA, QB, AB, QM, divfrac);

          fprintf(ci.fp_uchimealns, "Diffs Left %d: N %d, A %d, Y %d (%.1f%%); "
                  "Right %d: N %d, A %d, Y %d (%.1f%%), Score %.4f\n",
                  sumL, best_left_n, best_left_a, best_left_y,
                  100.0 * best_left_y / sumL,
                  sumR, best_right_n, best_right_a, best_right_y,
                  100.0 * best_right_y / sumR,
                  best_h);
        }

      if (opt_uchimeout)
        {
          if (opt_uchimeout5)
            {
              fprintf(ci.fp_uchimeout,
                      "%.4f\t%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t"
                      "%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%c\n",
                      best_h,
                      ci.query_head,
                      db_getheader(seqno_a),
                      db_getheader(seqno_b),
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
                      ischimeric ? 'Y' : (isnonchimeric ? 'N' : '?'));
            }
          else
            {
              fprintf(ci.fp_uchimeout,
                      "%.4f\t%s\t%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t"
                      "%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%c\n",
                      best_h,
                      ci.query_head,
                      db_getheader(seqno_a),
                      db_getheader(seqno_b),
                      db_getheader(QA >= QB ? seqno_a : seqno_b),
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
                      ischimeric ? 'Y' : (isnonchimeric ? 'N' : '?'));
            }
        }
    }
  else
    {
      isnonchimeric = 1;

      if (opt_uchimeout)
        {
          if (opt_uchimeout5)
            fprintf(ci.fp_uchimeout,
                    "0.0000\t%s\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n",
                    ci.query_head);
          else
            fprintf(ci.fp_uchimeout, 
                    "0.0000\t%s\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n",
                    ci.query_head);
        }
    }


  if (ischimeric)
    return 1;
  else if (isnonchimeric)
    return -1;
  else
    return 0;
}

void query_init(struct searchinfo_s * si)
{
  si->qsequence = 0;
  si->kmers = 0;
  si->hits = (struct hit *) xmalloc(sizeof(struct hit) * ci.tophits);
  si->kmers = (count_t *) xmalloc(db_getsequencecount() * 
                                  sizeof(count_t) + 32);
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
  si->m = minheap_init(ci.tophits);
}

void query_exit(struct searchinfo_s * si)
{
  search16_exit(si->s);
  unique_exit(si->uh);
  minheap_exit(si->m);
  nw_exit(si->nw);
  
  if (si->qsequence)
    free(si->qsequence);
  if (si->hits)
    free(si->hits);
  if (si->kmers)
    free(si->kmers);
}

void partition_query(int parts)
{
  int rest = ci.query_len;
  char * p = ci.query_seq;
  for (int i=0; i<parts; i++)
    {
      int len = (rest+(parts-1-i))/(parts-i);

      struct searchinfo_s * si = ci.si + i;

      si->query_no = ci.query_no;
      si->strand = 0;
      si->qsize = ci.query_size;
      si->query_head_len = ci.query_head_len;
      si->query_head = ci.query_head;
      si->qseqlen = len;
      strncpy(si->qsequence, p, len);
      si->qsequence[len] = 0;

      rest -= len;
      p += len;
    }

}

void open_chimera_file(FILE * * f, char * name)
{
  if (name)
    {
      *f = fopen(name, "w");
      if (!*f)
        fatal("Unable to open file %s for writing", name);
    }
  else
    *f = 0;
}

void close_chimera_file(FILE * f)
{
  if (f)
    fclose(f);
}

inline int hit2_compare_typed(struct hit * x, struct hit * y)
{
  // high id, then low id
  // early target (high abundance), then late target (low abundance)

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

int hit2_compare(const void * a, const void * b)
{
  return hit2_compare_typed((struct hit *) a, (struct hit *) b);
}

int cand_compare(const void * a, const void * b)
{
  int x = *(int*)a;
  int y = *(int*)b;

  if (ci.nwmatches[x] > ci.nwmatches[y])
    return -1;
  else if (ci.nwmatches[x] < ci.nwmatches[y])
    return +1;
  else if (x < y)
    return -1;
  else if (x > y)
    return +1;
  else
    return 0;
}

void chimera()
{
  open_chimera_file(&ci.fp_chimeras, opt_chimeras);
  open_chimera_file(&ci.fp_nonchimeras, opt_nonchimeras);
  open_chimera_file(&ci.fp_uchimealns, opt_uchimealns);
  open_chimera_file(&ci.fp_uchimeout, opt_uchimeout);

  /* override any options the user might have set */
  opt_maxaccepts = few;
  opt_maxrejects = 1;
  opt_id = 0.0;
  opt_strand = 1;

  if (opt_uchime_denovo)
    opt_maxsizeratio = 1.0 / opt_abskew;

  ci.tophits = opt_maxaccepts + opt_maxrejects;

  int progress_total;

  int queries = 0;
  ci.chimeras = 0;
  
  ci.qalloc = 0;
  ci.maxi = 0;
  ci.maxsmooth = 0;
  ci.match = 0;
  ci.smooth = 0;
  ci.paln[0] = 0;
  ci.paln[1] = 0;
  ci.qaln = 0;
  ci.diffs = 0;
  ci.votes = 0;
  ci.model = 0;
  ci.ignore = 0;
  
  struct hit * allhits_list = (struct hit *) xmalloc(maxcandidates * sizeof(struct hit));
  
  /* prepare queries / database */

  if (opt_uchime_ref)
    {
      db_read(opt_db, 1);
      dbindex_prepare(1);
      dbindex_addallsequences();
      query_open(opt_uchime_ref);
      progress_total = query_getfilesize();
    }
  else
    {
      db_read(opt_uchime_denovo, 1);
      db_sortbyabundance();
      dbindex_prepare(1);
      progress_total = db_getnucleotidecount();
    }

  for(int i = 0; i < parts; i++)
    query_init(ci.si + i);

  int progress = 0;

  ci.s = search16_init(opt_match,
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
  
  progress_init("Detecting chimeras", progress_total);

  int seqno = 0;
  int seqcount = db_getsequencecount();

  while(1)
    {
      if (opt_uchime_ref)
        {
          if (!query_getnext(& ci.query_head, & ci.query_head_len,
                             & ci.query_seq, & ci.query_len, & ci.query_no,
                             & ci.query_size, 1))
            break;
        }
      else
        {
          if (seqno >= seqcount)
            break;

          ci.query_head = db_getheader(seqno);
          ci.query_head_len = db_getheaderlen(seqno);
          ci.query_seq = db_getsequence(seqno);
          ci.query_len = db_getsequencelen(seqno);
          ci.query_no = seqno;
          ci.query_size = db_getabundance(seqno);
        }

      /* if necessary, expand memory for arrays dependent on query length */
      realloc_arrays(ci.query_len);
      
      int status = 0;

      /* partition query */
      partition_query(parts);

      /* perform searches and collect candidate parents */
      ci.cand_count = 0;
      int allhits_count = 0;

      for (int i=0; i<parts; i++)
        {
          struct hit * hits;
          int hit_count;
          search_onequery(ci.si+i);
          search_joinhits(ci.si+i, 0, & hits, & hit_count);
          for(int j=0; j<hit_count; j++)
            if (hits[j].accepted)
              allhits_list[allhits_count++] = hits[j];
          free(hits);
        }

      for(int i=0; i < allhits_count; i++)
        {
          unsigned int target = allhits_list[i].target;

          /* skip duplicates */
          int k;
          for(k = 0; k < ci.cand_count; k++)
            if (ci.cand_list[k] == target)
              break;

          if (k == ci.cand_count)
            ci.cand_list[ci.cand_count++] = target;

          /* deallocate cigar */
          if (allhits_list[i].nwalignment)
            free(allhits_list[i].nwalignment);
        }


      /* align full query to each candidate */

      search16_qprep(ci.s, ci.query_seq, ci.query_len);

      search16(ci.s,
               ci.cand_count,
               ci.cand_list,
               ci.nwscore,
               ci.nwalignmentlength,
               ci.nwmatches,
               ci.nwmismatches,
               ci.nwgaps,
               ci.nwcigar);
      
      /* find the best pair of parents */

      if (find_best_parents())
        {
          /* compute score for parents */
          status = eval_parents();
        }
      else
        {
          status = -1;
          /* output no parents, no chimeras */
          if (opt_uchimeout)
            {
              if (opt_uchimeout5)
                fprintf(ci.fp_uchimeout,
                        "0.0000\t%s\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n",
                        ci.query_head);
              else
                fprintf(ci.fp_uchimeout,
                        "0.0000\t%s\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n",
                        ci.query_head);
            }
        }

      if (status > 0)
        {
          ci.chimeras++;

          if (opt_chimeras)
            {
              fprint_fasta_hdr_only(ci.fp_chimeras, ci.query_head);
              fprint_fasta_seq_only(ci.fp_chimeras, ci.query_seq,
                                    ci.query_len, opt_fasta_width);
            }
        }

      if (status < 0)
        {
          /* uchime_denovo: add non-chimeras to db */
          if (opt_uchime_denovo)
            dbindex_addsequence(seqno);

          if (opt_nonchimeras)
            {
              fprint_fasta_hdr_only(ci.fp_nonchimeras, ci.query_head);
              fprint_fasta_seq_only(ci.fp_nonchimeras, ci.query_seq,
                                    ci.query_len, opt_fasta_width);
            }
        }
      
      for (int i=0; i < ci.cand_count; i++)
        if (ci.nwcigar[i])
          free(ci.nwcigar[i]);
      
      seqno++;
      queries++;

      if (opt_uchime_ref)
        progress = query_getfilepos();
      else
        progress += db_getsequencelen(seqno);

      progress_update(progress);
    }
  progress_done();

  for(int i = 0; i < parts; i++)
    query_exit(ci.si + i);
  
  fprintf(stderr, "Found %d chimeras in %d sequences (%.1f%%)\n",
          ci.chimeras, queries, 100.0 * ci.chimeras / queries);


  close_chimera_file(ci.fp_uchimeout);
  close_chimera_file(ci.fp_uchimealns);
  close_chimera_file(ci.fp_nonchimeras);
  close_chimera_file(ci.fp_chimeras);

  search16_exit(ci.s);

  if (allhits_list)
    free(allhits_list);

  if (ci.maxsmooth)
    free(ci.maxsmooth);
  if (ci.match)
    free(ci.match);
  if (ci.smooth)
    free(ci.smooth);
  if (ci.diffs)
    free(ci.diffs);
  if (ci.votes)
    free(ci.votes);
  if (ci.model)
    free(ci.model);
  if (ci.ignore)
    free(ci.ignore);
  if (ci.maxi)
    free(ci.maxi);
  if (ci.qaln)
    free(ci.qaln);
  if (ci.paln[0])
    free(ci.paln[0]);
  if (ci.paln[1])
    free(ci.paln[1]);
  
  if (opt_uchime_ref)
    query_close();
  
  dbindex_free();
  db_free();

  show_rusage();
}
