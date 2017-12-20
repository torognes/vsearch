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

/* Compute consensus sequence and msa of clustered sequences */

#define PROFSIZE 6

static char * aln;
static int alnpos;
static int * profile;

void msa_add(char c)
{
  int * p = profile + PROFSIZE * alnpos;

  switch(toupper(c))
    {
    case 'A':
      p[0]++;
      break;
    case 'C':
      p[1]++;
      break;
    case 'G':
      p[2]++;
      break;
    case 'T':
    case 'U':
      p[3]++;
      break;
    case 'R':
    case 'Y':
    case 'S':
    case 'W':
    case 'K':
    case 'M':
    case 'B':
    case 'D':
    case 'H':
    case 'V':
    case 'N':
      p[4]++;
      break;
    case '-':
      p[5]++;
      break;
    }

  aln[alnpos++] = c;
}

void msa(FILE * fp_msaout, FILE * fp_consout, FILE * fp_profile,
         int cluster,
         int target_count, struct msa_target_s * target_list,
         int64_t totalabundance)
{
  int centroid_seqno = target_list[0].seqno;
  int centroid_len = db_getsequencelen(centroid_seqno);

  /* find max insertions in front of each position in the centroid sequence */
  int * maxi = (int *) xmalloc((centroid_len + 1) * sizeof(int));
  memset(maxi, 0, (centroid_len + 1) * sizeof(int));

  for(int j=1; j<target_count; j++)
    {
      char * p = target_list[j].cigar;
      char * e = p + strlen(p);
      int pos = 0;
      while (p < e)
        {
          int64_t run = 1;
          int scanlength = 0;
          sscanf(p, "%" PRId64 "%n", &run, &scanlength);
          p += scanlength;
          char op = *p++;
          switch (op)
            {
            case 'M':
            case 'I':
              pos += run;
              break;
            case 'D':
              if (run > maxi[pos])
                maxi[pos] = run;
              break;
            }
        }
    }  

  /* find total alignment length */
  int alnlen = 0;
  for(int i=0; i < centroid_len+1; i++)
    alnlen += maxi[i];
  alnlen += centroid_len;

  /* allocate memory for profile (for consensus) and aligned seq */
  profile = (int *) xmalloc(PROFSIZE * sizeof(int) * alnlen);
  memset(profile, 0, PROFSIZE * sizeof(int) * alnlen);
  aln = (char *) xmalloc(alnlen+1);
  char * cons = (char *) xmalloc(alnlen+1);
  
  /* Find longest target sequence on reverse strand and allocate buffer */
  int64_t longest_reversed = 0;
  for(int i=0; i < target_count; i++)
    if (target_list[i].strand)
      {
        int64_t len = db_getsequencelen(target_list[i].seqno);
        if (len > longest_reversed)
          longest_reversed = len;
      }
  char * rc_buffer = 0;
  if (longest_reversed > 0)
    rc_buffer = (char*) xmalloc(longest_reversed + 1);

  /* blank line before each msa */
  if (fp_msaout)
    fprintf(fp_msaout, "\n");
  
  for(int j=0; j<target_count; j++)
    {
      int target_seqno = target_list[j].seqno;
      char * target_seq = db_getsequence(target_seqno);
      
      if (target_list[j].strand)
        {
          reverse_complement(rc_buffer, target_seq,
                             db_getsequencelen(target_seqno));
          target_seq = rc_buffer;
        }

      int inserted = 0;
      int qpos = 0;
      int tpos = 0;
      alnpos = 0;

      if (!j)
        {
          for(int x=0; x < centroid_len; x++)
            {
              for(int y=0; y < maxi[qpos]; y++)
                msa_add('-');
              msa_add(target_seq[tpos++]);
              qpos++;
            }
        }
      else
        {
          char * p = target_list[j].cigar;
          char * e = p + strlen(p);
          while (p < e)
            {
              int64_t run = 1;
              int scanlength = 0;
              sscanf(p, "%" PRId64 "%n", &run, &scanlength);
              p += scanlength;
              char op = *p++;
              
              if (op == 'D')
                {
                  for(int x=0; x < maxi[qpos]; x++)
                    {
                      if (x < run)
                        msa_add(target_seq[tpos++]);
                      else
                        msa_add('-');
                    }
                  inserted = 1;
                }
              else
                {
                  for(int x=0; x < run; x++)
                    {
                      if (!inserted)
                        for(int y=0; y < maxi[qpos]; y++)
                          msa_add('-');
                      
                      if (op == 'M')
                        msa_add(target_seq[tpos++]);
                      else
                        msa_add('-');
                      
                      qpos++;
                      inserted = 0;
                    }
                }
            }
        }
      
      if (!inserted)
        for(int x=0; x < maxi[qpos]; x++)
          msa_add('-');
      
      /* end of sequence string */
      aln[alnpos] = 0;

      /* print header & sequence */
      if (fp_msaout)
        fasta_print_general(fp_msaout,
                            j ? "" : "*",
                            aln, 
                            alnlen,
                            db_getheader(target_seqno),
                            db_getheaderlen(target_seqno),
                            db_getabundance(target_seqno), 
                            0, -1, -1, 0, 0.0);
    }

  if (rc_buffer)
    xfree(rc_buffer);

  /* consensus */

  int conslen = 0;

  /* Censor part of the consensus sequence outside the centroid sequence */

  int left_censored = maxi[0];
  int right_censored = maxi[centroid_len];

  for(int i=0; i<alnlen; i++)
    {
      if ((i < left_censored) || (i >= alnlen - right_censored))
        {
          aln[i] = '+';
        }
      else
        {
          /* find most common symbol of A, C, G and T */
          char best_sym = 0;
          int best_count = 0;
          for(int c=0; c<4; c++)
            {
              int count = profile[PROFSIZE*i+c];
              if (count > best_count)
                {
                  best_count = count;
                  best_sym = c+1;
                }
            }

          /* if no A, C, G, or T, check if there are any N's */
          int n_count = profile[PROFSIZE*i+4];
          if ((best_count == 0) && (n_count > 0))
            {
              best_count = n_count;
              best_sym = 15; // N
            }

          /* compare to the number of gap symbols */
          int gap_count = profile[PROFSIZE*i+5];
          if (best_count >= gap_count)
            {
              char sym = sym_nt_4bit[(int)best_sym];
              aln[i] = sym;
              cons[conslen++] = sym;
            }
          else
            aln[i] = '-';
        }
    }

  aln[alnlen] = 0;
  cons[conslen] = 0;

  if (fp_msaout)
    fasta_print(fp_msaout, "consensus", aln, alnlen);

  if (fp_consout)
    {
      fasta_print_general(fp_consout,
                          "consensus=",
                          cons,
                          conslen,
                          db_getheader(centroid_seqno),
                          db_getheaderlen(centroid_seqno),
                          totalabundance,
                          cluster+1,
                          target_count,
                          opt_clusterout_id ? cluster : -1,
                          0, 0.0);
    }
  
  if (fp_profile)
    {
      fasta_print_general(fp_profile,
                          "consensus=",
                          0,
                          0,
                          db_getheader(centroid_seqno),
                          db_getheaderlen(centroid_seqno),
                          totalabundance,
                          cluster+1,
                          target_count,
                          opt_clusterout_id ? cluster : -1,
                          0, 0.0);

      for (int i=0; i<alnlen; i++)
        {
          fprintf(fp_profile, "%d\t%c", i, aln[i]);
          // A, C, G and T
          for (int c=0; c<4; c++)
            fprintf(fp_profile, "\t%d", profile[PROFSIZE*i+c]);
          // Gap symbol
          fprintf(fp_profile, "\t%d", profile[PROFSIZE*i+5]);
          // Ambiguous nucleotide (Ns and others)
          fprintf(fp_profile, "\t%d", profile[PROFSIZE*i+4]);
          fprintf(fp_profile, "\n");
        }
      fprintf(fp_profile, "\n");
    }

  xfree(maxi);
  xfree(aln);
  xfree(cons);
  xfree(profile);
}
