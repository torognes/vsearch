/*
    Copyright (C) 2014-2015 Torbjorn Rognes

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

/* Compute consensus sequence and msa of clustered sequences */

static char * aln;
static int alnpos;
static int * profile;

void msa_add(char c)
{
  int * p = profile + 4 * alnpos;

  switch(toupper(c))
    {
    case 'A':
      p[0] += 12;
      break;
    case 'C':
      p[1] += 12;
      break;
    case 'G':
      p[2] += 12;
      break;
    case 'T':
    case 'U':
      p[3] += 12;
      break;
    case 'R':
      p[0] += 6;
      p[2] += 6;
      break;
    case 'Y':
      p[1] += 6;
      p[3] += 6;
      break;
    case 'S':
      p[1] += 6;
      p[2] += 6;
      break;
    case 'W':
      p[0] += 6;
      p[3] += 6;
      break;
    case 'K':
      p[2] += 6;
      p[3] += 6;
      break;
    case 'M':
      p[0] += 6;
      p[1] += 6;
      break;
    case 'B':
      p[1] += 4;
      p[2] += 4;
      p[3] += 4;
      break;
    case 'D':
      p[0] += 4;
      p[2] += 4;
      p[3] += 4;
      break;
    case 'H':
      p[0] += 4;
      p[1] += 4;
      p[3] += 4;
      break;
    case 'V':
      p[0] += 4;
      p[1] += 4;
      p[2] += 4;
      break;
    case 'N':
      p[0] += 3;
      p[1] += 3;
      p[2] += 3;
      p[3] += 3;
      break;
    }   

  aln[alnpos++] = c;
}

void msa(FILE * fp_msaout, FILE * fp_consout, FILE * fp_profile,
         int cluster,
         int target_count, struct msa_target_s * target_list,
         long totalabundance, abundance_t * abundance_handle)
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
          long run = 1;
          int scanlength = 0;
          sscanf(p, "%ld%n", &run, &scanlength);
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
  profile = (int *) xmalloc(4 * sizeof(int) * alnlen);
  memset(profile, 0, 4 * sizeof(int) * alnlen);
  aln = (char *) xmalloc(alnlen+1);
  char * cons = (char *) xmalloc(alnlen+1);
  
  /* Find longest target sequence on reverse strand and allocate buffer */
  long longest_reversed = 0;
  for(int i=0; i < target_count; i++)
    if (target_list[i].strand)
      {
        long len = db_getsequencelen(target_list[i].seqno);
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
              long run = 1;
              int scanlength = 0;
              sscanf(p, "%ld%n", &run, &scanlength);
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
        {
          fprintf(fp_msaout, ">%s%s\n", j ? "" : "*", 
                  db_getheader(target_seqno));
          fprint_fasta_seq_only(fp_msaout, aln, alnlen, opt_fasta_width);
        }
    }  

  if (rc_buffer)
    free(rc_buffer);

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
          /* find most common symbol */
          char best_sym = 0;
          int best_count = -1;
          int nongap_count = 0;
          for(int c=0; c<4; c++)
            {
              int count = profile[4*i+c];
              if (count > best_count)
                {
                  best_count = count;
                  best_sym = c;
                }
              nongap_count += count;
            }

          int gap_count = 12 * target_count - nongap_count;

          if (best_count >= gap_count)
            {
              char sym = sym_nt_2bit[(int)best_sym];
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
    {
      fprint_fasta_hdr_only(fp_msaout, "consensus");
      fprint_fasta_seq_only(fp_msaout, aln, alnlen, opt_fasta_width);
    }

  if (fp_consout)
    {
      if (opt_sizeout)
        {
          /* must remove old size info first */
          char * header_wo_size 
            = abundance_strip_size(abundance_handle,
                                   db_getheader(centroid_seqno), 
                                   db_getheaderlen(centroid_seqno));
          fprintf(fp_consout,
                  ">centroid=%s;seqs=%d;size=%ld;",
                  header_wo_size,
                  target_count,
                  totalabundance);
          free(header_wo_size);
        }
      else
          fprintf(fp_consout,
                  ">centroid=%s;seqs=%d;",
                  db_getheader(centroid_seqno),
                  target_count);

      if (opt_clusterout_id)
        fprintf(fp_consout, "clusterid=%d;", cluster);

      fprintf(fp_consout, "\n");

      fprint_fasta_seq_only(fp_consout, cons, conslen, opt_fasta_width);
    }
  
  if (fp_profile)
    {
      /* must remove old size info first */
      if (opt_sizeout)
        {
          char * header_wo_size 
            = abundance_strip_size(abundance_handle,
                                   db_getheader(centroid_seqno), 
                                   db_getheaderlen(centroid_seqno));
          fprintf(fp_profile,
                  ">centroid=%s;seqs=%d;size=%ld;",
                  header_wo_size, 
                  target_count,
                  totalabundance);
          free(header_wo_size);
        }
      else
        fprintf(fp_profile,
                ">centroid=%s;seqs=%d;",
                db_getheader(centroid_seqno),
                target_count);

      if (opt_clusterout_id)
        fprintf(fp_profile, "clusterid=%d;", cluster);

      fprintf(fp_profile, "\n");

      for (int i=0; i<alnlen; i++)
        {
          fprintf(fp_profile, "%d\t%c", i, aln[i]);
          int nongap_count = 0;
          for (int c=0; c<4; c++)
            {
              int count = profile[4*i+c];
              nongap_count += count;
              if (count % 12 == 0)
                fprintf(fp_profile, "\t%d", count / 12);
              else
                fprintf(fp_profile, "\t%.2f", 1.0 * count / 12.0);
            }
          if (nongap_count % 12 == 0)
            fprintf(fp_profile, "\t%d", target_count - nongap_count / 12);
          else
            fprintf(fp_profile, "\t%.2f", 1.0 * target_count - 1.0 * nongap_count / 12.0);
          fprintf(fp_profile, "\n");
        }
      fprintf(fp_profile, "\n");
    }

  free(maxi);
  free(aln);
  free(cons);
  free(profile);
}
