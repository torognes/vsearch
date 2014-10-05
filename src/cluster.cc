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

//#define COMPARENONVECTORIZED

static int tophits; /* the maximum number of hits to keep */
static int seqcount; /* number of database sequences */

typedef struct clusterinfo_s
{
  int seqno;
  int clusterno;
} clusterinfo_t;

int compare_bylength(const void * a, const void * b)
{
  seqinfo_t * x = (seqinfo_t *) a;
  seqinfo_t * y = (seqinfo_t *) b;

  /* longest first, otherwise keep order */

  if (x->seqlen < y->seqlen)
    return +1;
  else if (x->seqlen > y->seqlen)
    return -1;
  else
    if (x < y)
      return -1;
    else if (x > y)
      return +1;
    else
      return 0;
}

int compare_byclusterno(const void * a, const void * b)
{
  clusterinfo_t * x = (clusterinfo_t *) a;
  clusterinfo_t * y = (clusterinfo_t *) b;
  if (x->clusterno < y->clusterno)
    return -1;
  else if (x->clusterno > y->clusterno)
    return +1;
  else if (x->seqno < y->seqno)
    return -1;
  else if (x->seqno > y->seqno)
    return +1;
  else
    return 0;
}


void cluster(char * dbname, 
             char * cmdline,
             char * progheader,
             int sortbylength)
{
  FILE * fp_centroids = 0;
  FILE * fp_uc = 0;
  FILE * fp_alnout = 0;
  FILE * fp_userout = 0;
  FILE * fp_blast6out = 0;
  FILE * fp_fastapairs = 0;
  FILE * fp_matched = 0;
  FILE * fp_notmatched = 0;
  
  if (opt_centroids)
    {
      fp_centroids = fopen(opt_centroids, "w");
      if (!fp_centroids)
        fatal("Unable to open centroids file for writing");
    }

  if (ucfilename)
    {
      fp_uc = fopen(ucfilename, "w");
      if (!fp_uc)
        fatal("Unable to open uc file for writing");
    }

  if (opt_alnout)
    {
      fp_alnout = fopen(opt_alnout, "w");
      if (! fp_alnout)
        fatal("Unable to open alignment output file for writing");

      fprintf(fp_alnout, "%s\n", cmdline);
      fprintf(fp_alnout, "%s\n", progheader);
    }

  if (useroutfilename)
    {
      fp_userout = fopen(useroutfilename, "w");
      if (! fp_userout)
        fatal("Unable to open user-defined output file for writing");
    }

  if (opt_blast6out)
    {
      fp_blast6out = fopen(opt_blast6out, "w");
      if (! fp_blast6out)
        fatal("Unable to open blast6-like output file for writing");
    }

  if (opt_fastapairs)
    {
      fp_fastapairs = fopen(opt_fastapairs, "w");
      if (! fp_fastapairs)
        fatal("Unable to open fastapairs output file for writing");
    }

  if (opt_matched)
    {
      fp_matched = fopen(opt_matched, "w");
      if (! fp_matched)
        fatal("Unable to open matched output file for writing");
    }

  if (opt_notmatched)
    {
      fp_notmatched = fopen(opt_notmatched, "w");
      if (! fp_notmatched)
        fatal("Unable to open notmatched output file for writing");
    }

  db_read(dbname, opt_dbmask != MASK_SOFT);

  if (opt_dbmask == MASK_DUST)
    dust_all();
  else if ((opt_dbmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();
  
  show_rusage();
  
  seqcount = db_getsequencecount();
  
  if (sortbylength)
    {
      progress_init("Sorting by length", 100);
      qsort(seqindex, seqcount, sizeof(seqinfo_t), compare_bylength);
      progress_done();
    }
  
  dbindex_prepare(1);
  
  if ((maxrejects == 0) || (maxrejects > seqcount))
    maxrejects = seqcount;
  
  if (maxaccepts > seqcount)
    maxaccepts = seqcount;
  
  tophits = maxrejects + maxaccepts + 8;
  if (tophits > seqcount)
    tophits = seqcount;

  struct searchinfo_s si[1];

  si->uh = unique_init();
  si->kmers = (count_t *) xmalloc(seqcount * sizeof(count_t));
  si->m = minheap_init(tophits);
  si->targetlist = (unsigned int*) xmalloc(sizeof(unsigned int)*seqcount);
  si->hits = (struct hit *) xmalloc(sizeof(struct hit) * (tophits) * opt_strand);

  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = 0;
  si->seq_alloc = 0;
  si->qsequence = 0;
  
  if (opt_strand > 1)
    {
      si->rc_seq_alloc = db_getlongestsequence() + 1;
      si->rc = (char *) xmalloc(si->rc_seq_alloc);
    }
  else
    {
      si->rc_seq_alloc = 0;
      si->rc = 0;
    }

#ifdef COMPARENONVECTORIZED
  si->nw = nw_init();
#else
  si->nw = 0;
#endif

  si->s = search16_init(match_score,
                        mismatch_score,
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


  clusterinfo_t * clusterinfo = (clusterinfo_t *) xmalloc(seqcount * sizeof(clusterinfo_t));
  int clusters = 0;
  int lastlength = INT_MAX;
  
  progress_init("Clustering", seqcount);
  for (int i=0; i<seqcount; i++)
    {
      int seqno = i;
      int length = db_getsequencelen(seqno);

      if ((!opt_usersort) && (length > lastlength))
        fatal("Sequences not sorted by length and --usersort not specified.");

      lastlength = length;

      /* search selected db sequences, as indicated in bitmap */
      si->query_no = seqno;
      si->qsize = db_getabundance(seqno);
      si->query_head_len = db_getheaderlen(seqno);
      si->query_head = db_getheader(seqno);
      si->qseqlen = db_getsequencelen(seqno);
      si->qsequence = db_getsequence(seqno);

      if (opt_strand > 1)
        reverse_complement(si->rc, si->qsequence, si->qseqlen);

      search_onequery(si);
      
      if (si->hit_count)
        {
          int target = si->hits[0].target;
          clusterinfo[seqno].seqno = seqno;
          clusterinfo[seqno].clusterno = clusterinfo[target].clusterno;

          if (ucfilename)
            {
              fprintf(fp_uc, "H\t%d\t%ld\t%.1f\t%c\t0\t0\t%s\t%s\t%s\n",
                      clusterinfo[target].clusterno,
                      si->qseqlen,
                      si->hits[0].internal_id,
                      si->hits[0].strand ? '-' : '+',
                      si->hits[0].nwalignment,
                      si->query_head,
                      db_getheader(si->hits[0].target));
            }
          
          if (fp_alnout)
            results_show_alnout(fp_alnout,
                                si->hits, 1, si->query_head,
                                si->qsequence, si->qseqlen, si->rc);
          
          if (fp_fastapairs)
            results_show_fastapairs_one(fp_fastapairs,
                                        si->hits,
                                        si->query_head,
                                        si->qsequence,
                                        si->qseqlen,
                                        si->rc);
          
          if (fp_userout)
            results_show_userout_one(fp_userout, si->hits, si->query_head, 
                                     si->qsequence, si->qseqlen, si->rc);
          
          if (fp_blast6out)
            results_show_blast6out_one(fp_blast6out, si->hits, si->query_head,
                                       si->qsequence, si->qseqlen, si->rc);

          if (opt_matched)
            {
              fprintf(fp_matched, ">%s\n", si->query_head);
              fprint_fasta_seq_only(fp_matched, si->qsequence, si->qseqlen,
                                    opt_fasta_width);
            }
        }
      else
        {
          clusterinfo[seqno].seqno = seqno;
          clusterinfo[seqno].clusterno = clusters;
          dbindex_addsequence(seqno);
  
          if (ucfilename)
            {
              fprintf(fp_uc, "S\t%d\t%ld\t*\t*\t*\t*\t*\t%s\t*\n",
                      clusters, si->qseqlen, si->query_head);
            }

          if (fp_centroids)
            {
              fprintf(fp_centroids, ">%s\n", db_getheader(i));
              fprint_fasta_seq_only(fp_centroids, db_getsequence(i),
                                    db_getsequencelen(i),
                                    opt_fasta_width);
            }
          
          if (opt_output_no_hits)
            {
              if (fp_userout)
                results_show_userout_one(fp_userout, 0, si->query_head, 
                                         si->qsequence, si->qseqlen, si->rc);
              
              if (fp_blast6out)
                results_show_blast6out_one(fp_blast6out, 0, si->query_head,
                                           si->qsequence, si->qseqlen, si->rc);
            }

          if (opt_notmatched)
            {
              fprintf(fp_notmatched, ">%s\n", si->query_head);
              fprint_fasta_seq_only(fp_notmatched, si->qsequence, si->qseqlen,
                                    opt_fasta_width);
            }

          clusters++;
        }

      /* free memory for alignment strings */
      for(int j=0; j<si->hit_count; j++)
        free(si->hits[j].nwalignment);
      
      progress_update(i);
    }

  progress_done();

  qsort(clusterinfo, seqcount, sizeof(clusterinfo_t), compare_byclusterno);

  int size_min = INT_MAX;
  int size_max = 0;
  int lastcluster = 0;
  int size = 0;
  int singletons = 0;
  int centroid = 0;

  progress_init("Writing clusters", seqcount);

  FILE * fp_clusters = 0;
  char * fn_clusters = 0;
  if (opt_clusters)
    {
      fn_clusters = (char *) xmalloc(strlen(opt_clusters) + 25);
      sprintf(fn_clusters, "%s%d", opt_clusters, 0);
      fp_clusters = fopen(fn_clusters, "w");
      if (!fp_clusters)
        fatal("Unable to open clusters file for writing");
    };

  for(int i=0; i<seqcount; i++)
    {
      int clusterno = clusterinfo[i].clusterno;
      int seqno = clusterinfo[i].seqno;

      if (clusterno == lastcluster)
        {
          if (opt_clusters)
            {
              fprintf(fp_clusters, ">%s\n", db_getheader(seqno));
              fprint_fasta_seq_only(fp_clusters,
                                    db_getsequence(seqno),
                                    db_getsequencelen(seqno),
                                    opt_fasta_width);
            }
          size++;
        }
      else
        {
          if (size < size_min)
            size_min = size;
          if (size > size_max)
            size_max = size;
          if (size == 1)
            singletons++;

          if (opt_clusters)
            {
              fclose(fp_clusters);

              sprintf(fn_clusters, "%s%d", opt_clusters, clusterno);
              fp_clusters = fopen(fn_clusters, "w");
              if (!fp_clusters)
                fatal("Unable to open clusters file for writing");
              
              fprintf(fp_clusters, ">%s\n", db_getheader(seqno));
              fprint_fasta_seq_only(fp_clusters,
                                    db_getsequence(seqno),
                                    db_getsequencelen(seqno),
                                    opt_fasta_width);
            }              

          if (ucfilename)
            {
              fprintf(fp_uc, "C\t%d\t%d\t*\t*\t*\t*\t*\t%s\t*\n",
                      lastcluster, size, db_getheader(centroid));
              
            }
          
          centroid = clusterinfo[i].seqno;
          size = 1;
          lastcluster = clusterno;
        }
      progress_update(i);
    }
  if (size < size_min)
    size_min = size;
  if (size > size_max)
    size_max = size;
  if (size == 1)
    singletons++;

  if (opt_clusters)
    {
      fclose(fp_clusters);
      free(fn_clusters);
    }

  if (ucfilename)
    {
      fprintf(fp_uc, "C\t%d\t%d\t*\t*\t*\t*\t*\t%s\t*\n",
              lastcluster, size, db_getheader(centroid));
    }

  progress_done();


  fprintf(stderr, "Clusters: %d Size min %d, max %d, avg %.1f\n",
          clusters, size_min, size_max, 1.0 * seqcount / clusters);
  fprintf(stderr, "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
          singletons,
          100.0 * singletons / seqcount, 
          100.0 * singletons / clusters);

  /* thread specific clean up */
  search16_exit(si->s);
  unique_exit(si->uh);
  minheap_exit(si->m);

#ifdef COMPARENONVECTORIZED
  nw_exit(si->nw);
#endif

  if (si->targetlist)
    free(si->targetlist);
  if (si->hits)
    free(si->hits);
  if (si->kmers)
    free(si->kmers);
  if (si->rc)
    free(si->rc);

  if (opt_matched)
    fclose(fp_matched);
  if (opt_notmatched)
    fclose(fp_notmatched);
  if (opt_fastapairs)
    fclose(fp_fastapairs);
  if (fp_blast6out)
    fclose(fp_blast6out);
  if (fp_userout)
    fclose(fp_userout);
  if (fp_alnout)
    fclose(fp_alnout);

  if (fp_uc)
    fclose(fp_uc);
  if (fp_centroids)
    fclose(fp_centroids);
  
  dbindex_free();
  db_free();
  show_rusage();
}

void cluster_fast(char * cmdline, char * progheader)
{
  cluster(opt_cluster_fast, cmdline, progheader, 1);
}

void cluster_smallmem(char * cmdline, char * progheader)
{
  cluster(opt_cluster_smallmem, cmdline, progheader, 0);
}
