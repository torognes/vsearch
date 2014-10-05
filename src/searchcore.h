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

struct hit
{
  long target;
  int strand;
  long count; /* number of word matches */

  long matches;
  long mismatches;

  /* info for global alignment, entire sequences */

  long nwscore;
  long nwdiff; /* indels and mismatches in global alignment */
  long nwgaps; /* gaps in global alignment */
  long nwindels; /* indels in global alignment */
  long nwalignmentlength; /* length of global alignment */
  double nwid; /* percent identity of global alignment */
  char * nwalignment; /* alignment string (cigar) of global alignment */
  
  /* info for semi-global alignment, excluding gaps at ends */

  long internal_alignmentlength;
  long internal_gaps;
  long internal_indels;
  double internal_id;

  long trim_q_left;
  long trim_q_right;
  long trim_t_left;
  long trim_t_right;
  long trim_aln_left;
  long trim_aln_right;
};

/* type of kmer hit counter element remember possibility of overflow */
typedef unsigned short count_t;

struct searchinfo_s
{
  pthread_t pthread;            /* pthread info */
  long query_no;                /* query number, zero-based */
  long qsize;                   /* query abundance */
  long query_head_len;          /* query header length */
  long query_head_alloc;        /* bytes allocated for the header */
  char * query_head;            /* query header */
  long qseqlen;                 /* query length */
  long seq_alloc;               /* bytes allocated for the query sequence */
  long rc_seq_alloc;            /* bytes allocated for reverse complement q */
  char * qsequence;             /* query sequence */
  char * rc;                    /* query sequence, reverse complement */
  unsigned int kmersamplecount; /* number of kmer samples from query */
  unsigned int * kmersample;    /* list of kmers sampled from query */
  unsigned int * targetlist;    /* list of db seqs with >0 kmer match */
  count_t * kmers;              /* list of kmer counts for each db seq */
  struct hit * hits;            /* list of hits */
  int hit_count;                /* number of hits in the above list */
  struct uhandle_s * uh;        /* unique kmer finder instance */
  struct s16info_s * s;         /* SIMD aligner instance */
  struct nwinfo_s * nw;         /* NW aligner instance */
  int strand;                   /* strand of query being analysed */
  int accepts;                  /* number of accepts */
  int rejects;                  /* number of rejects */
  int candidatecount;           /* number of candidate seqs for SIMD align */
  int candidate_target[8];      /* sequence nos of top candidates */
  int candidate_kmercount[8];   /* kmer count in top candidates */
  minheap_t * m;                /* min heap with the top kmer db seqs */
};

void search_onequery(struct searchinfo_s * si);
