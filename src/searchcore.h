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

//#define COMPARENONVECTORIZED

/* the number of alignments that can be delayed */
#define MAXDELAYED 8

/* Default minimum number of word matches for word lengths 3-15 */
const int minwordmatches_defaults[] =
  { -1, -1, -1, 18, 17, 16, 15, 14, 12, 11, 10,  9,  8,  7,  5,  3 };

struct hit
{
  int target;
  int strand;

  /* candidate info */
  unsigned int count;     /* number of unique kmers shared with query */

  bool accepted;          /* is it accepted? */
  bool rejected;          /* is it rejected? */
  bool aligned;           /* has this hit been aligned */
  bool weak;              /* weak hits are aligned with id > weak_id */

  /* info about global alignment, including terminal gaps */

  int nwscore;           /* alignment score */
  int nwdiff;            /* indels and mismatches in global alignment */
  int nwgaps;            /* gaps in global alignment */
  int nwindels;          /* indels in global alignment */
  int nwalignmentlength; /* length of global alignment */
  double nwid;           /* percent identity of global alignment */
  char * nwalignment;    /* alignment string (cigar) of global alignment */
  int matches;
  int mismatches;
  
  /* info about alignment excluding terminal gaps */

  int internal_alignmentlength;
  int internal_gaps;
  int internal_indels;
  int trim_q_left;
  int trim_q_right;
  int trim_t_left;
  int trim_t_right;
  int trim_aln_left;
  int trim_aln_right;

  /* more info */

  double id;             /* identity used for ranking */
  double id0, id1, id2, id3, id4;

  int shortest;          /* length of shortest of query and target */
  int longest;           /* length of longest of query and target */
};

/* type of kmer hit counter element remember possibility of overflow */
typedef unsigned short count_t;

struct searchinfo_s
{
  int query_no;                 /* query number, zero-based */
  int strand;                   /* strand of query being analysed */
  int qsize;                    /* query abundance */
  int query_head_len;           /* query header length */
  int query_head_alloc;         /* bytes allocated for the header */
  char * query_head;            /* query header */
  int qseqlen;                  /* query length */
  int seq_alloc;                /* bytes allocated for the query sequence */
  char * qsequence;             /* query sequence */
  unsigned int kmersamplecount; /* number of kmer samples from query */
  unsigned int * kmersample;    /* list of kmers sampled from query */
  count_t * kmers;              /* list of kmer counts for each db seq */
  struct hit * hits;            /* list of hits */
  int hit_count;                /* number of hits in the above list */
  struct uhandle_s * uh;        /* unique kmer finder instance */
  struct s16info_s * s;         /* SIMD aligner instance */
  struct nwinfo_s * nw;         /* NW aligner instance */
  LinearMemoryAligner * lma;    /* Linear memory aligner instance pointer */
  int accepts;                  /* number of accepts */
  int rejects;                  /* number of rejects */
  minheap_t * m;                /* min heap with the top kmer db seqs */
  int finalized;
};

void search_onequery(struct searchinfo_s * si, int seqmask);

struct hit * search_findbest2_byid(struct searchinfo_s * si_p,
                                   struct searchinfo_s * si_m);

struct hit * search_findbest2_bysize(struct searchinfo_s * si_p,
                                     struct searchinfo_s * si_m);

int search_acceptable_unaligned(struct searchinfo_s * si, int target);

int search_acceptable_aligned(struct searchinfo_s * si,
                              struct hit * hit);

void align_trim(struct hit * hit);

void search_joinhits(struct searchinfo_s * si_p,
                     struct searchinfo_s * si_m,
                     struct hit * * hits,
                     int * hit_count);

bool search_enough_kmers(struct searchinfo_s * si,
                         unsigned int count);
