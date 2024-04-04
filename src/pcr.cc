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

const double prob_duplication = 0.95;
const double prob_chimera = 0.05;
const double prob_base_subst = 0.0004;
const int min_nwscore = 50;
const long cycles = 25;
const int border_left = 10;
const int border_right = 10;

void create_chimera(char * seq1,
		    char * seq2,
		    char * nwalignment,
		    int64_t nwmatches,
		    int64_t nwmismatches,
		    char ** chimera,
		    int64_t * chimera_length)
{
  printf("Cigar: %s\n", nwalignment);
  
  int breakpos1 = 0;
  int breakpos2 = 0;
  
  if (nwmatches + nwmismatches >= border_left + border_right)
    {
      int64_t breakpoint = border_left +
	random_int(nwmatches + nwmismatches + 1 - border_left - border_right);
      printf("Breakpoint: %lld\n", breakpoint);

      char * p = nwalignment;
      char op;
      int64_t run;
      int pos1 = 0;
      int pos2 = 0;
      int m = 0;

      while (*p)
	{
	  run = 1;
	  int scanlength = 0;
	  sscanf(p, "%" PRId64 "%n", &run, &scanlength);
	  op = *(p+scanlength);
	  switch(op)
	    {
	    case 'M':
	      for (int z = 0; z < run; z++)
		{
		  pos1++;
		  pos2++;
		  m++;
		  if (m >= breakpoint)
		    {
		      breakpos1 = pos1;
		      breakpos2 = pos2;
		      uint64_t chim_len = pos1 + strlen(seq2) - pos2;
		      char * chim = (char*) xmalloc(chim_len + 1);
		      strncpy(chim, seq1, pos1);
		      strcpy(chim + pos1, seq2 + pos2);
		      printf("Breakpoint: %d %d\n", breakpos1, breakpos2);
		      printf("Seq1: %s\n", seq1);
		      printf("Seq2: %s\n", seq2);
		      printf("SeqC: %s\n", chim);
		      return;
		    }
		}
	      break;
	    case 'D':
	      pos1 += run;
	      break;
	    case 'I':
	      pos2 += run;
	      break;
	    }
	  p += scanlength + 1;
	}
    }
  else
    {
      // too short
    }
}



void pcr()
{
  (void) prob_duplication;
  (void) prob_chimera;
  (void) prob_base_subst;

/*
    Simulate PCR with chimera formation and single nucleotide
    substitutions.
  */

  /* Input:
     - File with input sequences, FASTA. Header with sequence identifiers.
     - Abundance distribution type (log)
     - Duplication prob, e.g. 0.95
     - chimera formation prob, e.g. 0.05 ???
     - Base error frequency, e.g. 0.00004
     - Cycles, e.g. 25

     Output:
     - File with output sequences, FASTA. Header with sequence identifiers, history (parents, breakpoints), abundance.
  */

  /* pseudocode:

     read database from given input file
     apply abundance distribution, e.g. log, with abundance 10+
     keep track of total abundance
     for each cycle (1..25):
     for each sequence A in the database:
     for each copy (abundance) of the sequence:
     if random < chimera_formation_prob:
     pick a random other sequence B from database, scaled by abundance
     align A and B
     if A and B are sufficiently similar:
     choose breakpoint at a random position within aligned region
     create chimeric sequence from A and B at breakpoint

     if random < duplication_prob:
     make a duplicate sequence C
     for each base:
     if random < base_error_freq:
     substitute base randomly in C
     flag mutation
     add C to database (increase abundance if identical)
  */

  if (!opt_output)
    fatal("Output file for PCR simulation must be specified with --output");

  FILE * fp_output = fopen_output(opt_output);
  if (!fp_output)
    {
      fatal("Unable to open PCR simulation output file for writing");
    }

  db_read(opt_pcr, 0);

  long dbsequencecount = db_getsequencecount();

  /* apply abundance profile */
  
  uint64_t * abundance = (uint64_t *) xmalloc(dbsequencecount * sizeof(uint64_t));
  uint64_t total_abundance = 0;
  (void) total_abundance;
  for(long i = 0; i < dbsequencecount; i++)
    {
      /* Should be log normal distributed */
      abundance[i] = 10;
      total_abundance += abundance[i];
    }

  struct nwinfo_s * nwi = nw_init();


  LinearMemoryAligner * lma = new LinearMemoryAligner;

  int64_t * scorematrix = lma->scorematrix_create(opt_match, opt_mismatch);

  lma->set_parameters(scorematrix,
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
  
  
  progress_init("Simulating PCR", cycles);
  for (long cycle = 1; cycle <= cycles; cycle++)
    {
      for(long i = 0; i < dbsequencecount; i++)
	{
	  uint64_t a = abundance[i];
	  for(uint64_t x = 0; x < a ; x++)
	    {
	      long r = random_int(1000);
	      if (r < int(1000 * prob_chimera))
		{
		  // prob should be scaled by abundance
		  long j = random_int(dbsequencecount);
		  if (i != j)
		    {
		      printf("Aligning %ld and %ld!\n", i, j);
		      
		      char * seq1 = db_getsequence(i);
		      long seq1len = db_getsequencelen(i);
		      char * seq2 = db_getsequence(j);
		      long seq2len = db_getsequencelen(j);
		      
		      char * nwcigar = lma->align(seq1,
						  seq2,
						  seq1len,
						  seq2len);
		      
		      int64_t nwscore;
		      int64_t nwalignmentlength;
		      int64_t nwmatches, nwmismatches;
		      int64_t nwgaps;
		      
		      lma->alignstats(nwcigar,
                                      seq1,
                                      seq2,
                                      & nwscore,
                                      & nwalignmentlength,
                                      & nwmatches,
                                      & nwmismatches,
                                      & nwgaps);

		      printf("Score: %lld Matches: %lld Mismatches: %lld Length: %lld gaps: %lld\n", nwscore, nwmatches, nwmismatches, nwalignmentlength, nwgaps);

		      if (nwscore > min_nwscore)
			{
			  char * chimera;
			  int64_t chimera_length;
			  create_chimera(seq1,
					 seq2,
					 nwcigar,
					 nwmatches,
					 nwmismatches,
					 & chimera,
					 & chimera_length);
			  // xfree(chimera);
			}
		    }
		}
	    }
	}
      progress_update(cycle);
    }
  progress_done();

  nw_exit(nwi);

  progress_init("Writing output", dbsequencecount);
  for(long i = 0; i < dbsequencecount; i++)
    {
      fasta_print_db_relabel(fp_output, i, i+1);
      progress_update(i);
    }
  progress_done();

  xfree(abundance);
  db_free();
  fclose(fp_output);
}


/*

  make ; ../bin/vsearch --pcr pcr/pcr.fasta --output x.fasta --sizeout --fasta_width 0 --relabel seq ; cat x.fasta | cut -c1-79

*/
