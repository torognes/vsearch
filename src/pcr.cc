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

/*
  Simulate PCR with chimera formation and single nucleotide
  substitutions.

  Relevant paper:
  Potapov V, Ong JL (2017)
  Examining Sources of Error in PCR by Single-Molecule Sequencing
  PLOS ONE 12(1): e0169774.
  https://doi.org/10.1371/journal.pone.0169774

  Command:
  --pcr_sim input.fasta

  Required option:
  --output output.fasta

  Options:
  --pcr_cycles 20
  --pcr_chimera_p 0.01
  --pcr_subst_p 0.00015

  Instructions:
  First create the start.fasta file with the intital sequences, one
  for each copy. The "--rereplicate" command may be used to create
  multiple copies to obtain a skewed distribution. Then run "pcr_sim"
  as above. After the simulation, run "derep_id" to dereplicate normal
  and chimera sequences separately. Run "derep_full" to dereplicate all.

  Usage:
    vsearch --pcr_sim start.fasta --output mix.fasta
    vsearch --derep_id mix.fasta --output mix.derep.fasta --sizeout
    vsearch --derep_full mix.fasta --output all.derep.fasta --sizeout

  Input:
  - File with input sequences, FASTA. Headers are ignored, also abundances.

  Output:
  - File with output sequences, FASTA. Header with name normal or chimera.

  Pseudocode:
  read database from given input file
  for each cycle (1..25):
    for each sequence A in the database:
      if random < chimera_formation_prob:
        pick another random sequence B from the database
        align A and B
        if A and B are sufficiently similar:
          choose random breakpoint within aligned region with 10 bp border
          create chimeric sequence C from A and B at breakpoint
      else:
        make a duplicate sequence C from A
      for each base in C:
        if random < base_error_freq:
          substitute base randomly in C with another base
      add C to database

*/

const int big_int = 1000000000;
const int min_nwscore = 1;
const int border_left = 10;
const int border_right = 10;
const char * header_normal  = "normal";
const char * header_chimera = "chimera";

void mutate_sequence(char * seq,
                     int64_t len)
{
  /* apply random substitutions at random positions within sequence */
  /* use a uniform distribution of alternative bases (not correct) */
  for (int i = 0; i < len; i++)
    {
      if (random_int(big_int) < big_int * opt_pcr_subst_p)
        {
          int64_t x = seq[i];
          int64_t b = random_int(4);
          while (b == x)
            b = random_int(4);
          seq[i] = sym_nt_2bit[b];
        }
    }
}

void create_chimera(char * seq1,
                    char * seq2,
                    char * nwalignment,
                    int64_t nwmatches,
                    int64_t nwmismatches,
                    char ** chimera,
                    int64_t * chimera_length)
{
  if (nwmatches + nwmismatches >= border_left + border_right)
    {
      char * p = nwalignment;
      char op;
      int64_t run;
      int pos1 = 0;
      int pos2 = 0;
      int m = 0;
      int64_t breakpoint = border_left +
        random_int(nwmatches + nwmismatches + 1 - border_left - border_right);
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
                      uint64_t chim_len = pos1 + strlen(seq2) - pos2;
                      char * chim = (char*) xmalloc(chim_len + 1);
                      strncpy(chim, seq1, pos1);
                      strncpy(chim + pos1, seq2 + pos2, chim_len - pos1);
                      chim[chim_len] = 0;
                      * chimera_length = chim_len;
                      * chimera = chim;
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
  * chimera_length = 0;
  * chimera = nullptr;
}


void pcr()
{
  /* ignore abundance */
  opt_sizein = false;

  if ((opt_pcr_cycles < 0) || (opt_pcr_cycles > 100))
    fatal("The PCR cycles option argument must be between 0 and 100\n");

  if ((opt_pcr_chimera_p < 0.0) || (opt_pcr_chimera_p > 1.0))
    fatal("The PCR chimera formation probability must be between 0.0 and 1.0\n");

  if ((opt_pcr_subst_p < 0.0) || (opt_pcr_subst_p > 1.0))
    fatal("The PCR base substitution probability must be between 0.0 and 1.0\n");

  if (!opt_output)
    fatal("Output file for PCR simulation must be specified with --output");

  FILE * fp_output = fopen_output(opt_output);
  if (!fp_output)
    fatal("Unable to open PCR simulation output file for writing");

  db_read(opt_pcr_sim, 0);

  long dbsequencecount = db_getsequencecount();

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

  if (opt_log)
    fprintf(fp_log,
            "PCR with %" PRId64 " cycles, chimera prob. %3.10lg, substitution prob. %3.10lg\n",
            opt_pcr_cycles,
            opt_pcr_chimera_p,
            opt_pcr_subst_p);
  if (! opt_quiet)
    fprintf(stderr,
            "PCR with %" PRId64 " cycles, chimera prob. %3.10lg, substitution prob. %3.10lg\n",
            opt_pcr_cycles,
            opt_pcr_chimera_p,
            opt_pcr_subst_p);

  progress_init("Simulating PCR", opt_pcr_cycles);
  for (long cycle = 1; cycle <= opt_pcr_cycles; cycle++)
    {
      long count = dbsequencecount; /* ignore new sequences this cycle */

      for(long i = 0; i < count; i++)
        {
          char * seq1 = db_getsequence(i);
          long seq1len = db_getsequencelen(i);
          if (random_int(big_int) < int(big_int * opt_pcr_chimera_p))
            {
              long j = random_int(count);
              if (i != j)
                {
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

                      if (chimera_length)
                        {
                          mutate_sequence(chimera, chimera_length);
                          db_add(false,
                                 (char*) header_chimera,
                                 chimera,
                                 nullptr,
                                 strlen(header_chimera),
                                 chimera_length,
                                 1);
                          dbsequencecount++;
                          xfree(chimera);
                        }
                    }
                }
            }
          else /* duplication */
            {
              char * dup = strdup(seq1);
              long dup_length = seq1len;
              bool is_chim = strcmp(db_getheader(i), header_chimera) == 0;
              mutate_sequence(dup, dup_length);
              db_add(false,
                     (char*)(is_chim ? header_chimera : header_normal),
                     dup,
                     nullptr,
                     strlen(is_chim ? header_chimera : header_normal),
                     dup_length,
                     1);
              dbsequencecount++;
              xfree(dup);
            }
        }
      progress_update(cycle);
    }
  progress_done();

  uint64_t chimeric = 0;
  uint64_t non_chimeric = 0;

  progress_init("Writing output", dbsequencecount);
  for(long i = 0; i < dbsequencecount; i++)
    {
      char * header = (char*) header_normal;
      if (strcmp(db_getheader(i), header_chimera) == 0)
        {
          header = (char*) header_chimera;
          chimeric++;
        }
      else
        {
          non_chimeric++;
        }
      int headerlen = strlen(header);
      fasta_print_general(fp_output,
                          nullptr,
                          db_getsequence(i),
                          db_getsequencelen(i),
                          header,
                          headerlen,
                          0,
                          -1,
                          -1.0,
                          -1,
                          -1,
                          nullptr,
                          0.0);
      progress_update(i);
    }
  progress_done();

  if (opt_log)
    fprintf(fp_log,
            "Written %" PRIu64 " chimeric and %" PRIu64 " non-chimeric sequences\n",
            chimeric,
            non_chimeric);
  if (! opt_quiet)
    fprintf(stderr,
            "Written %" PRIu64 " chimeric and %" PRIu64 " non-chimeric sequences\n",
            chimeric,
            non_chimeric);

  db_free();
  fclose(fp_output);
}
