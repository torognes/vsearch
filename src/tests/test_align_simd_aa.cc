/*
 Copyright (C) 2015 Jakob Frielingsdorf

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

 Contact: Jakob Frielingsdorf <jfrielingsdorf@gmail.com>
 */

#include "tests.h"

#include "helper_functions.h"

#include "../vsearch.h"
#include "../align_simd_aa.h"
#include "../util.h"

static struct s16info_s * s16;

/* Is run once before each unit test */
static void setup()
{
  opt_maxseqlength = 2000;

  db_read("../data/uniprot_first_two_sequences.fasta", 0, DB_MODE_AA); // TODO what does upcase do?

  CELL match = 5;
  opt_match = match;
  CELL mismatch = -4;
  opt_mismatch = mismatch;
  CELL gap_open = 2;
  CELL gap_extension = 1;

  s16 = search16_aa_init(match, mismatch,
      gap_open, gap_open, gap_open, gap_open, gap_open, gap_open,
      gap_extension, gap_extension, gap_extension, gap_extension, gap_extension, gap_extension);
}

/* Is run once after each unit test */
static void teardown()
{
  search16_aa_exit(s16);

  db_free();
}

START_TEST (test_align_simd_simple)
    {
      char * query = (char *)"MSIIGATRLQNDKRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKTSQPLCAR";
      search16_aa_qprep(s16, query, strlen(query));

      unsigned int seq_count = 1;
      unsigned int seqnos[] = { 0 };
      CELL pscores[1];
      unsigned short paligned[1];
      unsigned short pmatches[1];
      unsigned short pmismatches[1];
      unsigned short pgaps[1];
      char * pcigar[1];
      pcigar[0] = 0;
      search16_aa(s16, seq_count, seqnos, pscores, paligned, pmatches, pmismatches, pgaps, pcigar);

      ck_assert_ptr_ne(0, pcigar[0]);
      ck_assert_str_eq("MD4MI3MI3MD37MI3M3I2M2DM", pcigar[0]);

      check_cigar_matches(pmatches[0], pmismatches[0], pcigar[0]);
    }END_TEST

START_TEST (test_align_simd_all)
    {
      char * query = (char *)"MSIIGATRLQNDKRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKTSQPLCAR";

      search16_aa_qprep(s16, query, strlen(query));

      unsigned long seq_count = 2;

      unsigned int seqnos[seq_count];
      for (unsigned int i = 0; i<seq_count; ++i)
        {
          seqnos[i] = i;
        }

      CELL pscores[seq_count];
      unsigned short paligned[seq_count];
      unsigned short pmatches[seq_count];
      unsigned short pmismatches[seq_count];
      unsigned short pgaps[seq_count];
      char * pcigar[seq_count];
      memset(pcigar, 0, seq_count);

      search16_aa(s16, seq_count, seqnos, pscores, paligned, pmatches, pmismatches, pgaps, pcigar);

      for (unsigned int i = 0; i<seq_count; ++i)
        {
          ck_assert_ptr_ne(0, pcigar[i]);
          check_cigar_matches(pmatches[i], pmismatches[i], pcigar[i]);
        }
    }END_TEST

void add_align_simd_aa_TC(Suite *s)
{
  TCase *tc_core = tcase_create("align simd for proteins");

  tcase_add_checked_fixture(tc_core, &setup, &teardown);

  tcase_add_test(tc_core, test_align_simd_simple);
  tcase_add_test(tc_core, test_align_simd_all);

  suite_add_tcase(s, tc_core);
}
