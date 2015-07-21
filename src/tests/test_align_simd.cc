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

#include "../vsearch.h"
#include "../align_simd.h"
#include "../util.h"

static struct s16info_s * s16;

/* Is run once before each unit test */
static void setup() {
	opt_maxseqlength = 5000;

	db_read("../data/AF091148.fsa", 0); // TODO what does upcase do?

	CELL match = 5;
	CELL mismatch = -4;
	CELL gap_open = 2;
	CELL gap_extension = 3;

	s16 = search16_init(match, mismatch,
			gap_open, gap_open, gap_open, gap_open, gap_open, gap_open,
			gap_extension, gap_extension, gap_extension, gap_extension, gap_extension, gap_extension);
}

/* Is run once after each unit test */
static void teardown() {
	search16_exit(s16);

	db_free();
}

START_TEST (test_align_simd_simple)
    {
      char * query = (char *)"ACAT";
      search16_qprep(s16, query, 4);

      unsigned int seq_count = 1;
      unsigned int seqnos[] = { 0 };
      CELL pscores[1];
      unsigned short paligned[1];
      unsigned short pmatches[1];
      unsigned short pmismatches[1];
      unsigned short pgaps[1];
      char * pcigar[1];
      pcigar[0] = 0;

      search16(s16, seq_count, seqnos, pscores, paligned, pmatches,
          pmismatches, pgaps, pcigar);

      ck_assert_ptr_ne(0, pcigar[0]);

    }END_TEST

START_TEST (test_align_simd_all)
    {
      char * query = (char *)"ACAT";
      search16_qprep(s16, query, 4);

      unsigned long seq_count = db_getsequencecount();

      unsigned int seqnos[seq_count];
      for (unsigned int i = 0; i < seq_count; ++i)
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

      search16(s16, seq_count, seqnos, pscores, paligned, pmatches,
          pmismatches, pgaps, pcigar);

      for (unsigned int i = 0; i<seq_count; ++i)
        {
          ck_assert_ptr_ne(0, pcigar[i]);
        }

    }END_TEST

void add_align_simd_TC( Suite *s ) {
    TCase *tc_core = tcase_create( "align simd" );

    tcase_add_checked_fixture(tc_core, &setup, &teardown);

    tcase_add_test(tc_core, test_align_simd_simple);
    tcase_add_test(tc_core, test_align_simd_all);

    suite_add_tcase(s, tc_core);
}
