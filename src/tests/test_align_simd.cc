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
	char * query = "ACAT";
	search16_qprep(s16, query, 4);

	unsigned int sequences = 1;
	unsigned int seqnos[] = { 0 };
	CELL pscores[1];
	unsigned short paligned[1];
	unsigned short pmatches[1];
	unsigned short pmismatches[1];
	unsigned short pgaps[1];
	char * pcigar = (char *)xmalloc(5);

	search16(s16, sequences, seqnos, pscores, paligned, pmatches,
				pmismatches, pgaps, &pcigar);

	printf("cigar: %s\n", pcigar);

	}END_TEST

void add_align_simd_TC( Suite *s ) {
    TCase *tc_core = tcase_create( "align simd" );

    tcase_add_checked_fixture(tc_core, &setup, &teardown);

    tcase_add_test(tc_core, test_align_simd_simple);

    suite_add_tcase(s, tc_core);
}
