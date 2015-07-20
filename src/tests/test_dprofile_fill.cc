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

#define CDEPTH 4
#define CHANNELS 8

static struct s16info_s * s16;

extern void dprofile_fill16(CELL * dprofile_word,
        CELL * score_matrix_word,
        BYTE * dseq);

/* Is run once before each unit test */
static void setup() {
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
}

START_TEST (test_dprofile_fill_nucleotide_simple)
	{

	BYTE * d_begin[CHANNELS];
	BYTE * d_end[CHANNELS];

	for (int c=0; c<CHANNELS; c++)
	{
	  d_begin[c] = 0;
	  d_end[c] = d_begin[c];
	}

	char * dseq_orig = (char*) xmalloc(4);
	strcpy(dseq_orig, "ACAT");

	d_begin[0] = (BYTE*) dseq_orig;

	__m128i dseqalloc[CDEPTH];
	BYTE * dseq = (BYTE*) &dseqalloc;

	for (int c = 0; c < CHANNELS; c++) {
		for (int j = 0; j < CDEPTH; j++) {
			if (d_begin[c] < d_end[c])
				dseq[CHANNELS * j + c] = chrmap_4bit[*(d_begin[c]++)];
			else
				dseq[CHANNELS * j + c] = 0;
		}
	}

	__m128i matrix[32];

	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++) {
			CELL value;
			if (i == j)
				value = opt_match;
			else if ((i == 0) || (j == 0) || (i > 4) || (j > 4))
				value = 0;
			else
				value = opt_mismatch;
			((CELL*) (matrix))[16 * i + j] = value;
		}

	__m128i * dprofile = (__m128i *) xmalloc(2*4*8*16);

	dprofile_fill16((CELL*) dprofile, (CELL*) matrix, dseq);

	}END_TEST

void add_dprofile_fill_TC( Suite *s ) {
    TCase *tc_core = tcase_create( "dprofile fill" );

    tcase_add_checked_fixture(tc_core, &setup, &teardown);

    tcase_add_test(tc_core, test_dprofile_fill_nucleotide_simple);

    suite_add_tcase(s, tc_core);
}
