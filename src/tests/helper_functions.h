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

#ifndef VSEARCH_TEST_HELPER_H_
#define VSEARCH_TEST_HELPER_H_

#include <check.h>

#include "../align_simd.h"

#define CDEPTH 4
#define CHANNELS 8
#define SCORE_MATRIX_DIM 16

void check_cigar_matches(unsigned short pmatches, unsigned short pmismatches, char* pcigar);

void print_profile(CELL * dprofile);
void print_matrix(CELL matrix[SCORE_MATRIX_DIM*SCORE_MATRIX_DIM]);
void print_search_window(BYTE * dseq);

#endif /* VSEARCH_TEST_HELPER_H_ */

