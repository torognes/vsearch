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

START_TEST (test_align_simd_simple)
	{
		fail( "not yet implemented" );
	}END_TEST

void add_align_simd_TC( Suite *s ) {
    TCase *tc_core = tcase_create( "align simd" );
    tcase_add_test( tc_core, test_align_simd_simple );

    suite_add_tcase( s, tc_core );
}
