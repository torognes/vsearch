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

Suite* vsearch_suite(void)
{
  Suite *s = suite_create("vsearch");

  add_maps_TC(s);
  add_matrices_TC(s);
//  add_db_TC(s);
  add_dprofile_fill_nuc_TC(s);
  add_align_simd_nuc_TC(s);
  add_align_simd_aa_TC(s);

  return s;
}

int main(void)
{
  printf("Using Check unit testing framework version %d.%d.%d\n", CHECK_MAJOR_VERSION, CHECK_MINOR_VERSION,
  CHECK_MICRO_VERSION);

  int number_failed;
  Suite *s = vsearch_suite();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed==0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
