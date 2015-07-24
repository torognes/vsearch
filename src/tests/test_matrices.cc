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

#include "../score_matrix.h"

// Checks the first two values for the first line of each matrix

static void test_first_two_values(ScoreMatrix * matrix, int val1, int val2)
{
  ck_assert_int_eq(val1, matrix->get64(1, 1));
  ck_assert_int_eq(val2, matrix->get64(1, 2));

  ck_assert_int_eq(val1, matrix->get16(1, 1));
  ck_assert_int_eq(val2, matrix->get16(1, 2));

  delete matrix;
}

START_TEST (test_matrices_buildin)
    {
      test_first_two_values(new ScoreMatrix(BLOSUM45), 5, -1);
      test_first_two_values(new ScoreMatrix(BLOSUM50), 5, -2);
      test_first_two_values(new ScoreMatrix(BLOSUM62), 4, -2);
      test_first_two_values(new ScoreMatrix(BLOSUM80), 5, -2);
      test_first_two_values(new ScoreMatrix(BLOSUM90), 5, -2);

      test_first_two_values(new ScoreMatrix(PAM30), 6, -3);
      test_first_two_values(new ScoreMatrix(PAM70), 5, -1);
      test_first_two_values(new ScoreMatrix(PAM250), 2, -0);
    }END_TEST

START_TEST (test_matrices_aa)
    {
      ScoreMatrix * matrix = new ScoreMatrix(BLOSUM62);

      ck_assert_int_eq(32, matrix->matrix_dim);
      ck_assert_int_eq(0, matrix->is_constant_scoring());

      delete matrix;

    }END_TEST

static void check_constant_score_matrix(ScoreMatrix * matrix, int match, int mismatch)
{
  for (int i = 0; i<matrix->matrix_dim; ++i)
    {
      for (int j = 0; j<matrix->matrix_dim; ++j)
        {
          if (i>0&&j>0)
            if (i==j)
              {
                ck_assert_int_eq(match, matrix->get16(i, j));
                ck_assert_int_eq(match, matrix->get64(i, j));
              }
            else
              {
                ck_assert_int_eq(mismatch, matrix->get16(i, j));
                ck_assert_int_eq(mismatch, matrix->get64(i, j));
              }
          else
            {
              ck_assert_int_eq(-1, matrix->get16(i, j));
              ck_assert_int_eq(-1, matrix->get64(i, j));
            }
        }
    }
}

START_TEST (test_matrices_constant_scoring)
    {
      // nucleotides
      ScoreMatrix * matrix = new ScoreMatrix(5, -4, MATRIX_MODE_NUC);

      ck_assert_int_eq(16, matrix->matrix_dim);
      ck_assert_int_eq(1, matrix->is_constant_scoring());

      check_constant_score_matrix(matrix, 5, -4);

      delete matrix;

      // amino acids
      matrix = new ScoreMatrix(6, -7, MATRIX_MODE_AA);

      ck_assert_int_eq(32, matrix->matrix_dim);
      ck_assert_int_eq(1, matrix->is_constant_scoring());

      check_constant_score_matrix(matrix, 6, -7);

      delete matrix;

    }END_TEST

void add_matrices_TC(Suite *s)
{
  TCase *tc_core = tcase_create("score matrices");

  tcase_add_test(tc_core, test_matrices_aa);
  tcase_add_test(tc_core, test_matrices_buildin);
  tcase_add_test(tc_core, test_matrices_constant_scoring);

  suite_add_tcase(s, tc_core);
}
