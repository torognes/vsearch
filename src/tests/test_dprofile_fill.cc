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

#include <string.h>

#include "../vsearch.h"
#include "../align_simd.h"
#include "../util.h"

#define CDEPTH 4
#define CHANNELS 8
#define SCORE_MATRIX_DIM 16

extern void dprofile_fill16(CELL * dprofile_word,
                            CELL * score_matrix_word,
                            BYTE * dseq);

/* Is run once before each unit test */
static void setup()
{
  opt_match = 5;
  opt_mismatch = -4;
}

/* Is run once after each unit test */
static void teardown()
{
}

void print_profile(CELL * dprofile)
  {
    for (int i = 0; i<SCORE_MATRIX_DIM; ++i)
      {
        for (int j = 0; j<CDEPTH; ++j)
          {
            for (int k = 0; k<CHANNELS; k++)
              {
                printf("%2d ", dprofile[CHANNELS*CDEPTH*i+CHANNELS*j+k]);
              }
            printf(" | ");
          }
        printf("\n");
      }
  }

void print_matrix(CELL matrix[SCORE_MATRIX_DIM*SCORE_MATRIX_DIM])
{
  // end copy
  for (int i = 0; i<SCORE_MATRIX_DIM; i++)
    {
      for (int j = 0; j<SCORE_MATRIX_DIM; j++)
        {
          printf("%2d, ", matrix[SCORE_MATRIX_DIM*i+j]);
        }
      printf("\n");
    }
  printf("\n");
}

void print_search_window(BYTE * dseq)
{
  for (int i = 0; i<CDEPTH; i++)
    {
      for (int j = 0; j<CHANNELS; j++)
        {
          printf("%2d ", dseq[CHANNELS*i+j]);
        }
      printf("\n");
    }
}

void check_profile(CELL matrix[SCORE_MATRIX_DIM*SCORE_MATRIX_DIM], CELL* dprofile, BYTE* dseq)
{
  for (int i = 0; i<SCORE_MATRIX_DIM; ++i)
    {
      for (int j = 0; j<CDEPTH; ++j)
        {
          for (int k = 0; k<CHANNELS; k++)
            {
              CELL pval = dprofile[CHANNELS*CDEPTH*i+CHANNELS*j+k];
              CELL dval = dseq[CHANNELS*j+k];

              ck_assert_int_eq(matrix[SCORE_MATRIX_DIM*dval+i], pval);
            }
        }
    }
}

/* copied from align_simd.c -> search16_init() */
void fill_matrix(CELL matrix[SCORE_MATRIX_DIM*SCORE_MATRIX_DIM])
{
  for (int i = 0; i<SCORE_MATRIX_DIM; i++)
    for (int j = 0; j<SCORE_MATRIX_DIM; j++)
      {
        CELL value;
        if (i==j)
          value = opt_match;
        else if ((i==0)||(j==0)||(i>4)||(j>4))
          value = 0;
        else
          value = opt_mismatch;
        matrix[SCORE_MATRIX_DIM*i+j] = value;
      }
}

START_TEST (test_dprofile_fill_nucleotide_simple)
    {
      BYTE * dseq = (BYTE *)calloc(1, CDEPTH * CHANNELS);

      int k = 1;

      int dseq_orig[][4] = {{ 'A', 'C', 'A', 'T' }, }; // ACAT
      for (int j = 0; j < k; ++j)
          for (int i = 0; i < CDEPTH; ++i) {
              dseq[CHANNELS*i+k] = chrmap_4bit[dseq_orig[j][i]];
          }
//      print_search_window(dseq);

      CELL matrix[SCORE_MATRIX_DIM*SCORE_MATRIX_DIM];

      fill_matrix(matrix);
//      print_matrix(matrix);

      CELL * dprofile = (CELL*)xmalloc(2*4*8*16);

      dprofile_fill16(dprofile, matrix, dseq);
//      print_profile(dprofile);

      check_profile(matrix, dprofile, dseq);
    }END_TEST

void add_dprofile_fill_TC(Suite *s)
{
  TCase *tc_core = tcase_create("dprofile fill");

  tcase_add_checked_fixture(tc_core, &setup, &teardown);

  tcase_add_test(tc_core, test_dprofile_fill_nucleotide_simple);

  suite_add_tcase(s, tc_core);
}
