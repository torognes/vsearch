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
#include <sys/time.h>
#include <time.h>

#include "../vsearch.h"
#include "../align_simd.h"
#include "../util.h"

#define CDEPTH 4
#define CHANNELS 8
#define SCORE_MATRIX_DIM 16

extern void dprofile_fill16(CELL * dprofile_word,
                            CELL * score_matrix_word,
                            BYTE * dseq);

static BYTE dseq[CDEPTH * CHANNELS];
static CELL matrix[SCORE_MATRIX_DIM*SCORE_MATRIX_DIM];

/* Is run once before each unit test */
static void setup()
{
  opt_match = 5;
  opt_mismatch = -4;

  memset(matrix, 0, SCORE_MATRIX_DIM*SCORE_MATRIX_DIM);
  memset(dseq, 0, CDEPTH * CHANNELS);
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
  //      print_profile(dprofile);
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

/* copied from search16_init() in align_simd.c */
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
  //      print_matrix(matrix);
}

void fill_search_window(int dseq_count, int db_sequences[][CDEPTH])
{
  for (int j = 0; j<dseq_count; ++j)
    for (int i = 0; i<CDEPTH; ++i)
      {
        dseq[CHANNELS*i+j] = chrmap_4bit[db_sequences[j][i]];
      }
//        print_search_window(dseq);
}

START_TEST (test_dprofile_fill_nucleotide_simple)
    {
      int dseq_count = 1;

      int db_sequences[][CDEPTH] = {{ 'A', 'C', 'A', 'T' }, };
      fill_search_window(dseq_count, db_sequences);

      fill_matrix(matrix);

      CELL * dprofile = (CELL*)xmalloc(2*4*8*16);

      dprofile_fill16(dprofile, matrix, dseq);

      check_profile(matrix, dprofile, dseq);
    }END_TEST

START_TEST (test_dprofile_fill_nucleotide_more)
    {
      int dseq_count = CHANNELS;

      int db_sequences[][CDEPTH] = {
          { 'A', 'C', 'A', 'T' },
          { 'A', 'T', 'C', 'C' },
          { 'T', 'T', 'T', 'T' },
          { 'A', 'A', 0, 0 },
          { 'C', 'T', 'C', 'C' },
          { 'A', 'C', 'T', 'C' },
          { 'A', 'T', 'A', 0 },
          { 'C', 'A', 'C', 'C' }};
      fill_search_window(dseq_count, db_sequences);

      fill_matrix(matrix);

      CELL * dprofile = (CELL*)xmalloc(2*4*8*16);

      dprofile_fill16(dprofile, matrix, dseq);

      check_profile(matrix, dprofile, dseq);
    }END_TEST

START_TEST (test_dprofile_fill_nucleotide_perf)
    {
      int dseq_count = 1;

      int db_sequences[][CDEPTH] = { { 'A', 'C', 'A', 'T' }, };
      fill_search_window(dseq_count, db_sequences);

      fill_matrix(matrix);

      CELL * dprofile = (CELL*)xmalloc(2*4*8*16);

      struct timeval start;
      struct timeval finish;

      gettimeofday(&start, NULL);

      long rounds = 10000000;
      for (int i = 0; i<rounds; ++i)
        {
          dprofile_fill16(dprofile, matrix, dseq);
        }

      gettimeofday(&finish, NULL);

      double elapsed = (finish.tv_sec-start.tv_sec);
      elapsed += (finish.tv_usec-start.tv_usec)/1000000.0;

      // TODO before: ca. 0.703s
      printf("\nRuntime used for %ld runs of dprofile_fill16: %lf sec\n\n", rounds, elapsed);

      check_profile(matrix, dprofile, dseq);
    }END_TEST

void add_dprofile_fill_TC(Suite *s)
{
  TCase *tc_core = tcase_create("dprofile fill nucleotides");

  tcase_add_checked_fixture(tc_core, &setup, &teardown);

  tcase_add_test(tc_core, test_dprofile_fill_nucleotide_simple);
  tcase_add_test(tc_core, test_dprofile_fill_nucleotide_more);
  tcase_add_test(tc_core, test_dprofile_fill_nucleotide_perf);

  suite_add_tcase(s, tc_core);
}
