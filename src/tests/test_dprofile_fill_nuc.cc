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

#include "helper_functions.h"

#include "../vsearch.h"
#include "../align_simd_dprofile.h"
#include "../util.h"
#include "../score_matrix.h"

extern void dprofile_fill16(CELL * dprofile_word,
                            CELL * score_matrix_word,
                            BYTE * dseq);

extern void dprofile_fill16_aa(CELL * dprofile_word,
                               CELL * score_matrix_word,
                               BYTE * dseq);

static BYTE dseq[CDEPTH*CHANNELS];

/* Is run once before each unit test */
static void setup()
{
  opt_match = 5;
  opt_mismatch = -4;

  ScoreMatrix::instance.init(opt_match, opt_mismatch, MATRIX_MODE_NUC);

  memset(dseq, 0, CDEPTH*CHANNELS);
}

/* Is run once after each unit test */
static void teardown()
{
}

static void check_profile(CELL * matrix, CELL* dprofile, BYTE* dseq)
{
//        print_profile(dprofile);
  for (int i = 0; i<ScoreMatrix::instance.get_dimension(); ++i)
    {
      for (int j = 0; j<CDEPTH; ++j)
        {
          for (int k = 0; k<CHANNELS; k++)
            {
              CELL pval = dprofile[CHANNELS*CDEPTH*i+CHANNELS*j+k];
              CELL dval = dseq[CHANNELS*j+k];

              ck_assert_int_eq(matrix[ScoreMatrix::instance.get_dimension()*dval+i], pval);
            }
        }
    }
}

static void fill_search_window(int dseq_count, int db_sequences[][CDEPTH])
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

      int db_sequences[][CDEPTH] = { { 'A', 'C', 'A', 'T' }, };
      fill_search_window(dseq_count, db_sequences);

      CELL * dprofile = (CELL*)xmalloc(2*4*8*16);

      dprofile_fill16(dprofile, ScoreMatrix::instance.score_matrix_16, dseq);

      check_profile(ScoreMatrix::instance.score_matrix_16, dprofile, dseq);
    }END_TEST

START_TEST (test_dprofile_fill_amino_acids_simple)
    {
      int dseq_count = 1;

      int db_sequences[][CDEPTH] = { { 'Q', 'R', 'S', 'T' }, };
      fill_search_window(dseq_count, db_sequences);

      CELL * dprofile = (CELL*)xmalloc(2*4*8*32);

      dprofile_fill16(dprofile, ScoreMatrix::instance.score_matrix_16, dseq);

      check_profile(ScoreMatrix::instance.score_matrix_16, dprofile, dseq);
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
          { 'C', 'A', 'C', 'C' } };
      fill_search_window(dseq_count, db_sequences);

      CELL * dprofile = (CELL*)xmalloc(sizeof(CELL)*CDEPTH*CHANNELS*ScoreMatrix::instance.get_dimension());

      dprofile_fill16(dprofile, ScoreMatrix::instance.score_matrix_16, dseq);

      check_profile(ScoreMatrix::instance.score_matrix_16, dprofile, dseq);
    }END_TEST

START_TEST (test_dprofile_fill_amino_acids_more)
    {
      int dseq_count = CHANNELS;

      int db_sequences[][CDEPTH] = {
          { 'A', 'Z', 'W', 'T' },
          { 'K', 'R', 'L', 'C' },
          { 'T', 'T', 'N', 'T' },
          { 'A', 'A', 0, 0 },
          { 'Q', 'T', 'C', 'N' },
          { 'T', 'C', 'T', 'C' },
          { 'G', 'U', 'U', 0 },
          { 'M', 'B', 'V', 'C' } };
      fill_search_window(dseq_count, db_sequences);

      CELL * dprofile = (CELL*)xmalloc(sizeof(CELL)*CDEPTH*CHANNELS*ScoreMatrix::instance.get_dimension());

      dprofile_fill16(dprofile, ScoreMatrix::instance.score_matrix_16, dseq);

      check_profile(ScoreMatrix::instance.score_matrix_16, dprofile, dseq);
    }END_TEST

static void run_perf_test(void (*dprofile_fill_func)(CELL *, CELL *, BYTE *), const char * desc)
{
  int dseq_count = 1;

  int db_sequences[][CDEPTH] = { { 'A', 'C', 'A', 'T' }, };
  fill_search_window(dseq_count, db_sequences);

  CELL * dprofile = (CELL*)xmalloc(sizeof(CELL)*CDEPTH*CHANNELS*ScoreMatrix::instance.get_dimension());

  struct timeval start;
  struct timeval finish;

  gettimeofday(&start, NULL);

  long rounds = 10000000;
  for (int i = 0; i<rounds; ++i)
    {
      dprofile_fill_func(dprofile, ScoreMatrix::instance.score_matrix_16, dseq);
    }

  gettimeofday(&finish, NULL);

  double elapsed = (finish.tv_sec-start.tv_sec);
  elapsed += (finish.tv_usec-start.tv_usec)/1000000.0;

  printf("\nRuntime used for %ld runs of dprofile_fill16 for %s: %lf sec\n\n", rounds, desc, elapsed);

  check_profile(ScoreMatrix::instance.score_matrix_16, dprofile, dseq);
}

START_TEST (test_dprofile_fill_nucleotide_perf)
    {
      run_perf_test(&dprofile_fill16, "nucleotides");
    }END_TEST

START_TEST (test_dprofile_fill_amino_acids_perf)
    {
      ScoreMatrix::instance.init(opt_match, opt_mismatch, MATRIX_MODE_AA);

      run_perf_test(&dprofile_fill16_aa, "amino acids");
    }END_TEST

void add_dprofile_fill_nuc_TC(Suite *s)
{
  TCase *tc_core = tcase_create("dprofile fill for nucleotide sequences");

  tcase_add_checked_fixture(tc_core, &setup, &teardown);

  tcase_add_test(tc_core, test_dprofile_fill_nucleotide_simple);
  tcase_add_test(tc_core, test_dprofile_fill_amino_acids_simple);
  tcase_add_test(tc_core, test_dprofile_fill_nucleotide_more);
  tcase_add_test(tc_core, test_dprofile_fill_amino_acids_more);
  tcase_add_test(tc_core, test_dprofile_fill_nucleotide_perf);
  tcase_add_test(tc_core, test_dprofile_fill_amino_acids_perf);

  suite_add_tcase(s, tc_core);
}
