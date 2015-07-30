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

#include "helper_functions.h"

#include "../vsearch.h"
#include "../align_simd.h"
#include "../score_matrix.h"
#include "../util.h"

static struct s16info_s * s16;

/* Is run once before each unit test */
static void setup()
{
  opt_maxseqlength = 500;

  db_read("../data/test_nucleotide_db.fasta", 0); // TODO what does upcase do?

  CELL match = 5;
  CELL mismatch = -4;
  CELL gap_open = 5;
  CELL gap_extension = 1;

  ScoreMatrix::instance.init(match, mismatch, MATRIX_MODE_NUC);

  s16 = search16_init_2(gap_open, gap_open, gap_open, gap_open, gap_open, gap_open,
      gap_extension, gap_extension, gap_extension, gap_extension, gap_extension, gap_extension);
}

/* Is run once after each unit test */
static void teardown()
{
  search16_exit(s16);

  db_free();
}

START_TEST (test_align_simd_simple)
    {
      char * query = (char *)"ACAC";
      search16_qprep(s16, query, strlen(query));

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

      ck_assert_int_eq(11, pscores[0]);
      ck_assert_str_eq("4M", pcigar[0]);

      check_cigar_matches(pmatches[0], pmismatches[0], pcigar[0]);
    }END_TEST

START_TEST (test_align_simd_all)
    {
      char * query = (char *)"ACAT";
      search16_qprep(s16, query, strlen(query));

      unsigned long seq_count = db_getsequencecount();

      unsigned int seqnos[seq_count];
      for (unsigned int i = 0; i<seq_count; ++i)
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

      CELL exp_scores[] = { 20, 4, 2, -7, 2, -7, -16, -4, -5, 0 };
      const char * exp_cigars[] = { "4M", "2M2I2M", "4M", "4M", "4M", "4M", "4M", "M2IM2D", "M2I3M", "M4I2MD" };

      for (unsigned int i = 0; i<seq_count; ++i)
        {
          ck_assert_ptr_ne(0, pcigar[i]);
          check_cigar_matches(pmatches[i], pmismatches[i], pcigar[i]);

          ck_assert_int_eq(exp_scores[i], pscores[i]);
          ck_assert_str_eq(exp_cigars[i], pcigar[i]);
        }
    }END_TEST

START_TEST (test_align_simd_AF091148)
    {
      db_read("../data/AF091148.fsa", 0);

      CELL match = 2;
      CELL mismatch = -2;
      CELL gap_open = 4;
      CELL gap_extension = 2;

      s16 = search16_init(match, mismatch,
          gap_open, gap_open, gap_open, gap_open, gap_open, gap_open,
          gap_extension, gap_extension, gap_extension, gap_extension, gap_extension, gap_extension);

      char * query = (char *)"ATGCCCAAGCTGAATAGCGTAGAGGGGTTTTCATCATTTGAGGACGATGTATAA";
      search16_qprep(s16, query, strlen(query));

      unsigned long seq_count = db_getsequencecount();

      unsigned int seqnos[seq_count];
      for (unsigned int i = 0; i<seq_count; ++i)
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

      /*
       * checks scores and cigar strings of the top matching sequences (score >= -112).
       */
      CELL exp_scores[] = { 8, -112, 23, -104, 75, -102, 378, -98, 612, -110, 908, -88, 938, -94, 1016, -112, 1050, -92,
          1069, -112, 1146, -112, 1148, -110, 1229, -106 };
      const char * exp_cigars[] = { "2MI2MI6M10I5M7I4M16I3M3I6M5I3M8I8MI4M14I11M",
          "2MI2MI6M3I6MD5M16I12M24I6M6I2M8I12M3I", "2MI2MI6M3I18M6I2M17I10M21I3M10I11M",
          "2MI2MI6M3I5M7I4M16I3M3I6M5I3M8I8MI4M14I11M", "2MI2MI6M3I6M4I12M6I2M17I10M21I3M10I11M",
          "2MI2MI6M3I6MD5M16I13M22I7M7I12M", "2MI2MI6M3I6MD5M16I12M24I11M2I3MD5M10I",
          "2MI2MI6M10I5M7I4M16I3M3I6M5I3M8I8MI4M14I11M", "2MI2MI6M3I6MD3M23I8M5I3M8I8MI4M14I11M",
          "2MI2MI6M3I6MD5M16I13M30I4MI4M14I11M", "2MI2MI6M10I5M7I4M16I3M3I6M5I3M8I8MI4M14I11M",
          "2MI2MI6M3I6M5I3M15I6M10I3M5I3M8I8MI4M14I11M", "2MI2MI6M3I6M16I12M6I2M17I10M6I2M8I12M3I" };

      int hit_counter = 0;

      for (unsigned int i = 0; i<seq_count; ++i)
        {
          ck_assert_ptr_ne(0, pcigar[i]);
          check_cigar_matches(pmatches[i], pmismatches[i], pcigar[i]);

          if( pscores[i] >= -112)
            {
              ck_assert_int_eq(exp_scores[2*hit_counter], i);
              ck_assert_int_eq(exp_scores[2*hit_counter+1], pscores[i]);
              ck_assert_str_eq(exp_cigars[hit_counter], pcigar[i]);

              hit_counter++;
            }
        }
    }END_TEST

void add_align_simd_nuc_TC(Suite *s)
{
  TCase *tc_core = tcase_create("align simd for nucleotide sequences");

  tcase_add_checked_fixture(tc_core, &setup, &teardown);

  tcase_add_test(tc_core, test_align_simd_simple);
  tcase_add_test(tc_core, test_align_simd_all);
  tcase_add_test(tc_core, test_align_simd_AF091148);

  suite_add_tcase(s, tc_core);
}
