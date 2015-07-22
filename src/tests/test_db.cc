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

/* Is run once before each unit test */
static void setup()
{
  opt_maxseqlength = 2000;
}

/* Is run once after each unit test */
static void teardown()
{
  db_free();
}

START_TEST (test_db_nuc)
    {
      db_read("../data/AF091148.fsa", 0); // TODO what does upcase do?

      ck_assert_str_eq("97485665bcded44c4d86c131ca714848", db_getheader(0));
      ck_assert_str_eq("443ddf5898dde8ad55a9abca6acb246a", db_getheader(1));

      ck_assert_str_eq("gtcgctcctaccgattgaatacgttggtgattgaattggataaagagatatcatcttaaatgatagcaaagcggtaaacatttgtaaactagattatttagaggaaggagaagtcgtaacaaggtttcc",
                       db_getsequence(0));
      ck_assert_str_eq("gtcgctcctaccgattgaatacattggtgattggattggataaagagatatcttcttaaatgataacaaaacggtaaacatttgtaaactagattatttagaggaaggagaagtcgtaacaaggtttcc",
                       db_getsequence(1));

      ck_assert_int_eq(1, db_getabundance(0));
      ck_assert_int_eq(129, db_getsequencelen(0));
      ck_assert_int_eq(32, db_getheaderlen(0));

      ck_assert_int_eq(32, db_getlongestheader());
      ck_assert_int_eq(137, db_getlongestsequence());
      ck_assert_int_eq(103, db_getshortestsequence());

      ck_assert_int_eq(1403, db_getsequencecount());
      ck_assert_int_eq(180704, db_getnucleotidecount());

    }END_TEST

START_TEST (test_db_aa)
    {
      opt_maxseqlength = 30000;

      db_read("../data/uniprot_sprot.fasta", 0); // TODO what does upcase do?

      ck_assert_int_eq(547964, db_getsequencecount());
      ck_assert_int_eq(195174196, db_getnucleotidecount());

      printf("db_getlongestheader: %ld\n", db_getlongestheader());
      printf("db_getlongestsequence: %ld\n", db_getlongestsequence());
      printf("db_getshortestsequence: %ld\n", db_getshortestsequence());

      printf("abundance: %ld\n", db_getabundance(0));
      printf("seq len: %ld\n", db_getsequencelen(0));
      printf("header len: %ld\n", db_getheaderlen(0));

    }END_TEST

void add_db_TC(Suite *s)
{
  TCase *tc_core = tcase_create("db");

  tcase_add_checked_fixture(tc_core, &setup, &teardown);

  tcase_add_test(tc_core, test_db_nuc);
  tcase_add_test(tc_core, test_db_aa);

  suite_add_tcase(s, tc_core);
}
