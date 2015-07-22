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

      ck_assert_str_eq(
          "gtcgctcctaccgattgaatacgttggtgattgaattggataaagagatatcatcttaaatgatagcaaagcggtaaacatttgtaaactagattatttagaggaaggagaagtcgtaacaaggtttcc",
          db_getsequence(0));
      ck_assert_str_eq(
          "gtcgctcctaccgattgaatacattggtgattggattggataaagagatatcttcttaaatgataacaaaacggtaaacatttgtaaactagattatttagaggaaggagaagtcgtaacaaggtttcc",
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
      opt_maxseqlength = 40000;
      opt_notrunclabels = 1;

      db_read("../data/uniprot_sprot.fasta", 0, DB_MODE_AA);

      ck_assert_int_eq(547964, db_getsequencecount());
      ck_assert_int_eq(195174196, db_getnucleotidecount());

      ck_assert_int_eq(280, db_getlongestheader());
      ck_assert_int_eq(35213, db_getlongestsequence());
      ck_assert_int_eq(2, db_getshortestsequence());

      ck_assert_int_eq(110, db_getheaderlen(0));
      ck_assert_int_eq(256, db_getsequencelen(0));
      ck_assert_int_eq(1, db_getabundance(0));

      ck_assert_str_eq(
          "sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) GN=FV3-001R PE=4 SV=1",
          db_getheader(0));
      ck_assert_str_eq(
          "sp|Q6GZX3|002L_FRG3G Uncharacterized protein 002L OS=Frog virus 3 (isolate Goorha) GN=FV3-002L PE=4 SV=1",
          db_getheader(1));

      ck_assert_str_eq(
          "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL",
          db_getsequence(0));
      ck_assert_str_eq(
          "MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASGFCTSQPLCARIKKTQVCGLRYSSKGKDPLVSAEWDSRGAPYVRCTYDADLIDTQAQVDQFVSMFGESPSLAERYCMRGVKNTAGELVSRVSSDADPAGGWCRKWYSAHRGPDQDAALGSFCIKNPGAADCKCINRASDPVYQKVKTLHAYPDQCWYVPCAADVGELKMGTQRDTPTNCPTQVCQIVFNMLDDGSVTMDDVKNTINCDFSKYVPPPPPPKPTPPTPPTPPTPPTPPTPPTPPTPRPVHNRKVMFFVAGAVLVAILISTVRW",
          db_getsequence(1));

    }END_TEST

START_TEST (test_db_aa_trunc_header)
    {
      opt_maxseqlength = 40000;
      opt_notrunclabels = 0;

      db_read("../data/uniprot_sprot.fasta", 0, DB_MODE_AA);

      ck_assert_int_eq(547964, db_getsequencecount());
      ck_assert_int_eq(195174196, db_getnucleotidecount());

      ck_assert_int_eq(25, db_getlongestheader());
      ck_assert_int_eq(20, db_getheaderlen(0));

      ck_assert_str_eq("sp|Q6GZX4|001R_FRG3G",
          db_getheader(0));
      ck_assert_str_eq("sp|Q6GZX3|002L_FRG3G",
          db_getheader(1));
    }END_TEST

void add_db_TC(Suite *s)
{
  TCase *tc_core = tcase_create("db");

  tcase_add_checked_fixture(tc_core, &setup, &teardown);

  tcase_add_test(tc_core, test_db_nuc);
  tcase_add_test(tc_core, test_db_aa);
  tcase_add_test(tc_core, test_db_aa_trunc_header);

  suite_add_tcase(s, tc_core);
}
