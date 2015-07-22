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

#include "../maps.h"

static const char * legal_nuc_symbols = "abcdghkmnrstuvwyABCDGHKMNRSTUVWY";
static const char * legal_aa_symbols = "abcdefghiklmnopqrstuvwxyzABCDEFGHIKLMNOPQRSTUVWXYZ";

static void check_chrstatus(const char * legal_symbols, unsigned int chr_status_map[256])
{
  for (int i = 0; i<256; ++i)
    {
      int status = chr_status_map[i];

      if (i<32)
        {
          if (i>=9&&i<=13)
            ck_assert_int_eq(3, status);
          else
            ck_assert_int_eq(2, status);
        }
      else if (i==45||i==46)
        ck_assert_int_eq(2, status);
      else
        {

          int symbol_is_legal = 0;

          for (unsigned int j = 0; j<strlen(legal_symbols); ++j)
            {
              if (i==legal_symbols[j])
                {
                  symbol_is_legal = 1;
                  break;
                }
            }
          ck_assert_int_eq(symbol_is_legal, status);
        }
    }
}

START_TEST (test_chrstatus)
    {
      check_chrstatus(legal_nuc_symbols, chrstatus);
    }END_TEST

START_TEST (test_chrstatus_aa)
    {
      check_chrstatus(legal_aa_symbols, chrstatus_aa);
    }END_TEST

void add_maps_TC(Suite *s)
{
  TCase *tc_core = tcase_create("maps");

  tcase_add_test(tc_core, test_chrstatus);
  tcase_add_test(tc_core, test_chrstatus_aa);

  suite_add_tcase(s, tc_core);
}
