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

#include "helper_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <regex.h>
#include <string.h>

void check_cigar_matches(unsigned short pmatches, unsigned short pmismatches, char* pcigar)
{
  regex_t re;
  regmatch_t rm;

  const char * pattern = "[0-9]*M";

  re_set_syntax (RE_SYNTAX_GREP);

  int status = regcomp(&re, pattern, REG_EXTENDED);
  if (status!=REG_NOERROR)
    {
      fail("Could no compile pattern");
    }

  int count = 0;

  int offset = 0;
  while (regexec(&re, pcigar+offset, 1, &rm, REG_EXTENDED)==0)
    {
      char otherString[] = { 0, 0, 0, 0, 0, 0 };
      strncpy(otherString, pcigar+offset+rm.rm_so, rm.rm_eo-rm.rm_so-1);

      if ((rm.rm_eo-rm.rm_so)>1)
        count += atoi(otherString);
      else
        count++;

      offset += rm.rm_eo;
    }

  ck_assert_int_eq(count, pmatches+pmismatches);
}

void print_profile(CELL * dprofile)
{
  for (int i = 0; i<ScoreMatrix::instance.get_dimension(); ++i)
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

void print_matrix(CELL * matrix)
{
  // end copy
  for (int i = 0; i<ScoreMatrix::instance.get_dimension(); i++)
    {
      for (int j = 0; j<ScoreMatrix::instance.get_dimension(); j++)
        {
          printf("%2d, ", matrix[ScoreMatrix::instance.get_dimension()*i+j]);
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

