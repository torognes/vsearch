/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.


  GNU General Public License version 3

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The BSD 2-Clause License

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "vsearch.h"

const char * tax_letters = "dkpcofgs";

bool tax_parse(const char * header,
               int header_length,
               int * tax_start,
               int * tax_end)
{
  /*
    Identify the first occurence of the pattern (^|;)tax=([^;]*)(;|$)
  */

  if (! header)
    {
      return false;
    }

  const char * attribute = "tax=";

  int hlen = header_length;
  int alen = strlen(attribute);

  int i = 0;

  while (i < hlen - alen)
    {
      char * r = (char *) strstr(header + i, attribute);

      /* no match */
      if (r == nullptr)
        {
          break;
        }

      i = r - header;

      /* check for ';' in front */
      if ((i > 0) && (header[i-1] != ';'))
        {
          i += alen + 1;
          continue;
        }

      * tax_start = i;

      /* find end (semicolon or end of header) */
      const char * s = strchr(header+i+alen, ';');
      if (s == nullptr)
        {
          * tax_end = hlen;
        }
      else
        {
          * tax_end = s - header;
        }

      return true;
    }
  return false;
}

void tax_split(int seqno, int * level_start, int * level_len)
{
  /* Parse taxonomy string into the following parts
     d domain
     k kingdom
     p phylum
     c class
     o order
     f family
     g genus
     s species
  */

  for (int i = 0; i < tax_levels; i++)
    {
      level_start[i] = 0;
      level_len[i] = 0;
    }

  int tax_start;
  int tax_end;
  char * h = db_getheader(seqno);
  int hlen = db_getheaderlen(seqno);
  if (tax_parse(h, hlen, & tax_start, & tax_end))
    {
      int t = tax_start + 4;

      while (t < tax_end)
        {
          /* Is the next char a recogized tax level letter? */
          const char * r = strchr(tax_letters, tolower(h[t]));
          if (r)
            {
              int level = r - tax_letters;

              /* Is there a colon after it? */
              if (h[t + 1] == ':')
                {
                  level_start[level] = t + 2;

                  char * z = strchr(h + t + 2, ',');
                  if (z)
                    {
                      level_len[level] = z - h - t - 2;
                    }
                  else
                    {
                      level_len[level] = tax_end - t - 2;
                    }
                }
            }

          /* skip past next comma */
          char * x = strchr(h + t, ',');
          if (x)
            {
              t = x - h + 1;
            }
          else
            {
              t = tax_end;
            }
        }
    }
}
