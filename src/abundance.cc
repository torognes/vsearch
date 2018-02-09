/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

bool header_find_attribute(const char * header,
                           int header_length,
                           const char * attribute,
                           int * start,
                           int * end)
{
  /*
    Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
    in the header string, where "size=" is the specified attribute.
  */

  const char * digit_chars = "0123456789";

  if ((! header) || (! attribute))
    return false;

  int hlen = header_length;
  int alen = strlen(attribute);

  int i = 0;

  while (i < hlen - alen)
    {
      char * r = (char *) strstr(header + i, attribute);

      /* no match */
      if (r == NULL)
        break;

      i = r - header;

      /* check for ';' in front */
      if ((i > 0) && (header[i-1] != ';'))
        {
          i += alen + 1;
          continue;
        }

      int digits = (int) strspn(header + i + alen, digit_chars);

      /* check for at least one digit */
      if (digits == 0)
        {
          i += alen + 1;
          continue;
        }

      /* check for ';' after */
      if ((i + alen + digits < hlen) && (header[i + alen + digits] != ';'))
        {
          i += alen + digits + 2;
          continue;
        }

      /* ok */
      * start = i;
      * end = i + alen + digits;
      return true;
    }
  return false;
}

int64_t abundance_get(char * header, int header_length)
{
  /* read size/abundance annotation */
  int64_t abundance = 1;
  int start = 0;
  int end = 0;
  if (header_find_attribute(header, header_length, "size=", & start, & end))
    {
      int64_t number = atol(header + start + 5);
      if (number > 0)
        abundance = number;
      else
        fatal("Invalid (zero) abundance annotation in FASTA file header");
    }
  return abundance;
}

void abundance_fprint_header_strip_size(FILE * fp,
                                        char * header,
                                        int header_length)
{
  int start = 0;
  int end = 0;
  if (header_find_attribute(header, header_length, "size=", & start, & end))
    {
      if (start <= 1)
        {
          if (end < header_length)
            fprintf(fp, "%s", header + end + 1);
        }
      else
        {
          if (end == header_length)
            fprintf(fp, "%.*s", start - 1, header);
          else
            fprintf(fp, "%.*s;%.*s",
                    start - 1, header,
                    header_length - end - 1, header + end + 1);
        }
    }
  else
    fprintf(fp, "%s", header);
}
