/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

abundance_t * abundance_init(void)
{
 abundance_t * a = (abundance_t *) xmalloc(sizeof(abundance_t));
 if (regcomp(&a->regex, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
   fatal("Compilation of regular expression for abundance annotation failed");
 return a;
}

void abundance_exit(abundance_t * a)
{
  regfree(&a->regex);
  xfree(a);
}

int64_t abundance_get(abundance_t * a, char * header)
{
  /* read size/abundance annotation */
  
  int64_t abundance = 1;
  regmatch_t pmatch[4];
  
  if (!regexec(&a->regex, header, 4, pmatch, 0))
    {
      int64_t number = atol(header + pmatch[2].rm_so);
      if (number > 0)
        abundance = number;
      else
        fatal("Invalid (zero) abundance annotation in fasta header");
    }
  return abundance;
}

void abundance_fprint_header_with_size(abundance_t * a,
                                       FILE * fp,
                                       char * header,
                                       int header_length,
                                       uint64_t size)
{
  /* remove any previous size annotation */
  /* regexp search for "(^|;)(\d+)(;|$)" */
  /* replace by ';' if not at either end */

  regmatch_t pmatch[1];

  if (!regexec(&a->regex, header, 1, pmatch, 0))
    {
      int pat_start = pmatch[0].rm_so;
      int pat_end = pmatch[0].rm_eo;

      fprintf(fp,
              "%.*s%s%.*s%ssize=%" PRIu64 ";",
              pat_start, header,
              (pat_start > 0 ? ";" : ""),
              header_length - pat_end, header + pat_end,
              (((pat_end < header_length) &&
                (header[header_length - 1] != ';')) ? ";" : ""),
              size);
    }
  else
    {
      fprintf(fp,
              "%s%ssize=%" PRIu64 ";",
              header,
              (((header_length == 0) || 
                (header[header_length - 1] != ';')) ? ";" : ""),
              size);
    }
}

void abundance_fprint_header_strip_size(abundance_t * a,
                                        FILE * fp,
                                        char * header,
                                        int header_length)
{
  regmatch_t pmatch[1];

  if (!regexec(&a->regex, header, 1, pmatch, 0))
    {
      int pat_start = pmatch[0].rm_so;
      int pat_end = pmatch[0].rm_eo;

      fprintf(fp,
              "%.*s%s%.*s",
              pat_start, header,
              ((pat_start > 0) && (pat_end < header_length)) ? ";" : "",
              header_length - pat_end, header + pat_end);
    }
  else
    fprintf(fp, "%s", header);
}

char * abundance_strip_size(abundance_t * a,
                            char * header,
                            int header_length)
{
  int ret;
  char * temp = 0;
  regmatch_t pmatch[1];
  

  if (!regexec(&a->regex, header, 1, pmatch, 0))
    {
      int pat_start = pmatch[0].rm_so;
      int pat_end = pmatch[0].rm_eo;

      ret = xsprintf(&temp,
                     "%.*s%s%.*s",
                     pat_start, header,
                     ((pat_start > 0) && (pat_end < header_length)) ? ";" : "",
                     header_length - pat_end, header + pat_end);
    }
  else
    ret = xsprintf(&temp, "%s", header);
  
  if (ret == -1)
    fatal("Out of memory");
  
  return temp;
}
