/*
    Copyright (C) 2014-2015 Torbjorn Rognes and Tomas Flouri

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

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

typedef struct abundance_s
{
  regex_t regex;
} abundance_t;

abundance_t * abundance_init(void);

void abundance_exit(abundance_t * a);

long abundance_get(abundance_t * a, char * header);

void abundance_fprint_header_with_size(abundance_t * a,
                                       FILE * fp,
                                       char * header,
                                       int header_length,
                                       unsigned long size);

void abundance_fprint_header_strip_size(abundance_t * a,
                                        FILE * fp,
                                        char * header,
                                        int header_length);

char * abundance_strip_size(abundance_t * a,
                            char * header,
                            int header_length);
