/*
    Copyright (C) 2014 Torbjorn Rognes

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

struct bucket_s;
struct uhandle_s;

struct uhandle_s * unique_init();

void unique_exit(struct uhandle_s * u);

void unique_count(struct uhandle_s * uh, 
                  int k,
                  int seqlen,
                  char * seq,
                  unsigned int * listlen,
                  unsigned int * * list);

int unique_count_shared(struct uhandle_s * uh,
                        int k,
                        int listlen,
                        unsigned int * list);
