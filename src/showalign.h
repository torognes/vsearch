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

char * align_getrow(char * seq, char * cigar, int alignlen, int origin);

void align_fprint_uncompressed_alignment(FILE * f, char * cigar);

void align_show(FILE * f,
                char * seq1,
                long seq1len,
                long seq1off,
                const char * seq1name,
                char * seq2,
                long seq2len,
                long seq2off,
                const char * seq2name,
                char * cigar,
                long cigarlen,
                int numwidth,
                int namewidth,
                int alignwidth,
                int strand);
