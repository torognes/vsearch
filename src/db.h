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

struct seqinfo_s
{
  char * header;
  char * seq;
  unsigned long headerlen;
  unsigned long headeridlen;
  unsigned long seqlen;
  unsigned long size;
};

typedef struct seqinfo_s seqinfo_t;

extern seqinfo_t * seqindex;
extern regex_t db_regexp;

inline char * db_getheader(unsigned long seqno)
{
  return seqindex[seqno].header;
}

inline char * db_getsequence(unsigned long seqno)
{
  return seqindex[seqno].seq;
}

inline unsigned long db_getabundance(unsigned long seqno)
{
  return seqindex[seqno].size;
}

inline unsigned long db_getsequencelen(unsigned long seqno)
{
  return seqindex[seqno].seqlen;
}

inline unsigned long db_getheaderlen(unsigned long seqno)
{
  return seqindex[seqno].headerlen;
}

void db_read(const char * filename, int upcase);
void db_free();

unsigned long db_getsequencecount();
unsigned long db_getnucleotidecount();
unsigned long db_getlongestheader();
unsigned long db_getlongestsequence();
unsigned long db_getshortestsequence();
