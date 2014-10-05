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

extern unsigned int * kmercount; /* number of matching seqnos for each kmer */
extern unsigned int * kmerhash;  /* index into the list below for each kmer */
extern unsigned int * kmerindex; /* the list of matching seqnos for kmers */
extern bitmap_t * * kmerbitmap;
extern unsigned int * dbindex_map;
extern unsigned int dbindex_count;

void fprint_kmer(FILE * f, unsigned int k, unsigned long kmer);

void dbindex_prepare(int use_bitmap);
void dbindex_addallsequences();
void dbindex_addsequence(unsigned int seqno);
void dbindex_free();

inline unsigned char * dbindex_getbitmap(unsigned int kmer)
{
  if (kmerbitmap[kmer])
    return kmerbitmap[kmer]->bitmap;
  else
    return 0;
}

inline unsigned int dbindex_getmatchcount(unsigned int kmer)
{
  return kmercount[kmer];
}

inline unsigned int * dbindex_getmatchlist(unsigned int kmer)
{
  return kmerindex + kmerhash[kmer];
}

inline unsigned int dbindex_getmapping(unsigned int index)
{
  return dbindex_map[index];
}

inline unsigned int dbindex_getcount()
{
  return dbindex_count;
}
