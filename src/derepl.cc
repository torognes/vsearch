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

#include "vsearch.h"

/* replicate sequences */

/* 

   Procedure:

   1. read all sequences into memory
   2. compute a hash for each sequence (sha1 or fnv_1a_64 or whatever)
   3. use a hash table to count the abundance of each sequence
   4. in each bucket of the hash table, store: (1) the hash value,
      (2) the abundance, and (3) a pointer to the first sequence
   5. use an additional table with one entry for each sequence with pointers
      to the "next" identical sequence
   6. In the case two different sequences have the same hash (collision),
      we need to check that the sequences are actually identical after identical
      hashes has been found. Different sequences should have separate entries in
      the hash table.
   7. Use linear probing and a hash table twice the size of the data.
   8. Sort the hash table by abundance (decreasing), hash value, sequence
   9. Output

*/

struct bucket
{
  unsigned long hash;
  long abundance;
  unsigned long seqno;
};



void db_dereplicate()
{
  db_read(databasefilename);
  long dbsequencecount = db_getsequencecount();
  
  hashtablesize = dbsequencecount * 2;

  struct bucket * hashtable =
    (struct bucket *) xmemalloc(sizeof(bucket) * hashtablesize);

  memset(hashtable, 0, sizeof(bucket) * hashtablesize);

  unsigned long nextseqno =
    (unsigned long *) xmemalloc(sizeof(unsigned long) * dbsequencecount);

  memset(nextseqno, 0, sizeof(unsigned long) * dbsequenccount);

  
  for(long i=0; i<dbsequencecount; i++)
    {
      char * seq;
      long seqlen;

      db_getsequenceandlength(i, seq, & seqlen);

      hash = hashfunc(seq, seqlen);

      j = hash % hashtablesize;

      ...
      
    }
  

  free(hashtable);
  
  db_free();

}
