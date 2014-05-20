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

void printkmers(unsigned char * kmervector)
{
  /* print kmervector */
  //  fprintf(stderr, "kmer vector: ");
  for(int i = 0; i < KMERVECTORBYTES; i++)
  {
    fprintf(stderr, "%02x", kmervector[i]);
    if ((i % 32) == 31)
      fprintf(stderr, "\n");
  }
}

void findkmers(unsigned char * seq, unsigned long seqlen, 
               unsigned char * kmervector)
{
  /* set kmer bit vector by xoring occurrences of kmers in sequence */
  /* 64 bit version */

  unsigned long * kmerv = (unsigned long *) kmervector;

  memset(kmerv, 0, KMERVECTORBYTES);
  
  unsigned long kmer = 0;
  unsigned long i = 0;

  while((i < KMERLENGTH-1) && (i<seqlen))
  {
    kmer = (kmer << 2) | (seq[i]-1);
    i++;
  }

  while(i < seqlen)
  {
    kmer = (kmer << 2) | (seq[i]-1);
    kmerv[(kmer >> 6) & ((KMERVECTORBYTES/8)-1)] ^= (((unsigned long)1) << (kmer & 63)); 
    i++;
  }
}

void findkmers_8(unsigned char * seq, unsigned long seqlen, 
                   unsigned char * kmervector)
{
  /* set kmer bit vector by xoring occurrences of kmers in sequence */

  memset(kmervector, 0, KMERVECTORBYTES);
  
  unsigned long kmer = 0;
  unsigned long i = 0;

  while((i < KMERLENGTH-1) && (i<seqlen))
  {
    kmer = (kmer << 2) | (seq[i]-1);
    i++;
  }

  while(i < seqlen)
  {
    kmer = (kmer << 2) | (seq[i]-1);
    kmervector[(kmer >> 3) & (KMERVECTORBYTES-1)] ^= (1 << (kmer & 7)); 
    i++;
  }
}

unsigned long comparekmervectors(unsigned char * a, unsigned char * b)
{
  /* count number of different bits */

  unsigned long count = 0;

#if 1

  /* 64 bit version */

  unsigned long *ap = (unsigned long*)a;
  unsigned long *bp = (unsigned long*)b;

#if 0

  /* with ordinary loop */

  while ( (unsigned char*) ap < a + KMERVECTORBYTES)
    count += popcount(*ap++ ^ *bp++);

#else

  /* unrolled loop */

  count += popcount(ap[ 0] ^ bp[ 0]);
  count += popcount(ap[ 1] ^ bp[ 1]);
  count += popcount(ap[ 2] ^ bp[ 2]);
  count += popcount(ap[ 3] ^ bp[ 3]);

#endif

#else

  /* 128 bit version with SSE2 */

  /* input MUST be 16-byte aligned */

  __m128i * ap = (__m128i *) a;
  __m128i * bp = (__m128i *) b;
  __m128i z;

#if 0

  /* with ordinary loop */

  while ( (unsigned char*) ap < a + KMERVECTORBYTES)
    {
      z = _mm_xor_si128(*ap++, *bp++);

#if 0

      /* with 64-bit popcount */

      count += popcount(_mm_extract_epi64(z, 0));
      count += popcount(_mm_extract_epi64(z, 1));

#else

      /* with 128-bit popcount */

      count += popcount_128(z);

#endif

    }

#else

  /* with unrolled loop */

#if 0

  /* with 64-bit popcount */

  z = _mm_xor_si128(*ap++, *bp++);
  count += popcount(_mm_extract_epi64(z, 0));
  count += popcount(_mm_extract_epi64(z, 1));
  z = _mm_xor_si128(*ap++, *bp++);
  count += popcount(_mm_extract_epi64(z, 0));
  count += popcount(_mm_extract_epi64(z, 1));

#else

  /* with 128-bit popcount */

  z = _mm_xor_si128(*ap++, *bp++);
  count += popcount_128(z);
  z = _mm_xor_si128(*ap++, *bp++);
  count += popcount_128(z);

#endif

#endif

#endif

  return count;
}

unsigned long kmer_diff(unsigned long a, unsigned long b)
{
  unsigned long diffkmers = comparekmervectors(db_getkmervector(a),
                                                 db_getkmervector(b));
  unsigned long mindiff = (diffkmers + 2*KMERLENGTH - 1)/(2*KMERLENGTH);
  return mindiff;
}
