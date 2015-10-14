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

#define MEMCHUNK 16777216

static fasta_handle h;
static unsigned long sequences = 0;
static unsigned long nucleotides = 0;
static unsigned long longest;
static unsigned long shortest;
static unsigned long longestheader;

seqinfo_t * seqindex;
char * datap;

void db_read(const char * filename, int upcase)
{
  /* compile regexp for abundance pattern */

  h = fasta_open(filename);

  long filesize = fasta_get_size(h);
  
  char * prompt;
  if (asprintf(& prompt, "Reading file %s", filename) == -1)
    fatal("Out of memory");

  progress_init(prompt, filesize);

  longest = 0;
  shortest = LONG_MAX;
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;

  long discarded_short = 0;
  long discarded_long = 0;

  /* allocate space for data */
  unsigned long dataalloc = 0;
  datap = 0;
  unsigned long datalen = 0;

  /* allocate space for index */
  size_t seqindex_alloc = 0;
  seqindex = 0;

  while(fasta_next(h,
                   ! opt_notrunclabels,
                   upcase ? chrmap_upcase : chrmap_no_change))
    {
      size_t headerlength = fasta_get_header_length(h);
      size_t sequencelength = fasta_get_sequence_length(h);

      unsigned int abundance = abundance_get(global_abundance,
                                             fasta_get_header(h));
      
      if (sequencelength < (size_t)opt_minseqlength)
        {
          discarded_short++;
        }
      else if (sequencelength > (size_t)opt_maxseqlength)
        {
          discarded_long++;
        }
      else
        {
          /* grow space for data, if necessary */
          size_t dataalloc_old = dataalloc;
          while (datalen + headerlength + sequencelength + 2 > dataalloc)
            dataalloc += MEMCHUNK;
          if (dataalloc > dataalloc_old)
            datap = (char *) xrealloc(datap, dataalloc);

          /* store the header */
          size_t header_p = datalen;
          memcpy(datap + header_p, fasta_get_header(h), headerlength + 1);
          datalen += headerlength + 1;
          
          /* store sequence */
          size_t sequence_p = datalen;
          memcpy(datap+sequence_p, fasta_get_sequence(h), sequencelength + 1);
          datalen += sequencelength + 1;
          
          /* grow space for index, if necessary */
          size_t seqindex_alloc_old = seqindex_alloc;
          while ((sequences + 1) * sizeof(seqinfo_t) > seqindex_alloc)
            seqindex_alloc += MEMCHUNK;
          if (seqindex_alloc > seqindex_alloc_old)
            seqindex = (seqinfo_t *) xrealloc(seqindex, seqindex_alloc);
          
          /* update index */
          seqinfo_t * seqindex_p = seqindex + sequences;
          seqindex_p->headerlen = headerlength;
          seqindex_p->header_p = header_p;
          seqindex_p->seqlen = sequencelength;
          seqindex_p->seq_p = sequence_p;
          seqindex_p->size = abundance;

          /* update statistics */
          sequences++;
          nucleotides += sequencelength;
          if (sequencelength > longest)
            longest = sequencelength;
          if (sequencelength < shortest)
            shortest = sequencelength;
          if (headerlength > longestheader)
            longestheader = headerlength;
        }
      progress_update(fasta_get_position(h));
    }

  progress_done();
  free(prompt);
  fasta_close(h);

  if (!opt_quiet)
    {
      if (sequences > 0)
        fprintf(stderr,
                "%'lu nt in %'lu seqs, min %'lu, max %'lu, avg %'.0f\n", 
                db_getnucleotidecount(),
                db_getsequencecount(),
                db_getshortestsequence(),
                db_getlongestsequence(),
                db_getnucleotidecount() * 1.0 / db_getsequencecount());
      else
        fprintf(stderr,
                "%'lu nt in %'lu seqs\n", 
                db_getnucleotidecount(),
                db_getsequencecount());
    }

  if (opt_log)
    {
      if (sequences > 0)
        fprintf(fp_log,
                "%'lu nt in %'lu seqs, min %'lu, max %'lu, avg %'.0f\n\n", 
                db_getnucleotidecount(),
                db_getsequencecount(),
                db_getshortestsequence(),
                db_getlongestsequence(),
                db_getnucleotidecount() * 1.0 / db_getsequencecount());
      else
        fprintf(fp_log,
                "%'lu nt in %'lu seqs\n\n", 
                db_getnucleotidecount(),
                db_getsequencecount());
    }

  /* Warn about discarded sequences */

  if (discarded_short)
    {
      fprintf(stderr,
              "WARNING: %ld sequences shorter than %ld nucleotides discarded.\n",
              discarded_short, opt_minseqlength);

      if (opt_log)
        fprintf(fp_log,
                "WARNING: %ld sequences shorter than %ld nucleotides discarded.\n\n",
                discarded_short, opt_minseqlength);
    }
  
  if (discarded_long)
    {
      fprintf(stderr,
              "WARNING: %ld sequences longer than %ld nucleotides discarded.\n",
              discarded_long, opt_maxseqlength);

      if (opt_log)
        fprintf(fp_log,
                "WARNING: %ld sequences longer than %ld nucleotides discarded.\n\n",
                discarded_long, opt_maxseqlength);
    }

  show_rusage();

}

unsigned long db_getsequencecount()
{
  return sequences;
}

unsigned long db_getnucleotidecount()
{
  return nucleotides;
}

unsigned long db_getlongestheader()
{
  return longestheader;
}

unsigned long db_getlongestsequence()
{
  return longest;
}

unsigned long db_getshortestsequence()
{
  return shortest;
}

void db_free()
{
  if (datap)
    free(datap);
  if (seqindex)
    free(seqindex);
}

int compare_bylength(const void * a, const void * b)
{
  seqinfo_t * x = (seqinfo_t *) a;
  seqinfo_t * y = (seqinfo_t *) b;

  /* longest first, then by abundance, then by label, otherwise keep order */

  if (x->seqlen < y->seqlen)
    return +1;
  else if (x->seqlen > y->seqlen)
    return -1;
  else
    {
      if (x->size < y->size)
        return +1;
      else if (x->size > y->size)
        return -1;
      else
        {
          int r = strcmp(datap + x->header_p, datap + y->header_p);
          if (r != 0)
            return r;
          else
            {
              if (x < y)
                return -1;
              else if (x > y)
                return +1;
              else
                return 0;
            }
        }
    }
}

int compare_bylength_shortest_first(const void * a, const void * b)
{
  seqinfo_t * x = (seqinfo_t *) a;
  seqinfo_t * y = (seqinfo_t *) b;

  /* shortest first, then by abundance, then by label, otherwise keep order */

  if (x->seqlen < y->seqlen)
    return -1;
  else if (x->seqlen > y->seqlen)
    return +1;
  else
    {
      if (x->size < y->size)
        return +1;
      else if (x->size > y->size)
        return -1;
      else
        {
          int r = strcmp(datap + x->header_p, datap + y->header_p);
          if (r != 0)
            return r;
          else
            {
              if (x < y)
                return -1;
              else if (x > y)
                return +1;
              else
                return 0;
            }
        }
    }
}

inline int compare_byabundance(const void * a, const void * b)
{
  seqinfo_t * x = (seqinfo_t *) a;
  seqinfo_t * y = (seqinfo_t *) b;

  /* most abundant first, then by label, otherwise keep order */

  if (x->size > y->size)
    return -1;
  else if (x->size < y->size)
    return +1;
  else
    {
      int r = strcmp(datap + x->header_p, datap + y->header_p);
      if (r != 0)
        return r;
      else
        {
          if (x < y)
            return -1;
          else if (x > y)
            return +1;
          else
            return 0;
        }
    }
}

void db_sortbylength()
{
  progress_init("Sorting by length", 100);
  qsort(seqindex, sequences, sizeof(seqinfo_t), compare_bylength);
  progress_done();
}

void db_sortbylength_shortest_first()
{
  progress_init("Sorting by length", 100);
  qsort(seqindex, sequences, sizeof(seqinfo_t), compare_bylength_shortest_first);
  progress_done();
}

void db_sortbyabundance()
{
  progress_init("Sorting by abundance", 100);
  qsort(seqindex, sequences, sizeof(seqinfo_t), compare_byabundance);
  progress_done();
}


