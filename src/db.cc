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

#define MEMCHUNK 16777216

static fastx_handle h = 0;
static bool is_fastq = 0;
static uint64_t sequences = 0;
static uint64_t nucleotides = 0;
static uint64_t longest = 0;
static uint64_t shortest = 0;
static uint64_t longestheader = 0;

seqinfo_t * seqindex = 0;
char * datap = 0;

void db_setinfo(bool new_is_fastq,
                uint64_t new_sequences,
                uint64_t new_nucleotides,
                uint64_t new_longest,
                uint64_t new_shortest,
                uint64_t new_longestheader)
{
  is_fastq = new_is_fastq;
  sequences = new_sequences;
  nucleotides = new_nucleotides;
  longest = new_longest;
  shortest = new_shortest;
  longestheader = new_longestheader;
}

bool db_is_fastq()
{
  return is_fastq;
}

char * db_getquality(uint64_t seqno)
{
  if (is_fastq)
    return datap + seqindex[seqno].qual_p;
  else
    return 0;
}

void db_read(const char * filename, int upcase)
{
  h = fastx_open(filename);

  if (!h)
    fatal("Unrecognized file type (not proper FASTA or FASTQ format)");

  is_fastq = fastx_is_fastq(h);

  int64_t filesize = fastx_get_size(h);
  
  char * prompt = 0;
  if (xsprintf(& prompt, "Reading file %s", filename) == -1)
    fatal("Out of memory");

  progress_init(prompt, filesize);

  longest = 0;
  shortest = LONG_MAX;
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;

  int64_t discarded_short = 0;
  int64_t discarded_long = 0;
  int64_t discarded_unoise = 0;

  /* allocate space for data */
  uint64_t dataalloc = 0;
  datap = 0;
  uint64_t datalen = 0;

  /* allocate space for index */
  size_t seqindex_alloc = 0;
  seqindex = 0;

  while(fastx_next(h,
                   ! opt_notrunclabels,
                   upcase ? chrmap_upcase : chrmap_no_change))
    {
      size_t headerlength = fastx_get_header_length(h);
      size_t sequencelength = fastx_get_sequence_length(h);
      int64_t abundance = fastx_get_abundance(h);
      
      if (sequencelength < (size_t)opt_minseqlength)
        {
          discarded_short++;
        }
      else if (sequencelength > (size_t)opt_maxseqlength)
        {
          discarded_long++;
        }
      else if (opt_cluster_unoise && (abundance < (int64_t)opt_minsize))
        {
          discarded_unoise++;
        }
      else
        {
          /* grow space for data, if necessary */
          size_t dataalloc_old = dataalloc;
          size_t needed = datalen + headerlength + 1 + sequencelength + 1;
          if (is_fastq)
            needed += sequencelength + 1;
          while (dataalloc < needed)
            dataalloc += MEMCHUNK;
          if (dataalloc > dataalloc_old)
            datap = (char *) xrealloc(datap, dataalloc);

          /* store the header */
          size_t header_p = datalen;
          memcpy(datap + header_p,
                 fastx_get_header(h),
                 headerlength + 1);
          datalen += headerlength + 1;
          
          /* store sequence */
          size_t sequence_p = datalen;
          memcpy(datap + sequence_p,
                 fastx_get_sequence(h),
                 sequencelength + 1);
          datalen += sequencelength + 1;
          
          size_t quality_p = datalen;
          if (is_fastq)
            {
              /* store quality */
              memcpy(datap+quality_p,
                     fastx_get_quality(h),
                     sequencelength + 1);
              datalen += sequencelength + 1;
            }

          /* grow space for index, if necessary */
          size_t seqindex_alloc_old = seqindex_alloc;
          while ((sequences + 1) * sizeof(seqinfo_t) > seqindex_alloc)
            seqindex_alloc += MEMCHUNK;
          if (seqindex_alloc > seqindex_alloc_old)
            seqindex = (seqinfo_t *) xrealloc(seqindex, seqindex_alloc);
          
          /* update index */
          seqinfo_t * seqindex_p = seqindex + sequences;
          seqindex_p->headerlen = headerlength;
          seqindex_p->seqlen = sequencelength;
          seqindex_p->header_p = header_p;
          seqindex_p->seq_p = sequence_p;
          seqindex_p->qual_p = quality_p;
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
      progress_update(fastx_get_position(h));
    }

  progress_done();
  xfree(prompt);
  fastx_close(h);

  if (!opt_quiet)
    {
      if (sequences > 0)
        fprintf(stderr,
                "%'" PRIu64 " nt in %'" PRIu64 " seqs, min %'" PRIu64 ", max %'" PRIu64 ", avg %'.0f\n", 
                db_getnucleotidecount(),
                db_getsequencecount(),
                db_getshortestsequence(),
                db_getlongestsequence(),
                db_getnucleotidecount() * 1.0 / db_getsequencecount());
      else
        fprintf(stderr,
                "%'" PRIu64 " nt in %'" PRIu64 " seqs\n", 
                db_getnucleotidecount(),
                db_getsequencecount());
    }

  if (opt_log)
    {
      if (sequences > 0)
        fprintf(fp_log,
                "%'" PRIu64 " nt in %'" PRIu64 " seqs, min %'" PRIu64 ", max %'" PRIu64 ", avg %'.0f\n\n", 
                db_getnucleotidecount(),
                db_getsequencecount(),
                db_getshortestsequence(),
                db_getlongestsequence(),
                db_getnucleotidecount() * 1.0 / db_getsequencecount());
      else
        fprintf(fp_log,
                "%'" PRIu64 " nt in %'" PRIu64 " seqs\n\n", 
                db_getnucleotidecount(),
                db_getsequencecount());
    }

  /* Warn about discarded sequences */

  if (discarded_short)
    {
      fprintf(stderr,
              "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_minseqlength, discarded_short, (discarded_short == 1 ? "sequence" : "sequences"));

      if (opt_log)
        fprintf(fp_log,
              "minseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
              opt_minseqlength, discarded_short, (discarded_short == 1 ? "sequence" : "sequences"));
    }
  
  if (discarded_long)
    {
      fprintf(stderr,
              "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_maxseqlength, discarded_long, (discarded_long == 1 ? "sequence" : "sequences"));

      if (opt_log)
        fprintf(fp_log,
              "maxseqlength %" PRId64 ": %" PRId64 " %s discarded.\n\n",
              opt_maxseqlength, discarded_long, (discarded_long == 1 ? "sequence" : "sequences"));
    }

    if (discarded_unoise)
    {
      fprintf(stderr,
              "minsize %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_minsize, discarded_unoise, (discarded_unoise == 1 ? "sequence" : "sequences"));

      if (opt_log)
        fprintf(fp_log,
              "minsize %" PRId64 ": %" PRId64 " %s discarded.\n",
              opt_minsize, discarded_unoise, (discarded_unoise == 1 ? "sequence" : "sequences"));
    }
            
  show_rusage();
}

uint64_t db_getsequencecount()
{
  return sequences;
}

uint64_t db_getnucleotidecount()
{
  return nucleotides;
}

uint64_t db_getlongestheader()
{
  return longestheader;
}

uint64_t db_getlongestsequence()
{
  return longest;
}

uint64_t db_getshortestsequence()
{
  return shortest;
}

void db_free()
{
  if (datap)
    xfree(datap);
  if (seqindex)
    xfree(seqindex);
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


