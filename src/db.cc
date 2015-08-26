/*
    Copyright (C) 2014-2015 Torbjorn Rognes & Tomas Flouri

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

#define MEMCHUNK 10485760

static fasta_handle h;
static unsigned long sequences = 0;
static unsigned long nucleotides = 0;
static unsigned long longest;
static unsigned long shortest;
static unsigned long longestheader;
static char * datap;
seqinfo_t * seqindex;
regex_t db_regexp;

void db_read(const char * filename, int upcase)
{
  h = fasta_open(filename);

  long filesize = fasta_get_size(h);
  
  char * prompt;
  (void) asprintf(& prompt, "Reading file %s", filename);
  progress_init(prompt, filesize);

  /* allocate space */

  unsigned long dataalloc = MEMCHUNK;
  datap = (char *) xmalloc(dataalloc);
  unsigned long datalen = 0;
  seqindex = 0;
  longest = 0;
  shortest = LONG_MAX;
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;

  long discarded_short = 0;
  long discarded_long = 0;

  while(fasta_next(h,
                   ! opt_notrunclabels,
                   char_action_std,
                   upcase ? chrmap_upcase : chrmap_no_change))
    {
      size_t sequencelength = fasta_get_sequence_length(h);

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
          size_t headerlength = fasta_get_header_length(h);
          
          /*
            make space for:

            + header length
            + header
            + trailing zero
            + sequence length
            + sequence
            + trailing zero
            + abundance
          */
          
          size_t dataalloc_old = dataalloc;
          while (datalen + headerlength + sequencelength + 2
                 + 3 * sizeof(unsigned int)
                 > dataalloc)
            dataalloc += MEMCHUNK;

          if (dataalloc > dataalloc_old)
            datap = (char *) xrealloc(datap, dataalloc);
          
          /* store the header length */
          *(unsigned int*)(datap + datalen) = headerlength;
          datalen += 4;
          
          /* store the header */
          memcpy(datap + datalen, fasta_get_header(h), headerlength + 1);
          datalen += headerlength + 1;
          
          /* store sequence length */
          * (unsigned int*)(datap + datalen) = sequencelength;
          datalen += 4;
          
          /* store sequence */
          memcpy(datap + datalen, fasta_get_sequence(h), sequencelength + 1);
          datalen += sequencelength + 1;
          
          /* store abundance */
          * (unsigned int*)(datap + datalen) = fasta_get_abundance(h);
          datalen += 4;
          
          /* update statistics */
          nucleotides += sequencelength;
          sequences++;
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
                "%'ld nt in %'ld seqs, min %'ld, max %'ld, avg %'.0f\n", 
                db_getnucleotidecount(),
                db_getsequencecount(),
                db_getshortestsequence(),
                db_getlongestsequence(),
                db_getnucleotidecount() * 1.0 / db_getsequencecount());
      else
        fprintf(stderr,
                "%'ld nt in %'ld seqs\n", 
                db_getnucleotidecount(),
                db_getsequencecount());
    }

  if (opt_log)
    {
      if (sequences > 0)
        fprintf(fp_log,
                "%'ld nt in %'ld seqs, min %'ld, max %'ld, avg %'.0f\n\n", 
                db_getnucleotidecount(),
                db_getsequencecount(),
                db_getshortestsequence(),
                db_getlongestsequence(),
                db_getnucleotidecount() * 1.0 / db_getsequencecount());
      else
        fprintf(fp_log,
                "%'ld nt in %'ld seqs\n\n", 
                db_getnucleotidecount(),
                db_getsequencecount());
    }

  /* Warn about discarded sequences */

  if (discarded_short)
    {
      fprintf(stderr,
              "WARNING: %lu sequences shorter than %lu nucleotides discarded.\n",
              discarded_short, opt_minseqlength);

      if (opt_log)
        fprintf(fp_log,
                "WARNING: %lu sequences shorter than %lu nucleotides discarded.\n\n",
                discarded_short, opt_minseqlength);
    }
  
  if (discarded_long)
    {
      fprintf(stderr,
              "WARNING: %lu sequences longer than %lu nucleotides discarded.\n",
              discarded_long, opt_maxseqlength);

      if (opt_log)
        fprintf(fp_log,
                "WARNING: %lu sequences longer than %lu nucleotides discarded.\n\n",
                discarded_long, opt_maxseqlength);
    }

  show_rusage();

  progress_init("Indexing sequences", sequences);

  /* Create index and parse abundance info, if specified */

  seqindex = (seqinfo_t *) xmalloc(sequences * sizeof(seqinfo_t));
  seqinfo_t * seqindex_p = seqindex;

  char * p = datap;

  for(unsigned long i=0; i<sequences; i++)
  {
    seqindex_p->headerlen = * (unsigned int *) p;
    p += 4;

    seqindex_p->header = p;
    p += seqindex_p->headerlen + 1;

    seqindex_p->seqlen = *(unsigned int*) p;
    p += 4;

    seqindex_p->seq = p;
    p += seqindex_p->seqlen + 1;

    seqindex_p->size = *(unsigned int*)p;
    p += 4;

    seqindex_p++;
    progress_update(i);
  }

  progress_done();

  show_rusage();

  if (regcomp(& db_regexp, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
    fatal("Compilation of regular expression for abundance annotation failed");
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
  regfree(&db_regexp);
  if (datap)
    free(datap);
  if (seqindex)
    free(seqindex);
}

void db_fprint_fasta_seq_only(FILE * fp, unsigned long seqno)
{
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  fprint_fasta_seq_only(fp, seq, seqlen, opt_fasta_width);
}

void db_fprint_fasta(FILE * fp, unsigned long seqno)
{
  char * hdr = db_getheader(seqno);
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  
  fprint_fasta_hdr_only(fp, hdr);
  fprint_fasta_seq_only(fp, seq, seqlen, opt_fasta_width);
}

void db_fprint_fasta_with_size(FILE * fp, unsigned long seqno, unsigned long size)
{
  char * hdr = db_getheader(seqno);
  int hdrlen = db_getheaderlen(seqno);
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  
  /* remove any previous size annotation */
  /* regexp search for "(^|;)(\d+)(;|$)" */
  /* replace by ';' if not at either end */

  regmatch_t pmatch[1];

  if (!regexec(&db_regexp, hdr, 1, pmatch, 0))
    {
      int pat_start = pmatch[0].rm_so;
      int pat_end = pmatch[0].rm_eo;

      fprintf(fp,
              ">%.*s%s%.*s%ssize=%lu;\n",
              pat_start, hdr,
              pat_start > 0 ? ";" : "",
              hdrlen - pat_end, hdr + pat_end,
              ((pat_end < hdrlen) && (hdr[hdrlen - 1] != ';')) ? ";" : "",
              size);
    }
  else
    {
      fprintf(fp,
              ">%s%ssize=%lu;\n", 
              hdr,
              ((hdrlen == 0) || (hdr[hdrlen - 1] != ';')) ? ";" : "", 
              size);
    }

  fprint_fasta_seq_only(fp, seq, seqlen, opt_fasta_width);
}

void db_fprint_fasta_strip_size(FILE * fp, unsigned long seqno)
{
  /* write FASTA but remove abundance information, as with --xsize option */

  char * hdr = db_getheader(seqno);
  int hdrlen = db_getheaderlen(seqno);
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  
  /* remove any previous size annotation */
  /* regexp search for "(^|;)(\d+)(;|$)" */
  /* replace by ';' if not at either end */

  regmatch_t pmatch[1];

  if (!regexec(&db_regexp, hdr, 1, pmatch, 0))
    {
      int pat_start = pmatch[0].rm_so;
      int pat_end = pmatch[0].rm_eo;

      fprintf(fp,
              ">%.*s%s%.*s\n",
              pat_start, hdr,
              ((pat_start > 0) && (pat_end < hdrlen)) ? ";" : "",
              hdrlen - pat_end, hdr + pat_end);
    }
  else
    {
      fprintf(fp, ">%s\n", hdr);
    }

  fprint_fasta_seq_only(fp, seq, seqlen, opt_fasta_width);
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
          int r = strcmp(x->header, y->header);
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
          int r = strcmp(x->header, y->header);
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
      int r = strcmp(x->header, y->header);
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
