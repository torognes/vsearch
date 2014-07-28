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

long stripped[256];

static unsigned long sequences = 0;
static unsigned long nucleotides = 0;
static unsigned long headerchars = 0;
static int longest = 0;
static long shortest = LONG_MAX;
static int longestheader = 0;

seqinfo_t * seqindex = 0;
static char * datap = 0;

#define MEMCHUNK 1048576
#define LINEALLOC 1048576

regex_t db_regexp;

void db_read(const char * filename)
{
  if (regcomp(&db_regexp, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
    fatal("Regular expression compilation failed");
  
  FILE * fp = NULL;
  if (filename)
    {
      fp = fopen(filename, "r");
      if (!fp)
        fatal("Error: Unable to open database file (%s)", filename);
    }
  else
    fp = stdin;


  /* get file size */

  if (fseek(fp, 0, SEEK_END))
    fatal("Error: Unable to seek in database file (%s)", filename);

  long filesize = ftell(fp);
  
  rewind(fp);

  progress_init("Reading database file", filesize);

  /* allocate space */

  unsigned long dataalloc = MEMCHUNK;
  datap = (char *) xmalloc(dataalloc);
  unsigned long datalen = 0;

  longest = 0;
  shortest = LONG_MAX;
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;
  headerchars = 0;

  char line[LINEALLOC];
  line[0] = 0;
  fgets(line, LINEALLOC, fp);

  long lineno = 1;

  long stripped_count = 0;
  for(int i=0; i<256; i++)
    stripped[i] = 0;

  long discarded_short = 0;
  long discarded_long = 0;

  while(line[0])
    {
      /* read header */

      unsigned long hdrbegin = datalen;

      if (line[0] != '>')
        fatal("Illegal header line in fasta file.");
      
      /* terminate header at first space or end of line */

      char * z0 = line + 1;
      char * z = z0;
      while (*z)
	{
	  if ((!opt_notrunclabels) && (*z == ' '))
	    break;
	  if (*z == '\n')
	    break;
	  z++;
	}
      long headerlen = z - z0;

      headerchars += headerlen;

      if (headerlen > longestheader)
        longestheader = headerlen;

      /* store the header */

      while (datalen + headerlen + 1 + 4 > dataalloc)
      {
        dataalloc += MEMCHUNK;
        datap = (char *) xrealloc(datap, dataalloc);
      }

      *(unsigned int*)(datap + datalen) = headerlen;
      datalen += 4;

      memcpy(datap + datalen, line + 1, headerlen);
      *(datap + datalen + headerlen) = 0;
      datalen += headerlen + 1;

      /* get next line */

      line[0] = 0;
      fgets(line, LINEALLOC, fp);
      lineno++;

      /* read sequence */

      unsigned long seqbegin = datalen;

      /* make space for sequence length to be filled in later */

      while (datalen + 4 > dataalloc)
      {
        dataalloc += MEMCHUNK;
        datap = (char *) xrealloc(datap, dataalloc);
      }

      * (unsigned int*)(datap + datalen) = 0;
      datalen += 4;

      while (line[0] && (line[0] != '>'))
        {
          char m;
          char c;
          char * p = line;
          while((c = *p++))
	    {
	      m = chrstatus[(int)c];

	      switch(m)
		{
		case 0:
		  /* character to be stripped */
		  stripped_count++;
		  stripped[(int)c]++;
		  break;

		case 1:
		  /* legal character */
		  while (datalen >= dataalloc)
		    {
		      dataalloc += MEMCHUNK;
		      datap = (char *) xrealloc(datap, dataalloc);
		    }
		  *(datap+datalen) = c;
		  datalen++;
		  break;

		case 2:
		  /* fatal character */
		  char msg[200];
		  if (c>=32)
		    snprintf(msg, 200, "illegal character '%c' on line %ld in the database file", c, lineno);
		  else
		    snprintf(msg, 200, "illegal unprintable character %#.2x (hexadecimal) on line %ld in the database file", c, lineno);
		  fatal(msg);
		  break;

		case 3:
		  /* character to be stripped, silently */
		  break;
		}
	    }

          line[0] = 0;
          fgets(line, LINEALLOC, fp);
	  lineno++;
        }
      

      long length = datalen - seqbegin - 4;
      

      /* discard sequence if too short or long */

      if (length < opt_minseqlength)
	{
	  discarded_short++;
	  datalen = hdrbegin;
	}
      else if (length > opt_maxseqlength)
	{
	  discarded_long++;
	  datalen = hdrbegin;
	}
      else
	{
	  /* store the length in its designated space */
	  
	  *(unsigned int*)(datap + seqbegin) = length;
	  
	  
	  /* add a zero after the sequence */
	  
	  while (datalen >= dataalloc)
	    {
	      dataalloc += MEMCHUNK;
	      datap = (char *) xrealloc(datap, dataalloc);
	    }
	  *(datap+datalen) = 0;
	  datalen++;
	  

	  /* update statistics */

	  nucleotides += length;

	  if (length > longest)
	    longest = length;
	  
	  if (length < shortest)
	    shortest = length;

	  sequences++;
	}

      progress_update(ftell(fp));
    }

  fclose(fp);

  progress_done();

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

  /* Warn about stripped chars */

  if (stripped_count)
    {
      fprintf(stderr, "Warning: invalid characters stripped from sequence:");
      for (int i=0; i<256;i++)
	if (stripped[i])
	  fprintf(stderr, " %c(%ld)", i, stripped[i]);
      fprintf(stderr, "\n");
    }


  /* Warn about discarded sequences */

  if (discarded_short)
    fprintf(stderr,
	    "Warning: %lu sequences shorter than %lu nucleotides discarded.\n",
	    discarded_short, opt_minseqlength);
  
  if (discarded_long)
    fprintf(stderr,
	    "Warning: %lu sequences longer than %lu nucleotides discarded.\n",
	    discarded_long, opt_maxseqlength);

  show_rusage();
  

  progress_init("Indexing sequences", sequences);


  /* Create index and parse abundance info, if specified */

  seqindex = (seqinfo_t *) xmalloc(sequences * sizeof(seqinfo_t));
  seqinfo_t * seqindex_p = seqindex;

  char * p = datap;

  regmatch_t pmatch[4];
  
  for(unsigned long i=0; i<sequences; i++)
  {
    seqindex_p->headerlen = * (unsigned int *) p;
    p += 4;

    seqindex_p->header = p;
    p += seqindex_p->headerlen + 1;

    seqindex_p->headeridlen = xstrchrnul(seqindex_p->header, ' ') 
      - seqindex_p->header;

    seqindex_p->seqlen = *(unsigned int*) p;
    p += 4;

    seqindex_p->seq = p;
    p += seqindex_p->seqlen + 1;

    seqindex_p->size = 1;

    if (opt_sizein && !regexec(&db_regexp, seqindex_p->header, 4, pmatch, 0))
      {
	unsigned long size = atol(seqindex_p->header + pmatch[2].rm_so);
	if (size > 0)
	  seqindex_p->size = size;
	else
	  fatal("size annotation zero");
      }

    seqindex_p++;
    progress_update(i);
  }

  progress_done();
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
