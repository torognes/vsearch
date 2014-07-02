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

/*

char map_nt[256] =
  {
    // N = A

    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,  0, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,  0, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

*/

char map_nt[256] =
  {
    // ARWMDHVN = 0, CYB = 1, GS = 2, TU = 3

    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0,  1,  1,  0, -1, -1,  2,  0, -1, -1,  0, -1,  0,  0, -1,
    -1, -1,  0,  2,  3,  3,  0,  0, -1,  1, -1, -1, -1, -1, -1, -1,
    -1,  0,  1,  1,  0, -1, -1,  2,  0, -1, -1,  0, -1,  0,  0, -1,
    -1, -1,  0,  2,  3,  3,  0,  0, -1,  1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

char sym_nt[] = "acgt";

unsigned long sequences = 0;
unsigned long nucleotides = 0;
unsigned long headerchars = 0;
int longest = 0;
long shortest = LONG_MAX;
int longestheader = 0;

seqinfo_t * seqindex = 0;
char * datap = 0;

#define MEMCHUNK 1048576
#define LINEALLOC 1048576

void showseq(char * seq)
{
  char * p = seq;
  while (char c = *p++)
  {
    putchar(sym_nt[(unsigned int)c]);
  }
}


void db_read(const char * filename)
{
  show_rusage();
  fprintf(stderr, "Reading database: ");

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

  FILE * fp = NULL;
  if (filename)
    {
      fp = fopen(filename, "r");
      if (!fp)
        fatal("Error: Unable to open database file (%s)", filename);
    }
  else
    fp = stdin;

  char line[LINEALLOC];
  line[0] = 0;
  fgets(line, LINEALLOC, fp);

  while(line[0])
    {
      /* read header */

      if (line[0] != '>')
        fatal("Illegal header line in fasta file.");
      
      /* terminate header at first space or end of line */

      char * z0 = line + 1;
      char * z = z0;
      while (1)
	{
	  if ((*z == ' ') || (*z == '\n') || (!*z))
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
            if ((m = map_nt[(int)c]) >= 0)
            {
              while (datalen >= dataalloc)
              {
                dataalloc += MEMCHUNK;
                datap = (char *) xrealloc(datap, dataalloc);
              }

              *(datap+datalen) = m;
              datalen++;
            }
            else if (c != '\n')
              fatal("Illegal character in sequence.");
          line[0] = 0;
          fgets(line, LINEALLOC, fp);
        }
      
      long length = datalen - seqbegin - 4;

      /* store the length in designated space */

      *(unsigned int*)(datap + seqbegin) = length;

      nucleotides += length;

      if (length > longest)
        longest = length;

      if (length < shortest)
        shortest = length;

      /*
      *(datap+datalen) = 0;
      datalen++;
      */
      
      sequences++;
    }

  fclose(fp);

  /* create indices */

  seqindex = (seqinfo_t *) xmalloc(sequences * sizeof(seqinfo_t));
  seqinfo_t * seqindex_p = seqindex;

  char * p = datap;
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
    p += seqindex_p->seqlen;

    seqindex_p++;
  }

  fprintf(stderr,
          "%'ld nt in %'ld seqs, from %'ld to %'ld nt (avg %'.0f)\n", 
          db_getnucleotidecount(),
          db_getsequencecount(),
          db_getshortestsequence(),
          db_getlongestsequence(),
          db_getnucleotidecount() * 1.0 / db_getsequencecount());

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

seqinfo_t * db_getseqinfo(unsigned long seqno)
{
  return seqindex+seqno;
}

char * db_getsequence(unsigned long seqno)
{
  return seqindex[seqno].seq;
}

void db_getsequenceandlength(unsigned long seqno,
                             char ** address,
                             long * length)
{
  *address = seqindex[seqno].seq;
  *length = (long)(seqindex[seqno].seqlen);
}

unsigned long db_getsequencelen(unsigned long seqno)
{
  return seqindex[seqno].seqlen;
}

char * db_getheader(unsigned long seqno)
{
  return seqindex[seqno].header;
}

unsigned long db_getheaderlen(unsigned long seqno)
{
  return seqindex[seqno].headerlen;
}

void db_putseq(long seqno)
{
  char * seq;
  long len;
  db_getsequenceandlength(seqno, & seq, & len);
  for(int i=0; i<len; i++)
    putchar(sym_nt[(int)(seq[i])]);
}

void db_free()
{
  if (datap)
    free(datap);
  if (seqindex)
    free(seqindex);
}
