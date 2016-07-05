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

fastx_handle fasta_open(const char * filename)
{
  fastx_handle h = fastx_open(filename);

  if (fastx_is_fastq(h))
    fatal("FASTA file expected, FASTQ file found (%s)", filename);

  return h;
}

void fasta_close(fastx_handle h)
{
  fastx_close(h);
}

void fasta_truncate_header(fastx_handle h, bool truncateatspace)
{
  /* Truncate the zero-terminated header string by inserting a new
     terminator (zero byte) at the first space/tab character
     (if truncateatspace) or first linefeed. */
  
  if (truncateatspace)
    h->header_buffer.length = strcspn(h->header_buffer.data, " \t\n");
  else
    h->header_buffer.length = strcspn(h->header_buffer.data, "\n");
  
  h->header_buffer.data[h->header_buffer.length] = 0;
}

void fasta_filter_sequence(fastx_handle h,
                           unsigned int * char_action,
                           const unsigned char * char_mapping)
{
  /* Strip unwanted characters from the sequence and raise warnings or
     errors on certain characters. */

  char * p = h->sequence_buffer.data;
  char * q = p;
  char c;
  char msg[200];

  while ((c = *p++))
    {
      char m = char_action[(int)c];

      switch(m)
        {
        case 0:
          /* stripped */
          h->stripped_all++;
          h->stripped[(int)c]++;
          break;
              
        case 1:
          /* legal character */
          *q++ = char_mapping[(int)(c)];
          break;
          
        case 2:
          /* fatal character */
          if ((c>=32) && (c<127))
            snprintf(msg,
                     200,
                     "illegal character '%c' on line %lu in fasta file",
                     c,
                     h->lineno);
          else
            snprintf(msg,
                     200,
                     "illegal unprintable character %#.2x (hexadecimal) on line %lu in fasta file",
                     c,
                     h->lineno);
          fatal(msg);
          break;
              
        case 3:
          /* silently stripped chars (whitespace) */
          break;

        case 4:
          /* newline (silently stripped) */
          h->lineno++;
          break;
        }
    }

  /* add zero after sequence */
  *q = 0;
  h->sequence_buffer.length = q - h->sequence_buffer.data;
}

bool fasta_next(fastx_handle h,
                bool truncateatspace,
                const unsigned char * char_mapping)
{
  h->lineno_start = h->lineno;

  h->header_buffer.length = 0;
  h->sequence_buffer.length = 0;

  unsigned long rest = fastx_file_fill_buffer(h);

  if (rest == 0)
    return 0;

  /* read header */

  /* check initial > character */
  
  if (h->file_buffer.data[h->file_buffer.position] != '>')
    {
      fprintf(stderr, "Found character %02x\n", (unsigned char)(h->file_buffer.data[h->file_buffer.position]));
      fatal("Invalid FASTA - header must start with > character");
    }
  h->file_buffer.position++;
  rest--;

  char * lf = 0;
  while (lf == 0)
    {
      /* get more data if buffer empty*/
      rest = fastx_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTA - header must be terminated with newline");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);

      /* copy to header buffer */
      unsigned long len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }
      buffer_extend(& h->header_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  /* read one or more sequence lines */

  while (1)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(h);

      /* end if no more data */
      if (rest == 0)
        break;

      /* end if new sequence starts */
      if (lf && (h->file_buffer.data[h->file_buffer.position] == '>'))
        break;

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);

      unsigned long len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
        }
      buffer_extend(& h->sequence_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  h->seqno++;

  fasta_truncate_header(h, truncateatspace);
  fasta_filter_sequence(h, char_fasta_action, char_mapping);

  return 1;
}

long fasta_get_abundance(fastx_handle h)
{
  return abundance_get(global_abundance, h->header_buffer.data);
}

unsigned long fasta_get_position(fastx_handle h)
{
  return h->file_position;
}

unsigned long fasta_get_size(fastx_handle h)
{
  return h->file_size;
}

unsigned long fasta_get_lineno(fastx_handle h)
{
  return h->lineno_start;
}

unsigned long fasta_get_seqno(fastx_handle h)
{
  return h->seqno;
}

unsigned long fasta_get_header_length(fastx_handle h)
{
  return h->header_buffer.length;
}

unsigned long fasta_get_sequence_length(fastx_handle h)
{
  return h->sequence_buffer.length;
}

char * fasta_get_header(fastx_handle h)
{
  return h->header_buffer.data;
}

char * fasta_get_sequence(fastx_handle h)
{
  return h->sequence_buffer.data;
}


/* fasta output */

void fasta_print_header(FILE * fp, const char * hdr)
{
  fprintf(fp, ">%s\n", hdr);
}

void fasta_print_sequence(FILE * fp, char * seq, unsigned long len, int width)
{
  /*
     The actual length of the sequence may be longer than "len", but only
     "len" characters are printed.

     Specify width of lines - zero (or <1)  means linearize (all on one line).
  */

  if (width < 1)
    fprintf(fp, "%.*s\n", (int)(len), seq);
  else
    {
      long rest = len;
      for(unsigned long i=0; i<len; i += width)
        {
          fprintf(fp, "%.*s\n", (int)(MIN(rest,width)), seq+i);
          rest -= width;
        }
    }
}

void fasta_print(FILE * fp, const char * hdr,
                 char * seq, unsigned long len)
{
  fasta_print_header(fp, hdr);
  fasta_print_sequence(fp, seq, len, opt_fasta_width);
}

void fasta_print_ee(FILE * fp, const char * hdr,
                    char * seq, unsigned long len,
                    double ee)
{
  fprintf(fp, ">%s;%6.4lf;\n", hdr, ee);
  fasta_print_sequence(fp, seq, len, opt_fasta_width);
}

void fasta_print_general(FILE * fp,
                         char * seq,
                         int len,
                         char * header,
                         int header_len,
                         int abundance,
                         int ordinal,
                         int clustersize,
                         bool showclusterid,
                         int clusterid,
                         double ee)
{
  fprintf(fp, ">");

  if (clustersize > 0)
    fprintf(fp, "centroid=");

  if (opt_relabel || opt_relabel_sha1 || opt_relabel_md5)
    {
      if (opt_relabel_sha1)
        fprint_seq_digest_sha1(fp, seq, len);
      else if (opt_relabel_md5)
        fprint_seq_digest_md5(fp, seq, len);
      else
        fprintf(fp, "%s%d", opt_relabel, ordinal);
    }
  else
    {
      if (opt_sizeout || opt_xsize)
        abundance_fprint_header_strip_size(global_abundance,
                                           fp,
                                           header,
                                           header_len);
      else
        fprintf(fp, "%s", header);
    }

  if ((clustersize > 0) || opt_sizeout || showclusterid || (ee >= 0.0))
    fprintf(fp, ";");

  if (clustersize > 0)
    fprintf(fp, "seqs=%d;", clustersize);

  if (opt_sizeout)
    fprintf(fp, "size=%u;", abundance);

  if (showclusterid)
    fprintf(fp, "clusterid=%d;", clusterid);

  if (ee >= 0.0)
    fprintf(fp, "ee=%6.4lf;", ee);

  if (opt_relabel_keep &&
      ((opt_relabel || opt_relabel_sha1 || opt_relabel_md5)))
    fprintf(fp, " %s", header);

  fprintf(fp, "\n");

  if (seq)
    fasta_print_sequence(fp, seq, len, opt_fasta_width);
}


void fasta_print_relabel_cluster(FILE * fp,
                                 char * seq,
                                 int len,
                                 char * header,
                                 int header_len,
                                 int abundance,
                                 int ordinal,
                                 int clustersize,
                                 bool showclusterid,
                                 int clusterid)
{
  fasta_print_general(fp, seq, len, header, header_len,
                      abundance, ordinal,
                      clustersize, showclusterid, clusterid,
                      -1.0);
}

void fasta_print_relabel_ee(FILE * fp,
                            char * seq,
                            int len,
                            char * header,
                            int header_len,
                            int abundance,
                            int ordinal,
                            double ee)
{
  fasta_print_general(fp, seq, len, header, header_len,
                      abundance, ordinal,
                      0, 0, 0,
                      ee);
}

void fasta_print_relabel(FILE * fp,
                         char * seq,
                         int len,
                         char * header,
                         int header_len,
                         int abundance,
                         int ordinal)
{
  fasta_print_general(fp, seq, len, header, header_len,
                      abundance, ordinal,
                      0, 0, 0,
                      -1.0);
}

void fasta_print_db_relabel(FILE * fp,
                            unsigned long seqno,
                            int ordinal)
{
  fasta_print_relabel(fp,
                      db_getsequence(seqno),
                      db_getsequencelen(seqno),
                      db_getheader(seqno),
                      db_getheaderlen(seqno),
                      db_getabundance(seqno),
                      ordinal);
}

void fasta_print_db_sequence(FILE * fp, unsigned long seqno)
{
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

void fasta_print_db(FILE * fp, unsigned long seqno)
{
  char * hdr = db_getheader(seqno);
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);

  fasta_print_header(fp, hdr);
  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

void fasta_print_db_size(FILE * fp, unsigned long seqno, unsigned long size)
{
  char * hdr = db_getheader(seqno);
  int hdrlen = db_getheaderlen(seqno);
  
  fprintf(fp, ">");
  abundance_fprint_header_with_size(global_abundance,
                                    fp,
                                    hdr,
                                    hdrlen,
                                    size);
  fprintf(fp, "\n");

  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);

  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

void fasta_print_db_strip_size(FILE * fp, unsigned long seqno)
{
  /* write FASTA but remove abundance information, as with --xsize option */

  char * hdr = db_getheader(seqno);
  int hdrlen = db_getheaderlen(seqno);

  fprintf(fp, ">");
  abundance_fprint_header_strip_size(global_abundance,
                                     fp,
                                     hdr,
                                     hdrlen);
  fprintf(fp, "\n");

  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);

  fasta_print_sequence(fp, seq, seqlen, opt_fasta_width);
}

