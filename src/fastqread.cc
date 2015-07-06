/*
  Copyright (C) 2015 Torbjorn Rognes

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

typedef struct fastq_entry
{
  char * header;
  char * sequence;
  char * dummy;
  char * quality;
  
  size_t header_length;
  size_t sequence_length;
  size_t dummy_length;
  size_t quality_length;

  size_t header_alloc;
  size_t sequence_alloc;
  size_t dummy_alloc;
  size_t quality_alloc;

  size_t sequence_chars[256];
  size_t quality_chars[256];
} fastq_entry_t;

void fastq_read_init(fastq_entry_t * fqe)
{
  fqe->header_length = 0;
  fqe->sequence_length = 0;
  fqe->dummy_length = 0;
  fqe->quality_length = 0;
  fqe->header_alloc = 0;
  fqe->sequence_alloc = 0;
  fqe->dummy_alloc = 0;
  fqe->quality_alloc = 0;
  fqe->header = 0;
  fqe->sequence = 0;
  fqe->dummy = 0;
  fqe->quality = 0;

  for(int c=0; c<256; c++)
    {
      fqe->sequence_chars[c] = 0;
      fqe->quality_chars[c] = 0;
    }
}

void fastq_read_one(FILE * fp, fastq_entry_t * fqe)
{
  /* init entry */
  fqe->header_length = 0;
  fqe->sequence_length = 0;
  fqe->header_length = 0;

  /* read four lines from the fastq file */
  ssize_t len;
  len = getline(& fqe->header, & fqe->header_alloc, fp);
  if (len < 1)
    fatal("Invalid FASTQ file header line: too short");
  if (fqe->header[0] != '@')
    fatal("Invalid FASTQ file header line: does not start with '@'");
  if (fqe->header[len-1] != 0x0a)
    fatal("Invalid FASTQ file header line: does not end with newline");
  fqe->header_length = len;

  len = getline(& fqe->sequence, & fqe->sequence_alloc, fp);
  if (len < 1)
    fatal("Invalid FASTQ file sequence line");
  if (fqe->sequence[len-1] != 0x0a)
    fatal("Invalid FASTQ file sequence line: does not end with newline");
  fqe->sequence_length = len;

  for(size_t i=0; i<fqe->sequence_length-1; i++)
    fqe->sequence_chars[(int)(fqe->sequence[i])]++;

  len = getline(& fqe->dummy, & fqe->dummy_alloc, fp);
  if (len < 1)
    fatal("Invalid FASTQ file header line");
  if (fqe->dummy[0] != '+')
    fatal("Invalid FASTQ file third line: does not start with '+'");
  if (fqe->dummy[len-1] != 0x0a)
    fatal("Invalid FASTQ file third line: does not end with newline");
  
  len = getline(& fqe->quality, & fqe->quality_alloc, fp);
  if (len < 1)
    fatal("Invalid FASTQ file quality line");
  if (fqe->quality[len-1] != 0x0a)
    fatal("Invalid FASTQ file quality line: does not end with newline");
  fqe->quality_length = len;

  for(size_t i=0; i<fqe->quality_length-1; i++)
    fqe->quality_chars[(int)(fqe->quality[i])]++;

  /* show lines */
  if (fp_log)
    {
      fprintf(fp_log, "\n");
      fprintf(fp_log, "Header:   %.*s\n", (int)(fqe->header_length-2),
              fqe->header+1);
      fprintf(fp_log, "Sequence: %.*s\n", (int)(fqe->sequence_length-1),
              fqe->sequence);
      fprintf(fp_log, "Quality:  %.*s\n", (int)(fqe->quality_length-1),
              fqe->quality);
    }
}

void fastq_chars()
{
  FILE * fp_fastq = 0;
  
  fp_fastq = fopen(opt_fastq_chars, "r");
  if (! fp_fastq)
    fatal("Unable to open FASTQ input file (%s)", opt_fastq_chars);
  
  if (fseek(fp_fastq, 0, SEEK_END))
    fatal("Error: Unable to seek in FASTQ input file (%s)", opt_fastq_chars);

  size_t fastq_filesize = ftell(fp_fastq);
  rewind(fp_fastq);

  long seq_count = 0;

  fastq_entry_t fqe;
  fastq_read_init(&fqe);

  size_t fastq_filepos = ftell(fp_fastq);

  progress_init("Reading FASTQ file", fastq_filesize);

  while (fastq_filepos < fastq_filesize)
    {
      fastq_read_one(fp_fastq, &fqe);
      seq_count++;
      fastq_filepos = ftell(fp_fastq);
      progress_update(fastq_filepos);
    }

  progress_done();

  fprintf(stderr, "Read %ld sequences.\n", seq_count);

  char qmin = 0;
  char qmax = 0;

  for(int c=0; c<256; c++)
    {
      if (fqe.quality_chars[c] > 0)
        {
          qmin = c;
          break;
        }
    }

  for(int c=255; c>=0; c--)
    {
      if (fqe.quality_chars[c] > 0)
        {
          qmax = c;
          break;
        }
    }

  size_t total_chars = 0;
  for(int c=255; c>=0; c--)
    total_chars += fqe.sequence_chars[c];

  char fastq_ascii, fastq_qmin, fastq_qmax;

  if ((qmin >= 59) && (qmax > 74))
    fastq_ascii = 64;
  else
    fastq_ascii = 33;
  fastq_qmax = qmax - fastq_ascii;
  fastq_qmin = qmin - fastq_ascii;

  fprintf(stderr, "Qmin %d, QMax %d, Range %d\n",
          qmin, qmax, qmax-qmin+1);
  fprintf(stderr, "Guess: -fastq_qmin %d -fastq_qmax %d -fastq_ascii %d\n",
          fastq_qmin, fastq_qmax, fastq_ascii);
  if (fastq_ascii == 64)
    {
      if (qmin < 64)
        fprintf(stderr, "Guess: Solexa format\n");
      else if (qmin < 66)
        fprintf(stderr, "Guess: Illumina 1.3+ format\n");
      else
        fprintf(stderr, "Guess: Illumina 1.5+ format\n");
    }
  else
    {
      if (qmax == 74)
        fprintf(stderr, "Guess: Illumina 1.8+ format\n");
      else
        fprintf(stderr, "Guess: Sanger format\n");
    }
  fprintf(stderr, "\n");
  fprintf(stderr, "Letter          N   Freq MaxRun\n");
  fprintf(stderr, "------ ---------- ------ ------\n");
  
  for(int c=0; c<256; c++)
    {
      if (fqe.sequence_chars[c] > 0)
        {
          fprintf(stderr, "     %c %10ld %5.1f%% %6d\n",
                  c,
                  fqe.sequence_chars[c],
                  100.0 * fqe.sequence_chars[c] / total_chars,
                  0);
        }
    }
}
