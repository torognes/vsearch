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

#define BLOCKSIZE (4096 * 4096)

static unsigned int udb_dbaccel = 0;

typedef struct wordfreq
{
  unsigned int kmer;
  unsigned int count;
} wordfreq_t;

int wc_compare(const void * a, const void * b)
{
  wordfreq_t * x = (wordfreq_t *) a;
  wordfreq_t * y = (wordfreq_t *) b;
  if (x->count < y->count)
    return -1;
  else if (x->count > y->count)
    return +1;
  else
    {
      if (x->kmer < y->kmer)
        return +1;
      else if (x->kmer > y->kmer)
        return -1;
      else
        return 0;
    }
}

uint64_t largeread(int fd, void * buf, uint64_t nbyte, uint64_t offset)
{
  /* call pread multiple times and update progress */

  uint64_t progress = offset;
  for(uint64_t i = 0; i < nbyte; i += BLOCKSIZE)
    {
      uint64_t res = xlseek(fd, offset + i, SEEK_SET);
      if (res != offset + i)
        fatal("Unable to seek in UDB file or invalid UDB file");

      uint64 rem = MIN(BLOCKSIZE, nbyte - i);
      uint64_t bytesread = read(fd, ((char*)buf) + i, rem);
      if (bytesread != rem)
        fatal("Unable to read from UDB file or invalid UDB file");

      progress += rem;
      progress_update(progress);
    }
  return nbyte;
}

uint64_t largewrite(int fd, void * buf, uint64_t nbyte, uint64_t offset)
{
  /* call write multiple times and update progress */

  uint64_t progress = offset;
  for(uint64_t i = 0; i < nbyte; i += BLOCKSIZE)
    {
      uint64_t res = xlseek(fd, offset + i, SEEK_SET);
      if (res != offset + i)
        fatal("Unable to seek in UDB file or invalid UDB file");

      uint64 rem = MIN(BLOCKSIZE, nbyte - i);
      uint64_t byteswritten = write(fd, ((char*)buf) + i, rem);
      if (byteswritten != rem)
        fatal("Unable to write to UDB file");

      progress += rem;
      progress_update(progress);
    }
  return nbyte;
}

bool udb_detect_isudb(const char * filename)
{
  /*
    Detect whether the given filename seems to refer to an UDB file.
    It must be an uncompressed regular file, not a pipe.
  */

  xstat_t fs;

  if (xstat(filename, & fs))
    fatal("Unable to get status for input file (%s)", filename);

  bool is_pipe = S_ISFIFO(fs.st_mode);
  if (is_pipe)
    return 0;

  int fd = 0;
  fd = xopen_read(filename);
  if (!fd)
    fatal("Unable to open input file for reading (%s)", filename);

  unsigned int magic = 0;
  uint64_t bytesread = read(fd, & magic, 4);
  close(fd);

  if ((bytesread == 4) && (magic == 0x55444246))
    return 1;

  return 0;
}

void udb_info()
{
  /* Read UDB header and show basic info */

  unsigned int buffer[50];

  int fd_udbinfo = 0;

  fd_udbinfo = xopen_read(opt_udbinfo);
  if (! fd_udbinfo)
    fatal("Unable to open UDB file for reading");

  uint64_t bytesread = read(fd_udbinfo, buffer, 4 * 50);
  if (bytesread != 4 * 50)
    fatal("Unable to read from UDB file or invalid UDB file");

  if ((buffer[0]  != 0x55444246) ||
      (buffer[2] != 32) ||
      (buffer[4] < 3) ||
      (buffer[4] > 15) ||
      (buffer[13] == 0) ||
      (buffer[17] != 0x0000746e) ||
      (buffer[49] != 0x55444266))
    fatal("Invalid UDB file");

  if (!opt_quiet)
    {
      fprintf(stderr, "           Seqs  %u\n", buffer[13]);
      fprintf(stderr, "     SeqIx bits  %u\n", buffer[2]);
      fprintf(stderr, "          Alpha  nt (4)\n");
      fprintf(stderr, "     Word width  %u\n", buffer[4]);
      fprintf(stderr, "          Slots  %u\n", buffer[11]);
      fprintf(stderr, "      Dict size  %u (%.1fk)\n",
              (1 << (2 * buffer[4])),
              (1 << (2 * buffer[4])) * 1.0 / 1000.0);
      fprintf(stderr, "         DBstep  %u\n", buffer[5]);
      fprintf(stderr, "        DBAccel  %u%%\n", buffer[6]);
    }

  if (opt_log)
    {
      fprintf(fp_log, "           Seqs  %u\n", buffer[13]);
      fprintf(fp_log, "     SeqIx bits  %u\n", buffer[2]);
      fprintf(fp_log, "          Alpha  nt (4)\n");
      fprintf(fp_log, "     Word width  %u\n", buffer[4]);
      fprintf(fp_log, "          Slots  %u\n", buffer[11]);
      fprintf(fp_log, "      Dict size  %u (%.1fk)\n",
              (1 << (2 * buffer[4])),
              (1 << (2 * buffer[4])) * 1.0 / 1000.0);
      fprintf(fp_log, "         DBstep  %u\n", buffer[5]);
      fprintf(fp_log, "        DBAccel  %u%%\n", buffer[6]);
    }

  close(fd_udbinfo);
}

void udb_read(const char * filename,
              bool create_bitmaps,
              bool parse_abundances)
{
  /* read UDB as indexed database */

  unsigned int seqcount = 0;
  unsigned int udb_wordlength = 0;
  uint64 nucleotides = 0;

  xstat_t fs;
  if (xstat(filename, & fs))
    fatal("Unable to get status for input file (%s)", filename);

  bool is_pipe = S_ISFIFO(fs.st_mode);
  if (is_pipe)
    fatal("Cannot read UDB file from a pipe");

  /* get file size */

  uint64_t filesize = fs.st_size;

  /* open UDB file */

  int fd_udb = 0;

  fd_udb = xopen_read(filename);
  if (! fd_udb)
    fatal("Unable to open UDB file for reading");

  char * prompt = 0;
  if (xsprintf(& prompt, "Reading UDB file %s", filename) == -1)
    fatal("Out of memory");

  progress_init(prompt, filesize);

  /* header */

  unsigned int buffer[50];
  uint64_t pos = 0;

  pos += largeread(fd_udb, buffer, 4 * 50, pos);

  if ((buffer[0]  != 0x55444246) ||
      (buffer[2] != 32) ||
      (buffer[4] < 3) ||
      (buffer[4] > 15) ||
      (buffer[13] == 0) ||
      (buffer[17] != 0x0000746e) ||
      (buffer[49] != 0x55444266))
    fatal("Invalid UDB file");

  udb_wordlength = buffer[4];
  seqcount = buffer[13];
  udb_dbaccel = buffer[6];

  if (udb_wordlength != opt_wordlength)
    {
      fprintf(stderr, "\nWARNING: Wordlength adjusted to %u as indicated in UDB file\n", udb_wordlength);
      opt_wordlength = udb_wordlength;
    }

  /* word match counts */

  kmerhashsize = 1 << (2 * udb_wordlength);
  kmercount = (unsigned int*) xmalloc(kmerhashsize * sizeof(unsigned int));
  kmerhash = (uint64_t *) xmalloc(kmerhashsize * sizeof(uint64_t));
  kmerbitmap = (bitmap_t * *) xmalloc(kmerhashsize * sizeof(bitmap_t**));

  memset(kmerbitmap, 0, kmerhashsize * sizeof(bitmap_t**));

  pos += largeread(fd_udb, kmercount, 4 * kmerhashsize, pos);

  kmerindexsize = 0;
  for(uint64_t i = 0; i < kmerhashsize; i++)
    {
      kmerhash[i] = kmerindexsize;
      kmerindexsize += kmercount[i];
    }

  /* signature */

  pos += largeread(fd_udb, buffer, 4, pos);

  if (buffer[0] != 0x55444233)
    fatal("Invalid UDB file");

  /* sequence numbers for word matches */

  kmerindex = (unsigned int *) xmalloc(kmerindexsize * 4);

  pos += largeread(fd_udb, kmerindex, 4 * kmerindexsize, pos);

  /* new header */

  pos += largeread(fd_udb, buffer, 4 * 8, pos);

  if ((buffer[0] != 0x55444234) ||
      (buffer[1] != 0x005e0db3) ||
      (buffer[2] != seqcount) ||
      (buffer[7] != 0x005e0db4))
    fatal("Invalid UDB file");

  nucleotides = (((uint64_t) buffer[4]) << 32) | buffer[3];
  uint64_t udb_headerchars = (((uint64_t) buffer[6]) << 32) | buffer[5];

  /* header index */

  seqindex = (seqinfo_t *) xmalloc(seqcount * sizeof(seqinfo_t));

  int * header_index = (int *) xmalloc(4 * (seqcount+1));

  pos += largeread(fd_udb, header_index, 4 * seqcount, pos);

  header_index[seqcount] = udb_headerchars;

  unsigned last = 0;
  for(unsigned int i = 0; i < seqcount; i++)
    {
      unsigned int x = header_index[i];
      if ((x < last) || (x >= udb_headerchars))
        fatal("Invalid UDB file");
      seqindex[i].header_p = x;
      seqindex[i].headerlen = header_index[i+1] - x - 1;
      seqindex[i].size = 1;
      last = x;
    }

  xfree(header_index);

  /* headers */

  datap = (char *) xmalloc(udb_headerchars + nucleotides + seqcount);

  pos += largeread(fd_udb, datap, udb_headerchars, pos);

  uint64_t longestheader = 0;
  for(unsigned int i = 0; i < seqcount; i++)
    {
      if (seqindex[i].headerlen > longestheader)
        longestheader = seqindex[i].headerlen;
    }

  /* sequence lengths */

  int * sequence_lengths = (int *) xmalloc(4 * seqcount);

  pos += largeread(fd_udb, sequence_lengths, 4 * seqcount, pos);

  uint64_t sum = 0;
  unsigned int shortest = UINT_MAX;
  unsigned int longest = 0;

  for(unsigned int i = 0; i < seqcount; i++)
    {
      unsigned int x = sequence_lengths[i];

      seqindex[i].seq_p = udb_headerchars + sum;
      seqindex[i].seqlen = x;
      seqindex[i].qual_p = 0;

      if (x < shortest)
        shortest = x;

      if (x > longest)
        longest = x;

      sum += x;

      if (sum > nucleotides)
        fatal("Invalid UDB file");
    }

  xfree(sequence_lengths);

  if (sum != nucleotides)
    fatal("Invalid UDB file");

  /* sequences */

  pos += largeread(fd_udb, datap + udb_headerchars, nucleotides, pos);

  if (pos != filesize)
    fatal("Incorrect UDB file size");

  /* close UDB file */

  close(fd_udb);

  progress_done();
  xfree(prompt);

  /* move sequences and insert zero at end of each sequence */

  progress_init("Reorganizing data in memory", seqcount);
  for(unsigned int i = seqcount-1; i > 0; i--)
    {
      size_t old_p = seqindex[i].seq_p;
      size_t new_p = seqindex[i].seq_p + i;
      size_t len   = seqindex[i].seqlen;
      memmove(datap + new_p, datap + old_p, len);
      *(datap + new_p + len) = 0;
      seqindex[i].seq_p = new_p;
      progress_update(seqcount - i);
    }
  *(datap + seqindex[0].seq_p + seqindex[0].seqlen) = 0;
  progress_done();

  /* Create bitmaps for the most frequent words */

  if (create_bitmaps)
    {
      progress_init("Creating bitmaps", kmerhashsize);
      unsigned int bitmap_mincount = seqcount / 8;
      for(unsigned int i = 0; i < kmerhashsize; i++)
        {
          if (kmercount[i] >= bitmap_mincount)
            {
              kmerbitmap[i] = bitmap_init(seqcount+127); // pad for xmm
              bitmap_reset_all(kmerbitmap[i]);
              for(unsigned j = 0; j < kmercount[i]; j++)
                bitmap_set(kmerbitmap[i], kmerindex[kmerhash[i]+j]);
            }
          progress_update(i+1);
        }
      progress_done();
    }

  /* get abundances and longest header */

  if (parse_abundances)
    {
      progress_init("Parsing abundances", seqcount);
      for(unsigned int i = 0; i < seqcount; i++)
        {
          seqindex[i].size = abundance_get(datap + seqindex[i].header_p,
                                           seqindex[i].headerlen);
          progress_update(i+1);
        }
      progress_done();
    }

  /* set database info */

  dbindex_uh = unique_init();

  db_setinfo(0,
             seqcount,
             nucleotides,
             longest,
             shortest,
             longestheader);

  /* make mapping from indexno to seqno */

  dbindex_map = (unsigned int *) xmalloc(seqcount * sizeof(unsigned int));
  dbindex_count = seqcount;

  for (unsigned int i = 0; i < seqcount; i++)
    dbindex_map[i] = i;

  /* done */

  /* some stats */

  if (!opt_quiet)
    {
      if (seqcount > 0)
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
      if (seqcount > 0)
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
}

void udb_fasta()
{
  /* open FASTA file for writing */

  FILE * fp_output = fopen_output(opt_output);
  if (!fp_output)
    fatal("Unable to open FASTA output file for writing");

  /* read UDB file */

  udb_read(opt_udb2fasta, 0, 0);

  /* dump fasta */

  unsigned int seqcount = db_getsequencecount();
  progress_init("Writing FASTA file", seqcount);
  for(unsigned int i = 0; i < seqcount; i++)
    {
      fasta_print_db_relabel(fp_output, i, i+1);
      progress_update(i+1);
    }
  progress_done();
  fclose(fp_output);

  dbindex_free();
  db_free();
}

void udb_stats()
{
  /* show word statistics for an UDB file */

  /* read UDB file */

  udb_read(opt_udbstats, 0, 0);

  /* analyze word counts */

  wordfreq_t * freqtable = (wordfreq_t *) xmalloc
    (sizeof(wordfreq_t) * kmerhashsize);

  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      freqtable[i].kmer = i;
      freqtable[i].count = kmercount[i];
    }

  qsort(freqtable, kmerhashsize, sizeof(wordfreq_t), wc_compare);

  unsigned int wcmax = freqtable[kmerhashsize-1].count;
  unsigned int wcmedian = ( freqtable[(kmerhashsize / 2) - 1].count +
                            freqtable[kmerhashsize / 2].count ) / 2;

  unsigned int seqcount = db_getsequencecount();
  uint64_t nt = db_getnucleotidecount();

  /* show stats */

  if (opt_log)
    {
      fprintf(fp_log, "      Alphabet  nt\n");
      fprintf(fp_log, "    Word width  %" PRIu64 "\n", opt_wordlength);
      fprintf(fp_log, "     Word ones  %" PRIu64 "\n", opt_wordlength);
      fprintf(fp_log, "        Spaced  No\n");
      fprintf(fp_log, "        Hashed  No\n");
      fprintf(fp_log, "         Coded  No\n");
      fprintf(fp_log, "       Stepped  No\n");
      fprintf(fp_log,
              "         Slots  %u (%.1fk)\n",
              kmerhashsize,
              1.0 * kmerhashsize / 1000.0);
      fprintf(fp_log, "       DBAccel  %u%%\n", udb_dbaccel);
      fprintf(fp_log, "\n");

      fprintf(fp_log,
              "%10" PRIu64 "  DB size (%.1fk)\n",
              nt,
              1.0 * nt / 1000.0);
      fprintf(fp_log, "%10" PRIu64 "  Words\n", kmerindexsize);
      fprintf(fp_log, "%10u  Median size\n", wcmedian);
      fprintf(fp_log,
              "%10.1f  Mean size\n",
              1.0 * kmerindexsize / kmerhashsize);
      fprintf(fp_log, "\n");

      fprintf(fp_log,
              "     iWord         sWord         Cap        Size  Row\n");
      fprintf(fp_log,
              "----------  ------------  ----------  ----------  ---\n");

      for(unsigned int i = 0; i < kmerhashsize; i++)
        {
          fprintf(fp_log,
                  "%10u  ",
                  freqtable[kmerhashsize-1-i].kmer);

          fprintf(fp_log,
                  "%.*s", MAX(12 - (int)(opt_wordlength), 0), "            ");

          fprint_kmer(fp_log, opt_wordlength, freqtable[kmerhashsize-1-i].kmer);

          fprintf(fp_log,
                  "  %10u  %10u",
                  0,
                  freqtable[kmerhashsize-1-i].count);

          fprintf(fp_log, " ");

          for(unsigned j = 0; j < freqtable[kmerhashsize-1-i].count; j++)
            {
              fprintf(fp_log,
                      " %u", kmerindex[kmerhash[freqtable[kmerhashsize-1-i].kmer]+j]);

              if (j == 7)
                break;
            }


          if (freqtable[kmerhashsize-1-i].count > 8)
            fprintf(fp_log, "...");

          fprintf(fp_log, "\n");

          if (i == 10)
            break;
        }

      fprintf(fp_log, "\n\n");

      fprintf(fp_log, "Word width  %" PRIu64 "\n", opt_wordlength);
      fprintf(fp_log, "Slots       %u\n", kmerhashsize);
      fprintf(fp_log, "Words       %" PRIu64 "\n", kmerindexsize);
      fprintf(fp_log, "Max size    %u (", wcmax);
      fprint_kmer(fp_log, opt_wordlength, freqtable[kmerhashsize-1].kmer);
      fprintf(fp_log, ")\n\n");

      fprintf(fp_log, "   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct\n");
      fprintf(fp_log, "----------  ----------  ----------  ----------  ------  ------\n");


      unsigned int size_lo = 0;
      unsigned int size_hi = 0;
      unsigned int x = 0;
      double totpct = 0.0;

      while (size_lo < seqcount)
        {

          int count = 0;
          int size = 0;
          while((x < kmerhashsize) && (freqtable[x].count <= size_hi))
            {
              count++;
              size += freqtable[x].count;
              x++;
            }

          double pct = 100.0 * count / kmerhashsize;
          totpct += pct;

          if (size_lo < size_hi)
            fprintf(fp_log, "%10u", size_lo);
          else
            fprintf(fp_log, "          ");

          fprintf(fp_log, "  %10u", size_hi);

          if (size >= 10000)
            fprintf(fp_log, "  %9.1fk", size * 0.001);
          else
            fprintf(fp_log, "  %10.1f", size * 1.0);

          if (count >= 10000)
            fprintf(fp_log, "  %9.1fk", count * 0.001);
          else
            fprintf(fp_log, "  %10.1f", count * 1.0);

          fprintf(fp_log, "  %5.1f%%  %5.1f%%", pct, totpct);

          int dots = int (pct / 3.0 + 0.5);

          if (dots > 0)
            fprintf(fp_log, "  ");

          for (int i = 0; i < dots ; i++)
            fprintf(fp_log, "*");

          fprintf(fp_log, "\n");

          size_lo = size_hi + 1;
          if (size_hi > 0)
            size_hi *= 2;
          else
            size_hi = 1;
          if (size_hi > seqcount)
            size_hi = seqcount;
        }

      fprintf(fp_log, "----------  ----------  ----------  ----------\n");
      fprintf(fp_log, "                      ");

      if (kmerindexsize >= 10000)
        fprintf(fp_log, "  %9.1fk", kmerindexsize * 0.001);
      else
        fprintf(fp_log, "  %10.1f", kmerindexsize * 1.0);

      if (kmerhashsize >= 10000)
        fprintf(fp_log, "  %9.1fk", kmerhashsize * 0.001);
      else
        fprintf(fp_log, "  %10.1f", kmerhashsize * 1.0);

      fprintf(fp_log, "\n\n");

      fprintf(fp_log, "%10" PRIu64 "  Upper\n", nt);
      fprintf(fp_log, "%10u  Lower (%.1f%%)\n", 0, 0.0);
      fprintf(fp_log, "%10" PRIu64 "  Total\n", nt);
      fprintf(fp_log, "%10" PRIu64 "  Indexed words\n", kmerindexsize);
  }

  xfree(freqtable);
  dbindex_free();
  db_free();
}

void udb_make()
{
  int fd_output = 0;

  fd_output = xopen_write(opt_output);
  if (!fd_output)
    fatal("Unable to open output file for writing");

  db_read(opt_makeudb_usearch, 1);

  if (opt_dbmask == MASK_DUST)
    dust_all();
  else if ((opt_dbmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();

  dbindex_prepare(1, opt_dbmask);
  dbindex_addallsequences(opt_dbmask);

  unsigned int seqcount = db_getsequencecount();
  uint64_t ntcount = db_getnucleotidecount();

  uint64_t header_characters = 0;
  for (unsigned int i=0; i<seqcount; i++)
    header_characters += db_getheaderlen(i) + 1;

  uint64_t kmerhashsize = 1 << (2 * opt_wordlength);

  /* count word matches */
  uint64_t wordmatches = 0;
  for(unsigned int i = 0; i < kmerhashsize; i++)
    wordmatches += kmercount[i];

  uint64_t pos = 0;
  uint64_t progress_all =
    4 * 50 +
    4 * kmerhashsize +
    4 * 1 +
    4 * wordmatches +
    4 * 8 +
    4 * seqcount +
    header_characters +
    4 * seqcount +
    ntcount;

  progress_init("Writing UDB file", progress_all);

  uint64_t buffersize = 4 * MAX(50, seqcount);
  unsigned int * buffer = (unsigned int *) xmalloc(buffersize);
  memset(buffer, 0, buffersize);

  /* Header */
  buffer[0]  = 0x55444246; /* FBDU UDBF */
  buffer[2]  = 32; /* bits */
  buffer[4]  = opt_wordlength; /* default 8 */
  buffer[5]  = 1; /* dbstep */
  buffer[6]  = 100; /* dbaccelpct % */
  buffer[11] = 0; /* slots */
  buffer[13] = (unsigned int) seqcount; /* number of sequences */
  buffer[17] = 0x0000746e; /* alphabet: "nt" */
  buffer[49] = 0x55444266; /* fBDU UDBf */
  pos += largewrite(fd_output, buffer, 50 * 4, 0);

  /* write 4^wordlength uint32's with word match counts */
  pos += largewrite(fd_output, kmercount, 4 * kmerhashsize, pos);

  /* 3BDU */
  buffer[0] = 0x55444233; /* 3BDU UDB3 */
  pos += largewrite(fd_output, buffer, 1 * 4, pos);

  /* lists of sequence no's with matches for all words */
  for(unsigned int i = 0; i < kmerhashsize; i++)
    {
      if (kmerbitmap[i])
        {
          memset(buffer, 0, 4 * kmercount[i]);
          unsigned int elements = 0;
          for (unsigned int j = 0; j < seqcount; j++)
            if (bitmap_get(kmerbitmap[i], j))
              buffer[elements++] = j;
          pos += largewrite(fd_output, buffer, 4 * elements, pos);
        }
      else
        {
          if (kmercount[i] > 0)
            {
              pos += largewrite(fd_output,
                                kmerindex + kmerhash[i],
                                4 * kmercount[i],
                                pos);
            }
        }
    }

  /* New header */
  buffer[0] = 0x55444234; /* 4BDU UDB4 */
  /* 0x005e0db3 */
  buffer[1] = 0x005e0db3;
  /* number of sequences, uint32 */
  buffer[2] = (unsigned int) seqcount;
  /* total number of nucleotides, uint64 */
  buffer[3] = (unsigned int)(ntcount & 0xffffffff);
  buffer[4] = (unsigned int)(ntcount >> 32);
  /* total number of header characters, incl zero-terminator, uint64 */
  buffer[5] = (unsigned int)(header_characters & 0xffffffff);
  buffer[6] = (unsigned int)(header_characters >> 32);
  /* 0x005e0db4 */
  buffer[7] = 0x005e0db4;
  pos += largewrite(fd_output, buffer, 4 * 8, pos);

  /* indices to headers (uint32) */
  unsigned int sum = 0;
  for (unsigned int i = 0; i < seqcount; i++)
    {
      buffer[i] = sum;
      sum += db_getheaderlen(i) + 1;
    }
  pos += largewrite(fd_output, buffer, 4 * seqcount, pos);

  /* headers (ascii, zero terminated, not padded) */
  for (unsigned int i = 0; i < seqcount; i++)
    {
      unsigned int len = db_getheaderlen(i);
      pos += largewrite(fd_output, db_getheader(i), len + 1, pos);
    }

  /* sequence lengths (uint32) */
  for (unsigned int i = 0; i < seqcount; i++)
    buffer[i] = db_getsequencelen(i);
  pos += largewrite(fd_output, buffer, 4 * seqcount, pos);

  /* sequences (ascii, no term, no pad) */
  for (unsigned int i = 0; i < seqcount; i++)
    {
      unsigned int len = db_getsequencelen(i);
      pos += largewrite(fd_output, db_getsequence(i), len, pos);
    }

  if (close(fd_output) != 0)
    fatal("Unable to close UDB file");

  progress_done();
  dbindex_free();
  db_free();
  xfree(buffer);
}
