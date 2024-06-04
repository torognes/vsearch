/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include <array>


constexpr uint32_t sff_magic = 0x2e736666;

struct sff_header_s
{
  uint32_t magic_number; /* .sff */
  uint32_t version;
  uint64_t index_offset;
  uint32_t index_length;
  uint32_t number_of_reads;
  uint16_t header_length;
  uint16_t key_length;
  uint16_t flows_per_read;
  uint8_t  flowgram_format_code;
} sff_header;

struct sff_read_header_s
{
  uint16_t read_header_length;
  uint16_t name_length;
  uint32_t number_of_bases;
  uint16_t clip_qual_left;
  uint16_t clip_qual_right;
  uint16_t clip_adapter_left;
  uint16_t clip_adapter_right;
} read_header;


auto fskip(FILE * fp, uint64_t length) -> uint64_t
{
  /* read given amount of data from a stream and ignore it */
  /* used instead of seeking in order to work with pipes   */
  static constexpr auto blocksize = uint64_t{4096};
  std::array<char, blocksize> buffer;

  uint64_t skipped = 0;
  uint64_t rest = length;

  while (rest > 0)
    {
      uint64_t const want = (rest > blocksize) ? blocksize : rest;
      uint64_t const got = fread(buffer.data(), 1, want, fp);
      skipped += got;
      rest -= got;
      if (got < want)
        {
          break;
        }
    }
  return skipped;
}


auto sff_convert() -> void
{
  if (! opt_fastqout)
    {
      fatal("No output file for sff_convert specified with --fastqout.");
    }

  FILE * fp_fastqout = fopen_output(opt_fastqout);
  if (! fp_fastqout)
    {
      fatal("Unable to open FASTQ output file for writing.");
    }

  FILE * fp_sff = fopen_input(opt_sff_convert);
  if (! fp_sff)
    {
      fatal("Unable to open SFF input file for reading.");
    }

  /* read and check header */

  uint64_t filepos = 0;

  if (fread(&sff_header, 1, 31, fp_sff) < 31)
    {
      fatal("Unable to read from SFF file. File may be truncated.");
    }
  filepos += 31;

  sff_header.magic_number    = bswap_32(sff_header.magic_number);
  sff_header.version         = bswap_32(sff_header.version);
  sff_header.index_offset    = bswap_64(sff_header.index_offset);
  sff_header.index_length    = bswap_32(sff_header.index_length);
  sff_header.number_of_reads = bswap_32(sff_header.number_of_reads);
  sff_header.header_length   = bswap_16(sff_header.header_length);
  sff_header.key_length      = bswap_16(sff_header.key_length);
  sff_header.flows_per_read  = bswap_16(sff_header.flows_per_read);

  if (sff_header.magic_number != sff_magic)
    {
      fatal("Invalid SFF file. Incorrect magic number. Must be 0x2e736666 (.sff).");
    }

  if (sff_header.version != 1)
    {
      fatal("Invalid SFF file. Incorrect version. Must be 1.");
    }

  if (sff_header.flowgram_format_code != 1)
    {
      fatal("Invalid SFF file. Incorrect flowgram format code. Must be 1.");
    }

  if (sff_header.header_length != 8 * ((31 + sff_header.flows_per_read + sff_header.key_length + 7) / 8))
    {
      fatal("Invalid SFF file. Incorrect header length.");
    }

  if (sff_header.key_length != 4)
    {
      fatal("Invalid SFF file. Incorrect key length. Must be 4.");
    }

  if ((sff_header.index_length > 0) && (sff_header.index_length < 8))
    {
      fatal("Invalid SFF file. Incorrect index size. Must be at least 8.");
    }

  /* read and check flow chars, key and padding */

  if (fskip(fp_sff, sff_header.flows_per_read) < sff_header.flows_per_read)
    {
      fatal("Invalid SFF file. Unable to read flow characters. File may be truncated.");
    }
  filepos += sff_header.flows_per_read;

  char * key_sequence = (char *) xmalloc(sff_header.key_length + 1);
  if (fread(key_sequence, 1, sff_header.key_length, fp_sff) < sff_header.key_length)
    {
      fatal("Invalid SFF file. Unable to read key sequence. File may be truncated.");
    }
  key_sequence[sff_header.key_length] = 0;
  filepos += sff_header.key_length;

  uint32_t padding_length = sff_header.header_length - sff_header.flows_per_read - sff_header.key_length - 31;
  if (fskip(fp_sff, padding_length) < padding_length)
    {
      fatal("Invalid SFF file. Unable to read padding. File may be truncated.");
    }
  filepos += padding_length;

  double totallength = 0.0;
  uint32_t minimum = UINT_MAX;
  uint32_t maximum = 0;

  bool index_done = (sff_header.index_offset == 0) || (sff_header.index_length == 0);
  bool index_odd = false;
  char index_kind[9];

  uint32_t index_padding = 0;
  if ((sff_header.index_length & 7) > 0)
    {
      index_padding = 8 - (sff_header.index_length & 7);
    }

  if (! opt_quiet)
    {
      fprintf(stderr, "Number of reads: %d\n", sff_header.number_of_reads);
      fprintf(stderr, "Flows per read:  %d\n", sff_header.flows_per_read);
      fprintf(stderr, "Key sequence:    %s\n", key_sequence);
    }

  if (opt_log)
    {
      fprintf(fp_log, "Number of reads: %d\n", sff_header.number_of_reads);
      fprintf(fp_log, "Flows per read:  %d\n", sff_header.flows_per_read);
      fprintf(fp_log, "Key sequence:    %s\n", key_sequence);
    }

  progress_init("Converting SFF: ", sff_header.number_of_reads);

  for (uint32_t read_no = 0; read_no < sff_header.number_of_reads; read_no++)
    {
      /* check if the index block is here */

      if (! index_done)
        {
          if (filepos == sff_header.index_offset)
            {
              if (fread(index_kind, 1, 8, fp_sff) < 8)
                {
                  fatal("Invalid SFF file. Unable to read index header. File may be truncated.");
                }
              filepos += 8;
              index_kind[8] = 0;

              uint64 index_size = sff_header.index_length - 8 + index_padding;
              if (fskip(fp_sff, index_size) != index_size)
                {
                  fatal("Invalid SFF file. Unable to read entire index. File may be truncated.");
                }

              filepos += index_size;
              index_done = true;
              index_odd = true;
            }
        }

      /* read and check each read header */

      if (fread(&read_header, 1, 16, fp_sff) < 16)
        {
          fatal("Invalid SFF file. Unable to read read header. File may be truncated.");
        }
      filepos += 16;

      read_header.read_header_length = bswap_16(read_header.read_header_length);
      read_header.name_length = bswap_16(read_header.name_length);
      read_header.number_of_bases = bswap_32(read_header.number_of_bases);
      read_header.clip_qual_left = bswap_16(read_header.clip_qual_left);
      read_header.clip_qual_right = bswap_16(read_header.clip_qual_right);
      read_header.clip_adapter_left = bswap_16(read_header.clip_adapter_left);
      read_header.clip_adapter_right = bswap_16(read_header.clip_adapter_right);

      if (read_header.read_header_length != 8 * ((16 + read_header.name_length + 7) / 8))
        {
          fatal("Invalid SFF file. Incorrect read header length.");
        }
      if (read_header.clip_qual_left > read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Incorrect clip_qual_left value.");
        }
      if (read_header.clip_adapter_left > read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Incorrect clip_adapter_left value.");
        }
      if (read_header.clip_qual_right > read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Incorrect clip_qual_right value.");
        }
      if (read_header.clip_adapter_right > read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Incorrect clip_adapter_right value.");
        }

      char * read_name = (char *) xmalloc(read_header.name_length + 1);
      if (fread(read_name, 1, read_header.name_length, fp_sff) < read_header.name_length)
        {
          fatal("Invalid SFF file. Unable to read read name. File may be truncated.");
        }
      filepos += read_header.name_length;
      read_name[read_header.name_length] = 0;

      uint32_t read_header_padding_length = read_header.read_header_length - read_header.name_length - 16;
      if (fskip(fp_sff, read_header_padding_length) < read_header_padding_length)
        {
          fatal("Invalid SFF file. Unable to read read header padding. File may be truncated.");
        }
      filepos += read_header_padding_length;

      /* read and check the flowgram and sequence */

      if (fskip(fp_sff, 2 * sff_header.flows_per_read) < sff_header.flows_per_read)
        {
          fatal("Invalid SFF file. Unable to read flowgram values. File may be truncated.");
        }
      filepos += 2 * sff_header.flows_per_read;

      if (fskip(fp_sff, read_header.number_of_bases) < read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Unable to read flow indices. File may be truncated.");
        }
      filepos += read_header.number_of_bases;

      char * bases = (char *) xmalloc(read_header.number_of_bases + 1);
      if (fread(bases, 1, read_header.number_of_bases, fp_sff) < read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Unable to read read length. File may be truncated.");
        }
      bases[read_header.number_of_bases] = 0;
      filepos += read_header.number_of_bases;

      char * qual = (char *) xmalloc(read_header.number_of_bases + 1);
      if (fread(qual, 1, read_header.number_of_bases, fp_sff) < read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Unable to read quality scores. File may be truncated.");
        }
      filepos += read_header.number_of_bases;

      /* convert quality scores to ascii characters */
      for (uint32_t base_no = 0; base_no < read_header.number_of_bases; base_no++)
        {
          int q = qual[base_no];
          if (q < opt_fastq_qminout)
            {
              q = opt_fastq_qminout;
            }
          if (q > opt_fastq_qmaxout)
            {
              q = opt_fastq_qmaxout;
            }
          qual[base_no] = opt_fastq_asciiout + q;
        }
      qual[read_header.number_of_bases] = 0;

      uint32_t read_data_length = ((2 * sff_header.flows_per_read) + (3 * read_header.number_of_bases));
      uint32_t read_data_padded_length = 8 * ((read_data_length + 7) / 8);
      uint32_t read_data_padding_length = read_data_padded_length - read_data_length;

      if (fskip(fp_sff, read_data_padding_length) < read_data_padding_length)
        {
          fatal("Invalid SFF file. Unable to read read data padding. File may be truncated.");
        }
      filepos += read_data_padding_length;

      uint32_t clip_start = 0;
      clip_start = MAX(1, MAX(read_header.clip_qual_left, read_header.clip_adapter_left)) - 1;

      uint32_t clip_end = read_header.number_of_bases;
      clip_end = MIN((read_header.clip_qual_right == 0 ? read_header.number_of_bases : read_header.clip_qual_right), (read_header.clip_adapter_right == 0 ? read_header.number_of_bases : read_header.clip_adapter_right));

      /* make the clipped bases lowercase and the rest uppercase */
      for (uint32_t i = 0; i < read_header.number_of_bases; i++)
        {
          if ((i < clip_start) || (i >= clip_end))
            {
              bases[i] = tolower(bases[i]);
            }
          else
            {
              bases[i] = toupper(bases[i]);
            }
        }

      if (opt_sff_clip)
        {
          bases[clip_end] = 0;
          qual[clip_end] = 0;
        }
      else
        {
          clip_start = 0;
          clip_end = read_header.number_of_bases;
        }

      uint32_t length = clip_end - clip_start;

      fastq_print_general(fp_fastqout,
                          bases + clip_start,
                          length,
                          read_name,
                          strlen(read_name),
                          qual + clip_start,
                          1, read_no + 1, -1.0);

      xfree(read_name);
      xfree(bases);
      xfree(qual);

      totallength += length;
      if (length < minimum)
        {
          minimum = length;
        }
      if (length > maximum)
        {
          maximum = length;
        }

      progress_update(read_no + 1);
    }
  progress_done();

  /* check if the index block is here */

  if (! index_done)
    {
      if (filepos == sff_header.index_offset)
        {
          if (fread(index_kind, 1, 8, fp_sff) < 8)
            {
              fatal("Invalid SFF file. Unable to read index header. File may be truncated.");
            }
          filepos += 8;
          index_kind[8] = 0;

          uint64 index_size = sff_header.index_length - 8;
          if (fskip(fp_sff, index_size) != index_size)
            {
              fatal("Invalid SFF file. Unable to read entire index. File may be truncated.");
            }

          filepos += index_size;
          index_done = true;

          /* try to skip padding, if any */

          if (index_padding > 0)
            {
              uint64_t got = fskip(fp_sff, index_padding);
              if ((got < index_padding) && (got != 0))
                {
                  fprintf(stderr, "WARNING: Additional data at end of SFF file ignored\n");
                }
            }
        }
    }

  if (! index_done)
    {
      fprintf(stderr, "WARNING: SFF index missing\n");
      if (opt_log)
        {
          fprintf(fp_log, "WARNING: SFF index missing\n");
        }
    }

  if (index_odd)
    {
      fprintf(stderr, "WARNING: Index at unusual position in file\n");
      if (opt_log)
        {
          fprintf(fp_log, "WARNING: Index at unusual position in file\n");
        }
    }

  /* ignore the rest of file */

  /* try reading just another byte */

  if (fskip(fp_sff, 1) > 0)
    {
      fprintf(stderr, "WARNING: Additional data at end of SFF file ignored\n");
      if (opt_log)
        {
          fprintf(fp_log, "WARNING: Additional data at end of SFF file ignored\n");
        }
    }

  fclose(fp_sff);
  fclose(fp_fastqout);

  double average = totallength / sff_header.number_of_reads;

  if (! opt_quiet)
    {
      if (sff_header.index_length > 0)
        {
          fprintf(stderr, "Index type:      %s\n", index_kind);
        }
      fprintf(stderr, "\nSFF file read successfully.\n");
      if (sff_header.number_of_reads > 0)
        {
          fprintf(stderr, "Sequence length: minimum %d, average %.1f, maximum %d\n",
                  minimum,
                  average,
                  maximum);
        }
    }

  if (opt_log)
    {
      if (sff_header.index_length > 0)
        {
          fprintf(fp_log, "Index type:      %s\n", index_kind);
        }
      fprintf(fp_log, "\nSFF file read successfully.\n");
      if (sff_header.number_of_reads > 0)
        {
          fprintf(fp_log, "Sequence length: minimum %d, average %.1f, maximum %d\n",
                  minimum,
                  average,
                  maximum);
        }
    }

  xfree(key_sequence);
}
