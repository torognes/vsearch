/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2018, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

uint32_t sff_magic = 0x2e736666;

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

void sff_convert()
{
  if (! opt_fastqout)
    fatal("No output file for sff_convert specified with --fastqout");

  FILE * fp_fastqout = fopen_output(opt_fastqout);
  if (!fp_fastqout)
    fatal("Unable to open fastq output file for writing");

  FILE * fp_sff = fopen_input(opt_sff_convert);
  if (!fp_sff)
    fatal("Unable to open sff input file for reading");

  /* read and check header */

  if (fread(&sff_header, 1, 31, fp_sff) < 31)
    fatal("Unable to read from sff file, or file too small");

  sff_header.magic_number    = bswap_32(sff_header.magic_number);
  sff_header.version         = bswap_32(sff_header.version);
  sff_header.index_offset    = bswap_64(sff_header.index_offset);
  sff_header.index_length    = bswap_32(sff_header.index_length);
  sff_header.number_of_reads = bswap_32(sff_header.number_of_reads);
  sff_header.header_length   = bswap_16(sff_header.header_length);
  sff_header.key_length      = bswap_16(sff_header.key_length);
  sff_header.flows_per_read  = bswap_16(sff_header.flows_per_read);

#if 0
  fprintf(stderr, "Magic number:    %d\n", sff_header.magic_number);
  fprintf(stderr, "Version:         %d\n", sff_header.version);
  fprintf(stderr, "Index offset:    %lld\n", sff_header.index_offset);
  fprintf(stderr, "Index length:    %d\n", sff_header.index_length);
  fprintf(stderr, "Number of reads: %d\n", sff_header.number_of_reads);
  fprintf(stderr, "Header length:   %d\n", sff_header.header_length);
  fprintf(stderr, "Key length:      %d\n", sff_header.key_length);
  fprintf(stderr, "Flows per read:  %d\n", sff_header.flows_per_read);
  fprintf(stderr, "Flowgram format: %d\n", sff_header.flowgram_format_code);
#endif

  if (sff_header.magic_number != sff_magic)
    fatal("Invalid SFF magic number. Must be 0x2e736666 (.sff).");

  if (sff_header.version != 1)
    fatal("Invalid SFF version. Must be 1.");

  if (sff_header.flowgram_format_code != 1)
    fatal("Invalid SFF flowgram format code. Must be 1.");

  if (sff_header.header_length != 8 * ((31 + sff_header.flows_per_read + sff_header.key_length + 7) / 8))
    fatal("Invalid SFF header length");

  if (sff_header.key_length != 4)
    fatal("Invalid SFF key length. Must be 4.");

  /* read and check flow chars, key and padding */

  char * flow_chars = (char *) xmalloc(sff_header.flows_per_read + 1);

  if (fread(flow_chars, 1, sff_header.flows_per_read, fp_sff) < sff_header.flows_per_read)
    fatal("Unable to read flow charactersfrom sff file, or file too small");

  flow_chars[sff_header.flows_per_read] = 0;

  char * key_sequence = (char *) xmalloc(sff_header.key_length + 1);

  if (fread(key_sequence, 1, sff_header.key_length, fp_sff) < sff_header.key_length)
    fatal("Unable to read key sequence from sff file, or file too small");

  key_sequence[sff_header.key_length] = 0;

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

  char padding[8];

  uint32_t padding_length = sff_header.header_length - sff_header.flows_per_read - sff_header.key_length - 31;

#if 0
  fprintf(stderr, "Padding length:  %d\n", padding_length);
  fprintf(stderr, "\n");
#endif

  if ((padding_length > 0) &&
      (fread(padding, 1, padding_length, fp_sff) < padding_length))
    fatal("Unable to read padding from sff file, or file too small");

  double totallength = 0.0;
  uint32_t minimum = UINT_MAX;
  uint32_t maximum = 0;

  progress_init("Converting SFF: ", sff_header.number_of_reads);

  for(uint32_t read_no = 0; read_no < sff_header.number_of_reads; read_no++)
    {
      /* read and check each read header */

      if (fread(&read_header, 1, 16, fp_sff) < 16)
        fatal("Unable to read read header from sff file, or file too small");

      read_header.read_header_length = bswap_16(read_header.read_header_length);
      read_header.name_length = bswap_16(read_header.name_length);
      read_header.number_of_bases = bswap_32(read_header.number_of_bases);
      read_header.clip_qual_left = bswap_16(read_header.clip_qual_left);
      read_header.clip_qual_right = bswap_16(read_header.clip_qual_right);
      read_header.clip_adapter_left = bswap_16(read_header.clip_adapter_left);
      read_header.clip_adapter_right = bswap_16(read_header.clip_adapter_right);

      if (read_header.read_header_length != 8 * ((16 + read_header.name_length + 7) / 8))
        fatal("Invalid SFF read header length");

      char * read_name = (char *) xmalloc(read_header.name_length + 1);

      if (fread(read_name, 1, read_header.name_length, fp_sff) < read_header.name_length)
        fatal("Unable to read read name from sff file, or file too small");

      read_name[read_header.name_length] = 0;

      uint32_t read_header_padding_length = read_header.read_header_length - read_header.name_length - 16;

      char read_header_padding[8];

      if (fread(read_header_padding, 1, read_header_padding_length, fp_sff) < read_header_padding_length)
        fatal("Unable to read read header padding from sff file, or file too small");

      if (read_header.clip_qual_left > read_header.number_of_bases)
        fatal("Invalid clip_qual_left value");
      if (read_header.clip_adapter_left > read_header.number_of_bases)
        fatal("Invalid clip_adapter_left value");
      if (read_header.clip_qual_right > read_header.number_of_bases)
        fatal("Invalid clip_qual_right value");
      if (read_header.clip_adapter_right > read_header.number_of_bases)
        fatal("Invalid clip_adapter_right value");

      /* read and check the flowgram and sequence */

      uint16_t * flowgram_values = (uint16_t *) xmalloc(sff_header.flows_per_read * 2);
      if (fread(flowgram_values, 1, 2 * sff_header.flows_per_read, fp_sff) < sff_header.flows_per_read)
        fatal("Unable to read flowgram values from sff file, or file too small");

      uint8_t * flow_index_per_base = (uint8_t *) xmalloc(read_header.number_of_bases);
      if (fread(flow_index_per_base, 1, read_header.number_of_bases, fp_sff) < read_header.number_of_bases)
        fatal("Unable to read flow indices from sff file, or file too small");

      char * bases = (char *) xmalloc(read_header.number_of_bases + 1);
      if (fread(bases, 1, read_header.number_of_bases, fp_sff) < read_header.number_of_bases)
        fatal("Unable to read read length from sff file, or file too small");
      bases[read_header.number_of_bases] = 0;

      uint8_t * quality_scores = (uint8_t *) xmalloc(read_header.number_of_bases);
      if (fread(quality_scores, 1, read_header.number_of_bases, fp_sff) < read_header.number_of_bases)
        fatal("Unable to read quality scores from sff file, or file too small");

      char * qual = (char *) xmalloc(read_header.number_of_bases + 1);

      for(uint32_t base_no = 0; base_no < read_header.number_of_bases; base_no++)
        {
          uint8_t q = quality_scores[base_no];
          if (q < opt_fastq_qminout)
            q = opt_fastq_qminout;
          if (q > opt_fastq_qmaxout)
            q = opt_fastq_qmaxout;
          char c = opt_fastq_asciiout + q;
          qual[base_no] = c;
        }
      qual[read_header.number_of_bases] = 0;

      uint32_t read_data_length = (2 * sff_header.flows_per_read + 3 * read_header.number_of_bases);
      uint32_t read_data_padded_length = 8 * ((read_data_length + 7) / 8);
      uint32_t read_data_padding_length = read_data_padded_length - read_data_length;
      char read_data_padding[8];
      if (fread(read_data_padding, 1, read_data_padding_length, fp_sff) < read_data_padding_length)
        fatal("Unable to read read data padding from sff file, or file too small");

      uint32_t clip_start = 0;
      uint32_t clip_end = read_header.number_of_bases;

      clip_start = MAX(1, MAX(read_header.clip_qual_left, read_header.clip_adapter_left)) - 1;
      clip_end = MIN((read_header.clip_qual_right == 0 ? read_header.number_of_bases : read_header.clip_qual_right), (read_header.clip_adapter_right == 0 ? read_header.number_of_bases : read_header.clip_adapter_right));

      /* make the clipped bases lowercase and the rest uppercase */
      for (uint32_t i = 0; i < read_header.number_of_bases; i++)
        {
          if ((i < clip_start) || (i >= clip_end))
            bases[i] = tolower(bases[i]);
          else
            bases[i] = toupper(bases[i]);
        }

      if (! opt_sff_clip)
        {
          clip_start = 0;
          clip_end = read_header.number_of_bases;
        }

      uint32_t length = clip_end - clip_start;

      bases[clip_end] = 0;
      qual[clip_end] = 0;

      fastq_print_general(fp_fastqout,
                          bases + clip_start,
                          length,
                          read_name,
                          strlen(read_name),
                          qual + clip_start,
                          0, 0, 0, 0);

      xfree(read_name);
      xfree(flowgram_values);
      xfree(flow_index_per_base);
      xfree(bases);
      xfree(quality_scores);
      xfree(qual);

      totallength += length;
      if (length < minimum)
        minimum = length;
      if (length > maximum)
        maximum = length;

      progress_update(read_no + 1);
    }
  progress_done();

  double average = totallength / sff_header.number_of_reads;

  if (! opt_quiet)
    {
      if (sff_header.number_of_reads > 0)
        fprintf(stderr, "Sequence length: minimum %d, average %.1f, maximum %d\n",
                minimum,
                average,
                maximum);
    }

  if (opt_log)
    {
      if (sff_header.number_of_reads > 0)
        fprintf(fp_log, "Sequence length: minimum %d, average %.1f, maximum %d\n",
                minimum,
                average,
                maximum);
    }

  if (sff_header.index_offset > 0)
    {
      size_t filepos = xftello(fp_sff);

      if (filepos != sff_header.index_offset)
        fatal("Invalid sff file (wrong index offset)");

      if (fseek(fp_sff, sff_header.index_length, SEEK_CUR) != 0)
        fatal("Unable to seek past index in SFF file");
    }

  /* ignore padding and rest of file */

  fclose(fp_sff);
  fclose(fp_fastqout);
}
