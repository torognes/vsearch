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
#include <algorithm>  // std::min, std::max
#include <array>
#include <cassert>
#include <cctype>  // std::tolower, std::toupper
#include <cstdint>  // uint64_t, uint32_t, uint16_t, uint8_t
#include <cstdio>  // std::fprintf, std::FILE, std:fclose, std::fread, std::size_t
#include <cstring>  // std::strlen
#include <limits>
#include <vector>

// sff file map:
// - common header
// - index (optional, skipped when present)
// - reads
// - index (optional, can be at the end) (index_offset tells us where it is)


constexpr uint32_t sff_magic = 0x2e736666;  // encoding the string ".sff"
constexpr std::size_t n_bytes_in_header = 31;  // first part of the header is 31 bytes in total
constexpr auto byte_size = sizeof(uint8_t);
constexpr auto memory_alignment = 8U;
constexpr auto max_padding_length = 7U;
constexpr auto expected_version_number = 1U;
constexpr auto expected_flowgram_format_code = 1U;
constexpr auto expected_key_length = 4U;  // key sequences always have 4 nucleotides?
constexpr auto minimal_index_length = 8U;

// SFF format expects the following to be true:
static_assert(sizeof(uint8_t) == 1, "sff expects a uint8_t of size 1");
static_assert(sizeof(uint16_t) == 2, "sff expects a uint16_t of size 2");
static_assert(sizeof(uint32_t) == 4, "sff expects a uint32_t of size 4");
static_assert(sizeof(uint64_t) == memory_alignment, "sff expects a uint64_t of size 8");

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
};

struct sff_read_header_s
{
  uint16_t read_header_length;
  uint16_t name_length;
  uint32_t number_of_bases;
  uint16_t clip_qual_left;
  uint16_t clip_qual_right;
  uint16_t clip_adapter_left;
  uint16_t clip_adapter_right;
};


auto fskip(std::FILE * file_handle, uint64_t length) -> uint64_t
{
  /* read given amount of data from a stream and ignore it */
  /* used instead of seeking in order to work with pipes   */
  static constexpr uint64_t blocksize = 4096;
  std::array<char, blocksize> buffer {{}};

  uint64_t skipped = 0;
  uint64_t rest = length;

  while (rest > 0)
    {
      auto const want = (rest > blocksize) ? blocksize : rest;
      uint64_t const got = std::fread(buffer.data(), byte_size, want, file_handle);
      skipped += got;
      rest -= got;
      if (got < want)
        {
          break;
        }
    }
  return skipped;
}


auto read_sff_header(std::FILE * sff_handle) -> struct sff_header_s {
  assert(sff_handle != nullptr);

  struct sff_header_s sff_header;
  auto const n_bytes_read = std::fread(&sff_header, byte_size, n_bytes_in_header, sff_handle);
  if (n_bytes_read < n_bytes_in_header) {
    fatal("Unable to read from SFF file. File may be truncated.");
  }

  // SFF multi-byte numeric values are stored using a big-endian byte order
  // vsearch expects little-endian, so we need to swap bytes
  // refactoring: C++23 std::byteswap()
  sff_header.magic_number    = bswap_32(sff_header.magic_number);
  sff_header.version         = bswap_32(sff_header.version);
  sff_header.index_offset    = bswap_64(sff_header.index_offset);
  sff_header.index_length    = bswap_32(sff_header.index_length);
  sff_header.number_of_reads = bswap_32(sff_header.number_of_reads);
  sff_header.header_length   = bswap_16(sff_header.header_length);
  sff_header.key_length      = bswap_16(sff_header.key_length);
  sff_header.flows_per_read  = bswap_16(sff_header.flows_per_read);

  return sff_header;
};


auto check_sff_header(struct sff_header_s const &sff_header) -> void {
  if (sff_header.magic_number != sff_magic)
    {
      fatal("Invalid SFF file. Incorrect magic number. Must be 0x2e736666 (.sff).");
    }

  if (sff_header.version != expected_version_number)
    {
      fatal("Invalid SFF file. Incorrect version. Must be 1.");
    }

  if (sff_header.flowgram_format_code != expected_flowgram_format_code)
    {
      fatal("Invalid SFF file. Incorrect flowgram format code. Must be 1.");
    }

  if (sff_header.header_length != 8 * ((n_bytes_in_header + sff_header.flows_per_read + sff_header.key_length + max_padding_length) / 8))
    {
      fatal("Invalid SFF file. Incorrect header length.");
    }

  if (sff_header.key_length != expected_key_length)
    {
      fatal("Invalid SFF file. Incorrect key length. Must be 4.");
    }

  if ((sff_header.index_length > 0) and (sff_header.index_length < minimal_index_length))
    {
      fatal("Invalid SFF file. Incorrect index size. Must be at least 8.");
    }
  // index_length includes the bytes of index_magic_number (uint32_t),
  // index_version (char[4]), the 8n bytes for the indexing method,
  // and padding bytes. And the length of the index section is
  // divisible by 8. So, index_length modulo 8 must be null
  assert((sff_header.index_length & memory_alignment) == 0);
};


auto skip_sff_flow_chars(std::FILE * sff_handle, uint64_t n_bytes_to_skip) -> void {
  auto const n_bytes_skipped = fskip(sff_handle, n_bytes_to_skip);
  if (n_bytes_skipped < n_bytes_to_skip) {
    fatal("Invalid SFF file. Unable to read flow characters. File may be truncated.");
  }
};


auto skip_sff_padding_length(std::FILE * sff_handle, uint64_t n_bytes_to_skip) -> void {
  auto const n_bytes_skipped = fskip(sff_handle, n_bytes_to_skip);
  if (n_bytes_skipped < n_bytes_to_skip) {
    fatal("Invalid SFF file. Unable to read padding. File may be truncated.");
  }
};


auto read_key_sequence(std::FILE * sff_handle, uint16_t n_bytes_to_read) -> std::vector<char> {
  std::vector<char> key_sequence(n_bytes_to_read + 1);
  auto const n_bytes_read = std::fread(key_sequence.data(), byte_size, n_bytes_to_read, sff_handle);
  if (n_bytes_read < n_bytes_to_read) {
    fatal("Invalid SFF file. Unable to read key sequence. File may be truncated.");
  }
  assert(key_sequence.back() == '\0');  // C-string should be null-terminated
  return key_sequence;
}


auto check_for_additional_tail_data(std::FILE * sff_handle) -> void {
  // try to read another byte
  auto const n_bytes_read = fskip(sff_handle, byte_size);
  if (n_bytes_read == 0) {
    return;  // no trailing data
  }
  auto const * const message = "WARNING: Additional data at end of SFF file ignored\n";
  std::fprintf(stderr, message);
  if (opt_log != nullptr) {
    std::fprintf(fp_log, message);
  }
};


auto sff_convert() -> void
{
  if (opt_fastqout == nullptr)
    {
      fatal("No output file for sff_convert specified with --fastqout.");
    }

  auto * fp_fastqout = fopen_output(opt_fastqout);
  if (fp_fastqout == nullptr)
    {
      fatal("Unable to open FASTQ output file for writing.");
    }

  auto * fp_sff = fopen_input(opt_sff_convert);
  if (fp_sff == nullptr)
    {
      fatal("Unable to open SFF input file for reading.");
    }

  /* read and check header */

  uint64_t filepos = 0;

  auto const sff_header = read_sff_header(fp_sff);
  filepos += n_bytes_in_header;
  check_sff_header(sff_header);


  /* read and check flow chars, key and padding */

  skip_sff_flow_chars(fp_sff, sff_header.flows_per_read);
  filepos += sff_header.flows_per_read;

  auto const key_sequence = read_key_sequence(fp_sff, sff_header.key_length);
  filepos += sff_header.key_length;

  uint32_t const padding_length = sff_header.header_length - sff_header.flows_per_read - sff_header.key_length - n_bytes_in_header;
  skip_sff_padding_length(fp_sff, padding_length);
  filepos += padding_length;

  double totallength = 0.0;
  uint32_t minimum = std::numeric_limits<uint32_t>::max();
  uint32_t maximum = 0;

  bool index_done = (sff_header.index_offset == 0) or (sff_header.index_length == 0);
  bool index_odd = false;
  std::array<char, 9> index_kind;

  uint32_t index_padding = 0;
  if ((sff_header.index_length & max_padding_length) > 0)
    {
      index_padding = 8 - (sff_header.index_length & max_padding_length);
    }

  if (not opt_quiet)
    {
      fprintf(stderr, "Number of reads: %d\n", sff_header.number_of_reads);
      fprintf(stderr, "Flows per read:  %d\n", sff_header.flows_per_read);
      fprintf(stderr, "Key sequence:    %s\n", key_sequence.data());
    }

  if (opt_log != nullptr)
    {
      fprintf(fp_log, "Number of reads: %d\n", sff_header.number_of_reads);
      fprintf(fp_log, "Flows per read:  %d\n", sff_header.flows_per_read);
      fprintf(fp_log, "Key sequence:    %s\n", key_sequence.data());
    }

  progress_init("Converting SFF: ", sff_header.number_of_reads);

  auto const fastq_qminout = static_cast<int>(opt_fastq_qminout);
  auto const fastq_qmaxout = static_cast<int>(opt_fastq_qmaxout);

  for (uint32_t read_no = 0; read_no < sff_header.number_of_reads; read_no++)
    {
      /* check if the index block is here */

      if (not index_done)
        {
          if (filepos == sff_header.index_offset)
            {
              if (std::fread(index_kind.data(), byte_size, 8, fp_sff) < 8)
                {
                  fatal("Invalid SFF file. Unable to read index header. File may be truncated.");
                }
              filepos += 8;
              index_kind[8] = 0;

              uint64 const index_size = sff_header.index_length - 8 + index_padding;
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

      struct sff_read_header_s read_header;

      if (std::fread(&read_header, byte_size, 16, fp_sff) < 16)
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

      if (read_header.read_header_length != 8 * ((16 + read_header.name_length + max_padding_length) / 8))
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

      std::vector<char> read_name(read_header.name_length + 1);
      if (std::fread(read_name.data(), byte_size, read_header.name_length, fp_sff) < read_header.name_length)
        {
          fatal("Invalid SFF file. Unable to read read name. File may be truncated.");
        }
      filepos += read_header.name_length;
      read_name[read_header.name_length] = '\0';

      uint32_t const read_header_padding_length = read_header.read_header_length - read_header.name_length - 16;
      if (fskip(fp_sff, read_header_padding_length) < read_header_padding_length)
        {
          fatal("Invalid SFF file. Unable to read read header padding. File may be truncated.");
        }
      filepos += read_header_padding_length;

      /* read and check the flowgram and sequence */

      if (fskip(fp_sff, 2UL * sff_header.flows_per_read) < sff_header.flows_per_read)
        {
          fatal("Invalid SFF file. Unable to read flowgram values. File may be truncated.");
        }
      filepos += 2UL * sff_header.flows_per_read;

      if (fskip(fp_sff, read_header.number_of_bases) < read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Unable to read flow indices. File may be truncated.");
        }
      filepos += read_header.number_of_bases;

      std::vector<char> bases(read_header.number_of_bases + 1);
      if (std::fread(bases.data(), byte_size, read_header.number_of_bases, fp_sff) < read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Unable to read read length. File may be truncated.");
        }
      bases[read_header.number_of_bases] = '\0';
      filepos += read_header.number_of_bases;

      std::vector<char> quality_scores(read_header.number_of_bases + 1);
      if (std::fread(quality_scores.data(), byte_size, read_header.number_of_bases, fp_sff) < read_header.number_of_bases)
        {
          fatal("Invalid SFF file. Unable to read quality scores. File may be truncated.");
        }
      filepos += read_header.number_of_bases;

      /* convert quality scores to ascii characters */
      for (uint32_t base_no = 0; base_no < read_header.number_of_bases; base_no++)
        {
          int quality_score = quality_scores[base_no];
          quality_score = std::max(quality_score, fastq_qminout);
          quality_score = std::min(quality_score, fastq_qmaxout);
          quality_scores[base_no] = opt_fastq_asciiout + quality_score;  // refactoring C++17: std::clamp(q, min, max) + offset
        }
      quality_scores[read_header.number_of_bases] = '\0';

      uint32_t const read_data_length = ((2 * sff_header.flows_per_read) + (3 * read_header.number_of_bases));
      uint32_t const read_data_padded_length = 8 * ((read_data_length + max_padding_length) / 8);
      uint32_t const read_data_padding_length = read_data_padded_length - read_data_length;

      if (fskip(fp_sff, read_data_padding_length) < read_data_padding_length)
        {
          fatal("Invalid SFF file. Unable to read read data padding. File may be truncated.");
        }
      filepos += read_data_padding_length;

      uint32_t clip_start = std::max({uint16_t{1}, read_header.clip_qual_left, read_header.clip_adapter_left}) - 1 ;

      uint32_t clip_end = std::min((read_header.clip_qual_right == 0 ? read_header.number_of_bases : read_header.clip_qual_right), (read_header.clip_adapter_right == 0 ? read_header.number_of_bases : read_header.clip_adapter_right));

      /* make the clipped bases lowercase and the rest uppercase */
      for (uint32_t i = 0; i < read_header.number_of_bases; i++)
        {
          if ((i < clip_start) or (i >= clip_end))
            {
              bases[i] = std::tolower(bases[i]);
            }
          else
            {
              bases[i] = std::toupper(bases[i]);
            }
        }

      if (opt_sff_clip)
        {
          bases[clip_end] = '\0';
          quality_scores[clip_end] = '\0';
        }
      else
        {
          clip_start = 0;
          clip_end = read_header.number_of_bases;
        }

      uint32_t const length = clip_end - clip_start;

      fastq_print_general(fp_fastqout,
                          bases.data() + clip_start,
                          length,
                          read_name.data(),
                          std::strlen(read_name.data()),
                          quality_scores.data() + clip_start,
                          1, read_no + 1, -1.0);


      totallength += length;
      minimum = std::min(length, minimum);
      maximum = std::max(length, maximum);

      progress_update(read_no + 1);
    }
  progress_done();

  /* check if the index block is here */

  if (not index_done)
    {
      if (filepos == sff_header.index_offset)
        {
          if (std::fread(index_kind.data(), byte_size, 8, fp_sff) < 8)
            {
              fatal("Invalid SFF file. Unable to read index header. File may be truncated.");
            }
          filepos += 8;
          index_kind[8] = 0;

          uint64 const index_size = sff_header.index_length - 8;
          if (fskip(fp_sff, index_size) != index_size)
            {
              fatal("Invalid SFF file. Unable to read entire index. File may be truncated.");
            }

          filepos += index_size;
          index_done = true;

          /* try to skip padding, if any */

          if (index_padding > 0)
            {
              uint64_t const got = fskip(fp_sff, index_padding);
              if ((got < index_padding) and (got != 0))
                {
                  fprintf(stderr, "WARNING: Additional data at end of SFF file ignored\n");
                }
            }
        }
    }

  if (not index_done)
    {
      fprintf(stderr, "WARNING: SFF index missing\n");
      if (opt_log != nullptr)
        {
          fprintf(fp_log, "WARNING: SFF index missing\n");
        }
    }

  if (index_odd)
    {
      fprintf(stderr, "WARNING: Index at unusual position in file\n");
      if (opt_log != nullptr)
        {
          fprintf(fp_log, "WARNING: Index at unusual position in file\n");
        }
    }

  /* ignore the rest of file */

  check_for_additional_tail_data(fp_sff);

  std::fclose(fp_sff);
  std::fclose(fp_fastqout);

  double const average = totallength / sff_header.number_of_reads;

  if (not opt_quiet)
    {
      if (sff_header.index_length > 0)
        {
          fprintf(stderr, "Index type:      %s\n", index_kind.data());
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

  if (opt_log != nullptr)
    {
      if (sff_header.index_length > 0)
        {
          fprintf(fp_log, "Index type:      %s\n", index_kind.data());
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

}
