/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t


/* fasta input */

auto fasta_open_rest(fastx_handle input_handle) -> void;
auto fasta_open(const char * filename) -> fastx_handle;
auto fasta_close(fastx_handle input_handle) -> void;
auto fasta_next(fastx_handle input_handle,
                bool truncateatspace,
                const unsigned char * char_mapping) -> bool;
auto fasta_get_position(fastx_handle input_handle) -> uint64_t;
auto fasta_get_size(fastx_handle input_handle) -> uint64_t;
auto fasta_get_lineno(fastx_handle input_handle) -> uint64_t;
auto fasta_get_seqno(fastx_handle input_handle) -> uint64_t;
auto fasta_get_header(fastx_handle input_handle) -> char *;
auto fasta_get_sequence(fastx_handle input_handle) -> char *;
auto fasta_get_header_length(fastx_handle input_handle) -> uint64_t;
auto fasta_get_sequence_length(fastx_handle input_handle) -> uint64_t;
auto fasta_get_abundance(fastx_handle input_handle) -> int64_t;
auto fasta_get_abundance_and_presence(fastx_handle input_handle) -> int64_t;

/* fasta output */

auto fasta_print(std::FILE * output_handle,
                 const char * header,
                 char * seq,
                 uint64_t len) -> void;

auto fasta_print_general(std::FILE * output_handle,
                         const char * prefix,
                         char * seq,
                         int len,
                         char * header,
                         int header_length,
                         unsigned int abundance,
                         int ordinal,
                         double ee,
                         int clustersize,
                         int clusterid,
                         const char * score_name,
                         double score) -> void;

auto fasta_print_db(std::FILE * output_handle,
                    uint64_t seqno) -> void;

auto fasta_print_db_relabel(std::FILE * output_handle,
                            uint64_t seqno,
                            int ordinal) -> void;

auto fasta_print_db_relabel(std::FILE * output_handle,
                            uint64_t seqno,
                            std::size_t ordinal) -> void;
