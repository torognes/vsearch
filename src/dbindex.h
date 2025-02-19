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

#include "bitmap.h"
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t


struct uhandle_s;

extern unsigned int * kmercount; /* number of matching seqnos for each kmer */
extern uint64_t * kmerhash;  /* index into the list below for each kmer */
extern unsigned int * kmerindex; /* the list of matching seqnos for kmers */
extern struct bitmap_s * * kmerbitmap;
extern unsigned int * dbindex_map;
extern unsigned int dbindex_count;
extern unsigned int kmerhashsize;
extern uint64_t kmerindexsize;
extern uhandle_s * dbindex_uh;


auto fprint_kmer(std::FILE * output_handle, unsigned int kmer_length, uint64_t kmer) -> void;

auto dbindex_prepare(int use_bitmap, int seqmask) -> void;

auto dbindex_addallsequences(int seqmask) -> void;

auto dbindex_addsequence(unsigned int seqno, int seqmask) -> void;

auto dbindex_free() -> void;

auto dbindex_getbitmap(unsigned int kmer) -> unsigned char *;

auto dbindex_getmatchcount(unsigned int kmer) -> unsigned int;

auto dbindex_getmatchlist(unsigned int kmer) -> unsigned int *;

auto dbindex_getmapping(unsigned int index) -> unsigned int;

auto dbindex_getcount() -> unsigned int;
