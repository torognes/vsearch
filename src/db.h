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

struct seqinfo_s
{
  size_t header_p;
  size_t seq_p;
  size_t qual_p;
  unsigned int headerlen;
  unsigned int seqlen;
  unsigned int size;
};

typedef struct seqinfo_s seqinfo_t;

extern char * datap;
extern seqinfo_t * seqindex;

inline char * db_getheader(uint64_t seqno)
{
  return datap + seqindex[seqno].header_p;
}

inline char * db_getsequence(uint64_t seqno)
{
  return datap + seqindex[seqno].seq_p;
}

inline uint64_t db_getabundance(uint64_t seqno)
{
  return seqindex[seqno].size;
}

inline uint64_t db_getsequencelen(uint64_t seqno)
{
  return seqindex[seqno].seqlen;
}

inline uint64_t db_getheaderlen(uint64_t seqno)
{
  return seqindex[seqno].headerlen;
}

void db_read(const char * filename, int upcase);
void db_free();

uint64_t db_getsequencecount();
uint64_t db_getnucleotidecount();
uint64_t db_getlongestheader();
uint64_t db_getlongestsequence();
uint64_t db_getshortestsequence();

/* Note: the sorting functions below must be called after db_read,
   but before dbindex_prepare */

void db_sortbylength();
void db_sortbylength_shortest_first();

void db_sortbyabundance();

bool db_is_fastq();
char * db_getquality(uint64_t seqno);

void db_setinfo(bool new_is_fastq,
                uint64_t new_sequences,
                uint64_t new_nucleotides,
                uint64_t new_longest,
                uint64_t new_shortest,
                uint64_t new_longestheader);
