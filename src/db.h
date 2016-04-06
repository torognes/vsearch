/*

  VSEARCH5D: a modified version of VSEARCH

  Copyright (C) 2016, Akifumi S. Tanabe

  Contact: Akifumi S. Tanabe
  https://github.com/astanabe/vsearch5d

  Original version of VSEARCH
  Copyright (C) 2014-2015, Torbjorn Rognes, Frederic Mahe and Tomas Flouri

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

inline char * db_getheader(unsigned long seqno)
{
  return datap + seqindex[seqno].header_p;
}

inline char * db_getsequence(unsigned long seqno)
{
  return datap + seqindex[seqno].seq_p;
}

inline unsigned long db_getabundance(unsigned long seqno)
{
  return seqindex[seqno].size;
}

inline unsigned long db_getsequencelen(unsigned long seqno)
{
  return seqindex[seqno].seqlen;
}

inline unsigned long db_getheaderlen(unsigned long seqno)
{
  return seqindex[seqno].headerlen;
}

void db_read(const char * filename, int upcase);
void db_free();

unsigned long db_getsequencecount();
unsigned long db_getnucleotidecount();
unsigned long db_getlongestheader();
unsigned long db_getlongestsequence();
unsigned long db_getshortestsequence();

/* Note: the sorting functions below must be called after db_read,
   but before dbindex_prepare */

void db_sortbylength();
void db_sortbylength_shortest_first();

void db_sortbyabundance();

bool db_is_fastq();
char * db_getquality(unsigned long seqno);
