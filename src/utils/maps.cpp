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

#include "maps.hpp"
#include <vector>


const std::vector<unsigned char> chrmap_no_change_vector = {
    /*

      Map from ascii to ascii - no change

     @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
     P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',

    'N', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
    'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y',
    'Z', 'N', 'N', 'N', 'N', 'N', 'N', 'a', 'b', 'c', 'd', 'e', 'f',
    'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's',
    't', 'u', 'v', 'w', 'x', 'y', 'z', 'N', 'N', 'N', 'N', 'N',

    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',

    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};


const std::vector<unsigned char> chrmap_4bit_vector =
  {
    /*
      Map from ascii to 4-bit nucleotide code

      Aa:  1    0001
      Bb: 14    1110   ex: 'B' & 'A' == 0000 while 'B' & anyother != 0000
      Cc:  2    0010
      Dd: 13    1101
      Gg:  4    0100
      Hh: 11    1011
      Kk: 12    1100
      Mm:  3    0011
      Nn: 15    1111   ex: 'N' & any != 0000
      Rr:  5    0101
      Ss:  6    0110
      Tt:  8    1000
      Uu:  8    1000
      Vv:  7    0111
      Ww:  9    1001
      Yy: 10    1010   ex: 'Y' & 'C' or 'T' == 0000
      Others: 0

     @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
     P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  1, 14,  2, 13,  0,  0,  4, 11,  0,  0, 12,  0,  3, 15,  0,
     0,  0,  5,  6,  8,  8,  7,  9,  0, 10,  0,  0,  0,  0,  0,  0,
     0,  1, 14,  2, 13,  0,  0,  4, 11,  0,  0, 12,  0,  3, 15,  0,
     0,  0,  5,  6,  8,  8,  7,  9,  0, 10,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
  };


const std::vector<unsigned char> chrmap_complement_vector =
  {
    /*

      Map from ascii to ascii, complementary nucleotide

     @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
     P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','T','V','G','H','N','N','C','D','N','N','M','N','K','N','N',
    'N','N','Y','S','A','A','B','W','N','R','N','N','N','N','N','N',
    'N','t','v','g','h','N','N','c','d','N','N','m','N','k','n','N',
    'N','N','y','s','a','a','b','w','N','r','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
  };


const std::vector<unsigned char> chrmap_normalize_vector =
  {
    /*

      Map from ascii to ascii
      Convert to upper case nucleotide, and replace U by T

     @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
     P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','A','B','C','D','N','N','G','H','N','N','K','N','M','N','N',
    'N','N','R','S','T','T','V','W','N','Y','N','N','N','N','N','N',
    'N','A','B','C','D','N','N','G','H','N','N','K','N','M','N','N',
    'N','N','R','S','T','T','V','W','N','Y','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
  };


const std::vector<unsigned char> chrmap_upcase_vector =
  {
    /*

      Map from ascii to ascii
      Convert to upper case nucleotide

     @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
     P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
    */

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
    'P','Q','R','S','T','U','V','W','X','Y','Z','N','N','N','N','N',
    'N','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
    'P','Q','R','S','T','U','V','W','X','Y','Z','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
  };
