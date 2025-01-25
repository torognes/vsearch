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
#include "utils/taxonomic_fields.h"
#include <cctype>  // std::tolower
#include <cstring>  // std::strlen, std::strstr, std::strchr
#include <iterator>  // std::distance


auto const * tax_letters = "dkpcofgst";


// very similar to header_find_attribute() in attributes.cc
auto tax_parse(char const * header,
               int const header_length,
               int * tax_start,
               int * tax_end) -> bool
{
  /*
    Identify the first occurence of the pattern (^|;)tax=([^;]*)(;|$)
  */

  if (header == nullptr)
    {
      return false;
    }

  auto const * attribute = "tax=";

  auto const attribute_length = static_cast<int>(std::strlen(attribute));

  auto offset = 0;

  while (offset < header_length - attribute_length)
    {
      auto const * first_occurence = std::strstr(header + offset, attribute);

      /* no match */
      if (first_occurence == nullptr)
        {
          break;
        }

      offset = std::distance(header, first_occurence);

      /* check for ';' in front */
      if ((offset > 0) and (header[offset - 1] != ';'))
        {
          offset += attribute_length + 1;
          continue;
        }

      * tax_start = offset;

      /* find end (semicolon or end of header) */
      auto const * terminus = std::strchr(header + offset + attribute_length, ';');
      if (terminus == nullptr)
        {
          * tax_end = header_length;
        }
      else
        {
          * tax_end = std::distance(header, terminus);
        }

      return true;
    }
  return false;
}


auto tax_split(int seqno, int * level_start, int * level_len) -> void
{
  /* Parse taxonomy string into the following 9 parts
     d domain
     k kingdom
     p phylum
     c class
     o order
     f family
     g genus
     s species
     t strain
  */
  static constexpr auto length_of_attribute_name = 4;  // "tax=" -> 4 letters

  for (auto i = 0; i < tax_levels; ++i)
    {
      level_start[i] = 0;
      level_len[i] = 0;
    }

  auto tax_start = 0;
  auto tax_end = 0;
  auto const * const header = db_getheader(seqno);
  int const header_length = db_getheaderlen(seqno);
  auto const attribute_is_present = tax_parse(header, header_length, & tax_start, & tax_end);
  if (not attribute_is_present) { return; }
  auto offset = tax_start + length_of_attribute_name;

  while (offset < tax_end)
    {
      /* Is the next char a recognized tax level letter? */
      auto const * next_level = std::strchr(taxonomic_fields.data(), std::tolower(header[offset]));
      if (next_level != nullptr)
        {
          int const level = std::distance(taxonomic_fields.data(), next_level);

          /* Is there a colon after it? */
          if (header[offset + 1] == ':')
            {
              level_start[level] = offset + 2;

              auto const * next_comma = std::strchr(header + offset + 2, ',');
              if (next_comma != nullptr)
                {
                  level_len[level] = std::distance(header, next_comma) - offset - 2;
                }
              else
                {
                  level_len[level] = tax_end - offset - 2;
                }
            }
        }

      /* skip past next comma */
      auto const * next_comma_bis = std::strchr(header + offset, ',');
      if (next_comma_bis != nullptr)
        {
          offset = std::distance(header, next_comma_bis) + 1;
        }
      else
        {
          offset = tax_end;
        }
    }
}
