/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "core/db.hpp"  // Database
#include "core/tax.hpp"  // TaxLevel, tax_split
#include "utils/ascii_case.hpp"  // to_lower
#include "utils/header_attribute.hpp"  // Attribute, header_find_attribute
#include "utils/taxonomic_fields.h"
#include "utils/view.hpp"  // View
#include <algorithm>  // std::find
#include <array>  // std::array
#include <cstddef>  // std::size_t
#include <cstdint>
#include <iterator>  // std::distance


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* the (^|;)tax=([^;]*)(;|$) scan lives in utils/header_attribute.hpp, shared
     with core/attributes.cpp and core/otutable.cpp. A near-copy of it used to
     sit here, and the copy could not see a "tax=" whose '=' was the header's
     last byte. */
  constexpr vsearch::Attribute tax_attribute {"tax=", 4,
                                              vsearch::Value_chars::not_semicolon,
                                              true};

}  // anonymous namespace


auto tax_split(int const seqno, std::array<TaxLevel, tax_levels> & levels,
               struct Database const & db) -> void
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
  auto const header = db.header_view(static_cast<uint64_t>(seqno));
  auto const annotation = vsearch::header_find_attribute(header, tax_attribute);
  if (not annotation.present) { return; }

  /* absolute indices into the header, not the value view: levels[].start is an
     offset into the header the caller re-reads later (see TaxLevel in
     core/tax.hpp), so the span is what this function wants */
  auto const tax_start = static_cast<int>(annotation.start);
  auto const tax_end = static_cast<int>(annotation.end);
  auto offset = tax_start + tax_attribute.length;

  /* taxon names cannot contain commas or semicolons (see the manual), so
     the comma searches below must stop at the end of the tax= field: an
     unbounded search used to extend the last taxon name across the closing
     ';' when a later header field contained a comma */
  auto const * const field_end = header.begin() + tax_end;

  while (offset < tax_end)
    {
      /* Is the next char a recognized tax level letter? */
      auto const * next_level = std::find(taxonomic_fields.begin(), taxonomic_fields.end(),
                                          to_lower(header[static_cast<std::size_t>(offset)]));
      if (next_level != taxonomic_fields.end())
        {
          auto const level = static_cast<std::size_t>(std::distance(taxonomic_fields.begin(), next_level));

          /* Is there a colon after it? */
          /* the level letter can be the header's last byte (">s;tax=d"), so
             the position after it is outside the view. Reading it used to work
             only because Database NUL-terminates its headers, and a '\0' is
             not a ':' -- the bounds test makes that explicit and keeps the
             same answer. */
          auto const after_level = static_cast<std::size_t>(offset) + 1;
          if ((after_level < header.size()) and (header[after_level] == ':'))
            {
              levels[level].start = offset + 2;

              auto const * const next_comma = std::find(header.begin() + offset + 2, field_end, ',');
              if (next_comma != field_end)
                {
                  levels[level].length = static_cast<int>(std::distance(header.begin(), next_comma)) - offset - 2;
                }
              else
                {
                  levels[level].length = tax_end - offset - 2;
                }
            }
        }

      /* skip past next comma */
      auto const * const next_comma_bis = std::find(header.begin() + offset, field_end, ',');
      if (next_comma_bis != field_end)
        {
          offset = static_cast<int>(std::distance(header.begin(), next_comma_bis)) + 1;
        }
      else
        {
          offset = tax_end;
        }
    }
}
