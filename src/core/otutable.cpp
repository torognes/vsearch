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

#include "otutable.hpp"
#include "utils/ascii_case.hpp"  // is_alnum
#include "utils/header_attribute.hpp"  // Attribute, header_find_attribute
#include "utils/view.hpp"
#include "vsearch.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/timestamp.hpp"  // iso8601_local_timestamp
#include "utils/prog_id.hpp"  // PROG_NAME, PROG_VERSION
#include <algorithm>  // std::find, std::find_if
#include <array>
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <iterator>  // std::distance, std::next
#include <map>
#include <string>
#include <utility>  // std::pair


/*

  Identify sample and otu identifiers in headers, and count
  abundance of the samples in different OTUs.

  http://www.drive5.com/usearch/manual/upp_labels_sample.html
  http://www.drive5.com/usearch/manual/upp_labels_otus.html

  TODO:
  - add relabel @

*/

using string_no_map_t = std::map<std::string, uint64_t>;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* The four annotations this file reads, previously three extended regular
     expressions compiled once per OtuTable. All four take everything up to the
     next ';' and all four accept an empty value, which is what the regular
     expressions did and what the taxonomy column depends on (an OTU whose
     target header ends in ";tax=" enters otu_tax_map_ with an empty string,
     and it is otu_tax_map_.empty() that decides whether --otutabout emits a
     taxonomy column at all). */
  constexpr auto any_value = vsearch::Value_chars::not_semicolon;

  constexpr std::array<vsearch::Attribute, 2> sample_names {{
      {"sample=", 7, any_value, true},
      {"barcodelabel=", 13, any_value, true}}};

  constexpr vsearch::Attribute otu_name_attribute {"otu=", 4, any_value, true};
  constexpr vsearch::Attribute tax_attribute {"tax=", 4, any_value, true};

  auto sample_alternatives() -> View<vsearch::Attribute>
  {
    return View<vsearch::Attribute>{sample_names.data(), sample_names.size()};
  }

}  // end of anonymous namespace


auto OtuTable::add(View<char> const query_header, View<char> const target_header, int64_t const abundance) -> void
{
  /* read sample annotation in query */

  bool const has_sample = (query_header.data() != nullptr);
  std::string sample_name;

  if (has_sample)
    {
      std::size_t len_sample = 0;
      char const * start_sample = query_header.data();
      /* the leftmost of the two names, not the first one listed: the extended
         regular expression this replaces alternated over the subject, so a
         header carrying both answers with whichever comes first in it */
      auto const match = vsearch::header_find_first_attribute(query_header,
                                                              sample_alternatives());
      if (match.span.present)
        {
          /* match: use the matching sample name */
          auto const value = vsearch::attribute_value(query_header, match.span,
                                                      *match.attribute);
          len_sample = value.size();
          start_sample = value.data();
        }
      else
        {
          /* no match: use first name in header with A-Za-z0-9_ */
          auto const * const first_other = std::find_if(query_header.begin(), query_header.end(),
              [](char const chr) -> bool { return (not is_alnum(chr)) and (chr != '_'); });
          len_sample = static_cast<std::size_t>(std::distance(query_header.begin(), first_other));
        }

      sample_name.assign(start_sample, len_sample);
    }


  /* read OTU annotation in target */

  bool const has_otu = (target_header.data() != nullptr);
  std::string otu_name;

  if (has_otu)
    {
      std::size_t len_otu = 0;
      char const * start_otu = target_header.data();
      auto const otu_annotation = vsearch::header_find_attribute(target_header,
                                                                 otu_name_attribute);
      if (otu_annotation.present)
        {
          /* match: use the matching otu name */
          auto const value = vsearch::attribute_value(target_header, otu_annotation,
                                                      otu_name_attribute);
          len_otu = value.size();
          start_otu = value.data();
        }
      else
        {
          /* no match: use first name in header up to ; */
          auto const * const first_semicolon = std::find(target_header.begin(), target_header.end(), ';');
          len_otu = static_cast<std::size_t>(std::distance(target_header.begin(), first_semicolon));
        }

      otu_name.assign(start_otu, len_otu);

      /* read tax annotation in target */

      auto const tax_annotation = vsearch::header_find_attribute(target_header,
                                                                 tax_attribute);
      if (tax_annotation.present)
        {
          /* match: use the matching tax name. An empty value is a match and is
             stored as an empty string, which is what puts the OTU in the map
             and so decides the taxonomy column. */
          auto const value = vsearch::attribute_value(target_header, tax_annotation,
                                                      tax_attribute);
          otu_tax_map_[otu_name] = std::string(value.data(), value.size());
        }
    }

  /* store data */

  if (has_sample) {
    sample_set_.insert(sample_name);
  }

  if (has_otu) {
    otu_set_.insert(otu_name);
  }

  if (has_sample and has_otu and (abundance != 0))
    {
      sample_otu_count_[string_pair_t(sample_name, otu_name)] += static_cast<uint64_t>(abundance);
      otu_sample_count_[string_pair_t(otu_name, sample_name)] += static_cast<uint64_t>(abundance);
    }

}


auto OtuTable::print_otutabout(std::FILE * output_handle, struct Parameters const & parameters) const -> void
{
  int64_t progress = 0;
  Progress progress_bar("Writing OTU table (classic)", otu_set_.size(), parameters);

  fprint(output_handle, "#OTU ID");
  for (auto const & it_sample : sample_set_)
    {
      fprint(output_handle, '\t');
      fprint(output_handle, make_view(it_sample));
    }
  if (not otu_tax_map_.empty())
    {
      fprint(output_handle, "\ttaxonomy");
    }
  fprint(output_handle, '\n');

  auto it_map = otu_sample_count_.begin();
  for (auto it_otu = otu_set_.begin();
       it_otu != otu_set_.end();
       ++it_otu)
    {
      std::fputs(it_otu->c_str(), output_handle);

      for (auto it_sample = sample_set_.begin();
           it_sample != sample_set_.end();
           ++it_sample)
        {
          uint64_t a = 0;
          if ((it_map != otu_sample_count_.end()) and
              (it_map->first.first == *it_otu) and
              (it_map->first.second == *it_sample))
            {
              a = it_map->second;
              ++it_map;
            }
          fprint(output_handle, '\t');
          fprint_integer(output_handle, a);
        }
      if (not otu_tax_map_.empty())
        {
          fprint(output_handle, '\t');
          auto it
            = otu_tax_map_.find(*it_otu);
          if (it != otu_tax_map_.end())
            {
              fprint(output_handle, make_view(it->second));
            }
        }
      fprint(output_handle, '\n');
      progress_bar.update(static_cast<uint64_t>(++progress));
    }
}


auto OtuTable::print_mothur_shared_out(std::FILE * output_handle, struct Parameters const & parameters) const -> void
{
  int64_t progress = 0;
  Progress progress_bar("Writing OTU table (mothur)", sample_set_.size(), parameters);

  fprint(output_handle, "label\tGroup\tnumOtus");
  int64_t numotus = 0;
  for (auto const & it_otu : otu_set_)
    {
      fprint(output_handle, '\t');
      fprint(output_handle, make_view(it_otu));
      ++numotus;
    }
  fprint(output_handle, '\n');

  auto it_map = sample_otu_count_.begin();

  for (auto it_sample = sample_set_.begin();
       it_sample != sample_set_.end();
       ++it_sample)
    {
      fprint(output_handle, "vsearch\t");
      std::fputs(it_sample->c_str(), output_handle);
      fprint(output_handle, '\t');
      fprint_integer(output_handle, numotus);

      for (auto it_otu = otu_set_.begin();
           it_otu != otu_set_.end();
           ++it_otu)
        {
          uint64_t a = 0;
          if ((it_map != sample_otu_count_.end()) and
              (it_map->first.first == *it_sample) and
              (it_map->first.second == *it_otu))
            {
              a = it_map->second;
              ++it_map;
            }
          fprint(output_handle, '\t');
          fprint_integer(output_handle, a);
        }

      fprint(output_handle, '\n');
      progress_bar.update(static_cast<uint64_t>(++progress));
    }
}


auto OtuTable::print_biomout(std::FILE * output_handle, struct Parameters const & parameters) const -> void
{
  int64_t progress = 0;
  Progress progress_bar("Writing OTU table (biom 1.0)", otu_sample_count_.size(), parameters);

  auto const rows = static_cast<int64_t>(otu_set_.size());
  auto const columns = static_cast<int64_t>(sample_set_.size());

  static std::string const date = iso8601_local_timestamp();

  fprint(output_handle, "{\n\t\"id\":\"");
  std::fputs(parameters.opt_biomout, output_handle);
  fprint(output_handle, "\",\n\t\"format\": \"Biological Observation Matrix 1.0\",\n\t\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\",\n\t\"type\": \"OTU table\",\n\t\"generated_by\": \"");
  std::fputs(PROG_NAME, output_handle);
  fprint(output_handle, ' ');
  std::fputs(PROG_VERSION, output_handle);
  fprint(output_handle, "\",\n\t\"date\": \"");
  fprint(output_handle, make_view(date));
  fprint(output_handle, "\",\n\t\"matrix_type\": \"sparse\",\n\t\"matrix_element_type\": \"int\",\n\t\"shape\": [");
  fprint_integer(output_handle, rows);
  fprint(output_handle, ',');
  fprint_integer(output_handle, columns);
  fprint(output_handle, "],\n");

  string_no_map_t otu_no_map;
  uint64_t otu_no = 0;

  fprint(output_handle, "\t\"rows\":[");
  for (auto it_otu = otu_set_.begin();
       it_otu != otu_set_.end();
       ++it_otu)
    {
      if (it_otu != otu_set_.begin())
        {
          fprint(output_handle, ',');
        }
      char const * otu_name = it_otu->c_str();
      fprint(output_handle, "\n\t\t{\"id\":\"");
      std::fputs(otu_name, output_handle);
      fprint(output_handle, "\", \"metadata\":");
      if (otu_tax_map_.empty())
        {
          fprint(output_handle, "null");
        }
      else
        {
          fprint(output_handle, R"({"taxonomy":")");
          auto it = otu_tax_map_.find(otu_name);
          if (it != otu_tax_map_.end())
            {
              fprint(output_handle, make_view(it->second));
            }
          fprint(output_handle, "\"}");
        }
      fprint(output_handle, '}');
      otu_no_map[*it_otu] = otu_no;
      ++otu_no;
    }
  fprint(output_handle, '\n');
  fprint(output_handle, "\t],\n");

  string_no_map_t sample_no_map;
  uint64_t sample_no = 0;

  fprint(output_handle, "\t\"columns\":[");
  for (auto it_sample = sample_set_.begin();
       it_sample != sample_set_.end();
       ++it_sample)
    {
      if (it_sample != sample_set_.begin())
        {
          fprint(output_handle, ',');
        }
      fprint(output_handle, "\n\t\t{\"id\":\"");
      std::fputs(it_sample->c_str(), output_handle);
      fprint(output_handle, "\", \"metadata\":null}");
      sample_no_map[*it_sample] = sample_no++;
    }
  fprint(output_handle, "\n\t],\n");

  auto first = true;
  fprint(output_handle, "\t\"data\": [");

  for (auto const & it_map : otu_sample_count_)
    {
      if (not first)
        {
          fprint(output_handle, ',');
        }

      otu_no = otu_no_map[it_map.first.first];
      sample_no = sample_no_map[it_map.first.second];

      fprint(output_handle, "\n\t\t[");
      fprint_integer(output_handle, otu_no);
      fprint(output_handle, ',');
      fprint_integer(output_handle, sample_no);
      fprint(output_handle, ',');
      fprint_integer(output_handle, it_map.second);
      fprint(output_handle, ']');
      first = false;
      ++progress;
      progress_bar.update(static_cast<uint64_t>(progress));
    }
  fprint(output_handle, "\n\t]\n");

  fprint(output_handle, "}\n");
}
