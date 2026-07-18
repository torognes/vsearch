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
#include "vsearch.hpp"
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/timestamp.hpp"  // iso8601_local_timestamp
#include "utils/prog_id.hpp"  // PROG_NAME, PROG_VERSION
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstring>  // std::strcspn, std::strspn
#include <iterator>  // std::next
#include <map>
#include <set>
#include <string>
#include <utility>  // std::pair

// refactoring: is there a reason to prefer regex.h over <regex>?
#ifdef HAVE_REGEX_H
#include <regex.h>  // C: regcomp, regexec, regfree (POSIX functions)
#else
#include <regex>  // C++: std::regex_search
#endif


/*

  Identify sample and otu identifiers in headers, and count
  abundance of the samples in different OTUs.

  http://www.drive5.com/usearch/manual/upp_labels_sample.html
  http://www.drive5.com/usearch/manual/upp_labels_otus.html

  TODO:
  - add relabel @

*/

#ifndef HAVE_REGEX_H
const std::regex regex_sample("(^|;)(sample|barcodelabel)=([^;]*)($|;)",
                              std::regex::extended);
const std::regex regex_otu("(^|;)otu=([^;]*)($|;)",
                           std::regex::extended);
const std::regex regex_tax("(^|;)tax=([^;]*)($|;)",
                           std::regex::extended);
#endif

using string_no_map_t = std::map<std::string, uint64_t>;


OtuTable::OtuTable()
{
#ifdef HAVE_REGEX_H
  /* compile regular expression matchers */
  if (regcomp(&regex_sample_,
              "(^|;)(sample|barcodelabel)=([^;]*)($|;)",
              REG_EXTENDED) != 0)
    {
      fatal("Compilation of regular expression for sample annotation failed");
    }

  if (regcomp(&regex_otu_,
              "(^|;)otu=([^;]*)($|;)",
              REG_EXTENDED) != 0)
    {
      fatal("Compilation of regular expression for otu annotation failed");
    }

  if (regcomp(&regex_tax_,
              "(^|;)tax=([^;]*)($|;)",
              REG_EXTENDED) != 0)
    {
      fatal("Compilation of regular expression for taxonomy annotation failed");
    }
#endif
}


OtuTable::~OtuTable()
{
#ifdef HAVE_REGEX_H
  regfree(&regex_sample_);
  regfree(&regex_otu_);
  regfree(&regex_tax_);
#endif
}


auto OtuTable::add(char const * query_header, char const * target_header, int64_t const abundance) -> void
{
  /* read sample annotation in query */

  bool const has_sample = (query_header != nullptr);
  std::string sample_name;

  if (has_sample)
    {
      std::size_t len_sample = 0;
      char const * start_sample = query_header;
#ifdef HAVE_REGEX_H
      std::array<regmatch_t, 5> pmatch_sample {{}};
      if (regexec(&regex_sample_, query_header, 5, pmatch_sample.data(), 0) == 0)
        {
          /* match: use the matching sample name */
          len_sample = static_cast<std::size_t>(pmatch_sample[3].rm_eo - pmatch_sample[3].rm_so);
          start_sample = std::next(query_header, pmatch_sample[3].rm_so);
        }
#else
      std::cmatch cmatch_sample;
      if (std::regex_search(query_header, cmatch_sample, regex_sample))
        {
          len_sample = static_cast<std::size_t>(cmatch_sample.length(3));
          start_sample = std::next(query_header, cmatch_sample.position(3));
        }
#endif
      else
        {
          /* no match: use first name in header with A-Za-z0-9_ */
          len_sample = std::strspn(query_header,
                              "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                              "abcdefghijklmnopqrstuvwxyz"
                              "_"
                              "0123456789");
        }

      sample_name.assign(start_sample, len_sample);
    }


  /* read OTU annotation in target */

  bool const has_otu = (target_header != nullptr);
  std::string otu_name;

  if (has_otu)
    {
      std::size_t len_otu = 0;
      char const * start_otu = target_header;
#ifdef HAVE_REGEX_H
      std::array<regmatch_t, 4> pmatch_otu {{}};
      if (regexec(&regex_otu_, target_header, 4, pmatch_otu.data(), 0) == 0)
        {
          /* match: use the matching otu name */
          len_otu = static_cast<std::size_t>(pmatch_otu[2].rm_eo - pmatch_otu[2].rm_so);
          start_otu = std::next(target_header, pmatch_otu[2].rm_so);
        }
#else
      std::cmatch cmatch_otu;
      if (std::regex_search(target_header, cmatch_otu, regex_otu))
        {
          len_otu = static_cast<std::size_t>(cmatch_otu.length(2));
          start_otu = std::next(target_header, cmatch_otu.position(2));
        }
#endif
      else
        {
          /* no match: use first name in header up to ; */
          len_otu = std::strcspn(target_header, ";");
        }

      otu_name.assign(start_otu, len_otu);

      /* read tax annotation in target */

#ifdef HAVE_REGEX_H
      std::array<regmatch_t, 4> pmatch_tax {{}};
      if (regexec(&regex_tax_, target_header, 4, pmatch_tax.data(), 0) == 0)
        {
          /* match: use the matching tax name */
          std::size_t const len_tax = static_cast<std::size_t>(pmatch_tax[2].rm_eo - pmatch_tax[2].rm_so);
          char const * const start_tax = std::next(target_header, pmatch_tax[2].rm_so);
          otu_tax_map_[otu_name] = std::string(start_tax, len_tax);
        }
#else
      std::cmatch cmatch_tax;
      if (std::regex_search(target_header, cmatch_tax, regex_tax))
        {
          otu_tax_map_[otu_name] = cmatch_tax.str(2);
        }
#endif
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

  std::fprintf(output_handle, "#OTU ID");
  for (auto const & it_sample : sample_set_)
    {
      std::fprintf(output_handle, "\t%s", it_sample.c_str());
    }
  if (not otu_tax_map_.empty())
    {
      std::fprintf(output_handle, "\ttaxonomy");
    }
  std::fprintf(output_handle, "\n");

  auto it_map = otu_sample_count_.begin();
  for (auto it_otu = otu_set_.begin();
       it_otu != otu_set_.end();
       ++it_otu)
    {
      std::fprintf(output_handle, "%s", it_otu->c_str());

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
          std::fprintf(output_handle, "\t%" PRIu64, a);
        }
      if (not otu_tax_map_.empty())
        {
          std::fprintf(output_handle, "\t");
          auto it
            = otu_tax_map_.find(*it_otu);
          if (it != otu_tax_map_.end())
            {
              std::fprintf(output_handle, "%s", it->second.c_str());
            }
        }
      std::fprintf(output_handle, "\n");
      progress_bar.update(static_cast<uint64_t>(++progress));
    }
}


auto OtuTable::print_mothur_shared_out(std::FILE * output_handle, struct Parameters const & parameters) const -> void
{
  int64_t progress = 0;
  Progress progress_bar("Writing OTU table (mothur)", sample_set_.size(), parameters);

  std::fprintf(output_handle, "label\tGroup\tnumOtus");
  int64_t numotus = 0;
  for (auto const & it_otu : otu_set_)
    {
      char const * otu_name = it_otu.c_str();
      std::fprintf(output_handle, "\t%s", otu_name);
      ++numotus;
    }
  std::fprintf(output_handle, "\n");

  auto it_map = sample_otu_count_.begin();

  for (auto it_sample = sample_set_.begin();
       it_sample != sample_set_.end();
       ++it_sample)
    {
      std::fprintf(output_handle, "vsearch\t%s\t%" PRId64, it_sample->c_str(), numotus);

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
          std::fprintf(output_handle, "\t%" PRIu64, a);
        }

      std::fprintf(output_handle, "\n");
      progress_bar.update(static_cast<uint64_t>(++progress));
    }
}


auto OtuTable::print_biomout(std::FILE * output_handle, struct Parameters const & parameters) const -> void
{
  int64_t progress = 0;
  Progress progress_bar("Writing OTU table (biom 1.0)", otu_sample_count_.size(), parameters);

  int64_t const rows = static_cast<int64_t>(otu_set_.size());
  int64_t const columns = static_cast<int64_t>(sample_set_.size());

  static std::string const date = iso8601_local_timestamp();

  std::fprintf(output_handle,
          "{\n"
          "\t\"id\":\"%s\",\n"
          "\t\"format\": \"Biological Observation Matrix 1.0\",\n"
          "\t\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\",\n"
          "\t\"type\": \"OTU table\",\n"
          "\t\"generated_by\": \"%s %s\",\n"
          "\t\"date\": \"%s\",\n"
          "\t\"matrix_type\": \"sparse\",\n"
          "\t\"matrix_element_type\": \"int\",\n"
          "\t\"shape\": [%" PRId64 ",%" PRId64 "],\n",
          parameters.opt_biomout,
          PROG_NAME, PROG_VERSION,
          date.c_str(),
          rows,
          columns);

  string_no_map_t otu_no_map;
  uint64_t otu_no = 0;

  std::fprintf(output_handle, "\t\"rows\":[");
  for (auto it_otu = otu_set_.begin();
       it_otu != otu_set_.end();
       ++it_otu)
    {
      if (it_otu != otu_set_.begin())
        {
          std::fprintf(output_handle, ",");
        }
      char const * otu_name = it_otu->c_str();
      std::fprintf(output_handle, "\n\t\t{\"id\":\"%s\", \"metadata\":", otu_name);
      if (otu_tax_map_.empty())
        {
          std::fprintf(output_handle, "null");
        }
      else
        {
          std::fprintf(output_handle, R"({"taxonomy":")");
          auto it = otu_tax_map_.find(otu_name);
          if (it != otu_tax_map_.end())
            {
              fprintf(output_handle, "%s", it->second.c_str());
            }
          fprintf(output_handle, "\"}");
        }
      std::fprintf(output_handle, "}");
      otu_no_map[*it_otu] = otu_no;
      ++otu_no;
    }
  std::fprintf(output_handle, "\n");
  std::fprintf(output_handle, "\t],\n");

  string_no_map_t sample_no_map;
  uint64_t sample_no = 0;

  std::fprintf(output_handle, "\t\"columns\":[");
  for (auto it_sample = sample_set_.begin();
       it_sample != sample_set_.end();
       ++it_sample)
    {
      if (it_sample != sample_set_.begin())
        {
          std::fprintf(output_handle, ",");
        }
      std::fprintf(output_handle, "\n\t\t{\"id\":\"%s\", \"metadata\":null}", it_sample->c_str());
      sample_no_map[*it_sample] = sample_no++;
    }
  std::fprintf(output_handle, "\n\t],\n");

  auto first = true;
  std::fprintf(output_handle, "\t\"data\": [");

  for (auto const & it_map : otu_sample_count_)
    {
      if (not first)
        {
          std::fprintf(output_handle, ",");
        }

      otu_no = otu_no_map[it_map.first.first];
      sample_no = sample_no_map[it_map.first.second];

      std::fprintf(output_handle, "\n\t\t[%" PRIu64 ",%" PRIu64 ",%" PRIu64 "]", otu_no, sample_no, it_map.second);
      first = false;
      ++progress;
      progress_bar.update(static_cast<uint64_t>(progress));
    }
  std::fprintf(output_handle, "\n\t]\n");

  std::fprintf(output_handle, "}\n");
}
