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

#include "vsearch.h"
#include <array>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <ctime>  // std::strftime, std::localtime, std::time, std::time_t, std::tm
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstring>  // std::strncpy, std::strcspn, std::strspn
#include <map>
#include <set>
#include <string>
#include <utility>  // std::pair
#include <vector>

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

using string_set_t = std::set<std::string>;
using string_pair_t = std::pair<std::string, std::string>;
using string_pair_map_t = std::map<string_pair_t, uint64_t>;
using otu_tax_map_t = std::map<std::string, std::string>;
using string_no_map_t = std::map<std::string, uint64_t>;

struct otutable_s
{
#ifdef HAVE_REGEX_H
  regex_t regex_sample;
  regex_t regex_otu;
  regex_t regex_tax;
#endif

  string_set_t otu_set;
  string_set_t sample_set;
  string_pair_map_t sample_otu_count;
  string_pair_map_t otu_sample_count;
  otu_tax_map_t otu_tax_map;
};

static otutable_s * otutable;


auto otutable_init() -> void
{
  otutable = new otutable_s;

#ifdef HAVE_REGEX_H
  /* compile regular expression matchers */
  if (regcomp(&otutable->regex_sample,
              "(^|;)(sample|barcodelabel)=([^;]*)($|;)",
              REG_EXTENDED) != 0)
    {
      fatal("Compilation of regular expression for sample annotation failed");
    }

  if (regcomp(&otutable->regex_otu,
              "(^|;)otu=([^;]*)($|;)",
              REG_EXTENDED) != 0)
    {
      fatal("Compilation of regular expression for otu annotation failed");
    }

  if (regcomp(&otutable->regex_tax,
              "(^|;)tax=([^;]*)($|;)",
              REG_EXTENDED) != 0)
    {
      fatal("Compilation of regular expression for taxonomy annotation failed");
    }
#endif
}


auto otutable_done() -> void
{
#ifdef HAVE_REGEX_H
  regfree(&otutable->regex_sample);
  regfree(&otutable->regex_otu);
  regfree(&otutable->regex_tax);
#endif

  otutable->otu_set.clear();
  otutable->sample_set.clear();
  otutable->sample_otu_count.clear();
  otutable->otu_sample_count.clear();
  delete otutable;
}


auto otutable_add(char * query_header, char * target_header, int64_t abundance) -> void
{
  /* read sample annotation in query */

  int len_sample = 0;
  char * start_sample = query_header;
  char * sample_name = nullptr;

  if (query_header != nullptr)
    {
#ifdef HAVE_REGEX_H
      std::array<regmatch_t, 5> pmatch_sample {{}};
      if (regexec(&otutable->regex_sample, query_header, 5, pmatch_sample.data(), 0) == 0)
        {
          /* match: use the matching sample name */
          len_sample = pmatch_sample[3].rm_eo - pmatch_sample[3].rm_so;
          start_sample += pmatch_sample[3].rm_so;
        }
#else
      std::cmatch cmatch_sample;
      if (std::regex_search(query_header, cmatch_sample, regex_sample))
        {
          len_sample = cmatch_sample.length(3);
          start_sample += cmatch_sample.position(3);
        }
#endif
      else
        {
          /* no match: use first name in header with A-Za-z0-9_ */
          len_sample = strspn(query_header,
                              "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                              "abcdefghijklmnopqrstuvwxyz"
                              "_"
                              "0123456789");
        }

      sample_name = (char *) xmalloc(len_sample + 1);
      std::strncpy(sample_name, start_sample, len_sample);
      sample_name[len_sample] = '\0';
    }


  /* read OTU annotation in target */

  int len_otu = 0;
  char * start_otu = target_header;
  char * otu_name = nullptr;

  if (target_header != nullptr)
    {
#ifdef HAVE_REGEX_H
      std::array<regmatch_t, 4> pmatch_otu {{}};
      if (regexec(&otutable->regex_otu, target_header, 4, pmatch_otu.data(), 0) == 0)
        {
          /* match: use the matching otu name */
          len_otu = pmatch_otu[2].rm_eo - pmatch_otu[2].rm_so;
          start_otu += pmatch_otu[2].rm_so;
        }
#else
      std::cmatch cmatch_otu;
      if (std::regex_search(target_header, cmatch_otu, regex_otu))
        {
          len_otu = cmatch_otu.length(2);
          start_otu += cmatch_otu.position(2);
        }
#endif
      else
        {
          /* no match: use first name in header up to ; */
          len_otu = strcspn(target_header, ";");
        }

      otu_name = (char *) xmalloc(len_otu + 1);
      std::strncpy(otu_name, start_otu, len_otu);
      otu_name[len_otu] = 0;

      /* read tax annotation in target */

#ifdef HAVE_REGEX_H
      char * start_tax = target_header;

      std::array<regmatch_t, 4> pmatch_tax {{}};
      if (regexec(&otutable->regex_tax, target_header, 4, pmatch_tax.data(), 0) == 0)
        {
          /* match: use the matching tax name */
          int const len_tax = pmatch_tax[2].rm_eo - pmatch_tax[2].rm_so;
          start_tax += pmatch_tax[2].rm_so;

          std::vector<char> tax_name(len_tax + 1);
          std::strncpy(tax_name.data(), start_tax, len_tax);
          tax_name[len_tax] = '\0';
          otutable->otu_tax_map[otu_name] = tax_name.data();
        }
#else
      std::cmatch cmatch_tax;
      if (std::regex_search(target_header, cmatch_tax, regex_tax))
        {
          otutable->otu_tax_map[otu_name] = cmatch_tax.str(2);
        }
#endif
    }

  /* store data */

  if (sample_name != nullptr) {
    otutable->sample_set.insert(sample_name);
  }

  if (otu_name != nullptr) {
    otutable->otu_set.insert(otu_name);
  }

  if ((sample_name != nullptr) && (otu_name != nullptr) && (abundance != 0))
    {
      otutable->sample_otu_count[string_pair_t(sample_name,otu_name)] += abundance;
      otutable->otu_sample_count[string_pair_t(otu_name,sample_name)] += abundance;
    }

  if (otu_name != nullptr) {
    xfree(otu_name);
  }

  if (sample_name != nullptr) {
    xfree(sample_name);
  }
}


auto otutable_print_otutabout(std::FILE * output_handle) -> void
{
  int64_t progress = 0;
  progress_init("Writing OTU table (classic)", otutable->otu_set.size());

  fprintf(output_handle, "#OTU ID");
  for (auto const & it_sample : otutable->sample_set)
    {
      fprintf(output_handle, "\t%s", it_sample.c_str());
    }
  if (! otutable->otu_tax_map.empty())
    {
      fprintf(output_handle, "\ttaxonomy");
    }
  fprintf(output_handle, "\n");

  auto it_map = otutable->otu_sample_count.begin();
  for (auto it_otu = otutable->otu_set.begin();
       it_otu != otutable->otu_set.end();
       ++it_otu)
    {
      fprintf(output_handle, "%s", it_otu->c_str());

      for (auto it_sample = otutable->sample_set.begin();
           it_sample != otutable->sample_set.end();
           ++it_sample)
        {
          uint64_t a = 0;
          if ((it_map != otutable->otu_sample_count.end()) &&
              (it_map->first.first == *it_otu) &&
              (it_map->first.second == *it_sample))
            {
              a = it_map->second;
              ++it_map;
            }
          fprintf(output_handle, "\t%" PRIu64, a);
        }
      if (! otutable->otu_tax_map.empty())
        {
          fprintf(output_handle, "\t");
          auto it
            = otutable->otu_tax_map.find(*it_otu);
          if (it != otutable->otu_tax_map.end())
            {
              fprintf(output_handle, "%s", it->second.c_str());
            }
        }
      fprintf(output_handle, "\n");
      progress_update(++progress);
    }
  progress_done();
}


auto otutable_print_mothur_shared_out(std::FILE * output_handle) -> void
{
  int64_t progress = 0;
  progress_init("Writing OTU table (mothur)", otutable->sample_set.size());

  fprintf(output_handle, "label\tGroup\tnumOtus");
  int64_t numotus = 0;
  for (const auto & it_otu : otutable->otu_set)
    {
      const char * otu_name = it_otu.c_str();
      fprintf(output_handle, "\t%s", otu_name);
      ++numotus;
    }
  fprintf(output_handle, "\n");

  auto it_map = otutable->sample_otu_count.begin();

  for (auto it_sample = otutable->sample_set.begin();
       it_sample != otutable->sample_set.end();
       ++it_sample)
    {
      fprintf(output_handle, "vsearch\t%s\t%" PRId64, it_sample->c_str(), numotus);

      for (auto it_otu = otutable->otu_set.begin();
           it_otu != otutable->otu_set.end();
           ++it_otu)
        {
          uint64_t a = 0;
          if ((it_map != otutable->sample_otu_count.end()) &&
              (it_map->first.first == *it_sample) &&
              (it_map->first.second == *it_otu))
            {
              a = it_map->second;
              ++it_map;
            }
          fprintf(output_handle, "\t%" PRIu64, a);
        }

      fprintf(output_handle, "\n");
      progress_update(++progress);
    }
  progress_done();
}


auto otutable_print_biomout(std::FILE * output_handle) -> void
{
  int64_t progress = 0;
  progress_init("Writing OTU table (biom 1.0)", otutable->otu_sample_count.size());

  int64_t const rows = otutable->otu_set.size();
  int64_t const columns = otutable->sample_set.size();

  static const time_t time_now = time(nullptr);
  struct tm * tm_now = localtime(& time_now);
  std::array<char, 50> date {{}};
  strftime(date.data(), 50, "%Y-%m-%dT%H:%M:%S", tm_now);

  fprintf(output_handle,
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
          opt_biomout,
          PROG_NAME, PROG_VERSION,
          date.data(),
          rows,
          columns);

  string_no_map_t otu_no_map;
  uint64_t otu_no = 0;

  fprintf(output_handle, "\t\"rows\":[");
  for (auto it_otu = otutable->otu_set.begin();
       it_otu != otutable->otu_set.end();
       ++it_otu)
    {
      if (it_otu != otutable->otu_set.begin())
        {
          fprintf(output_handle, ",");
        }
      const char * otu_name = it_otu->c_str();
      fprintf(output_handle, "\n\t\t{\"id\":\"%s\", \"metadata\":", otu_name);
      if (otutable->otu_tax_map.empty())
        {
          fprintf(output_handle, "null");
        }
      else
        {
          fprintf(output_handle, R"({"taxonomy":")");
          auto it = otutable->otu_tax_map.find(otu_name);
          if (it != otutable->otu_tax_map.end())
            {
              fprintf(output_handle, "%s", it->second.c_str());
            }
          fprintf(output_handle, "\"}");
        }
      fprintf(output_handle, "}");
      otu_no_map[*it_otu] = otu_no++;
    }
  fprintf(output_handle, "\n");
  fprintf(output_handle, "\t],\n");

  string_no_map_t sample_no_map;
  uint64_t sample_no = 0;

  fprintf(output_handle, "\t\"columns\":[");
  for (auto it_sample = otutable->sample_set.begin();
       it_sample != otutable->sample_set.end();
       ++it_sample)
    {
      if (it_sample != otutable->sample_set.begin())
        {
          fprintf(output_handle, ",");
        }
      fprintf(output_handle, "\n\t\t{\"id\":\"%s\", \"metadata\":null}", it_sample->c_str());
      sample_no_map[*it_sample] = sample_no++;
    }
  fprintf(output_handle, "\n\t],\n");

  bool first = true;
  fprintf(output_handle, "\t\"data\": [");

  for (auto & it_map : otutable->otu_sample_count)
    {
      if (! first)
        {
          fprintf(output_handle, ",");
        }

      otu_no = otu_no_map[it_map.first.first];
      sample_no = sample_no_map[it_map.first.second];

      fprintf(output_handle, "\n\t\t[%" PRIu64 ",%" PRIu64 ",%" PRIu64 "]", otu_no, sample_no, it_map.second);
      first = false;
      progress_update(++progress);
    }
  fprintf(output_handle, "\n\t]\n");

  fprintf(output_handle, "}\n");
  progress_done();
}
