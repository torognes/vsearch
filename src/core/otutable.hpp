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

#pragma once

#ifdef HAVE_CONFIG_H
#include "config.h"  // HAVE_REGEX_H
#endif

#include "utils/view.hpp"  // View<char>
#include <cstdio>  // std::FILE
#include <cstdint>  // int64_t, uint64_t
#include <map>
#include <set>
#include <string>
#include <utility>  // std::pair

#ifdef HAVE_REGEX_H
#include <regex.h>  // regex_t
#endif


// Identify sample and OTU identifiers in sequence headers and accumulate the
// abundance of each sample in the different OTUs, then write the result as a
// classic OTU table, a mothur shared file, or a biom 1.0 document. A single
// instance owns the compiled matchers (RAII) and the accumulated counts. add()
// mutates shared state and is not thread-safe; its callers serialize access.

class OtuTable
{
public:
  OtuTable();
  ~OtuTable();
  OtuTable(OtuTable const &) = delete;
  auto operator=(OtuTable const &) -> OtuTable & = delete;

  auto add(View<char> query_header, View<char> target_header, int64_t abundance) -> void;
  auto print_otutabout(std::FILE * output_handle, struct Parameters const & parameters) const -> void;
  auto print_mothur_shared_out(std::FILE * output_handle, struct Parameters const & parameters) const -> void;
  auto print_biomout(std::FILE * output_handle, struct Parameters const & parameters) const -> void;

private:
  using string_set_t = std::set<std::string>;
  using string_pair_t = std::pair<std::string, std::string>;
  using string_pair_map_t = std::map<string_pair_t, uint64_t>;
  using otu_tax_map_t = std::map<std::string, std::string>;

#ifdef HAVE_REGEX_H
  regex_t regex_sample_ {};
  regex_t regex_otu_ {};
  regex_t regex_tax_ {};
#endif

  string_set_t otu_set_;
  string_set_t sample_set_;
  string_pair_map_t sample_otu_count_;
  string_pair_map_t otu_sample_count_;
  otu_tax_map_t otu_tax_map_;
};
