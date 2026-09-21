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


#include <algorithm>  // std::find_if, std::min
#include <string>  // std::string


/* The sample identifier --relabel @ derives from the name of the input file,
   for compatibility with usearch (torognes/vsearch#202).

   usearch documents the rule as "truncating at the first underscore (_) or
   period (.)", which reads as "whichever comes first". Measured against
   usearch 11.0.667 over twenty file names, that is not what it does: the
   underscore has absolute precedence and the period is only a fallback, so
   'a.b_c.fastq' yields 'a.b' and not 'a'. The two readings agree on every
   name in which no period precedes the first underscore, which is why the
   documented one is easy to believe.

   The base name is taken first. Cutting before that would turn
   '/home/user/my_data/sampleA.fastq' into '/home/user/my'.

   The last step is vsearch's own, not usearch's: the identifier is truncated
   again at the first ';' or blank, the same character set --sample truncates
   its argument at (cli.cc). usearch leaves those in, which yields '>my
   sample.1' for 'my sample.fastq' -- a header that splits at the space for
   any reader not given --notrunclabels -- and '>a;b.1' for 'a;b.fastq', whose
   ';' collides with the separator that introduces ;size=, ;sample= and
   ;length=. Since --fastq_mergepairs gained --tabbedout the identifier also
   reaches a tab-separated field, where a blank makes 'relabel=my sample.1'
   ambiguous to a parser. The step only ever fires on names for which
   usearch's own label is malformed.

   It must be applied AFTER the cut, not by widening the cut's fallback set:
   the underscore branch never consults that set, so widening it would leave
   'my sample_R1.fastq' as 'my sample'.

   Derived from the command-line string alone and never from the opened
   stream, so the result does not depend on what the path resolves to. A named
   pipe called 'runA_1.fifo' is a stream and still yields 'runA', which is the
   reason the rule is spelled on the name rather than on a stat() of the
   handle. */
namespace vsearch
{

  namespace sample_identifier_detail
  {
    /* Everything that ends a directory prefix. Windows accepts '/' as well as
       '\\', and ':' ends a drive letter ('C:file.fastq' names a file in the
       current directory of drive C), so all three end the prefix there. */
    constexpr auto is_path_separator(char const character) -> bool
    {
#ifdef _WIN32
      return (character == '/') or (character == '\\') or (character == ':');
#else
      return (character == '/');
#endif
    }

    /* The characters --sample truncates its argument at: ';' plus the six
       standard blanks. A blank splits the header for any reader not given
       --notrunclabels, and ';' introduces an annotation. */
    constexpr auto is_label_terminator(char const character) -> bool
    {
      return (character == ';') or (character == ' ') or (character == '\t') or
        (character == '\n') or (character == '\r') or (character == '\v') or
        (character == '\f');
    }
  }  // namespace sample_identifier_detail


  /* Not noexcept: it returns a std::string, whose allocation may throw. */
  inline auto sample_identifier(std::string const & path) -> std::string
  {
    auto const after_last_separator =
      std::find_if(path.rbegin(), path.rend(),
                   sample_identifier_detail::is_path_separator).base();
    std::string identifier {after_last_separator, path.end()};

    /* usearch's cut: the first underscore if there is one, the first period
       otherwise, and the whole base name when it has neither. */
    auto cut = identifier.find('_');
    if (cut == std::string::npos)
      {
        cut = identifier.find('.');
      }
    identifier.resize(std::min(cut, identifier.size()));

    /* vsearch's own hygiene, applied to the result of that cut */
    auto const terminator =
      std::find_if(identifier.begin(), identifier.end(),
                   sample_identifier_detail::is_label_terminator);
    identifier.erase(terminator, identifier.end());

    return identifier;
  }

}  // namespace vsearch
