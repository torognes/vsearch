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

#include "utils/view.hpp"
#include "core/seq_record.hpp"
#include "vsearch.hpp"
#include "core/attributes.hpp"
#include "core/db.hpp"
#include "core/fastx.hpp"
#include "core/fastx_char_class.hpp"  // vsearch::CharClass, class_of
#include "utils/fatal.hpp"
#include "utils/maps.hpp"  // Mapping, map_accepted_base, chrmap_*
#include "utils/print_view.hpp"  // fprint
#include <algorithm>  // std::min
#include <array>
#include <cassert>  // assert
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::size_t
#include <iterator>  // std::next
#include <memory>  // std::unique_ptr


// anonymous namespace: limit visibility and usage to this translation unit
namespace {


  enum struct Action : unsigned char {
    warn,    // (0) symbol is stripped, with a warning
    accept,  // (1)
    reject,  // (2) fatal printable symbol ('.', '-')
    show,    // (3) fatal non-printable symbol (0-8, 14-31, but not 127)
    skip,    // (4) symbol is stripped, silently (tab, VT, FF, CR)
    count,    // (5) track the number of lines
  };




  /* How to handle an input character in a FASTA sequence: one row per byte
     class (see core/fastx_char_class.hpp). The fallback carries 190 of the 256
     bytes -- a printable that is not IUPAC, a space and a high byte are all
     stripped with a warning.

     C++11 constexpr allows a single return statement -- no if, no switch -- so
     the chain below can only be nested conditional operators. */
  // NOLINTBEGIN(readability-avoid-nested-conditional-operator)
  constexpr auto sequence_policy(vsearch::CharClass const character_class) noexcept -> Action {
    return character_class == vsearch::CharClass::iupac           ? Action::accept
      : character_class == vsearch::CharClass::line_feed          ? Action::count
      : character_class == vsearch::CharClass::carriage_return    ? Action::skip
      : character_class == vsearch::CharClass::blank_control      ? Action::skip
      : character_class == vsearch::CharClass::control            ? Action::show
      : character_class == vsearch::CharClass::dot_dash           ? Action::reject
      : Action::warn;  // space, del_or_high, printable_other
  }
  // NOLINTEND(readability-avoid-nested-conditional-operator)


  constexpr std::array<Action, byte_range> char_actions =
    vsearch::expanded_policy<Action, sequence_policy>();


  auto map_action(char const nucleotide) -> Action {
    auto const current_char = static_cast<unsigned char>(nucleotide);
    return char_actions[current_char];
  }


  auto report_illegal_symbol_and_exit(fastx_handle input_handle, unsigned char const symbol, uint64_t const line_number) -> void {
    std::string const message =
      "Illegal character '" + std::string(1, static_cast<char>(symbol))
      + "' in sequence on line " + decimal::to_text(line_number)
      + " of FASTA file";
    /* deferred-error mode (see fastx.h): record and return instead of
       exiting, so a worker thread does not std::exit() with siblings live */
    if (input_handle->defers_errors()) {
      input_handle->set_deferred_error(message);
      return;
    }
    fatal(message);
  }


  auto report_unprintable_symbol_and_exit(fastx_handle input_handle, unsigned char const symbol, uint64_t const line_number) -> void {
    std::string const message =
      "Illegal unprintable ASCII character no " + decimal::to_text(symbol)
      + " in sequence on line " + decimal::to_text(line_number)
      + " of FASTA file";
    if (input_handle->defers_errors()) {
      input_handle->set_deferred_error(message);
      return;
    }
    fatal(message);
  }

}  // end of anonymous namespace


auto fasta_open(const char * filename, struct Parameters const & parameters) -> std::unique_ptr<fastx_s>
{
  // fastx_open hands back an owning unique_ptr; on the fatal() path below it
  // frees the handle as the stack unwinds (library session), otherwise it is
  // moved out to the caller, which owns the reader by RAII.
  auto input_handle = fastx_open(filename, parameters);

  if (input_handle->is_fastq_input() and not input_handle->is_empty_input())
    {
      fatal(std::string("FASTA file expected, FASTQ file found (")
            + std::string(filename)
            + ")");
    }

  assert(input_handle != nullptr);
  return input_handle;
}


template <Mapping mapping>
auto fasta_filter_sequence(fastx_handle input_handle) -> void
{
  /* Strip unwanted characters from the sequence and raise warnings or
     errors on certain characters. */

  auto * dest = input_handle->sequence_buffer.data();

  /* bounded by the recorded length, not by a NUL terminator: a NUL inside a
     sequence line used to end the scan, silently dropping the rest of the
     record (and leaving it unchecked for illegal characters). The length bound
     is what lets byte 0 reach the table, which reports it like the other
     unprintable symbols. Writes trail reads by construction, so the filtering
     stays in place. */
  for (auto const symbol : input_handle->sequence_buffer.view())
    {
      auto const current_char = static_cast<unsigned char>(symbol);
      auto const action = map_action(symbol);

      /* Fast path: nucleotides dominate the input, so test the common
         'accept' action with a single, well-predicted branch. Handling it
         inside the switch let the compiler lower the switch to a jump table
         whose per-character indirect branch mispredicts on every byte. */
      if (action == Action::accept)
        {
          /* legal character */
          *dest = map_accepted_base<mapping>(symbol);
          dest = std::next(dest);
        }
      else
        {
          switch (action)
            {
            case Action::warn:
              /* stripped */
              input_handle->record_stripped(current_char);
              break;

            case Action::reject:
              /* fatal character */
              report_illegal_symbol_and_exit(input_handle, current_char, input_handle->lineno);
              if (input_handle->error) { return; }
              break;

            case Action::show:
              /* fatal unprintable character */
              report_unprintable_symbol_and_exit(input_handle, current_char, input_handle->lineno);
              if (input_handle->error) { return; }
              break;

            case Action::skip:
              /* silently stripped chars (whitespace) */
              break;

            case Action::count:
              /* newline (silently stripped) */
              ++input_handle->lineno;
              break;

            case Action::accept:
              /* handled above on the fast path */
              break;
            }
        }
    }

  /* add nullchar after sequence */
  *dest = '\0';
  input_handle->sequence_buffer.length = static_cast<uint64_t>(dest - input_handle->sequence_buffer.data());
}


auto fasta_next(fastx_handle input_handle,
                bool const truncateatspace,
                const unsigned char * char_mapping) -> bool
{
  input_handle->lineno_start = input_handle->lineno;

  input_handle->header_buffer.length = 0;
  input_handle->header_buffer.data()[0] = 0;
  input_handle->sequence_buffer.length = 0;
  input_handle->sequence_buffer.data()[0] = 0;

  std::size_t rest = fastx_file_fill_buffer(input_handle);

  if (rest == 0)
    {
      return false;
    }

  /* read header */

  /* check initial > character */

  /* guaranteed by two invariants, so an assertion rather than a check:
     fastx_open() rejects any input whose first byte is neither '>' nor '@'
     ("File type not recognized."), and the sequence loop below stops only at
     end of file or at the '>' starting the next record -- so every call lands
     on '>' or on end of file (returned above) */
  assert(input_handle->file_buffer.peek() == '>');
  ++input_handle->file_buffer.position;

  bool header_complete = false;
  while (not header_complete)
    {
      /* get more data if buffer empty*/
      rest = fastx_file_fill_buffer(input_handle);
      if (rest == 0)
        {
          if (input_handle->defers_errors())
            {
              input_handle->set_deferred_error("Invalid FASTA - header must be terminated with newline");
              return false;
            }
          fatal("Invalid FASTA - header must be terminated with newline");
        }

      /* copy to header buffer */
      auto const fragment = scan_line_fragment(input_handle);
      input_handle->header_buffer.extend(fragment.view);
      consume_fragment(input_handle, fragment);
      if (fragment.has_newline)
        {
          ++input_handle->lineno;
        }
      header_complete = fragment.has_newline;
    }

  /* read one or more sequence lines */

  // the header loop above consumed a full LF-terminated line
  bool previous_line_complete = true;
  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(input_handle);

      /* end if no more data */
      if (rest == 0)
        {
          break;
        }

      /* end if new sequence starts */
      if (previous_line_complete and (input_handle->file_buffer.peek() == '>'))
        {
          break;
        }

      auto const fragment = scan_line_fragment(input_handle);
      input_handle->sequence_buffer.extend(fragment.view);
      consume_fragment(input_handle, fragment);
      previous_line_complete = fragment.has_newline;
    }

  ++input_handle->seqno;

  fastx_filter_header(input_handle, truncateatspace);
  /* The mapping is a compile-time fact at every caller (see maps.hpp), and
     these are the only two tables in the tree, so the dispatch here is
     exhaustive -- the assert says so, because a third table would otherwise
     take the pass-through path in silence. One pointer comparison per record
     replaces a table load per accepted byte. */
  assert((char_mapping == chrmap_no_change()) or (char_mapping == chrmap_upcase()));
  if (char_mapping == chrmap_upcase())
    {
      fasta_filter_sequence<Mapping::upcase>(input_handle);
    }
  else
    {
      fasta_filter_sequence<Mapping::none>(input_handle);
    }
  fastx_filter_sequence_length(input_handle);

  return true;
}


/* fasta output */

namespace {
auto fasta_print_sequence(std::FILE * output_handle, View<char> const seq, std::size_t const len, int const width) -> void
{
  /*
    The actual length of the sequence may be longer than "len", but only
    "len" characters are printed.

    Specify width of lines - zero (or <1) means linearize (all on one line).
  */

  if (width < 1)  // no sequence folding
    {
      fprint(output_handle, seq.first(len));
      fprint(output_handle, '\n');
    }
  else  // sequence folding every 'width'
    {
      auto const width_u = static_cast<std::size_t>(width);
      for (std::size_t i = 0; i < len; i += width_u)
        {
          fprint(output_handle, seq.subspan(i, std::min(len - i, width_u)));
          fprint(output_handle, '\n');
        }
    }
}
}  // anonymous namespace


auto fasta_print(std::FILE * output_handle, View<char> const header,
                 View<char> const seq,
                 struct Parameters const & parameters) -> void
{
  fprint(output_handle, '>');
  fprint(output_handle, header);
  fprint(output_handle, '\n');
  fasta_print_sequence(output_handle, seq, seq.size(),
                       static_cast<int>(parameters.opt_fasta_width));
}


auto fasta_print_general(std::FILE * output_handle,
                         char const * prefix,
                         View<char> const seq,
                         View<char> const header,
                         OutputAnnotations const & annotations,
                         struct Parameters const & parameters) -> void
{
  fprint(output_handle, '>');

  if (prefix != nullptr)
    {
      std::fputs(prefix, output_handle);
    }

  fprint_header_annotations(output_handle, seq, header, annotations, parameters);

  fprint(output_handle, '\n');

  fasta_print_sequence(output_handle, seq, seq.size(),
                       static_cast<int>(parameters.opt_fasta_width));
}


auto fasta_print_general(std::FILE * output_handle,
                         char const * prefix,
                         SeqRecord const & record,
                         OutputAnnotations const & annotations,
                         struct Parameters const & parameters) -> void
{
  fasta_print_general(output_handle,
                      prefix,
                      record.sequence,
                      record.header,
                      annotations,
                      parameters);
}


// A single uint64_t ordinal parameter: it is the widest unsigned type, so
// every caller (passing int, size_t or uint64_t, all non-negative 1-based
// counters) converts without narrowing or sign-change. Two overloads taking
// int and size_t were ambiguous for a uint64_t argument on platforms where
// uint64_t, size_t and int are all distinct types (e.g. macOS).
auto fasta_print_db_relabel(std::FILE * output_handle,
                            uint64_t const seqno,
                            uint64_t const ordinal,
                            struct Database const & db,
                            struct Parameters const & parameters) -> void
{
  fasta_print_general(output_handle,
                      nullptr,
                      db.record(seqno),
                      OutputAnnotations{db.getabundance(seqno), static_cast<int64_t>(ordinal)},
                      parameters);
}


auto fasta_print_db(std::FILE * output_handle, uint64_t const seqno,
                    struct Database const & db,
                    struct Parameters const & parameters) -> void
{
  fasta_print_general(output_handle,
                      nullptr,
                      db.record(seqno),
                      OutputAnnotations{db.getabundance(seqno), 0},
                      parameters);
}
