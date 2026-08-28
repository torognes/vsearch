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
#include "vsearch.hpp"
#include "core/attributes.hpp"
#include "core/seq_record.hpp"  // struct SeqRecord (for the record overload)
#include "core/fastx.hpp"
#include "core/fastx_char_class.hpp"  // vsearch::CharClass, class_of
#include "core/illegal_character.hpp"  // vsearch::illegal_character_message
#include "utils/fatal.hpp"
#include "utils/maps.hpp"  // Mapping, map_accepted_base, chrmap_*
#include "utils/print_view.hpp"  // fprint
#include "utils/quality_encoding.hpp"  // classify_encoding, is_phred64, offset_of
#include "utils/warn.hpp"  // vsearch::warn
#include <array>
#include <cassert>  // assert
#include <cstddef>  // std::ptrdiff_t
#include <cstdint> // uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::snprintf, std::size_t
#include <iterator>  // std::next
#include <memory>  // std::unique_ptr
#include <string>  // std::string, std::to_string


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // How to handle an input character, mirroring fasta.cc's Action. The tables
  // below map every byte to one of these; buffer_filter_extend() switches on
  // it. Unlike fasta, 'reject' is deferred (recorded and reported by the
  // caller), 'newline' is a no-op here (the caller keeps the line count), and
  // there is no 'warn': FASTQ strips nothing with a warning ("the rest is
  // fatal", see the policies below), so no byte maps to it and record_stripped()
  // is
  // reached from the FASTA parser only.
  enum struct Action : unsigned char {
    accept,   // (0) legal character
    reject,   // (1) fatal character (recorded, reported by the caller)
    skip,     // (2) silently stripped (e.g. CR)
    newline,   // (3) LF; silently stripped here
  };






  /* How to handle an input character in a FASTQ sequence: all IUPAC characters
     are valid, CR (^M) is silently stripped, LF is newline, the rest is fatal.
     One row per byte class (see core/fastx_char_class.hpp).

     C++11 constexpr allows a single return statement -- no if, no switch -- so
     the two chains below can only be nested conditional operators. */
  constexpr auto sequence_policy(vsearch::CharClass const character_class) noexcept -> Action {
    return character_class == vsearch::CharClass::iupac           ? Action::accept
      : character_class == vsearch::CharClass::line_feed          ? Action::newline
      : character_class == vsearch::CharClass::carriage_return    ? Action::skip
      : Action::reject;
  }


  /* How to handle a FASTQ quality character: any from 33 to 126 is valid
     (legal), CR (^M) is silently stripped, LF is newline, the rest is fatal. */
  constexpr auto quality_policy(vsearch::CharClass const character_class) noexcept -> Action {
    return character_class == vsearch::CharClass::line_feed       ? Action::newline
      : character_class == vsearch::CharClass::carriage_return    ? Action::skip
      : character_class == vsearch::CharClass::iupac              ? Action::accept
      : character_class == vsearch::CharClass::printable_other    ? Action::accept
      : character_class == vsearch::CharClass::dot_dash           ? Action::accept
      : Action::reject;  // space, del_or_high, control, blank_control
  }


  constexpr std::array<Action, byte_range> char_fq_action_seq =
    vsearch::expanded_policy<Action, sequence_policy>();

  constexpr std::array<Action, byte_range> char_fq_action_qual =
    vsearch::expanded_policy<Action, quality_policy>();

}  // end of anonymous namespace


namespace {
/* Report an already-complete message. The two illegal-character sites carry
   their own location clause (see core/illegal_character.hpp) and so cannot go
   through the "Invalid line N in FASTQ file: " prefix below; they share this
   dispatch with it instead of duplicating it.

   deferred-error mode (see fastx.h): record the message and return
   so the worker can stop cooperatively instead of std::exit()-ing
   from a worker thread. The caller (fastq_next) returns false right
   after this call. */
auto fastq_report(fastx_handle input_handle, std::string const & message) -> void
{
  if (input_handle->defers_errors())
    {
      input_handle->set_deferred_error(message);
      return;
    }
  fatal(message);
}


/* msg is a std::string rather than a char const *: the six remaining callers
   pass a literal and build one temporary each, on a path that has already
   given up on the file. (The two that used to assemble msg by concatenation,
   and hand over msg.c_str() only for the line below to concatenate it back,
   now build a complete message and call fastq_report() directly.) */
auto fastq_fatal(fastx_handle input_handle, uint64_t const lineno, std::string const & msg) -> void
{
  /* decimal::to_text, not std::to_string: on libstdc++ <= 10 the latter is a
     std::vsnprintf call with a format string (see decimal_digits.hpp). */
  std::string const message = "Invalid line " + decimal::to_text(lineno)
    + " in FASTQ file: " + msg;
  fastq_report(input_handle, message);
}


// no fastx_handle parameter: with no 'warn' action there is nothing to record
// on the handle, so this is a pure transform over the two buffers (the symbol
// range below is passed in by reference for the same reason)
//
// track_range is a compile-time switch, not a runtime flag: the quality line
// wants the min/max of the bytes it copies and the sequence line does not, and
// the two share this instantiation. At false the compiler drops the two
// comparisons entirely, so the sequence path pays nothing.
template <Mapping mapping, bool track_range>
auto buffer_filter_extend(FastxBuffer & dest_buffer,
                          View<char> const source,
                          Action const * char_action,
                          bool & ok,
                          char & illegal_char,
                          QualitySymbolRange & range) -> void
{
  dest_buffer.makespace(source.size() + 1);

  /* Strip unwanted characters from the string and record the first illegal
     one for the caller to report. */

  auto * d = std::next(dest_buffer.data(),
                       static_cast<std::ptrdiff_t>(dest_buffer.length));
  auto * q = d;
  ok = true;

  for (auto const symbol : source)
    {
      auto const action = char_action[static_cast<unsigned char>(symbol)];

      /* Fast path: legal characters dominate, so test 'accept' with a single
         predictable branch and keep it off the switch (which the compiler may
         lower to a mispredicting jump table). Same precaution as in fasta.cc. */
      if (action == Action::accept)
        {
          /* legal character */
          if (track_range)
            {
              range.observe(static_cast<unsigned char>(symbol));
            }
          *q++ = map_accepted_base<mapping>(symbol);
        }
      else
        {
          switch (action)
            {
            case Action::reject:
              /* fatal character */
              if (ok)
                {
                  illegal_char = symbol;
                }
              ok = false;
              break;

            case Action::skip:
              /* silently stripped chars (whitespace) */
              break;

            case Action::newline:
              /* newline (silently stripped) */
              break;

            case Action::accept:
              /* handled above on the fast path */
              break;
            }
        }
    }

  /* add zero after sequence */
  *q = 0;
  dest_buffer.length += static_cast<uint64_t>(q - d);
}
}  // anonymous namespace


auto fastq_open(const char * filename, struct Parameters const & parameters) -> std::unique_ptr<fastx_s>
{
  // fastx_open hands back an owning unique_ptr; on the fatal() path below it
  // frees the handle as the stack unwinds (library session), otherwise it is
  // moved out to the caller, which owns the reader by RAII.
  auto input_handle = fastx_open(filename, parameters);

  if (not input_handle->is_fastq_input())
    {
      fatal(std::string("FASTQ file expected, FASTA file found (")
            + std::string(filename)
            + ")");
    }

  assert(input_handle != nullptr);
  return input_handle;
}


/* Warn once per file when the observed quality symbols contradict the
   requested --fastq_ascii.

   Until 3.0 this case was caught by accident: a phred+64 file read at the
   default offset 33 decodes to Q71, --fastq_qmax defaulted to 41, and vsearch
   stopped. Now that the bound accepts everything the encoding can represent,
   Q71 is a legal Sanger score and the run goes through with every error
   probability understated by 31 Phred units -- a Q20 base is scored Q51, and
   an --fastq_maxee filter that should discard the read keeps it.

   The heuristic is classify_encoding(), the same one --fastq_chars prints its
   guess from, so there is one set of thresholds and not two. It is a warning
   and not a fatal because the two encodings genuinely overlap: a PacBio HiFi
   file whose scores all exceed Q30 is indistinguishable from a phred+64 file
   by its symbols alone. That ambiguity is exactly when the user should look,
   which is what the message says. */
auto fastx_s::warn_if_offset_looks_wrong() -> void
{
  if (not warn_on_suspicious_offset) { return; }
  if (not is_fastq) { return; }
  if (not quality_range.seen()) { return; }
  if (seqno < minimum_records_for_offset_guess) { return; }

  auto const lowest = static_cast<char>(quality_range.lowest);
  auto const highest = static_cast<char>(quality_range.highest);
  auto const encoding = classify_encoding(lowest, highest);
  auto const implied_offset = offset_of(encoding);
  if (implied_offset == quality_offset) { return; }

  // warn once: this runs on the end-of-file branch of fastq_next()
  warn_on_suspicious_offset = false;

  /* decimal::to_text, not std::to_string: see fastq_fatal() above */
  auto const shift = implied_offset - quality_offset;
  vsearch::warn(
    "quality symbols in this file range from '" + std::string(1, lowest)
    + "' to '" + std::string(1, highest) + "', which looks like "
    + (is_phred64(encoding) ? "phred+64" : "phred+33")
    + ", but --fastq_ascii " + decimal::to_text(quality_offset)
    + " was used.\nIf that offset is wrong, every quality score is decoded "
    + (shift > 0 ? decimal::to_text(shift) + " too high"
                 : decimal::to_text(-shift) + " too low")
    + " and the error probabilities are wrong.\n"
      "Run 'vsearch --fastq_chars' on this file to check the encoding.");
}


auto fastq_next(fastx_handle input_handle,
                bool const truncateatspace,
                const unsigned char * char_mapping) -> bool
{
  input_handle->header_buffer.length = 0;
  input_handle->header_buffer.data()[0] = 0;
  input_handle->sequence_buffer.length = 0;
  input_handle->sequence_buffer.data()[0] = 0;
  input_handle->plusline_buffer.length = 0;
  input_handle->plusline_buffer.data()[0] = 0;
  input_handle->quality_buffer.length = 0;
  input_handle->quality_buffer.data()[0] = 0;

  input_handle->lineno_start = input_handle->lineno;

  auto ok = true;
  char illegal_char = '\0';

  auto rest = fastx_file_fill_buffer(input_handle);

  /* check end of file */

  if (rest == 0)
    {
      /* the whole file has been read, so the symbol range is final */
      input_handle->warn_if_offset_looks_wrong();
      return false;
    }

  /* read header */

  /* check initial @ character */

  /* guaranteed, so an assertion rather than a check: fastx_open() sets
     is_fastq (which is what routes here) only when the first byte is '@', and
     the quality loop below returns true only after stopping on a '@' whose
     record had equal sequence and quality lengths -- every other ending is
     reported as "Sequence and quality lines must be equally long" first */
  assert(input_handle->file_buffer.peek() == '@');
  input_handle->file_buffer.position++;

  bool header_complete = false;
  while (not header_complete)
    {
      /* get more data if buffer empty */
      rest = fastx_file_fill_buffer(input_handle);
      if (rest == 0)
        {
          fastq_fatal(input_handle, input_handle->lineno, "Unexpected end of file");
          return false;
        }

      /* copy to header buffer */
      auto const fragment = scan_line_fragment(input_handle);
      input_handle->header_buffer.extend(fragment.view);
      consume_fragment(input_handle, fragment);
      if (fragment.has_newline)
        {
          input_handle->lineno++;
        }
      header_complete = fragment.has_newline;
    }

  /* read sequence line(s) */
  bool previous_line_complete = false;
  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(input_handle);

      /* cannot end here */
      if (rest == 0)
        {
          fastq_fatal(input_handle, input_handle->lineno, "Unexpected end of file");
          return false;
        }

      /* end when new line starting with + is seen */
      if (previous_line_complete and (input_handle->file_buffer.peek() == '+'))
        {
          break;
        }

      /* copy to sequence buffer */
      auto const fragment = scan_line_fragment(input_handle);
      /* The mapping is a compile-time fact at every caller (see maps.hpp), and
         these are the only two tables in the tree, so the dispatch is
         exhaustive -- the assert says so, because a third table would
         otherwise take the pass-through path in silence. One pointer
         comparison per line replaces a table load per accepted byte. */
      assert((char_mapping == chrmap_no_change()) or (char_mapping == chrmap_upcase()));
      if (char_mapping == chrmap_upcase())
        {
          buffer_filter_extend<Mapping::upcase, false>(input_handle->sequence_buffer,
                                                       fragment.view,
                                                       char_fq_action_seq.data(),
                                                       ok, illegal_char,
                                                       input_handle->quality_range);
        }
      else
        {
          buffer_filter_extend<Mapping::none, false>(input_handle->sequence_buffer,
                                                     fragment.view,
                                                     char_fq_action_seq.data(),
                                                     ok, illegal_char,
                                                     input_handle->quality_range);
        }
      consume_fragment(input_handle, fragment);
      if (fragment.has_newline)
        {
          input_handle->lineno++;
        }
      previous_line_complete = fragment.has_newline;

      if (not ok)
        {
          fastq_report(input_handle,
                       vsearch::illegal_character_message(
                         vsearch::IllegalField::sequence,
                         vsearch::IllegalFormat::fastq,
                         static_cast<unsigned char>(illegal_char),
                         input_handle->lineno - (fragment.has_newline ? 1 : 0)));
          return false;
        }
    }

  /* read + line */

  /* skip + character */
  input_handle->file_buffer.position++;

  bool plusline_complete = false;
  while (not plusline_complete)
    {
      /* get more data if buffer empty */
      rest = fastx_file_fill_buffer(input_handle);

      /* cannot end here */
      if (rest == 0)
        {
          fastq_fatal(input_handle, input_handle->lineno, "Unexpected end of file");
          return false;
        }

      /* copy to plusline buffer */
      auto const fragment = scan_line_fragment(input_handle);
      input_handle->plusline_buffer.extend(fragment.view);
      consume_fragment(input_handle, fragment);
      if (fragment.has_newline)
        {
          input_handle->lineno++;
        }
      plusline_complete = fragment.has_newline;
    }

  /* check that the plus line is empty or identical to @ line */

  auto const header = input_handle->header_buffer.view();
  auto const plusline = input_handle->plusline_buffer.view();

  auto plusline_invalid = false;
  if (header.size() == plusline.size())
    {
      if (header != plusline)
        {
          plusline_invalid = true;
        }
    }
  else
    {
      /* a plus line of a different length than the header is legal only when it
         carries nothing but the line terminator: the LF on its own, or CR LF.
         The loop above exits only on a newline, so the LF is always there. */
      static constexpr auto crlf_length = std::size_t{2};
      if ((plusline.size() > crlf_length) or
          ((plusline.size() == crlf_length) and (plusline.front() != '\r')))
        {
          plusline_invalid = true;
        }
    }
  if (plusline_invalid)
    {
      /* the '+' line above is always LF-terminated here (its loop only
         exits on a newline), so the offending line is lineno - 1 */
      fastq_fatal(input_handle, input_handle->lineno - 1,
                  "'+' line must be empty or identical to header");
      return false;
    }

  /* read quality line(s) */

  bool last_line_complete = false;
  while (true)
    {
      /* get more data, if necessary */
      rest = fastx_file_fill_buffer(input_handle);

      /* end if no more data */
      if (rest == 0)
        {
          break;
        }

      /* end if next entry starts : LF + '@' + correct length */
      if (last_line_complete and
          (input_handle->file_buffer.peek() == '@') and
          (input_handle->quality_buffer.length == input_handle->sequence_buffer.length))
        {
          break;
        }

      /* copy to quality buffer */
      auto const fragment = scan_line_fragment(input_handle);
      /* no mapping: the table this used to index was the identity on all 256
         values, so the quality symbols are written through verbatim */
      /* one branch per line fragment, not per byte: after the sample the
         tracking instantiation is not called at all (see fastx.hpp) */
      if (input_handle->samples_quality_range())
        {
          buffer_filter_extend<Mapping::none, true>(input_handle->quality_buffer,
                                                    fragment.view,
                                                    char_fq_action_qual.data(),
                                                    ok, illegal_char,
                                                    input_handle->quality_range);
        }
      else
        {
          buffer_filter_extend<Mapping::none, false>(input_handle->quality_buffer,
                                                     fragment.view,
                                                     char_fq_action_qual.data(),
                                                     ok, illegal_char,
                                                     input_handle->quality_range);
        }
      consume_fragment(input_handle, fragment);
      if (fragment.has_newline)
        {
          input_handle->lineno++;
        }
      last_line_complete = fragment.has_newline;

      /* break if quality line already too long */
      if (input_handle->quality_buffer.length > input_handle->sequence_buffer.length)
        {
          break;
        }

      if (not ok)
        {
          fastq_report(input_handle,
                       vsearch::illegal_character_message(
                         vsearch::IllegalField::quality,
                         vsearch::IllegalFormat::fastq,
                         static_cast<unsigned char>(illegal_char),
                         input_handle->lineno - (fragment.has_newline ? 1 : 0)));
          return false;
        }
    }

  if (input_handle->sequence_buffer.length != input_handle->quality_buffer.length)
    {
      fastq_fatal(input_handle, input_handle->lineno - (last_line_complete ? 1 : 0),
                  "Sequence and quality lines must be equally long");
      return false;
    }

  fastx_filter_header(input_handle, truncateatspace);
  fastx_filter_sequence_length(input_handle);

  input_handle->seqno++;

  return true;
}


auto fastq_print_general(FILE * output_handle,
                         View<char> const seq,
                         View<char> const header,
                         View<char> const quality,
                         OutputAnnotations const & annotations,
                         struct Parameters const & parameters) -> void
{
  fprint(output_handle, '@');

  /* the same chain as the FASTA writer, and the five fields FASTQ output has
     no caller for suppress themselves at their sentinel defaults */
  fprint_header_annotations(output_handle, seq, header, annotations, parameters);

  fprint(output_handle, '\n');
  fprint(output_handle, seq);
  fprint(output_handle, "\n+\n");
  fprint(output_handle, quality);
  fprint(output_handle, '\n');
}


auto fastq_print_general(std::FILE * output_handle,
                         SeqRecord const & record,
                         OutputAnnotations const & annotations,
                         struct Parameters const & parameters) -> void
{
  fastq_print_general(output_handle,
                      record.sequence,
                      record.header,
                      record.quality,
                      annotations,
                      parameters);
}

