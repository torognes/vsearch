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

#include "vsearch.hpp"
#include "commands/search_oligodb.hpp"
#include "core/db.hpp"
#include "core/fastx.hpp"
#include "core/mask.hpp"  // apply_masking
#include "core/match_counts.hpp"  // vsearch::MatchCounts, vsearch::print_match_counts
#include "core/oligo_scan.hpp"
#include "core/results.hpp"
#include "core/searchcore.hpp"  // struct hit, align_trim, difference_count
#include "utils/chunk_reorder.hpp"  // ChunkReorder
#include "utils/fatal.hpp"
#include "utils/maps/four_bit.hpp"  // vsearch::maps::four_bit::map
#include "utils/open_file.hpp"
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/reverse_complement.hpp"
#include "utils/span.hpp"
#include "utils/threads.hpp"
#include "utils/view.hpp"
#include "utils/worker_loop.hpp"
#include <algorithm>  // std::any_of, std::reverse
#include <cassert>  // assert
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE
#include <mutex>  // std::mutex, std::lock_guard
#include <string>  // std::string, std::to_string
#include <vector>


namespace {

  namespace oligo = vsearch::oligo;


  /* One query of a chunk: where its header and sequence sit in the chunk's
     buffers (offsets, not views, so that a chunk can be moved), and its hits,
     ordered by position. */
  struct oligo_query_s
  {
    int64_t qsize = 0;
    uint64_t position = 0;  /* input read so far, for the progress bar */
    std::size_t head_offset = 0;
    std::size_t head_length = 0;
    std::size_t seq_offset = 0;
    std::size_t seq_length = 0;
    std::vector<struct hit> hits;
  };


  /* The queries a worker claims at once. rc_sequences holds the reverse
     complements the minus-strand hits are reported against, at the same
     offsets as the sequences. */
  struct oligo_chunk_s
  {
    unsigned long rank = 0;  /* claim order, which is also the input order */
    std::vector<char> headers;
    std::vector<char> sequences;
    std::vector<char> rc_sequences;
    std::vector<struct oligo_query_s> queries;
    std::size_t count = 0;  /* queries[0 .. count) are live */
  };


  /* Queries claimed at once, as --search_exact does: parsing happens under
     the input lock, and large chunks keep the lock hand-offs rare. The
     second bound is on bases, because long reads range from a few
     nucleotides to hundreds of kilobases: 8192 reads of 100 kb would make
     one claim of 800 Mb, and leave the other threads idle. A chunk of 2 Mb
     costs 2 Mb x lanes x ~1 ns of scanning, a fraction of a second for a
     hundred oligos on both strands. */
  constexpr std::size_t queries_per_chunk = 8192;
  constexpr std::size_t bases_per_chunk = std::size_t{2} * 1024 * 1024;


  struct search_oligodb_state_s
  {
    struct Parameters const & parameters;
    struct Database db;  /* the oligos */
    oligo::LaneSet const * lanes = nullptr;
    fastx_handle query_fastx_h = nullptr;

    /* accessed by the worker threads; access serialized by the mutexes */
    std::mutex mutex_input;
    std::mutex mutex_output;
    unsigned long next_claim_rank = 0;  /* under mutex_input */
    ChunkReorder<struct oligo_chunk_s> waiting;  /* under mutex_output */
    int queries = 0;
    uint64_t queries_abundance = 0;
    int qmatches = 0;
    uint64_t qmatches_abundance = 0;
    uint64_t hit_count = 0;
    std::FILE * fp_alnout = nullptr;
    std::FILE * fp_userout = nullptr;
    std::FILE * fp_blast6out = nullptr;
    Progress * progress = nullptr;

    explicit search_oligodb_state_s(struct Parameters const & params) : parameters(params) {}
  };


  /* The oligos are read with no length filter, so that an empty oligo is an
     error rather than a warning, and each is then checked against what the
     scan supports. */
  auto read_oligos(struct search_oligodb_state_s & state) -> void
  {
    struct Parameters unfiltered = state.parameters;
    unfiltered.opt_minseqlength = 0;
    state.db.read(state.parameters.opt_db, 0, unfiltered);

    auto const count = state.db.getsequencecount();
    if (count == 0)
      {
        fatal("No oligo in the database file given with --db");
      }
    for (uint64_t index = 0; index < count; ++index)
      {
        auto const sequence = state.db.sequence_view(index);
        auto const label = std::string(state.db.header_view(index).cbegin(),
                                       state.db.header_view(index).cend());
        if (sequence.empty())
          {
            fatal("Oligo " + label + " is empty");
          }
        if (sequence.size() > oligo::max_oligo_length)
          {
            fatal("Oligo " + label + " is longer than " +
                  std::to_string(oligo::max_oligo_length) + " nucleotides");
          }
        auto const is_not_nucleotide = [](char const symbol) -> bool {
          return vsearch::maps::four_bit::map(symbol) == 0;
        };
        if (std::any_of(sequence.cbegin(), sequence.cend(), is_not_nucleotide))
          {
            fatal("Oligo " + label + " contains a symbol that is not a nucleotide");
          }
      }
  }


  auto append_run(std::string & cigar, int64_t const run, char const operation) -> void
  {
    if (run == 0) { return; }
    cigar += std::to_string(run);
    cigar += operation;
  }


  /* An occurrence as a vsearch hit: a global alignment of the oligo (the
     target) against the whole query, where the query before and after the
     occurrence, and the oligo bases hanging off a query end, are terminal
     gaps. align_trim() then derives the internal alignment the writers
     report. A minus-strand occurrence was found by scanning the reverse
     complement of the oligo along the query; the hit is the same alignment
     seen from the other strand, the oligo against the reverse-complemented
     query, so its columns are reversed. */
  auto make_hit(oligo::Occurrence const & occurrence,
                oligo::Lane const & lane,
                int64_t const query_length,
                struct Parameters const & parameters) -> struct hit
  {
    auto const pattern_length = static_cast<int64_t>(lane.pattern.size());
    auto const flank_left = occurrence.read_start;
    auto const flank_right = query_length - 1 - occurrence.read_end;
    auto const overhang_left = static_cast<int64_t>(occurrence.pattern_start);
    auto const overhang_right = pattern_length - occurrence.pattern_end;
    auto const minus = (lane.strand == oligo::Strand::minus);

    /* the columns in the lane's orientation, then in the hit's */
    std::string columns;
    columns.reserve(occurrence.columns.size() + 4);
    columns.append(static_cast<std::size_t>(flank_left), 'D');
    columns.append(static_cast<std::size_t>(overhang_left), 'I');
    columns += occurrence.columns;
    columns.append(static_cast<std::size_t>(overhang_right), 'I');
    columns.append(static_cast<std::size_t>(flank_right), 'D');
    if (minus)
      {
        std::reverse(columns.begin(), columns.end());
      }

    struct hit hit {};
    hit.target = static_cast<int>(lane.oligo);
    hit.strand = minus ? 1 : 0;
    hit.count = 0;
    hit.accepted = true;
    hit.rejected = false;
    hit.aligned = true;
    hit.weak = false;

    /* run-length encode, and count the gap runs, terminal ones included */
    auto gap_runs = 0;
    auto run = int64_t{0};
    auto previous = '\0';
    for (auto const operation : columns)
      {
        if ((operation != previous) and (run > 0))
          {
            append_run(hit.nwalignment, run, previous);
            if (previous != 'M') { ++gap_runs; }
            run = 0;
          }
        previous = operation;
        ++run;
      }
    append_run(hit.nwalignment, run, previous);
    if ((run > 0) and (previous != 'M')) { ++gap_runs; }

    auto const terminal_columns = flank_left + flank_right + overhang_left + overhang_right;
    hit.nwalignmentlength = static_cast<int>(columns.size());
    hit.matches = occurrence.matches;
    hit.mismatches = occurrence.mismatches;
    hit.nwindels = static_cast<int>(terminal_columns) + occurrence.gap_columns;
    hit.nwgaps = gap_runs;
    hit.nwdiff = hit.mismatches + hit.nwindels;
    hit.nwid = (hit.nwalignmentlength > 0) ?
      100.0 * hit.matches / hit.nwalignmentlength : 0.0;
    hit.shortest = static_cast<int>(std::min(query_length, pattern_length));
    hit.longest = static_cast<int>(std::max(query_length, pattern_length));
    /* raw: the edit distance of the reported alignment */
    hit.nwscore = occurrence.mismatches + occurrence.gap_columns;

    align_trim(hit, parameters);
    assert(difference_count(hit) == occurrence.mismatches + occurrence.gap_columns);
    assert(hit.internal_gaps == occurrence.gap_openings);
    return hit;
  }


  /* Write the hits of one query. Called with mutex_output held. */
  auto output_query(struct search_oligodb_state_s & state,
                    struct oligo_query_s const & query,
                    View<char> const header,
                    View<char> const sequence,
                    View<char> const rc_sequence) -> void
  {
    struct Parameters const & parameters = state.parameters;
    if (query.hits.empty()) { return; }

    if (state.fp_alnout != nullptr)
      {
        results_show_alnout(state.fp_alnout, make_view(query.hits),
                            header, sequence, state.db, parameters);
      }

    for (auto const & hit : query.hits)
      {
        if (state.fp_userout != nullptr)
          {
            results_show_userout_one(state.fp_userout, &hit, header, sequence,
                                     rc_sequence, state.db, parameters,
                                     HitSpan::aligned_region);
          }
        if (state.fp_blast6out != nullptr)
          {
            results_show_blast6out_one(state.fp_blast6out, &hit, header,
                                       static_cast<int64_t>(sequence.size()),
                                       state.db, HitSpan::aligned_region);
          }
      }
  }


  /* Write the queries of a chunk, in order, with the statistics and the
     progress they account for. Called with mutex_output held. */
  auto output_chunk(struct search_oligodb_state_s & state,
                    struct oligo_chunk_s const & chunk) -> void
  {
    for (std::size_t nth = 0; nth < chunk.count; ++nth)
      {
        auto const & query = chunk.queries[nth];
        auto const header = make_view(chunk.headers).subspan(query.head_offset, query.head_length);
        auto const sequence = make_view(chunk.sequences).subspan(query.seq_offset, query.seq_length);
        auto const rc_sequence = make_view(chunk.rc_sequences).subspan(query.seq_offset, query.seq_length);
        output_query(state, query, header, sequence, rc_sequence);

        ++state.queries;
        state.queries_abundance += static_cast<uint64_t>(query.qsize);
        state.hit_count += query.hits.size();
        if (not query.hits.empty())
          {
            ++state.qmatches;
            state.qmatches_abundance += static_cast<uint64_t>(query.qsize);
          }
      }
    /* one progress update per chunk, not per query */
    if (chunk.count > 0)
      {
        state.progress->update(chunk.queries[chunk.count - 1].position);
      }
  }


  auto search_oligodb_thread_run(struct search_oligodb_state_s & state) -> void
  {
    struct Parameters const & parameters = state.parameters;
    auto const & lanes = *state.lanes;
    oligo::Scanner scanner(lanes);
    std::vector<oligo::Occurrence> occurrences;
    struct oligo_chunk_s chunk;

    /* parse the next queries into the chunk, and give it its rank */
    auto const has_work_to_claim = [&]() -> bool {
      chunk.headers.clear();
      chunk.sequences.clear();
      chunk.count = 0;
      while ((chunk.count < queries_per_chunk) and
             (chunk.sequences.size() < bases_per_chunk))
        {
          if (not state.query_fastx_h->next(header_truncation(parameters.opt_notrunclabels), Mapping::none))
            {
              break;
            }
          if (chunk.queries.size() <= chunk.count)
            {
              chunk.queries.resize(chunk.count + 1);
            }
          auto & query = chunk.queries[chunk.count];
          auto const qhead = state.query_fastx_h->header_view();
          auto const qseq = state.query_fastx_h->sequence_view();
          query.qsize = state.query_fastx_h->get_abundance();
          query.head_offset = chunk.headers.size();
          query.head_length = qhead.size();
          chunk.headers.insert(chunk.headers.end(), qhead.cbegin(), qhead.cend());
          query.seq_offset = chunk.sequences.size();
          query.seq_length = qseq.size();
          chunk.sequences.insert(chunk.sequences.end(), qseq.cbegin(), qseq.cend());
          /* get progress as amount of input file read */
          query.position = state.query_fastx_h->get_position();
          ++chunk.count;
        }
      if (chunk.count == 0)
        {
          return false;
        }
      chunk.rank = state.next_claim_rank;
      ++state.next_claim_rank;
      return true;
    };

    auto const process_chunk = [&]() -> void {
      chunk.rc_sequences.resize(chunk.sequences.size());
      for (std::size_t nth = 0; nth < chunk.count; ++nth)
        {
          auto & query = chunk.queries[nth];
          auto const sequence = make_span(chunk.sequences).subspan(query.seq_offset, query.seq_length);
          /* the query is masked in place, in the chunk, where the output
             writers read it (by default --qmask none: nothing to do) */
          apply_masking(sequence, parameters.opt_qmask, parameters);
          scanner.scan(View<char>{sequence}, occurrences);
          query.hits.clear();
          for (auto const & occurrence : occurrences)
            {
              query.hits.push_back(make_hit(occurrence, lanes.lanes()[occurrence.lane],
                                            static_cast<int64_t>(query.seq_length),
                                            parameters));
            }
          auto const is_minus = [](struct hit const & hit) -> bool { return hit.strand != 0; };
          if (std::any_of(query.hits.cbegin(), query.hits.cend(), is_minus))
            {
              reverse_complement(make_span(chunk.rc_sequences).subspan(query.seq_offset, query.seq_length),
                                 View<char>{sequence});
            }
        }

      std::lock_guard<std::mutex> const output_lock(state.mutex_output);
      /* chunks are written in claim order, which is the input order: a chunk
         searched ahead of its turn waits, and the worker moves on */
      state.waiting.submit(chunk.rank, chunk,
                           [&state](struct oligo_chunk_s const & ready) -> void {
                             output_chunk(state, ready);
                           });
    };

    run_worker_loop(state.mutex_input, has_work_to_claim, process_chunk);
  }


  auto model_from(struct Parameters const & parameters) -> oligo::Model
  {
    return (parameters.opt_maxgaps == 0) ? oligo::Model::substitutions_only
      : oligo::Model::edit_distance;
  }

}  // anonymous namespace


auto search_oligodb(struct Parameters const & parameters) -> void
{
  search_oligodb_state_s state(parameters);

  /* open output files; the handles are owned here so they outlive the worker
     pool, which reads the non-owning state.fp_* under the output lock */
  OutputFileHandle alnout_handle = open_optional_output_file(parameters.opt_alnout, OutputOption{"--alnout"});
  state.fp_alnout = alnout_handle.get();
  if (state.fp_alnout != nullptr)
    {
      fprint(state.fp_alnout, make_view(parameters.runtime.command_line));
      fprint(state.fp_alnout, '\n');
      fprint(state.fp_alnout, make_view(parameters.runtime.prog_header));
      fprint(state.fp_alnout, '\n');
    }
  OutputFileHandle userout_handle = open_optional_output_file(parameters.opt_userout, OutputOption{"--userout"});
  state.fp_userout = userout_handle.get();
  OutputFileHandle blast6out_handle = open_optional_output_file(parameters.opt_blast6out, OutputOption{"--blast6out"});
  state.fp_blast6out = blast6out_handle.get();

  read_oligos(state);

  std::vector<View<char>> oligos;
  for (uint64_t index = 0; index < state.db.getsequencecount(); ++index)
    {
      oligos.push_back(state.db.sequence_view(index));
    }
  oligo::Limits limits;
  limits.max_diffs = parameters.opt_maxdiffs;
  limits.max_gap_openings = parameters.opt_maxgaps;
  limits.target_cov = parameters.opt_target_cov;
  limits.model = model_from(parameters);
  oligo::LaneSet const lanes(oligos,
                             parameters.opt_strand ? oligo::Strands::both : oligo::Strands::plus_only,
                             parameters.opt_n_mismatch ? oligo::NPositions::mismatch : oligo::NPositions::match,
                             limits);
  state.lanes = &lanes;

  auto const query_owner = fastx_open(parameters.input_filename, parameters);
  state.query_fastx_h = query_owner.get();  // workers borrow the raw handle

  /* The query file is parsed inside the worker threads. Defer parse errors
     so a malformed query stops the pool cooperatively instead of calling
     fatal() from a worker while siblings are writing output; reported below
     from the main thread after join. */
  state.query_fastx_h->enable_deferred_errors();

  {
    Progress progress_bar("Searching", state.query_fastx_h->get_size(), parameters);
    state.progress = &progress_bar;
    ThreadRunner threadrunner(static_cast<std::size_t>(parameters.opt_threads),
                              [&state](uint64_t const /* t */) -> void
                              { search_oligodb_thread_run(state); });
    threadrunner.run();
  }
  /* every claimed chunk was searched, so every rank up to the last was
     written */
  assert(state.waiting.empty());

  if (state.query_fastx_h->get_error())
    {
      fatal(state.query_fastx_h->get_errmsg());
    }

  query_owner->report_stripped_warning(parameters);

  auto const match_counts = vsearch::MatchCounts{state.qmatches, state.queries,
                                                 state.qmatches_abundance,
                                                 state.queries_abundance};
  if (not parameters.opt_quiet)
    {
      vsearch::print_match_counts(stderr, match_counts, parameters.opt_sizein);
      fprint(stderr, "Occurrences found: ");
      fprint_integer(stderr, state.hit_count);
      fprint(stderr, '\n');
    }
  if (parameters.fp_log != nullptr)
    {
      vsearch::print_match_counts(parameters.fp_log, match_counts, parameters.opt_sizein);
      fprint(parameters.fp_log, "Occurrences found: ");
      fprint_integer(parameters.fp_log, state.hit_count);
      fprint(parameters.fp_log, '\n');
    }

  state.db.clear();
}
