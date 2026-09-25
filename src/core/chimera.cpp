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
#include "utils/progress.hpp"
#include "core/align_simd.hpp"
#include "core/attributes.hpp"
#include "core/chimera.hpp"
#include "core/query_record.hpp"  // struct query_record_s
#include "core/chimera_internal.hpp"
#include "core/db.hpp"
#include "core/dbindex.hpp"
#include "core/fasta.hpp"
#include "core/fastx.hpp"
#include "core/linmemalign.hpp"
#include "core/mask.hpp"
#include "core/minheap.hpp"
#include "core/searchcore.hpp"
#include "core/udb.hpp"
#include "core/unique.hpp"
#include "utils/ascii_case.hpp"  // to_lower
#include "utils/base_mapping.hpp"
#include "utils/cigar.hpp"
#include "utils/fatal.hpp"
#include "utils/grow_to_fit.hpp"  // vsearch::grow_to_fit
#include "utils/make_unique.hpp"  // make_unique
#include "utils/maps/four_bit.hpp"
#include "utils/open_file.hpp"
#include "utils/span.hpp"
#include "utils/threads.hpp"
#include "utils/worker_loop.hpp"
#include "utils/print_view.hpp"  // fprint
#include <algorithm>  // std::copy, std::fill, std::fill_n, std::max, std::max_element, std::min, std::sort, std::transform
#include <array>
#include <cassert>
#include <cstddef> // std::ptrdiff_t, std::size_t
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fputs
#include <iterator>  // std::next
#include <limits>
#include <map>  // std::map
#include <memory>
#include <mutex>  // std::mutex, std::lock_guard
#include <numeric>  // std::accumulate
#include <string>  // std::string
#include <utility>  // std::move
#include <vector>
#include "utils/maps/upcase.hpp"

namespace upcase = vsearch::maps::upcase;

namespace four_bit = vsearch::maps::four_bit;


/*
  This code implements the method described in this paper:

  Robert C. Edgar, Brian J. Haas, Jose C. Clemente, Christopher Quince
  and Rob Knight (2011)
  UCHIME improves sensitivity and speed of chimera detection
  Bioinformatics, 27, 16, 2194-2200
  https://doi.org/10.1093/bioinformatics/btr381
*/

/* global constants/data, no need for synchronization */
constexpr auto maxparts = 100;
constexpr auto window = 32;
constexpr auto few = 4;
constexpr auto maxcandidates = few * maxparts;
constexpr auto rejects = 16;
constexpr auto chimera_id = 0.55;
/* mutex_output, fp_uchimealns and fp_uchimeout are no longer file-static: they
   live in chimera_cli_state_s and are used only by the CLI report writers
   (print_report/print_report_long, called from process_query); the detection
   core (eval_parents/eval_parents_long) writes no output on either path
   (E6). mutex_input is
   likewise not here — it only serializes input reading on the CLI path and is
   owned as a local by chimera_threads_run(). */

/* information for each query sequence to be checked */
/* The parent segments that tile the query: starts[nth] and lengths[nth]
   describe parent nth. Kept as one value so the two same-typed views cannot be
   passed in the wrong order. */
struct parent_tiling_s {
  View<int> starts;
  View<int> lengths;
};


/* The four per-query rows of an alignment block, all the same length. */
struct alignment_rows_s {
  Span<char> query;
  Span<char> model;
  Span<char> diffs;
  Span<char> votes;
};


/* The figures eval_parents() or eval_parents_long() computed for the current
   query, kept for the report writers (--uchimealns and --uchimeout, or
   --alnout and --tabbedout for --chimeras_denovo). Detection fills them and
   writes nothing; the writers run later, when the query's result is output.
   Valid when the detection status is low_score or higher. The output step
   reads it from the query's chimera_query_result_s, which also carries
   best_h and, when needed, the alignment rows. Named after the matching
   chimera_result_s fields. */
struct chimera_report_s {
  int parent_a = 0;  /* seqno of the parent printed as A */
  int parent_b = 0;
  int parent_c = -1;  /* third parent (eval_parents_long only), or -1 */
  bool parents_swapped = false;  /* paln[1] is parent A (eval_parents only) */
  double id_query_model = 0.0;
  double id_query_a = 0.0;
  double id_query_b = 0.0;
  double id_query_c = 0.0;  /* eval_parents_long only */
  double id_a_b = 0.0;  /* eval_parents only */
  double id_query_top = 0.0;
  double divergence = 0.0;  /* QM - QT, eval_parents only */
  double divergence_percent = 0.0;  /* 100 * (QM - QT) / QT */
  int left_yes = 0;  /* the six vote counts: eval_parents only */
  int left_no = 0;
  int left_abstain = 0;
  int right_yes = 0;
  int right_no = 0;
  int right_abstain = 0;
};


struct chimera_info_s
{
  /* run configuration, set by chimera_thread_init and read by the detection
     core instead of the opt_* globals (E1 shared-infra phase); a pointer so
     chimera_info_s stays default-constructible for the arrays that hold it.
     On the CLI path it points at chimera_cli_state_s::detection_parameters; on
     the library path it points at detection_parameters below. The detection
     knobs (maxaccepts/maxrejects/id/self/selfid/threads/maxsizeratio/weak_id)
     are carried in that copy, so neither path mutates the opt_* globals (E1). */
  struct Parameters const * parameters = nullptr;

  /* which of the five chimera algorithms this detection is running, set by
     chimera_thread_init beside parameters -- from the command on the CLI path,
     from the caller's argument on the library one. The detection core reads it
     where the five differ (the --chimeras_denovo partitioning, the
     uchime2/uchime3 scoring rule) instead of asking which opt_<command>
     pointer is non-null. The initializer matches the default argument of
     chimera_detect_thread_init(), which is the reference-based detection this
     API has always documented (LIBRARY_API.md). */
  ChimeraMode mode = ChimeraMode::uchime_ref;

  /* the chimera-detection configuration for the library path: a copy of the
     caller's Parameters with the detection knobs (maxaccepts/maxrejects/id, the
     weak_id clamp, and in denovo mode self/selfid/maxsizeratio) applied, built
     by chimera_detect_thread_init. The detection core reads these through
     `parameters` above. Unused on the CLI path, which threads the equivalent
     copy held in chimera_cli_state_s. */
  struct Parameters detection_parameters;

  /* the sequence database this thread queries, installed by chimera_thread_init
     and read by the detection core (realloc_arrays, find_matches, ...) instead
     of a process-wide global; a pointer so chimera_info_s stays
     default-constructible. Must outlive ci. */
  struct Database const * db = nullptr;

  int query_no = 0;
  std::vector<char> query_head_v;  /* owned header storage, grown monotonically
                                      to the longest header seen, so it is
                                      generally longer than the header */
  View<char> query_head {nullptr, 0};  /* the header itself: a view into
                                          query_head_v, cut to the current
                                          header's length */
  int64_t query_size = 0;
  std::vector<char> query_seq;
  int query_len = 0;

  std::array<struct searchinfo_s, maxparts> si {{}};

  std::array<unsigned int, maxcandidates> cand_list {{}};
  int cand_count = 0;

  std::unique_ptr<s16info_s, s16info_deleter> s;  /* SIMD aligner instance (owned) */
  std::array<CELL, maxcandidates> snwscore {{}};
  std::array<unsigned short, maxcandidates> snwalignmentlength {{}};
  std::array<unsigned short, maxcandidates> snwmatches {{}};
  std::array<unsigned short, maxcandidates> snwmismatches {{}};
  std::array<unsigned short, maxcandidates> snwgaps {{}};
  std::array<int64_t, maxcandidates> nwscore {{}};
  std::array<int64_t, maxcandidates> nwalignmentlength {{}};
  std::array<int64_t, maxcandidates> nwmatches {{}};
  std::array<int64_t, maxcandidates> nwmismatches {{}};
  std::array<int64_t, maxcandidates> nwgaps {{}};
  std::vector<std::string> nwcigar = std::vector<std::string>(maxcandidates);

  std::vector<int> match;
  std::vector<int> insert;
  std::vector<int> smooth;
  std::vector<int> maxsmooth;

  std::vector<double> scan_p;
  std::vector<double> scan_q;

  int parents_found = 0;
  std::array<int, maxparents> best_parents {{}};
  std::array<int, maxparents> best_start {{}};
  std::array<int, maxparents> best_len {{}};

  int best_target = 0;
  char * best_cigar = nullptr;

  std::vector<int> maxi;  // longest insertion per position
  std::vector<std::vector<char>> paln;
  std::vector<char> qaln;
  std::vector<char> diffs;
  std::vector<char> votes;
  std::vector<char> model;
  std::vector<bool> ignore;

  double best_h = 0;
  struct chimera_report_s report;  /* filled by eval_parents*, read by the report writers */

  int parts = 0;  /* number of query parts for chimera detection */

  /* si[0 .. parts_ready) have been through query_init(). The rest are built on
     demand by chimera_process_query, because each one owns an aligner and a
     hit buffer, and a uchime run uses four of the maxparts. (Each also used
     to own a k-mer counter array with one entry per database sequence; they
     now share kmer_counters below.) dbindex/tophits are the two query_init()
     arguments that are not already reachable from ci. */
  int parts_ready = 0;
  struct Dbindex const * dbindex = nullptr;
  int tophits = 0;

  /* The k-mer counters of this thread's part searches: one array, with one
     entry per database sequence, lent to each part's searchinfo_s while it
     is searched (see CounterLoan). The parts are searched one after the
     other, and search_topscores() clears the counters it reads, so they
     can share it: one array per thread instead of one per part, of which
     --chimeras_denovo uses up to 100. Sized by chimera_thread_init. */
  std::vector<count_t> kmer_counters;

  /* API result fields — populated by eval_parents when result_out is set */
  struct chimera_result_s * result_out = nullptr;

  /* API per-thread working state — initialized by chimera_detect_init,
     reused across chimera_detect_single calls */
  std::vector<struct hit> api_allhits_list;
  std::unique_ptr<LinearMemoryAligner> api_lma_ptr;

  /* Views cut to the current record. The buffers above are high-water marks
     reused across queries, so every reader has to cut them to the length that
     is live right now; these do it in one place. Past the cut they hold stale
     data from earlier records (or, for maxi, zeros), and an index beyond the
     cut now asserts rather than reading it -- which is the disagreement that
     cost a heap overflow in fill_in_alignment_string_for_query.

     Accessors rather than stored members on purpose: a view cached across a
     grow_to_fit dangles when the buffer reallocates, which realloc_arrays
     already has to work around for si[].qsequence. */
  auto query() const -> View<char> {
    return make_view(query_seq).first(static_cast<std::size_t>(query_len));
  }
  // one insertion count per query position, plus one for the terminal run
  auto insertions() const -> View<int> {
    return make_view(maxi).first(static_cast<std::size_t>(query_len) + 1);
  }
  auto candidates() const -> View<unsigned int> {
    return make_view(cand_list).first(static_cast<std::size_t>(cand_count));
  }
  auto parents() const -> View<int> {
    return make_view(best_parents).first(static_cast<std::size_t>(parents_found));
  }
  /* The alignment rows, cut to the alignment being built. alnlen is passed in
     rather than stored: find_total_alignment_length() derives it from maxi,
     and a cached copy would be exactly the second length this struct is trying
     to stop keeping. The rows are written to, hence Spans. */
  auto rows(int const alnlen) -> struct alignment_rows_s {
    auto const length = static_cast<std::size_t>(alnlen);
    return {make_span(qaln).first(length), make_span(model).first(length),
            make_span(diffs).first(length), make_span(votes).first(length)};
  }
  auto parent_row(int const nth, int const alnlen) -> Span<char> {
    return make_span(paln[static_cast<std::size_t>(nth)])
             .first(static_cast<std::size_t>(alnlen));
  }
  auto tiling() const -> struct parent_tiling_s {
    auto const count = static_cast<std::size_t>(parents_found);
    return {make_view(best_start).first(count), make_view(best_len).first(count)};
  }
};


/* Per-invocation CLI state for the chimera() command — the statics used only
   by the CLI path: the query file handle, progress, the six stats
   counters/abundances, the chimeras/nonchimeras/borderline output handles, the
   per-thread chimera_info array, and the report output handles
   (fp_uchimealns/fp_uchimeout) with the mutex serializing all CLI writes.
   Threaded through chimera_threads_run() and chimera_thread_core() so the CLI
   command is reentrant (E4). The detection core (eval_parents /
   eval_parents_long, reached from both the CLI and the library
   chimera_detect_single) never sees this state: it leaves its figures in
   chimera_info_s::report, and the CLI writes them from process_query (E6
   split of the core from its CLI output). */
struct chimera_cli_state_s
{
  /* the run configuration, threaded through the CLI-path helpers instead of the
     opt_* globals (E1/F3); set once at construction, read-only thereafter. */
  struct Parameters const & parameters;
  /* which of the five chimera commands this run is; the engine reads it
     wherever the five behave differently, instead of asking which
     opt_<command> pointer is non-null */
  ChimeraMode const mode;
  /* a copy of parameters with the chimera-detection knobs applied (maxaccepts/
     maxrejects/id and, in denovo mode, self/selfid/threads/maxsizeratio); built
     once in chimera(). The detection core reads these through si->parameters and
     chimera_threads_run sizes the pool from opt_threads here, so chimera() no
     longer mutates the opt_* globals (E1 trap-writer). */
  struct Parameters detection_parameters;
  /* the in-memory sequence database this run owns (RAII); declared before
     dbindex so it is destroyed after it, and si->db points here */
  struct Database db;
  struct Dbindex dbindex;  /* the k-mer index this run owns (RAII); si->dbindex points here */
  std::mutex mutex_output;  /* serializes all CLI output + stats updates */
  std::FILE * fp_uchimealns = nullptr;
  std::FILE * fp_uchimeout = nullptr;
  std::unique_ptr<fastx_s> query_fasta_h;
  unsigned int seqno = 0;
  uint64_t progress = 0;
  int chimera_count = 0;
  int nonchimera_count = 0;
  int borderline_count = 0;
  int total_count = 0;
  int64_t chimera_abundance = 0;
  int64_t nonchimera_abundance = 0;
  int64_t borderline_abundance = 0;
  int64_t total_abundance = 0;
  std::FILE * fp_chimeras = nullptr;
  std::FILE * fp_nonchimeras = nullptr;
  std::FILE * fp_borderline = nullptr;
  std::vector<struct chimera_info_s> cia;  /* one per worker thread */

  Progress * progress_bar = nullptr;  /* owner progress bar; worker updates it under output_lock (state.progress is the counter) */

  chimera_cli_state_s(struct Parameters const & params, ChimeraMode const chimera_mode)
    : parameters(params), mode(chimera_mode) {}
};


enum struct Status : unsigned char {
  no_parents,   // (0) non-chimeric
  no_alignment, // (1) score < 0, non-chimeric
  low_score,    // (2) score < minh, non-chimeric
  suspicious,   // (3) score >= minh, not available with uchime2_denovo and uchime3_denovo
  chimeric,      // (4) score >= minh && divdiff >= opt_mindiv && ...
};


/* The rows of one query's --uchimealns (--alnout) block, copied out of the
   worker's chimera_info_s, which the worker's next query overwrites. All
   rows have the alignment's length. */
struct chimera_alignment_s {
  std::vector<char> query;
  std::vector<std::vector<char>> parents;  /* in chimera_info_s::paln order */
  std::vector<int> parent_seqnos;  /* parallel to parents */
  std::vector<char> diffs;
  std::vector<char> votes;  /* eval_parents only; empty otherwise */
  std::vector<char> model;
};


/* One part search of a query, as the denovo batch driver needs it to check
   a result computed against an index that missed some earlier queries: the
   heap threshold the search left (see TopscoresThreshold) and the part's
   unique k-mers. */
struct chimera_part_search_s {
  struct TopscoresThreshold threshold;
  std::vector<unsigned int> kmers;
};


/* Everything the CLI writes or counts for one query, taken from the worker's
   chimera_info_s once detection is over, so that the output step reads
   this and nothing else. header and sequence are views: into the database
   in denovo mode, into the worker's chimera_info_s for --uchime_ref (valid
   until that worker claims its next query). The alignment rows are copied
   only for a chimera, and only when --uchimealns (--alnout) is set. */
struct chimera_query_result_s {
  Status status = Status::no_parents;
  unsigned int seqno = 0;  /* denovo: the query's database index */
  uint64_t query_position = 0;  /* --uchime_ref: progress in the query file */
  int64_t abundance = 0;
  View<char> header {nullptr, 0};
  View<char> sequence {nullptr, 0};
  double best_h = 0.0;
  struct chimera_report_s report;
  struct chimera_alignment_s alignment;
  /* denovo batch driver only: one per part search, none when the query was
     too short to be split and searched */
  std::vector<struct chimera_part_search_s> part_searches;
  /* owned copies of the header and sequence, taken when the result must
     outlive its worker's next claim (keep_query_text) */
  std::vector<char> header_copy;
  std::vector<char> sequence_copy;
};


/* Puts the results of the worker pool back into query order: a worker
   writes its query if it is the next one due, then the queries that were
   waiting for it; otherwise it leaves its result here and moves on. No
   worker waits for another. Guarded by chimera_cli_state_s::mutex_output. */
struct ordered_output_s {
  unsigned int next_rank = 0;  /* the claim rank of the next query to write */
  std::map<unsigned int, struct chimera_query_result_s> waiting;
};


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* Copy a stored header or sequence into a scratch buffer reused across
     queries. Nothing in the detection core reads a terminator any more -- every
     consumer of these copies is driven by a View's size -- so the copy is
     exactly source.size() bytes and the buffer needs no room past them. */
  // Returns a view over the copy, so the caller can store the bytes and the
  // length as one value instead of tracking a separate length field.
  auto copy_into_scratch(View<char> const source,
                         std::vector<char> & destination) -> View<char> {
    assert(destination.size() >= source.size());
    std::copy(source.cbegin(), source.cend(), destination.begin());
    return make_view(destination).first(source.size());
  }

  /* Store a label in one of the library result struct's fixed-size char
     arrays, truncating it if it does not fit and always terminating it -- the
     same contract snprintf("%.*s") had, expressed once here instead of at each
     of the seven sites, and without narrowing the length to the int that the
     "%.*s" precision requires. */
  template <std::size_t Capacity>
  auto copy_label(std::array<char, Capacity> & label, View<char> const text) -> void {
    static_assert(Capacity > 0, "a label buffer must have room for its terminator");
    auto const stored = std::min(text.size(), Capacity - 1);
    std::copy_n(text.cbegin(), stored, label.begin());
    label[stored] = '\0';
  }

  /* One sequence row of an alignment block in the --uchimealns report: a
     one-character label ('Q', 'A', 'B'), the row's start position, one block's
     worth of the row, then its end position. The row is a window into a longer
     buffer and is not terminated at the block's width, hence the View. */
  auto print_alignment_row(std::FILE * output_handle, char const label,
                           int const start, View<char> const row,
                           int const end) -> void {
    fprint(output_handle, static_cast<char>(label));
    fprint(output_handle, ' ');
    fprint_integer(output_handle, start, 5);
    fprint(output_handle, ' ');
    fprint(output_handle, row);
    fprint(output_handle, ' ');
    fprint_integer(output_handle, end);
    fprint(output_handle, '\n');
  }

  /* An annotation row of the same block (Diffs, Votes, Model): no positions,
     just the name padded to the width of the label column above. */
  auto print_annotation_row(std::FILE * output_handle, char const * padded_name,
                            View<char> const row) -> void {
    std::fputs(padded_name, output_handle);
    fprint(output_handle, row);
    fprint(output_handle, '\n');
  }

}  // end of anonymous namespace


namespace {
// header_length is passed in rather than read from chimera_info: it sizes
// query_head_v, so it has to be known before the buffer exists and therefore
// before chimera_info->query_head (a view into that buffer) can be set.
auto realloc_arrays(struct chimera_info_s * chimera_info, struct Database const & db,
                    std::size_t const header_length) -> void
{
  struct Parameters const & parameters = *chimera_info->parameters;
  if (chimera_info->mode == ChimeraMode::chimeras_denovo)
    {
      if (parameters.opt_chimeras_parts == 0) {
        chimera_info->parts = (chimera_info->query_len + 99) / 100;
      }
      else {
        chimera_info->parts = parameters.opt_chimeras_parts;
      }
      if (chimera_info->parts < 2) {
        chimera_info->parts = 2;
      }
      else if (chimera_info->parts > maxparts) {
        chimera_info->parts = maxparts;
      }
    }
  else
    {
      /* default for uchime, uchime2, and uchime3 */
      chimera_info->parts = 4;
    }

  int const maxhlen = std::max(static_cast<int>(header_length), 1);
  vsearch::grow_to_fit(chimera_info->query_head_v, static_cast<size_t>(maxhlen));

  /* realloc arrays based on query length */

  int const maxqlen = std::max(chimera_info->query_len, 1);
  auto const max_2x2_size = static_cast<size_t>(maxcandidates) * static_cast<size_t>(maxqlen);
  int64_t const maxalnlen = static_cast<int64_t>(maxqlen) + (2 * static_cast<int64_t>(db.getlongestsequence()));

  /* Each buffer states the size it needs. They all grow monotonically with
     maxqlen, so the single query_alloc guard they used to share was correct --
     but it hid the fact that they ask for four different sizes, and it kept a
     copy of a length the vectors already know. */
  vsearch::grow_to_fit(chimera_info->query_seq, static_cast<size_t>(maxqlen));
  vsearch::grow_to_fit(chimera_info->maxi, static_cast<size_t>(maxqlen) + 1);
  vsearch::grow_to_fit(chimera_info->maxsmooth, static_cast<size_t>(maxqlen));
  vsearch::grow_to_fit(chimera_info->match, max_2x2_size);
  vsearch::grow_to_fit(chimera_info->insert, max_2x2_size);
  vsearch::grow_to_fit(chimera_info->smooth, max_2x2_size);

  vsearch::grow_to_fit(chimera_info->scan_p, static_cast<size_t>(maxqlen) + 1);
  vsearch::grow_to_fit(chimera_info->scan_q, static_cast<size_t>(maxqlen) + 1);

  vsearch::grow_to_fit(chimera_info->paln, maxparents);
  for (auto & a_parent_alignment : chimera_info->paln) {
    vsearch::grow_to_fit(a_parent_alignment, static_cast<size_t>(maxalnlen) + 1);
  }
  vsearch::grow_to_fit(chimera_info->qaln, static_cast<size_t>(maxalnlen) + 1);
  vsearch::grow_to_fit(chimera_info->diffs, static_cast<size_t>(maxalnlen) + 1);
  vsearch::grow_to_fit(chimera_info->votes, static_cast<size_t>(maxalnlen) + 1);
  vsearch::grow_to_fit(chimera_info->model, static_cast<size_t>(maxalnlen) + 1);
  vsearch::grow_to_fit(chimera_info->ignore, static_cast<size_t>(maxalnlen) + 1);

  // resize query parts if longer than earlier, minimum 100
  int const maxpartlen =
    std::max((maxqlen + chimera_info->parts - 1) / chimera_info->parts, 100);
  for (auto & query_info: chimera_info->si)
    {
      /* the span is reset unconditionally now that the growth is: a grow may
         reallocate and leave it dangling, and partition_query -- the first
         thing chimera_process_query does -- sets it for real before any reader
         sees it */
      vsearch::grow_to_fit(query_info.qsequence_v, static_cast<size_t>(maxpartlen));
      query_info.qsequence = make_span(query_info.qsequence_v).first(0);
    }
}


auto reset_matches(struct chimera_info_s * a_chimera_info) -> void {
  // refactoring: initialization to zero? (useless), or reset to zero??
  /* match and insert are row-major, one row of query_len entries per
     candidate, grown once to the high-water mark maxcandidates *
     longest-query-seen. Only the first cand_count * query_len entries are ever
     written (find_matches) or read (find_best_parents, find_best_parents_long),
     so clearing the whole allocation is wasted: at the uchime defaults
     cand_count averages ~9 of the 400 rows. */
  auto const live = static_cast<std::size_t>(a_chimera_info->cand_count) *
                    static_cast<std::size_t>(a_chimera_info->query_len);
  assert(live <= a_chimera_info->match.size());
  assert(live <= a_chimera_info->insert.size());
  std::fill_n(a_chimera_info->match.begin(), live, 0);
  std::fill_n(a_chimera_info->insert.begin(), live, 0);
}


auto find_matches(struct chimera_info_s * chimera_info, struct Database const & db) -> void
{
  /* find the positions with matches for each potential parent */
  /* also note the positions with inserts in front */

  auto const & qseq = chimera_info->query_seq;

  for (auto i = 0; i < chimera_info->cand_count; ++i)
    {
      auto const tseq = db.sequence_view(chimera_info->cand_list[static_cast<size_t>(i)]);

      auto qpos = 0;
      auto tpos = 0;

      auto const & cigar = chimera_info->nwcigar[static_cast<size_t>(i)];
      auto const cigar_pairs = parse_cigar_string(make_view(cigar));

      for (auto const & a_pair: cigar_pairs) {
        auto const operation = a_pair.first;
        auto const runlength = a_pair.second;
        switch (operation) {
        case Operation::match:
          for (auto j = 0; j < runlength; ++j)
            {
              if ((four_bit::map(qseq[static_cast<size_t>(qpos)]) &
                   four_bit::map(tseq[static_cast<std::size_t>(tpos)])) != 0U)
                {
                  chimera_info->match[static_cast<size_t>((i * chimera_info->query_len) + qpos)] = 1;
                }
              ++qpos;
              ++tpos;
            }
          break;

        case Operation::insertion:
          chimera_info->insert[static_cast<size_t>((i * chimera_info->query_len) + qpos)] = static_cast<int>(runlength);
          tpos += static_cast<int>(runlength);
          break;

        case Operation::deletion:
          qpos += static_cast<int>(runlength);
          break;
        }
      }
    }
}
}  // anonymous namespace


struct parents_info_s
{
  int cand = -1;
  int start = -1;
  int len = 0;
};


namespace {
auto scan_matches(struct chimera_info_s * ci,
                  int const * matches,
                  int const len,
                  double const percentage,
                  int & best_start,
                  int & best_len) -> bool
{
  /*
    Scan matches array of zeros and ones, and find the longest subsequence
    having a match fraction above or equal to the given percentage (e.g. 2%).
    Based on an idea of finding the longest positive sum substring:
    https://stackoverflow.com/questions/28356453/longest-positive-sum-substring
    If the percentage is 2%, matches are given a score of 2 and mismatches -98.
  */

  auto const score_match = percentage;
  auto const score_mismatch = percentage - 100.0;

  auto & p = ci->scan_p;
  auto & q = ci->scan_q;

  p[0] = 0.0;
  for (auto i = 0; i < len; ++i) {
    p[static_cast<size_t>(i) + 1] = p[static_cast<size_t>(i)] + ((matches[i] != 0) ? score_match : score_mismatch);
  }

  q[static_cast<size_t>(len)] = p[static_cast<size_t>(len)];
  for (auto i = len - 1; i >= 0; --i) {
    q[static_cast<size_t>(i)] = std::max(q[static_cast<size_t>(i) + 1], p[static_cast<size_t>(i)]);
  }

  auto best_i = 0;
  auto best_d = -1;
  auto best_c = -1.0;
  auto i = 1;
  auto j = 1;
  while (j <= len)
    {
      auto const c = q[static_cast<size_t>(j)] - p[static_cast<size_t>(i - 1)];
      if (c >= 0.0)
        {
          auto const d = j - i + 1;
          if (d > best_d)
            {
              best_i = i;
              best_d = d;
              best_c = c;
            }
          j += 1;
        }
      else
        {
          i += 1;
        }
    }

  if (best_c >= 0.0)
    {
      best_start = best_i - 1;
      best_len = best_d;
      return true;
    }
  return false;
}


auto find_best_parents_long(struct chimera_info_s * ci) -> int
{
  struct Parameters const & parameters = *ci->parameters;
  /* Find parents with longest matching regions, without indels, allowing
     a given percentage of mismatches (specified with --chimeras_diff_pct),
     and excluding regions matched by previously identified parents. */

  reset_matches(ci);
  find_matches(ci, *ci->db);

  std::vector<struct parents_info_s> best_parents(maxparents);
  std::vector<bool> position_used(static_cast<size_t>(ci->query_len), false);

  int pos_remaining = ci->query_len;
  int parents_found = 0;

  /* Stop at maxparents even if opt_chimeras_parents_max is larger: best_parents
     (above) and ci->best_parents/best_start/best_len are all maxparents-sized,
     so f must stay below maxparents. parameters_validate() already rejects an
     out-of-range option; this bound is the guard at the write site. */
  for (int f = 0; (f < parameters.opt_chimeras_parents_max) and (f < maxparents); ++f)
    {
      /* scan each candidate and find longest matching region */

      int best_start = 0;
      int best_len = 0;
      int best_cand = -1;

      for (int i = 0; i < ci->cand_count; ++i)
        {
          int j = 0;
          while (j < ci->query_len)
            {
              int const start = j;
              int len = 0;
              while ((j < ci->query_len) &&
                     (not position_used[static_cast<size_t>(j)]) &&
                     ((len == 0) or (ci->insert[static_cast<size_t>((i * ci->query_len) + j)] == 0)))
                {
                  ++len;
                  ++j;
                }
              if (len > best_len)
                {
                  int scan_best_start = 0;
                  int scan_best_len = 0;
                  if (scan_matches(ci,
                                   &ci->match[static_cast<size_t>((i * ci->query_len) + start)],
                                   len,
                                   parameters.opt_chimeras_diff_pct,
                                   scan_best_start,
                                   scan_best_len) and (scan_best_len > best_len))
                    {
                      best_cand = i;
                      best_start = start + scan_best_start;
                      best_len = scan_best_len;
                    }
                }
              ++j;
            }
        }

      if (best_len >= parameters.opt_chimeras_length_min)
        {
          best_parents[static_cast<size_t>(f)].cand = best_cand;
          best_parents[static_cast<size_t>(f)].start = best_start;
          best_parents[static_cast<size_t>(f)].len = best_len;
          ++parents_found;

          /* mark positions used */
          for (int j = best_start; j < best_start + best_len; ++j)
            {
              position_used[static_cast<size_t>(j)] = true;
            }
          pos_remaining -= best_len;
        }
      else {
        break;
      }
    }

  /* sort parents by position */
  std::sort(best_parents.begin(),
            best_parents.begin() + parents_found,
            [](parents_info_s const & lhs, parents_info_s const & rhs) -> bool
            { return lhs.start < rhs.start; });

  ci->parents_found = parents_found;

  for (int f = 0; f < parents_found; ++f)
    {
      ci->best_parents[static_cast<size_t>(f)] = best_parents[static_cast<size_t>(f)].cand;
      ci->best_start[static_cast<size_t>(f)] = best_parents[static_cast<size_t>(f)].start;
      ci->best_len[static_cast<size_t>(f)] = best_parents[static_cast<size_t>(f)].len;
    }

  return static_cast<int>((parents_found > 1) and (pos_remaining == 0));
}


auto find_best_parents(struct chimera_info_s * ci) -> int
{
  reset_matches(ci);
  find_matches(ci, *ci->db);

  std::array<int, maxparents> best_parent_cand {{}};

  for (int f = 0; f < 2; ++f)
    {
      best_parent_cand[static_cast<size_t>(f)] = -1;
      ci->best_parents[static_cast<size_t>(f)] = -1;
    }

  std::vector<bool> cand_selected(static_cast<size_t>(ci->cand_count), false);

  for (int f = 0; f < 2; ++f)
    {
      if (f > 0)
        {
          /* for all parents except the first */

          /* wipe out matches for all candidates in positions
             covered by the previous parent */

          /* Every winning qpos clears the window ending at it, and qpos only
             increases, so consecutive winners ask for windows overlapping in
             all but one position. Clearing a position is idempotent, so
             skipping what an earlier window already cleared leaves exactly the
             same array: wiped_upto is the first position not yet cleared. */
          int wiped_upto = 0;
          for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
            {
              int const z = (best_parent_cand[static_cast<size_t>(f - 1)] * ci->query_len) + qpos;
              if (ci->smooth[static_cast<size_t>(z)] == ci->maxsmooth[static_cast<size_t>(qpos)])
                {
                  int const first = std::max(qpos + 1 - window, wiped_upto);
                  /* wiped_upto is a previous qpos plus one and qpos grows, so
                     the window is never entirely behind the cleared prefix */
                  assert(first <= qpos);
                  for (int j = 0; j < ci->cand_count; ++j)
                    {
                      auto const row = static_cast<std::ptrdiff_t>(j) *
                                       static_cast<std::ptrdiff_t>(ci->query_len);
                      std::fill(std::next(ci->match.begin(), row + first),
                                std::next(ci->match.begin(), row + qpos + 1),
                                0);
                    }
                  wiped_upto = qpos + 1;
                }
            }
        }


      /* Compute smoothed score in a 32bp window for each candidate. */
      /* Record max smoothed score for each position among candidates left. */

      /* a reset, not an initialization: maxsmooth is a per-thread buffer sized
         once to maxqlen (see the resize in chimera_thread_init) and reused
         across queries and across the rounds of this parent search, so it
         still holds the previous round's maxima here. The whole buffer is
         cleared rather than its first query_len entries: the tail past this
         query's length is never read, so either would do, and clearing all of
         it keeps the reset independent of the query at hand. */
      std::fill(ci->maxsmooth.begin(), ci->maxsmooth.end(), 0);

      for (int i = 0; i < ci->cand_count; ++i)
        {
          if (not cand_selected[static_cast<size_t>(i)])
            {
              int sum = 0;
              for (int qpos = 0; qpos < ci->query_len; ++qpos)
                {
                  size_t const z = (static_cast<size_t>(i) * static_cast<size_t>(ci->query_len)) + static_cast<size_t>(qpos);
                  sum += ci->match[z];
                  if (qpos >= window)
                    {
                      sum -= ci->match[z - static_cast<size_t>(window)];
                    }
                  if (qpos >= window - 1)
                    {
                      ci->smooth[z] = sum;
                      ci->maxsmooth[static_cast<size_t>(qpos)] = std::max(ci->smooth[z], ci->maxsmooth[static_cast<size_t>(qpos)]);
                    }
                }
            }
        }


      /* find parent with the most wins */

      std::vector<int> wins(static_cast<size_t>(ci->cand_count), 0);

      for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
        {
          if (ci->maxsmooth[static_cast<size_t>(qpos)] != 0)
            {
              for (int i = 0; i < ci->cand_count; ++i)
                {
                  if (not cand_selected[static_cast<size_t>(i)])
                    {
                      size_t const z = (static_cast<size_t>(i) * static_cast<size_t>(ci->query_len)) + static_cast<size_t>(qpos);
                      if (ci->smooth[z] == ci->maxsmooth[static_cast<size_t>(qpos)])
                        {
                          ++wins[static_cast<size_t>(i)];
                        }
                    }
                }
            }
        }

      /* select best parent based on most wins */

      int maxwins = 0;
      for (int i = 0; i < ci->cand_count; ++i)
        {
          int const w = wins[static_cast<size_t>(i)];
          if (w > maxwins)
            {
              maxwins = w;
              best_parent_cand[static_cast<size_t>(f)] = i;
            }
        }

      /* terminate loop if no parent found */

      if (best_parent_cand[static_cast<size_t>(f)] < 0) {
        break;
      }

      ci->best_parents[static_cast<size_t>(f)] = best_parent_cand[static_cast<size_t>(f)];
      cand_selected[static_cast<size_t>(best_parent_cand[static_cast<size_t>(f)])] = true;
    }

  /* Check if at least 2 candidates selected */

  return static_cast<int>((best_parent_cand[0] >= 0) and (best_parent_cand[1] >= 0));
}


auto find_total_alignment_length(struct chimera_info_s const * chimera_info) -> int {
  // query_len, plus the sum of the longest insertion runs (I) for each position
  return std::accumulate(chimera_info->maxi.begin(),
                         chimera_info->maxi.end(),
                         chimera_info->query_len);
}


auto fill_max_alignment_length(struct chimera_info_s * chimera_info) -> void
{
  /* find max insertions in front of each position in the query sequence */

  std::fill(chimera_info->maxi.begin(), chimera_info->maxi.end(), 0);

  for (auto const best_parent : chimera_info->parents()) {
    auto pos = 0LL;
    auto const & cigar = chimera_info->nwcigar[static_cast<size_t>(best_parent)];
    auto const cigar_pairs = parse_cigar_string(make_view(cigar));

    for (auto const & a_pair: cigar_pairs) {
      auto const operation = a_pair.first;
      auto const runlength = a_pair.second;
      switch (operation) {
      case Operation::match:
      case Operation::deletion:
        pos += runlength;
        break;

      case Operation::insertion:
        assert(runlength <= std::numeric_limits<int>::max());
        chimera_info->maxi[static_cast<size_t>(pos)] = std::max(static_cast<int>(runlength), chimera_info->maxi[static_cast<size_t>(pos)]);
        break;
      }
    }
  }
}


/* Write a run of filler characters into an alignment row and return the
   advanced cursor. Every alignment row in this file is built the same way --
   a filler run, then one character -- so the capacity check lives here rather
   than being repeated, and correctly-sized rows cannot drift apart from the
   checks that guard them. Only the run itself is checked here: whatever is
   written after it goes through Span::operator[], which checks itself. */
auto fill_run(Span<char> const row, std::size_t const cursor,
              int const run_length, char const filler) -> std::size_t {
  assert(run_length >= 0);
  auto const length = static_cast<std::size_t>(run_length);
  assert(cursor + length <= row.size());
  std::fill_n(std::next(row.begin(), static_cast<std::ptrdiff_t>(cursor)),
              length, filler);
  return cursor + length;
}


auto fill_alignment_parents(struct chimera_info_s * ci, struct Database const & db,
                            int const alnlen) -> void
{
  /* fill in alignment strings for the parents */

  for (int i = 0; i < ci->parents_found; ++i)
    {
      /* cut to alnlen, the same bound the readers use, so a row that does not
         fill exactly trips an assert here rather than being read back short */
      auto const alignment = ci->parent_row(i, alnlen);
      int const cand = ci->best_parents[static_cast<size_t>(i)];
      int const target_seqno = static_cast<int>(ci->cand_list[static_cast<size_t>(cand)]);
      auto const target_seq = db.sequence_view(static_cast<uint64_t>(target_seqno));

      auto is_inserted = false;
      int qpos = 0;
      int tpos = 0;
      std::size_t alnpos = 0;

      auto const & cigar = ci->nwcigar[static_cast<size_t>(cand)];
      auto const cigar_pairs = parse_cigar_string(make_view(cigar));
      for (auto const & a_pair: cigar_pairs) {
        auto const operation = a_pair.first;
        auto const runlength = a_pair.second;
        switch (operation) {
        case Operation::insertion:
          for (int j = 0; j < ci->maxi[static_cast<size_t>(qpos)]; ++j)
            {
              if (j < runlength)
                {
                  alignment[alnpos] = upcase::map(target_seq[static_cast<std::size_t>(tpos)]);
                  ++tpos;
                  ++alnpos;
                }
              else
                {
                  alignment[alnpos] = '-';
                  ++alnpos;
                }
            }
          is_inserted = true;
          break;

        case Operation::match:
        case Operation::deletion:
          for (int j = 0; j < runlength; ++j)
            {
              if (not is_inserted)
                {
                  alnpos = fill_run(alignment, alnpos, ci->maxi[static_cast<size_t>(qpos)], '-');
                }

              if (operation == Operation::match)
                {
                  alignment[alnpos] = upcase::map(target_seq[static_cast<std::size_t>(tpos)]);
                  ++tpos;
                  ++alnpos;
                }
              else
                {
                  alignment[alnpos] = '-';
                  ++alnpos;
                }

              ++qpos;
              is_inserted = false;
            }
        }
      }

      /* add any gaps at the end */

      if (not is_inserted)
        {
          alnpos = fill_run(alignment, alnpos, ci->maxi[static_cast<size_t>(qpos)], '-');
        }

      assert(alnpos == alignment.size());  // the row filled exactly
    }
}


/* Expand a query into its alignment row: emit insertions[qpos] filler
   characters before each query position's character, then the terminal filler
   run, then a terminator.

   The query and the insertion counts arrive as views, so the bound comes from
   the data itself rather than from a query_len field that has to be kept in
   agreement with the separately-sized high-water buffers it indexes into --
   the agreement that broke when a range-for walked the whole query_seq buffer
   and wrote past the end of the row. The row is written through a Span, so its
   capacity is available at the point of writing and an overrun trips an assert
   in a debug build instead of smashing the heap. */
auto fill_in_alignment_string_for_query(View<char> const query,
                                        View<int> const insertions,
                                        Span<char> const alignment) -> void {
  // one insertion count per query position, plus one for the terminal run
  assert(insertions.size() == query.size() + 1);

  std::size_t alnpos = 0;
  for (std::size_t qpos = 0; qpos < query.size(); ++qpos) {
    // add insertion (if any):
    alnpos = fill_run(alignment, alnpos, insertions[qpos], '-');

    // add (mis-)matching position:
    alignment[alnpos] = upcase::map(query[qpos]);
    ++alnpos;
  }
  // add terminal gap (if any):
  alnpos = fill_run(alignment, alnpos, insertions[query.size()], '-');
  assert(alnpos == alignment.size());  // the row filled exactly
}


/* Fill the model row in lockstep with the query row: same bound, same checked
   writes, one letter per parent. qpos stays an int because the tiling compares
   it as a value, not as an index. */
auto fill_in_model_string_for_query(View<int> const insertions,
                                    struct parent_tiling_s const tiling,
                                    Span<char> const model) -> void {
  assert(tiling.starts.size() == tiling.lengths.size());
  auto const query_length = static_cast<int>(insertions.size()) - 1;
  auto const parents_found = static_cast<int>(tiling.starts.size());
  int nth_parent = 0;
  std::size_t alnpos = 0;
  for (int qpos = 0; qpos < query_length; ++qpos)
    {
      /* Advance to the next parent only while one exists. The parent segments
         are expected to tile the query exactly, in which case nth_parent
         reaches parents_found - 1 at the last segment and stops naturally. But
         if the tiling leaves a tail uncovered, the default-zero tail slots make
         "qpos >= best_start + best_len" (0 + 0) true at every remaining
         position, so an unclamped ++nth_parent would run past parents_found —
         reading best_start[]/best_len[] beyond the valid entries and eventually
         past the maxparents-element array, and emitting model letters beyond
         the last parent. Clamping keeps nth_parent in [0, parents_found - 1]
         and attributes any uncovered tail to the last parent (S19). */
      if ((nth_parent + 1 < parents_found) and
          (qpos >= (tiling.starts[static_cast<std::size_t>(nth_parent)]
                    + tiling.lengths[static_cast<std::size_t>(nth_parent)]))) {
        ++nth_parent;
      }
      // add insertion (if any):
      auto const parent_letter = static_cast<char>('A' + nth_parent);
      alnpos = fill_run(model, alnpos, insertions[static_cast<std::size_t>(qpos)], parent_letter);

      // add (mis-)matching position:
      model[alnpos] = parent_letter;
      ++alnpos;
    }
  // add terminal gap (if any):
  alnpos = fill_run(model, alnpos, insertions[insertions.size() - 1],
                    static_cast<char>('A' + nth_parent));
  assert(alnpos == model.size());  // the row filled exactly
}


auto count_matches_with_parents(struct chimera_info_s const * chimera_info,
                                int const alignment_length) -> std::array<int, maxparents> {
  std::array<int, maxparents> matches {{}};

  for (auto i = 0; i < alignment_length; ++i)
    {
      auto const qsym = four_bit::map(chimera_info->qaln[static_cast<size_t>(i)]);

      for (auto f = 0; f < chimera_info->parents_found; ++f)
        {
          auto const psym = four_bit::map(chimera_info->paln[static_cast<size_t>(f)][static_cast<size_t>(i)]);
          if (qsym == psym) {
            ++matches[static_cast<size_t>(f)];
          }
        }
    }
  return matches;
}


auto compute_global_similarities_with_parents(
    std::array<int, maxparents> const & match_counts,
    int const alignment_length) -> std::array<double, maxparents> {
  std::array<double, maxparents> similarities {{}};
  auto compute_percentage = [alignment_length](int const match_count) -> double {
    return 100.0 * match_count / alignment_length;
  };
  std::transform(match_counts.begin(), match_counts.end(),
                 similarities.begin(), compute_percentage);
  return similarities;
}


auto compute_diffs(struct chimera_info_s const * ci,
                   std::vector<unsigned char> const & psym,
                   unsigned char const qsym) -> char {
  auto const all_defined = (qsym != 0U) and
    std::all_of(psym.begin(),
                psym.end(),
                [](unsigned char const symbol) -> bool{ return symbol != 0U; });

  char diff = ' ';

  if (not all_defined) { return diff; }

  auto z = 0;
  for (auto f = 0; f < ci->parents_found; ++f) {
    if (psym[static_cast<size_t>(f)] == qsym) {
      diff = static_cast<char>('A' + f);
      ++z;
    }
  }
  if (z > 1) {
    diff = ' ';
  }
  return diff;
}


auto eval_parents_long(struct chimera_info_s * ci, struct Database const & db) -> Status
{
  /* always chimeric if called */
  auto const status = Status::chimeric;

  fill_max_alignment_length(ci);
  auto const alnlen = find_total_alignment_length(ci);

  fill_alignment_parents(ci, db, alnlen);

  auto const rows = ci->rows(alnlen);
  fill_in_alignment_string_for_query(ci->query(), ci->insertions(), rows.query);
  fill_in_model_string_for_query(ci->insertions(), ci->tiling(), rows.model);

  std::vector<unsigned char> psym;
  psym.reserve(maxparents);

  for (int i = 0; i < alnlen; ++i)
    {
      auto const qsym = four_bit::map(ci->qaln[static_cast<size_t>(i)]);
      for (int f = 0; f < ci->parents_found; ++f) {
        psym.emplace_back(four_bit::map(ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i)]));
      }

      /* lower case parent symbols that differ from query */

      for (int f = 0; f < ci->parents_found; ++f) {
        if ((psym[static_cast<size_t>(f)] != 0U) and (psym[static_cast<size_t>(f)] != qsym)) {
          ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i)] = to_lower(ci->paln[static_cast<size_t>(f)][static_cast<size_t>(i)]);
        }
      }

      /* compute diffs */
      ci->diffs[static_cast<size_t>(i)] = compute_diffs(ci, psym, qsym);
      psym.clear();
    }


  auto const match_QP = count_matches_with_parents(ci, alnlen);

  int const seqno_a = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[0])]);
  int const seqno_b = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[1])]);
  int const seqno_c = ci->parents_found > 2 ? static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[2])]) : -1;

  auto const QP = compute_global_similarities_with_parents(match_QP, alnlen);
  auto const QT = *std::max_element(QP.begin(), QP.end());

  double const QA = QP[0];
  double const QB = QP[1];
  double const QC = ci->parents_found > 2 ? QP[2] : 0.00;
  double const QM = 100.00;
  double const divfrac = 100.00 * (QM - QT) / QT;  // divergence of the model with the closest parent

  /* Populate API result struct if requested */
  if (ci->result_out != nullptr)
    {
      auto * r = ci->result_out;
      r->score = 99.9999;  /* chimeras_denovo always reports chimeric */
      auto const parent_a_header = db.header_view(static_cast<uint64_t>(seqno_a));
      auto const parent_b_header = db.header_view(static_cast<uint64_t>(seqno_b));
      copy_label(r->query_label, ci->query_head);
      copy_label(r->parent_a_label, parent_a_header);
      copy_label(r->parent_b_label, parent_b_header);
      /* closest parent = max of QA, QB */
      copy_label(r->closest_parent_label, (QA >= QB) ? parent_a_header : parent_b_header);
      r->id_query_model = QM;
      r->id_query_a = QA;
      r->id_query_b = QB;
      r->id_a_b = QC;  /* AB not computed in long path; use QC as proxy */
      r->id_query_top = QT;
      r->left_yes = 0;
      r->left_no = 0;
      r->left_abstain = 0;
      r->right_yes = 0;
      r->right_no = 0;
      r->right_abstain = 0;
      r->divergence = divfrac;
      r->flag = 'Y';  /* eval_parents_long is always chimeric */
    }

  ci->report.parent_a = seqno_a;
  ci->report.parent_b = seqno_b;
  ci->report.parent_c = seqno_c;
  ci->report.id_query_model = QM;
  ci->report.id_query_a = QA;
  ci->report.id_query_b = QB;
  ci->report.id_query_c = QC;
  ci->report.id_query_top = QT;
  ci->report.divergence_percent = divfrac;

  return status;
}


auto eval_parents(struct chimera_info_s * ci, struct Database const & db) -> Status
{
  struct Parameters const & parameters = *ci->parameters;
  auto status = Status::no_alignment;
  ci->parents_found = 2;

  fill_max_alignment_length(ci);
  auto const alnlen = find_total_alignment_length(ci);

  fill_alignment_parents(ci, db, alnlen);

  /* fill in alignment string for query */

  fill_in_alignment_string_for_query(ci->query(), ci->insertions(),
                                     ci->rows(alnlen).query);

  /* mark positions to ignore in voting */
  std::fill(ci->ignore.begin(), ci->ignore.end(), false);

  for (int i = 0; i < alnlen; ++i)
    {
      auto const qsym  = four_bit::map(ci->qaln[static_cast<size_t>(i)]);
      auto const p1sym = four_bit::map(ci->paln[0][static_cast<size_t>(i)]);
      auto const p2sym = four_bit::map(ci->paln[1][static_cast<size_t>(i)]);

      /* ignore gap positions and those next to the gap */
      if ((qsym == 0U) or (p1sym == 0U) or (p2sym == 0U))
        {
          ci->ignore[static_cast<size_t>(i)] = true;
          if (i > 0)
            {
              ci->ignore[static_cast<size_t>(i - 1)] = true;
            }
          if (i < alnlen - 1)
            {
              ci->ignore[static_cast<size_t>(i + 1)] = true;
            }
        }

      /* ignore ambiguous symbols */
      if (four_bit::is_ambiguous(qsym) or
          four_bit::is_ambiguous(p1sym) or
          four_bit::is_ambiguous(p2sym))
        {
          ci->ignore[static_cast<size_t>(i)] = true;
        }

      /* lower case parent symbols that differ from query */

      if ((p1sym != 0U) and (p1sym != qsym))
        {
          ci->paln[0][static_cast<size_t>(i)] = to_lower(ci->paln[0][static_cast<size_t>(i)]);
        }

      if ((p2sym != 0U) and (p2sym != qsym))
        {
          ci->paln[1][static_cast<size_t>(i)] = to_lower(ci->paln[1][static_cast<size_t>(i)]);
        }

      /* compute diffs */

      char diff = '\0';

      if ((qsym != 0U) and (p1sym != 0U) and (p2sym != 0U))
        {
          if (p1sym == p2sym)
            {
              if (qsym == p1sym)
                {
                  diff = ' ';
                }
              else
                {
                  diff = 'N';
                }
            }
          else
            {
              if (qsym == p1sym)
                {
                  diff = 'A';
                }
              else if (qsym == p2sym)
                {
                  diff = 'B';
                }
              else
                {
                  diff = '?';
                }
            }
        }
      else
        {
          diff = ' ';
        }

      ci->diffs[static_cast<size_t>(i)] = diff;
    }


  /* compute score */

  int sumA = 0;
  int sumB = 0;
  int sumN = 0;

  // refactoring: extract to function, use a struct to pass results
  // std::transform(ci->diffs.begin(),
  //                std::next(ci->diffs.begin(), alnlen),
  //                ci->ignore.begin(),
  //                [&sumA, &sumB, &sumN](char const diff, bool const is_ignored) -> void {
  //                         if (is_ignored) { return; }
  //                         if (diff == 'A') {
  //                             ++sumA;
  //                            }
  //                         else if (diff == 'B') {
  //                             ++sumB;
  //                            }
  //                         else if (diff != ' ') {
  //                             ++sumN;
  //                         }
  //                         return;
  //                 }
  //                );

  for (auto i = 0; i < alnlen; ++i)
    {
      if (ci->ignore[static_cast<size_t>(i)]) { continue; }
      auto const diff = ci->diffs[static_cast<size_t>(i)];

      if (diff == 'A')
        {
          ++sumA;
        }
      else if (diff == 'B')
        {
          ++sumB;
        }
      else if (diff != ' ')
        {
          ++sumN;
        }
    }

  int left_n = 0;
  int left_a = 0;
  int left_y = 0;
  int right_n = sumA;
  int right_a = sumN;
  int right_y = sumB;

  double best_h = -1;
  int best_i = -1;
  auto best_is_reverse = false;

  int best_left_y = 0;
  int best_right_y = 0;
  int best_left_n = 0;
  int best_right_n = 0;
  int best_left_a = 0;
  int best_right_a = 0;

  for (int i = 0; i < alnlen; ++i)
    {
      if (not ci->ignore[static_cast<size_t>(i)])
        {
          char const diff = ci->diffs[static_cast<size_t>(i)];
          if (diff != ' ')
            {
              if (diff == 'A')
                {
                  ++left_y;
                  --right_n;
                }
              else if (diff == 'B')
                {
                  ++left_n;
                  --right_y;
                }
              else
                {
                  ++left_a;
                  --right_a;
                }

              double left_h = 0;
              double right_h = 0;
              double h = 0;

              if ((left_y > left_n) and (right_y > right_n))
                {
                  left_h = left_y / ((parameters.opt_xn * (left_n + parameters.opt_dn)) + left_a);
                  right_h = right_y / ((parameters.opt_xn * (right_n + parameters.opt_dn)) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_is_reverse = false;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_n;
                      best_left_y = left_y;
                      best_left_a = left_a;
                      best_right_n = right_n;
                      best_right_y = right_y;
                      best_right_a = right_a;
                    }
                }
              else if ((left_n > left_y) and (right_n > right_y))
                {
                  /* swap left/right and yes/no */

                  left_h = left_n / ((parameters.opt_xn * (left_y + parameters.opt_dn)) + left_a);
                  right_h = right_n / ((parameters.opt_xn * (right_y + parameters.opt_dn)) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_is_reverse = true;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_y;
                      best_left_y = left_n;
                      best_left_a = left_a;
                      best_right_n = right_y;
                      best_right_y = right_n;
                      best_right_a = right_a;
                    }
                }
            }
        }
    }

  ci->best_h = best_h > 0 ? best_h : 0.0;

  if (best_h >= 0.0)
    {
      status = Status::low_score;

      /* flip A and B if necessary */

      if (best_is_reverse)
        {
          for (auto & diff : make_span(ci->diffs).first(static_cast<std::size_t>(alnlen)))
            {
              if (diff == 'A')
                {
                  diff = 'B';
                }
              else if (diff == 'B')
                {
                  diff = 'A';
                }
            }
        }

      /* fill in votes and model */

      for (int i = 0; i < alnlen; ++i)
        {
          char const m = i <= best_i ? 'A' : 'B';
          ci->model[static_cast<size_t>(i)] = m;

          char v = ' ';
          if (not ci->ignore[static_cast<size_t>(i)])
            {
              char const d = ci->diffs[static_cast<size_t>(i)];

              if ((d == 'A') or (d == 'B'))
                {
                  if (d == m)
                    {
                      v = '+';
                    }
                  else
                    {
                      v = '!';
                    }
                }
              else if ((d == 'N') or (d == '?'))
                {
                  v = '0';
                }
            }
          ci->votes[static_cast<size_t>(i)] = v;

          /* lower case diffs for no votes */
          if (v == '!')
            {
              ci->diffs[static_cast<size_t>(i)] = to_lower(ci->diffs[static_cast<size_t>(i)]);
            }
        }

      /* fill in crossover region */

      for (int i = best_i + 1; i < alnlen; ++i)
        {
          if ((ci->diffs[static_cast<size_t>(i)] == ' ') or (ci->diffs[static_cast<size_t>(i)] == 'A'))
            {
              ci->model[static_cast<size_t>(i)] = 'x';
            }
          else
            {
              break;
            }
        }

      ci->votes[static_cast<size_t>(alnlen)] = 0;
      ci->model[static_cast<size_t>(alnlen)] = 0;

      /* count matches */

      auto const index_a = best_is_reverse ? 1U : 0U;
      auto const index_b = best_is_reverse ? 0U : 1U;

      int match_QA = 0;
      int match_QB = 0;
      int match_AB = 0;
      int match_QM = 0;
      int cols = 0;

      for (auto i = 0; i < alnlen; i++)
        {
          if (not ci->ignore[static_cast<size_t>(i)])
            {
              ++cols;

              auto const qsym = four_bit::map(ci->qaln[static_cast<size_t>(i)]);
              auto const asym = four_bit::map(ci->paln[index_a][static_cast<size_t>(i)]);
              auto const bsym = four_bit::map(ci->paln[index_b][static_cast<size_t>(i)]);
              auto const msym = (i <= best_i) ? asym : bsym;

              if (qsym == asym)
                {
                  ++match_QA;
                }

              if (qsym == bsym)
                {
                  ++match_QB;
                }

              if (asym == bsym)
                {
                  ++match_AB;
                }

              if (qsym == msym)
                {
                  ++match_QM;
                }
            }
        }

      int const seqno_a = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[index_a])]);
      int const seqno_b = static_cast<int>(ci->cand_list[static_cast<size_t>(ci->best_parents[index_b])]);

      double const QA = 100.0 * match_QA / cols;
      double const QB = 100.0 * match_QB / cols;
      double const AB = 100.0 * match_AB / cols;
      double const QT = std::max(QA, QB);
      double const QM = 100.0 * match_QM / cols;
      double const divdiff = QM - QT;
      double const divfrac = 100.0 * divdiff / QT;

      int const sumL = best_left_n + best_left_a + best_left_y;
      int const sumR = best_right_n + best_right_a + best_right_y;

      if ((ci->mode == ChimeraMode::uchime2_denovo) or (ci->mode == ChimeraMode::uchime3_denovo))
        {
          // fix -Wfloat-equal: if match_QM == cols, then QM == 100.0
          if ((match_QM == cols) and (QT < 100.0))
            {
              status = Status::chimeric;
            }
        }
      else
        if (best_h >= parameters.opt_minh)
          {
            status = Status::suspicious;
            if ((divdiff >= parameters.opt_mindiv) and
                (sumL >= parameters.opt_mindiffs) and
                (sumR >= parameters.opt_mindiffs))
              {
                status = Status::chimeric;
              }
          }

      /* Populate API result struct if requested */
      if (ci->result_out != nullptr)
        {
          auto * r = ci->result_out;
          r->score = best_h;
          auto const parent_a_header = db.header_view(static_cast<uint64_t>(seqno_a));
          auto const parent_b_header = db.header_view(static_cast<uint64_t>(seqno_b));
          copy_label(r->query_label, ci->query_head);
          copy_label(r->parent_a_label, parent_a_header);
          copy_label(r->parent_b_label, parent_b_header);
          copy_label(r->closest_parent_label, (QA >= QB) ? parent_a_header : parent_b_header);
          r->id_query_model = QM;
          r->id_query_a = QA;
          r->id_query_b = QB;
          r->id_a_b = AB;
          r->id_query_top = QT;
          r->left_yes = best_left_y;
          r->left_no = best_left_n;
          r->left_abstain = best_left_a;
          r->right_yes = best_right_y;
          r->right_no = best_right_n;
          r->right_abstain = best_right_a;
          r->divergence = divdiff;
          r->flag = (status == Status::chimeric) ? 'Y' :
                    (status == Status::low_score ? 'N' : '?');
        }

      ci->report.parent_a = seqno_a;
      ci->report.parent_b = seqno_b;
      ci->report.parents_swapped = best_is_reverse;
      ci->report.id_query_model = QM;
      ci->report.id_query_a = QA;
      ci->report.id_query_b = QB;
      ci->report.id_a_b = AB;
      ci->report.id_query_top = QT;
      ci->report.divergence = divdiff;
      ci->report.divergence_percent = divfrac;
      ci->report.left_yes = best_left_y;
      ci->report.left_no = best_left_n;
      ci->report.left_abstain = best_left_a;
      ci->report.right_yes = best_right_y;
      ci->report.right_no = best_right_n;
      ci->report.right_abstain = best_right_a;
    }

  return status;
}


/* The --alnout and --tabbedout records of one --chimeras_denovo query, from
   the figures and rows its result carries. Called with the output lock
   held, only when eval_parents_long() ran (it always reports a chimera). */
auto print_report_long(struct chimera_cli_state_s const & cli,
                       struct chimera_query_result_s const & result,
                       struct Database const & db) -> void
{
  struct Parameters const & parameters = cli.detection_parameters;
  auto const & report = result.report;
  auto const & alignment = result.alignment;
  auto const status = result.status;
  auto const alnlen = static_cast<int>(alignment.query.size());
  auto const parent_count = static_cast<int>(alignment.parents.size());

  if ((parameters.opt_alnout != nullptr) and (status == Status::chimeric))
    {
      fprint(cli.fp_uchimealns, '\n');
      fprint(cli.fp_uchimealns, "----------------------------------------"
                                 "--------------------------------\n");
      fprint(cli.fp_uchimealns, "Query   (");
      fprint_integer(cli.fp_uchimealns, result.sequence.size(), 5);
      fprint(cli.fp_uchimealns, " nt) ");
      header_fprint_strip(cli.fp_uchimealns,
                          result.header,
                          attributes_to_strip(parameters));

      if (parent_count > maxparents)  // 20 parents max ('A' to 'U')
        {
          fatal("Internal error: chimera parents_found exceeds maxparents");
        }
      for (int f = 0; f < parent_count; ++f)
        {
          int const parent_seqno = alignment.parent_seqnos[static_cast<size_t>(f)];
          fprint(cli.fp_uchimealns, "\nParent");
          fprint(cli.fp_uchimealns, static_cast<char>('A' + f));
          fprint(cli.fp_uchimealns, " (");
          fprint_integer(cli.fp_uchimealns, db.getsequencelen(static_cast<uint64_t>(parent_seqno)), 5);
          fprint(cli.fp_uchimealns, " nt) ");
          header_fprint_strip(cli.fp_uchimealns,
                              db.header_view(static_cast<uint64_t>(parent_seqno)),
                              attributes_to_strip(parameters));
        }

      fprint(cli.fp_uchimealns, "\n\n");


      int const width = parameters.opt_alignwidth > 0 ? parameters.opt_alignwidth : alnlen;
      int qpos = 0;
      std::array<int, maxparents> ppos {{}};
      int rest = alnlen;

      for (int i = 0; i < alnlen; i += width)
        {
          /* count non-gap symbols on current line */

          int qnt = 0;
          std::array<int, maxparents> pnt {{}};

          int const w = std::min(rest, width);

          for (int j = 0; j < w; ++j)
            {
              if (alignment.query[static_cast<size_t>(i + j)] != '-')
                {
                  ++qnt;
                }

              for (int f = 0; f < parent_count; ++f) {
                if (alignment.parents[static_cast<size_t>(f)][static_cast<size_t>(i + j)] != '-')
                  {
                    ++pnt[static_cast<size_t>(f)];
                  }
              }
            }

          print_alignment_row(cli.fp_uchimealns, 'Q', qpos + 1,
                  View<char>{&alignment.query[static_cast<size_t>(i)], static_cast<std::size_t>(w)}, qpos + qnt);

          for (int f = 0; f < parent_count; ++f)
            {
              print_alignment_row(cli.fp_uchimealns, static_cast<char>('A' + f),
                      ppos[static_cast<size_t>(f)] + 1,
                      View<char>{&alignment.parents[static_cast<size_t>(f)][static_cast<size_t>(i)], static_cast<std::size_t>(w)},
                      ppos[static_cast<size_t>(f)] + pnt[static_cast<size_t>(f)]);
            }

          print_annotation_row(cli.fp_uchimealns, "Diffs   ",
                  View<char>{&alignment.diffs[static_cast<size_t>(i)], static_cast<std::size_t>(w)});
          print_annotation_row(cli.fp_uchimealns, "Model   ",
                  View<char>{&alignment.model[static_cast<size_t>(i)], static_cast<std::size_t>(w)});
          fprint(cli.fp_uchimealns, '\n');

          rest -= width;
          qpos += qnt;
          for (int f = 0; f < parent_count; ++f) {
            ppos[static_cast<size_t>(f)] += pnt[static_cast<size_t>(f)];
          }
        }

      fprint(cli.fp_uchimealns, "Ids.  QA ");
      std::fprintf(cli.fp_uchimealns, "%.2f", report.id_query_a);
      fprint(cli.fp_uchimealns, "%, QB ");
      std::fprintf(cli.fp_uchimealns, "%.2f", report.id_query_b);
      fprint(cli.fp_uchimealns, "%, QC ");
      std::fprintf(cli.fp_uchimealns, "%.2f", report.id_query_c);
      fprint(cli.fp_uchimealns, "%, QT ");
      std::fprintf(cli.fp_uchimealns, "%.2f", report.id_query_top);
      fprint(cli.fp_uchimealns, "%, QModel ");
      std::fprintf(cli.fp_uchimealns, "%.2f", report.id_query_model);
      fprint(cli.fp_uchimealns, "%, Div. ");
      std::fprintf(cli.fp_uchimealns, "%+.2f", report.divergence_percent);
      fprint(cli.fp_uchimealns, "%\n");
    }

  if (parameters.opt_tabbedout != nullptr)
    {
      std::fprintf(cli.fp_uchimeout, "%.4f", 99.9999);
      fprint(cli.fp_uchimeout, '\t');

      header_fprint_strip(cli.fp_uchimeout,
                          result.header,
                          attributes_to_strip(parameters));
      fprint(cli.fp_uchimeout, '\t');
      header_fprint_strip(cli.fp_uchimeout,
                          db.header_view(static_cast<uint64_t>(report.parent_a)),
                          attributes_to_strip(parameters));
      fprint(cli.fp_uchimeout, '\t');
      header_fprint_strip(cli.fp_uchimeout,
                          db.header_view(static_cast<uint64_t>(report.parent_b)),
                          attributes_to_strip(parameters));
      fprint(cli.fp_uchimeout, '\t');
      if (report.parent_c >= 0)
        {
          header_fprint_strip(cli.fp_uchimeout,
                              db.header_view(static_cast<uint64_t>(report.parent_c)),
                              attributes_to_strip(parameters));
        }
      else
        {
          fprint(cli.fp_uchimeout, '*');
        }
      fprint(cli.fp_uchimeout, '\t');

      std::fprintf(cli.fp_uchimeout,
              "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t"
              "%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%c\n",
              report.id_query_model,
              report.id_query_a,
              report.id_query_b,
              report.id_query_c,
              report.id_query_top,
              0, /* ignore, left yes */
              0, /* ignore, left no */
              0, /* ignore, left abstain */
              0, /* ignore, right yes */
              0, /* ignore, right no */
              0, /* ignore, right abstain */
              0.00,
              status == Status::chimeric ? 'Y' : (status == Status::low_score ? 'N' : '?'));
    }

}


/* The --uchimealns and --uchimeout records of one uchime query, from the
   figures and rows its result carries. Called with the output lock held,
   only when eval_parents() scored the query (status low_score or higher); a
   query with no parents or no alignment gets its --uchimeout line from
   output_query_result instead. */
auto print_report(struct chimera_cli_state_s const & cli,
                  struct chimera_query_result_s const & result,
                  struct Database const & db) -> void
{
  struct Parameters const & parameters = cli.detection_parameters;
  auto const & report = result.report;
  auto const & alignment = result.alignment;
  auto const status = result.status;
  auto const alnlen = static_cast<int>(alignment.query.size());
  int const sumL = report.left_no + report.left_abstain + report.left_yes;
  int const sumR = report.right_no + report.right_abstain + report.right_yes;

  /* print alignment */

  if ((parameters.opt_uchimealns != nullptr) and (status == Status::chimeric))
    {
      fprint(cli.fp_uchimealns, '\n');
      fprint(cli.fp_uchimealns, "----------------------------------------"
                                 "--------------------------------\n");
      fprint(cli.fp_uchimealns, "Query   (");
      fprint_integer(cli.fp_uchimealns, result.sequence.size(), 5);
      fprint(cli.fp_uchimealns, " nt) ");

      header_fprint_strip(cli.fp_uchimealns,
                          result.header,
                          attributes_to_strip(parameters));

      fprint(cli.fp_uchimealns, "\nParentA (");
      fprint_integer(cli.fp_uchimealns, db.getsequencelen(static_cast<uint64_t>(report.parent_a)), 5);
      fprint(cli.fp_uchimealns, " nt) ");
      header_fprint_strip(cli.fp_uchimealns,
                          db.header_view(static_cast<uint64_t>(report.parent_a)),
                          attributes_to_strip(parameters));

      fprint(cli.fp_uchimealns, "\nParentB (");
      fprint_integer(cli.fp_uchimealns, db.getsequencelen(static_cast<uint64_t>(report.parent_b)), 5);
      fprint(cli.fp_uchimealns, " nt) ");
      header_fprint_strip(cli.fp_uchimealns,
                          db.header_view(static_cast<uint64_t>(report.parent_b)),
                          attributes_to_strip(parameters));
      fprint(cli.fp_uchimealns, "\n\n");

      auto const width = parameters.opt_alignwidth > 0 ? parameters.opt_alignwidth : alnlen;
      auto qpos = 0;
      auto p1pos = 0;
      auto p2pos = 0;
      auto rest = alnlen;

      for (auto i = 0; i < alnlen; i += width)
        {
          /* count non-gap symbols on current line */

          auto qnt = 0;
          auto p1nt = 0;
          auto p2nt = 0;

          auto const w = std::min(rest, width);

          for (auto j = 0; j < w; ++j)
            {
              if (alignment.query[static_cast<size_t>(i + j)] != '-')
                {
                  ++qnt;
                }
              if (alignment.parents[0][static_cast<size_t>(i + j)] != '-')
                {
                  ++p1nt;
                }
              if (alignment.parents[1][static_cast<size_t>(i + j)] != '-')
                {
                  ++p2nt;
                }
            }

          auto const parent_a_row = View<char>{&alignment.parents[0][static_cast<size_t>(i)], static_cast<std::size_t>(w)};
          auto const parent_b_row = View<char>{&alignment.parents[1][static_cast<size_t>(i)], static_cast<std::size_t>(w)};
          auto const query_row = View<char>{&alignment.query[static_cast<size_t>(i)], static_cast<std::size_t>(w)};

          if (not report.parents_swapped)
            {
              print_alignment_row(cli.fp_uchimealns, 'A', p1pos + 1, parent_a_row, p1pos + p1nt);
              print_alignment_row(cli.fp_uchimealns, 'Q', qpos + 1, query_row, qpos + qnt);
              print_alignment_row(cli.fp_uchimealns, 'B', p2pos + 1, parent_b_row, p2pos + p2nt);
            }
          else
            {
              print_alignment_row(cli.fp_uchimealns, 'A', p2pos + 1, parent_b_row, p2pos + p2nt);
              print_alignment_row(cli.fp_uchimealns, 'Q', qpos + 1, query_row, qpos + qnt);
              print_alignment_row(cli.fp_uchimealns, 'B', p1pos + 1, parent_a_row, p1pos + p1nt);
            }

          print_annotation_row(cli.fp_uchimealns, "Diffs   ",
                  View<char>{&alignment.diffs[static_cast<size_t>(i)], static_cast<std::size_t>(w)});
          print_annotation_row(cli.fp_uchimealns, "Votes   ",
                  View<char>{&alignment.votes[static_cast<size_t>(i)], static_cast<std::size_t>(w)});
          print_annotation_row(cli.fp_uchimealns, "Model   ",
                  View<char>{&alignment.model[static_cast<size_t>(i)], static_cast<std::size_t>(w)});
          fprint(cli.fp_uchimealns, '\n');

          qpos += qnt;
          p1pos += p1nt;
          p2pos += p2nt;
          rest -= width;
        }

      fprint(cli.fp_uchimealns, "Ids.  QA ");
      std::fprintf(cli.fp_uchimealns, "%.1f", report.id_query_a);
      fprint(cli.fp_uchimealns, "%, QB ");
      std::fprintf(cli.fp_uchimealns, "%.1f", report.id_query_b);
      fprint(cli.fp_uchimealns, "%, AB ");
      std::fprintf(cli.fp_uchimealns, "%.1f", report.id_a_b);
      fprint(cli.fp_uchimealns, "%, QModel ");
      std::fprintf(cli.fp_uchimealns, "%.1f", report.id_query_model);
      fprint(cli.fp_uchimealns, "%, Div. ");
      std::fprintf(cli.fp_uchimealns, "%+.1f", report.divergence_percent);
      fprint(cli.fp_uchimealns, "%\n");

      fprint(cli.fp_uchimealns, "Diffs Left ");
      fprint_integer(cli.fp_uchimealns, sumL);
      fprint(cli.fp_uchimealns, ": N ");
      fprint_integer(cli.fp_uchimealns, report.left_no);
      fprint(cli.fp_uchimealns, ", A ");
      fprint_integer(cli.fp_uchimealns, report.left_abstain);
      fprint(cli.fp_uchimealns, ", Y ");
      fprint_integer(cli.fp_uchimealns, report.left_yes);
      fprint(cli.fp_uchimealns, " (");
      std::fprintf(cli.fp_uchimealns, "%.1f", 100.0 * report.left_yes / sumL);
      fprint(cli.fp_uchimealns, "%); Right ");
      fprint_integer(cli.fp_uchimealns, sumR);
      fprint(cli.fp_uchimealns, ": N ");
      fprint_integer(cli.fp_uchimealns, report.right_no);
      fprint(cli.fp_uchimealns, ", A ");
      fprint_integer(cli.fp_uchimealns, report.right_abstain);
      fprint(cli.fp_uchimealns, ", Y ");
      fprint_integer(cli.fp_uchimealns, report.right_yes);
      fprint(cli.fp_uchimealns, " (");
      std::fprintf(cli.fp_uchimealns, "%.1f", 100.0 * report.right_yes / sumR);
      fprint(cli.fp_uchimealns, "%), Score ");
      std::fprintf(cli.fp_uchimealns, "%.4f", result.best_h);
      fprint(cli.fp_uchimealns, '\n');
    }

  if (parameters.opt_uchimeout != nullptr)
    {
      std::fprintf(cli.fp_uchimeout, "%.4f", result.best_h);
      fprint(cli.fp_uchimeout, '\t');

      header_fprint_strip(cli.fp_uchimeout,
                          result.header,
                          attributes_to_strip(parameters));
      fprint(cli.fp_uchimeout, '\t');
      header_fprint_strip(cli.fp_uchimeout,
                          db.header_view(static_cast<uint64_t>(report.parent_a)),
                          attributes_to_strip(parameters));
      fprint(cli.fp_uchimeout, '\t');
      header_fprint_strip(cli.fp_uchimeout,
                          db.header_view(static_cast<uint64_t>(report.parent_b)),
                          attributes_to_strip(parameters));
      fprint(cli.fp_uchimeout, '\t');

      if (parameters.opt_uchimeout5 == 0)
        {
          if (report.id_query_a >= report.id_query_b)
            {
              header_fprint_strip(cli.fp_uchimeout,
                                  db.header_view(static_cast<uint64_t>(report.parent_a)),
                                  attributes_to_strip(parameters));
            }
          else
            {
              header_fprint_strip(cli.fp_uchimeout,
                                  db.header_view(static_cast<uint64_t>(report.parent_b)),
                                  attributes_to_strip(parameters));
            }
          fprint(cli.fp_uchimeout, '\t');
        }

      std::fprintf(cli.fp_uchimeout, "%.1f", report.id_query_model);
      fprint(cli.fp_uchimeout, '\t');
      std::fprintf(cli.fp_uchimeout, "%.1f", report.id_query_a);
      fprint(cli.fp_uchimeout, '\t');
      std::fprintf(cli.fp_uchimeout, "%.1f", report.id_query_b);
      fprint(cli.fp_uchimeout, '\t');
      std::fprintf(cli.fp_uchimeout, "%.1f", report.id_a_b);
      fprint(cli.fp_uchimeout, '\t');
      std::fprintf(cli.fp_uchimeout, "%.1f", report.id_query_top);
      fprint(cli.fp_uchimeout, '\t');
      fprint_integer(cli.fp_uchimeout, report.left_yes);
      fprint(cli.fp_uchimeout, '\t');
      fprint_integer(cli.fp_uchimeout, report.left_no);
      fprint(cli.fp_uchimeout, '\t');
      fprint_integer(cli.fp_uchimeout, report.left_abstain);
      fprint(cli.fp_uchimeout, '\t');
      fprint_integer(cli.fp_uchimeout, report.right_yes);
      fprint(cli.fp_uchimeout, '\t');
      fprint_integer(cli.fp_uchimeout, report.right_no);
      fprint(cli.fp_uchimeout, '\t');
      fprint_integer(cli.fp_uchimeout, report.right_abstain);
      fprint(cli.fp_uchimeout, '\t');
      std::fprintf(cli.fp_uchimeout, "%.1f", report.divergence);
      fprint(cli.fp_uchimeout, '\t');
      fprint(cli.fp_uchimeout, status == Status::chimeric ? 'Y' : (status == Status::low_score ? 'N' : '?'));
      fprint(cli.fp_uchimeout, '\n');
    }
}
}  // anonymous namespace


static auto query_init(struct searchinfo_s * search_info, int const tophits,
                       struct Database const & db,
                       struct Parameters const & parameters,
                       struct Dbindex const & dbindex) -> void
{
  search_info->parameters = &parameters;  /* searchcore reads config through the si (E1) */
  search_info->dbindex = &dbindex;  /* searchcore reads the k-mer index through the si */
  search_info->db = &db;  /* searchcore reads the sequences through the si */
  search_info->hits_v.resize(static_cast<size_t>(tophits));
  /* no kmers_v: the part borrows ci->kmer_counters while it is searched */
  search_info->hit_count = 0;
  /* search_info->uh (a Uniquer value member) is ready to use as default-constructed */
  // scoring.n_mismatch is always false: no chimera command accepts --n_mismatch
  search_info->s.reset(search16_init(scoring_from_options(parameters)));
  search_info->m = Minheap(tophits);
}


namespace {
auto query_exit(struct searchinfo_s & search_info) -> void
{
  /* The handles are also freed by ~searchinfo_s if an exception unwinds
     before this runs. */
  search_info.s.reset();
  search_info.uh = Uniquer();
  search_info.m = Minheap();

  search_info.qsequence = Span<char>{};
  /* the hit and kmer-count buffers (hits_v/kmers_v) free their own storage */
}


auto partition_query(struct chimera_info_s * chimera_info) -> void
{
  auto rest = chimera_info->query_len;
  auto * cursor = chimera_info->query_seq.data();
  for (auto i = 0; i < chimera_info->parts; ++i)
    {
      auto const length =
        (rest + (chimera_info->parts - i - 1)) / (chimera_info->parts - i);

      auto & search_info = chimera_info->si[static_cast<size_t>(i)];

      search_info.query_no = chimera_info->query_no;
      search_info.strand = 0;
      search_info.qsize = chimera_info->query_size;
      search_info.query_head = chimera_info->query_head;
      search_info.full_qsequence =
        make_span(chimera_info->query_seq).first(static_cast<std::size_t>(chimera_info->query_len));
      assert(static_cast<std::size_t>(length) <= search_info.qsequence_v.size());
      std::copy(cursor, std::next(cursor, length), search_info.qsequence_v.begin());
      search_info.qsequence = make_span(search_info.qsequence_v).first(static_cast<std::size_t>(length));

      rest -= length;
      cursor = std::next(cursor, length);
    }
}


auto chimera_thread_init(struct chimera_info_s * ci, int const tophits,
                         struct Parameters const & parameters,
                         struct Dbindex const & dbindex,
                         struct Database const & db,
                         ChimeraMode const mode) -> void
{
  ci->parameters = &parameters;  /* detection core reads config through ci (E1) */
  ci->mode = mode;  /* detection core reads the command variant through ci */
  ci->db = &db;  /* detection core reads the sequences through ci */

  /* the per-part searchinfo_s are built on demand (see parts_ready) */
  ci->dbindex = &dbindex;
  ci->tophits = tophits;
  ci->parts_ready = 0;

  /* headroom past the logical end for the SIMD counter stores */
  static constexpr auto overflow_padding = 16U;  // 16 * sizeof(short) = 32 bytes
  ci->kmer_counters.reserve(db.getsequencecount() + overflow_padding);
  ci->kmer_counters.resize(db.getsequencecount());

  // scoring.n_mismatch is always false: no chimera command accepts --n_mismatch
  ci->s.reset(search16_init(scoring_from_options(parameters)));
}


auto chimera_thread_exit(struct chimera_info_s * ci) -> void
{
  ci->s.reset();

  for (auto & a_search_info : ci->si) {
    query_exit(a_search_info);
  }
}
}  // anonymous namespace


/* Lends a thread's k-mer counters to one part search, and takes them back
   when the search is over, also when fatal() throws (library sessions),
   so that the next query finds them where it expects them. A swap: it
   keeps the reserved headroom and costs nothing. */
class CounterLoan
{
public:
  CounterLoan(std::vector<count_t> & lender, struct searchinfo_s & borrower) noexcept
    : lender_(lender), borrower_(borrower)
  {
    borrower_.kmers_v.swap(lender_);
  }
  ~CounterLoan()
  {
    borrower_.kmers_v.swap(lender_);
  }
  CounterLoan(CounterLoan const &) = delete;
  CounterLoan(CounterLoan &&) = delete;
  auto operator=(CounterLoan const &) -> CounterLoan & = delete;
  auto operator=(CounterLoan &&) -> CounterLoan & = delete;

private:
  std::vector<count_t> & lender_;
  struct searchinfo_s & borrower_;
};


/* Process a single query that has already been loaded into ci.
   Shared by chimera_thread_core (CLI) and chimera_detect_single (API).
   ci->query_seq, query_head, query_len, query_size must
   be populated. allhits_list must be pre-allocated to maxcandidates.
   lma is the per-thread linear memory aligner (fallback for SIMD overflow). */
static auto chimera_process_query(struct chimera_info_s * ci,
                                  std::vector<struct hit> & allhits_list,
                                  LinearMemoryAligner & lma,
                                  struct Database const & db) -> Status
{
  struct Parameters const & parameters = *ci->parameters;

  /* build the per-part search state this query needs, once */
  assert(ci->dbindex != nullptr);
  assert(ci->parts <= maxparts);
  for (int i = ci->parts_ready; i < ci->parts; ++i)
    {
      query_init(&ci->si[static_cast<size_t>(i)], ci->tophits, db, parameters, *ci->dbindex);
    }
  ci->parts_ready = std::max(ci->parts_ready, ci->parts);

  /* partition query */
  partition_query(ci);

  /* perform searches and collect candidate parents */
  ci->cand_count = 0;
  ci->best_h = 0.0;
  auto allhits_count = 0;

  if (ci->query_len >= ci->parts)
    {
      std::vector<struct hit> hits;
      for (auto i = 0; i < ci->parts; ++i)
        {
          {
            CounterLoan const loan(ci->kmer_counters, ci->si[static_cast<size_t>(i)]);
            search_onequery(&ci->si[static_cast<size_t>(i)], parameters.opt_qmask);
          }
          search_joinhits(&ci->si[static_cast<size_t>(i)], nullptr, hits);
          for (auto & hit : hits) {
            if (hit.accepted and allhits_count < maxcandidates)
              {
                allhits_list[static_cast<size_t>(allhits_count)] = hit;
                ++allhits_count;
              }
            else
              {
                // Unallocate alignments for weak hits
                hit.nwalignment.clear();  // std::string; frees the weak-hit alignment
              }
          }
          hits.clear();
        }
    }

  for (auto i = 0; i < allhits_count; ++i)
    {
      auto const target = static_cast<unsigned int>(allhits_list[static_cast<size_t>(i)].target);

      /* skip duplicates */
      auto k = 0;
      for (k = 0; k < ci->cand_count; ++k)
        {
          if (ci->cand_list[static_cast<size_t>(k)] == target)
            {
              break;
            }
        }

      if (k == ci->cand_count)
        {
          ci->cand_list[static_cast<size_t>(ci->cand_count)] = target;
          ++ci->cand_count;
        }

      /* deallocate cigar */
      allhits_list[static_cast<size_t>(i)].nwalignment.clear();  // std::string; frees after use
    }


  /* align full query to each candidate */

  search16_qprep(*ci->s, ci->query());

  /* the candidates found above, not the whole maxcandidates buffers */
  auto const candidates = static_cast<std::size_t>(ci->cand_count);
  search16(*ci->s,
           make_view(ci->cand_list).first(candidates),
           make_span(ci->snwscore).first(candidates),
           make_span(ci->snwalignmentlength).first(candidates),
           make_span(ci->snwmatches).first(candidates),
           make_span(ci->snwmismatches).first(candidates),
           make_span(ci->snwgaps).first(candidates),
           make_span(ci->nwcigar).first(candidates),
           db);

  for (auto i = 0; i < ci->cand_count; ++i)
    {
      int64_t const target = ci->cand_list[static_cast<size_t>(i)];
      int64_t const nwscore = ci->snwscore[static_cast<size_t>(i)];

      if (nwscore == std::numeric_limits<short>::max())
        {
          /* In case the SIMD aligner cannot align,
             perform a new alignment with the
             linear memory aligner */

          auto const tseq = db.sequence_view(static_cast<uint64_t>(target));
          auto const qseq = ci->query();

          std::string nwcigar = lma.align(qseq, tseq);
          auto const stats = lma.alignstats(nwcigar.c_str(), qseq, tseq);

          ci->nwcigar[static_cast<size_t>(i)] = std::move(nwcigar);
          ci->nwscore[static_cast<size_t>(i)] = stats.score;
          ci->nwalignmentlength[static_cast<size_t>(i)] = stats.alignmentlength;
          ci->nwmatches[static_cast<size_t>(i)] = stats.matches;
          ci->nwmismatches[static_cast<size_t>(i)] = stats.mismatches;
          ci->nwgaps[static_cast<size_t>(i)] = stats.gaps;
        }
      else
        {
          ci->nwscore[static_cast<size_t>(i)] = ci->snwscore[static_cast<size_t>(i)];
          ci->nwalignmentlength[static_cast<size_t>(i)] = ci->snwalignmentlength[static_cast<size_t>(i)];
          ci->nwmatches[static_cast<size_t>(i)] = ci->snwmatches[static_cast<size_t>(i)];
          ci->nwmismatches[static_cast<size_t>(i)] = ci->snwmismatches[static_cast<size_t>(i)];
          ci->nwgaps[static_cast<size_t>(i)] = ci->snwgaps[static_cast<size_t>(i)];
        }
    }


  /* find the best pair of parents, then compute score for them */

  if (ci->mode == ChimeraMode::chimeras_denovo)
    {
      /* long high-quality reads */
      if (find_best_parents_long(ci) != 0)
        {
          return eval_parents_long(ci, db);
        }
      return Status::no_parents;
    }
  if (find_best_parents(ci) != 0)
    {
      return eval_parents(ci, db);
    }
  return Status::no_parents;
}


/* Take from the worker's chimera_info_s what the output step needs, once
   detection of the current query is over (see chimera_query_result_s). */
static auto collect_query_result(struct chimera_info_s const & ci,
                                 Status const status,
                                 uint64_t const query_position,
                                 struct chimera_query_result_s & result) -> void
{
  result.status = status;
  result.seqno = static_cast<unsigned int>(ci.query_no);
  result.query_position = query_position;
  result.abundance = ci.query_size;
  if (chimera_is_denovo(ci.mode))
    {
      /* the same bytes ci copied, in storage that outlives ci's next query */
      result.header = ci.db->header_view(result.seqno);
      result.sequence = ci.db->sequence_view(result.seqno);
    }
  else
    {
      result.header = ci.query_head;
      result.sequence = ci.query();
    }
  result.best_h = ci.best_h;
  result.report = ci.report;

  auto & alignment = result.alignment;
  alignment.query.clear();
  alignment.parents.clear();
  alignment.parent_seqnos.clear();
  alignment.diffs.clear();
  alignment.votes.clear();
  alignment.model.clear();

  auto const * const alignment_output =
    (ci.mode == ChimeraMode::chimeras_denovo) ?
    ci.parameters->opt_alnout : ci.parameters->opt_uchimealns;
  if ((status != Status::chimeric) or (alignment_output == nullptr))
    {
      return;
    }

  /* the rows are high-water buffers: keep the alignment's length of each */
  auto const alnlen = static_cast<std::size_t>(find_total_alignment_length(&ci));
  auto const copy_row = [alnlen](std::vector<char> const & row,
                                 std::vector<char> & destination) -> void {
    assert(row.size() >= alnlen);
    destination.assign(row.cbegin(),
                       std::next(row.cbegin(), static_cast<std::ptrdiff_t>(alnlen)));
  };
  copy_row(ci.qaln, alignment.query);
  auto const parent_count = static_cast<std::size_t>(ci.parents_found);
  alignment.parents.resize(parent_count);
  alignment.parent_seqnos.resize(parent_count);
  for (std::size_t nth = 0; nth < parent_count; ++nth)
    {
      copy_row(ci.paln[nth], alignment.parents[nth]);
      alignment.parent_seqnos[nth] = static_cast<int>(
        ci.cand_list[static_cast<std::size_t>(ci.best_parents[nth])]);
    }
  copy_row(ci.diffs, alignment.diffs);
  copy_row(ci.model, alignment.model);
  if (ci.mode != ChimeraMode::chimeras_denovo)
    {
      copy_row(ci.votes, alignment.votes);
    }
}


/* Write, count and index one query's result, under the output lock. Reads
   the result and the run state only: never a worker's chimera_info_s. */
static auto output_query_result(struct chimera_cli_state_s & state,
                                struct chimera_query_result_s const & result,
                                struct Database const & db) -> void
{
  auto const status = result.status;

  /* the detection core writes nothing: its alignment and tabbed records
     are written here, under the same lock as the rest of the query's
     output */
  if (status >= Status::low_score)
    {
      if (state.mode == ChimeraMode::chimeras_denovo)
        {
          print_report_long(state, result, db);
        }
      else
        {
          print_report(state, result, db);
        }
    }

  ++state.total_count;
  state.total_abundance += result.abundance;

  /* the three FASTA outputs below annotate the same query the same way and
     differ only in the per-status counter that supplies the ordinal */
  auto const query_annotations = [&](int64_t const ordinal) -> OutputAnnotations {
    OutputAnnotations annotations {static_cast<uint64_t>(result.abundance), ordinal};
    annotations.score_name = state.parameters.opt_fasta_score ?
      ( (state.mode == ChimeraMode::uchime_ref) ?
        "uchime_ref" : "uchime_denovo" ) : nullptr;
    annotations.score = result.best_h;
    return annotations;
  };

  if (status == Status::chimeric)
    {
      ++state.chimera_count;
      state.chimera_abundance += result.abundance;

      if (state.parameters.opt_chimeras != nullptr)
        {
          fasta_print_general(state.fp_chimeras,
                              result.sequence,
                              result.header,
                              query_annotations(state.chimera_count),
                              state.parameters);

        }
    }

  if (status == Status::suspicious)
    {
      ++state.borderline_count;
      state.borderline_abundance += result.abundance;

      if (state.parameters.opt_borderline != nullptr)
        {
          fasta_print_general(state.fp_borderline,
                              result.sequence,
                              result.header,
                              query_annotations(state.borderline_count),
                              state.parameters);

        }
    }

  if (status < Status::suspicious)
    {
      ++state.nonchimera_count;
      state.nonchimera_abundance += result.abundance;

      /* output no parents, no chimeras */
      if ((status < Status::low_score) and (state.parameters.opt_uchimeout != nullptr))
        {
          std::fprintf(state.fp_uchimeout, "%.4f", result.best_h);
          fprint(state.fp_uchimeout, '\t');

          header_fprint_strip(state.fp_uchimeout,
                              result.header,
                              attributes_to_strip(state.parameters));

          if (state.parameters.opt_uchimeout5 != 0)
            {
              fprint(state.fp_uchimeout, "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
            }
          else
            {
              fprint(state.fp_uchimeout, "\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
            }
        }

      if (state.parameters.opt_nonchimeras != nullptr)
        {
          fasta_print_general(state.fp_nonchimeras,
                              result.sequence,
                              result.header,
                              query_annotations(state.nonchimera_count),
                              state.parameters);
        }
    }

  if (status < Status::suspicious)
    {
      /* uchime_denovo: add non-chimeras to db */
      if (chimera_is_denovo(state.mode))
        {
          state.dbindex.add_sequence(result.seqno, state.parameters.opt_qmask, db);
        }
    }

  if (state.mode == ChimeraMode::uchime_ref)
    {
      state.progress = result.query_position;
    }
  else
    {
      state.progress += db.getsequencelen(result.seqno);
    }

  state.progress_bar->update(state.progress);
}


namespace {
/* Make the result's header and sequence its own: for --uchime_ref they are
   views into the worker's buffers, which its next claim overwrites. */
static auto keep_query_text(struct chimera_query_result_s & result) -> void
{
  result.header_copy.assign(result.header.cbegin(), result.header.cend());
  result.sequence_copy.assign(result.sequence.cbegin(), result.sequence.cend());
  result.header = make_view(result.header_copy);
  result.sequence = make_view(result.sequence_copy);
}


/* Write the result of the query claimed at `rank`, in claim order (see
   ordered_output_s). Called with the output lock held. */
static auto output_in_order(struct chimera_cli_state_s & state,
                            struct ordered_output_s & ordered,
                            unsigned int const rank,
                            struct chimera_query_result_s & result,
                            struct Database const & db) -> void
{
  if (rank != ordered.next_rank)
    {
      keep_query_text(result);
      ordered.waiting.emplace(rank, std::move(result));
      return;
    }
  output_query_result(state, result, db);
  ++ordered.next_rank;
  auto next = ordered.waiting.find(ordered.next_rank);
  while (next != ordered.waiting.end())
    {
      output_query_result(state, next->second, db);
      ordered.waiting.erase(next);
      ++ordered.next_rank;
      next = ordered.waiting.find(ordered.next_rank);
    }
}


/* Copy database sequence seqno into chimera_info as the query to process
   (denovo) */
auto load_denovo_query(struct chimera_info_s * chimera_info,
                       struct Database const & database,
                       unsigned int const seqno) -> void
{
  auto const query_record = database.record(seqno);
  chimera_info->query_no = static_cast<int>(seqno);
  chimera_info->query_len = static_cast<int>(query_record.sequence.size());
  chimera_info->query_size = static_cast<int64_t>(database.getabundance(seqno));

  /* if necessary expand memory for arrays based on query length */
  realloc_arrays(chimera_info, database, query_record.header.size());

  chimera_info->query_head = copy_into_scratch(query_record.header, chimera_info->query_head_v);
  copy_into_scratch(query_record.sequence, chimera_info->query_seq);
}
}  // anonymous namespace


static auto chimera_thread_core(struct chimera_cli_state_s & state,
                         struct chimera_info_s * ci,
                         std::mutex & mutex_input,
                         struct ordered_output_s & ordered,
                         struct Database const & db) -> uint64_t
{
  /* tophits sizes the per-part minheaps; it is maxaccepts + maxrejects from
     the chimera-detection copy chimera() built before spawning the pool.
     Computed here rather than read from a shared file-static (E4). */
  int const tophits = static_cast<int>(state.detection_parameters.opt_maxaccepts +
                                       state.detection_parameters.opt_maxrejects);
  chimera_thread_init(ci, tophits, state.detection_parameters, state.dbindex, db, state.mode);

  std::vector<struct hit> allhits_list(maxcandidates);

  struct Scoring const scoring = scoring_from_options(state.parameters);

  LinearMemoryAligner lma(scoring);

  uint64_t query_position = 0;

  /* this worker's current query, handed from detection to output */
  struct chimera_query_result_s result;
  unsigned int claimed_rank = 0;  /* its rank in claim order */

  auto const has_work_to_claim = [&]() -> bool {
    /* get next sequence */

    /* Query-file progress position. Read here, under input_lock, into a
       worker-local: fasta_get_position() reads the shared query handle,
       which another worker advances in fasta_next() under the same lock.
       Reading it later under mutex_output (as before) raced that writer. */
    query_position = 0;

    if (state.mode == ChimeraMode::uchime_ref)
      {
        if (state.query_fasta_h->next(header_truncation(state.parameters.opt_notrunclabels),
                       Mapping::none))
          {
            auto const query_record = state.query_fasta_h->record();
            ci->query_len = static_cast<int>(query_record.sequence.size());
            ci->query_no = static_cast<int>(state.query_fasta_h->get_seqno());
            ci->query_size = state.query_fasta_h->get_abundance();

            /* if necessary expand memory for arrays based on query length */
            realloc_arrays(ci, db, query_record.header.size());

            /* copy the data locally (query seq, head) */
            ci->query_head = copy_into_scratch(query_record.header, ci->query_head_v);
            copy_into_scratch(query_record.sequence, ci->query_seq);
            query_position = state.query_fasta_h->get_position();
          }
        else
          {
            return false;
          }
      }
    else
      {
        if (state.seqno < db.getsequencecount())
          {
            load_denovo_query(ci, db, state.seqno);
          }
        else
          {
            return false;
          }
      }

    /* One more query claimed, in either mode: denovo takes its next query
       number from state.seqno, and the --log summary counts the queries
       with it. Claim and advance in the same critical section
       (mutex_input), so two workers can never claim the same denovo query.
       The output step reads the query's own number from its result. This
       loop runs denovo detection with one worker only: the index has to
       grow in query order, which claiming alone does not ensure, so more
       threads go to chimera_denovo_batches instead (see chimera()). */
    claimed_rank = state.seqno;
    ++state.seqno;
    return true;
  };

  auto const process_query = [&]() -> void {
    auto const status = chimera_process_query(ci, allhits_list, lma, db);
    collect_query_result(*ci, status, query_position, result);

    /* output results */

    std::lock_guard<std::mutex> const output_lock(state.mutex_output);

    output_in_order(state, ordered, claimed_rank, result, db);
  };

  run_worker_loop(mutex_input, has_work_to_claim, process_query);

  chimera_thread_exit(ci);


  return 0;
}


static auto chimera_threads_run(struct chimera_cli_state_s & state) -> void
{
  /* mutex_input serializes input reading among the CLI workers; it is
     owned here rather than at file scope (the API path does not use it). */
  std::mutex mutex_input;

  /* results are written in query order, whatever the thread count */
  struct ordered_output_s ordered;

  /* run the worker pool; each worker processes queries until the input
     is exhausted. chimera_thread_core returns a value that the previous
     pthread_join already discarded, so it is ignored here too. */
  ThreadRunner threadrunner(static_cast<std::size_t>(state.detection_parameters.opt_threads),
                            [&state, &mutex_input, &ordered](uint64_t const nth_thread) -> void {
                              chimera_thread_core(state, &state.cia[nth_thread], mutex_input,
                                                  ordered, state.db);
                            });
  threadrunner.run();
  /* every claimed query was processed, so every rank up to the last was
     written */
  assert(ordered.waiting.empty());
}


namespace {
/* Keep what the batch driver needs to check a query's result later: the
   heap threshold and the unique k-mers of each part search. Called
   right after chimera_process_query, before ci's searches are reused. */
auto collect_part_searches(struct chimera_info_s const & chimera_info,
                           struct chimera_query_result_s & result) -> void
{
  /* chimera_process_query searches the parts only when the query is at
     least as long as their number */
  auto const searched = (chimera_info.query_len >= chimera_info.parts) ?
    static_cast<std::size_t>(chimera_info.parts) : 0U;
  result.part_searches.resize(searched);
  for (std::size_t nth = 0; nth < searched; ++nth)
    {
      auto const & search_info = chimera_info.si[nth];
      auto & part = result.part_searches[nth];
      part.threshold = search_info.topscores_threshold;
      part.kmers.assign(search_info.kmersample.cbegin(), search_info.kmersample.cend());
    }
}


/* No k-mer has this code: codes stay below 4^15 (--wordlength 15) */
constexpr auto empty_kmer_slot = std::numeric_limits<unsigned int>::max();


/* A query found not chimeric earlier in the current batch: the index the
   batch was detected against did not hold it yet. Its unique k-mers, as
   Dbindex::add_sequence counts them, are kept in a small open-addressing
   hash set (linear probing, at most half full), so that counting the
   k-mers a part search shares with it costs one probe per k-mer of the
   part. */
struct batch_target_s {
  unsigned int seqno = 0;
  unsigned int length = 0;
  std::size_t kmer_count = 0;
  std::vector<unsigned int> slots;  /* a power of two of them */

  auto set_kmers(View<unsigned int> const kmers) -> void
  {
    kmer_count = kmers.size();
    std::size_t size = minimum_slots;
    while (size < 2 * kmer_count)
      {
        size *= 2;
      }
    slots.assign(size, empty_kmer_slot);
    for (auto const kmer : kmers)
      {
        assert(kmer != empty_kmer_slot);
        slots[find_slot(kmer)] = kmer;  /* the k-mers are distinct */
      }
  }

  auto contains(unsigned int const kmer) const noexcept -> bool
  {
    return slots[find_slot(kmer)] == kmer;
  }

private:
  static constexpr std::size_t minimum_slots = 16;
  static constexpr unsigned int hash_multiplier = 2654435761U;  /* Knuth's */
  static constexpr unsigned int hash_fold = 16U;  /* mix the high half into the low */

  /* the slot holding kmer, or the empty slot where it would go */
  auto find_slot(unsigned int const kmer) const noexcept -> std::size_t
  {
    auto const mask = slots.size() - 1;
    auto mixed = kmer * hash_multiplier;
    mixed ^= mixed >> hash_fold;
    auto slot = static_cast<std::size_t>(mixed) & mask;
    while ((slots[slot] != empty_kmer_slot) and (slots[slot] != kmer))
      {
        slot = (slot + 1) & mask;
      }
    return slot;
  }
};


/* The number of the part's k-mers the target holds too, when it is at
   least `needed`; otherwise some smaller number, as counting stops once
   the k-mers left to probe could no longer make up the difference */
auto count_shared_kmers(std::vector<unsigned int> const & part_kmers,
                        struct batch_target_s const & target,
                        std::size_t const needed) noexcept -> std::size_t
{
  std::size_t shared = 0;
  auto remaining = part_kmers.size();
  for (auto const kmer : part_kmers)
    {
      if (shared + remaining < needed)
        {
          break;
        }
      --remaining;
      if (target.contains(kmer))
        {
          ++shared;
        }
    }
  return shared;
}


/* Would target have entered the heap of this part search, had it been
   indexed before the search ran? The two tests search_topscores() applies
   (see TopscoresThreshold): its k-mer count reaches its threshold, then
   the heap takes it. The target has a higher sequence number than every
   sequence the search saw, which topscore_ranks_below() accounts for. */
auto enters_part_heap(struct chimera_part_search_s const & part,
                      struct batch_target_s const & target,
                      unsigned int const minwordmatches) noexcept -> bool
{
  auto const & threshold = part.threshold;
  if (threshold.capacity == 0)
    {
      return false;
    }
  auto needed = static_cast<std::size_t>(threshold.minmatches);
  auto const target_kmers = target.kmer_count;
  if ((target_kmers != 0) and (target_kmers < minwordmatches))
    {
      /* a low-k-mer target: offered at its own, lower threshold */
      needed = std::min(needed, target_kmers);
    }
  auto const heap_full = (threshold.filled >= threshold.capacity);
  if (heap_full)
    {
      /* fewer k-mer hits than the weakest element never displace it */
      needed = std::max(needed, static_cast<std::size_t>(threshold.weakest.count));
    }
  /* the k-mer counters saturate at INT16_MAX (accumulate_slice_counts), so
     no count ever reaches more */
  constexpr auto counter_ceiling = static_cast<std::size_t>(INT16_MAX);
  if (needed > counter_ceiling)
    {
      return false;
    }
  auto const shared = std::min(count_shared_kmers(part.kmers, target, needed),
                               counter_ceiling);
  if (shared < needed)
    {
      return false;
    }
  if (not heap_full)
    {
      return true;
    }
  elem_t const candidate {static_cast<unsigned int>(shared), target.seqno, target.length};
  return topscore_ranks_below(threshold.weakest, candidate);
}


/* True when one of the targets would have changed one of the query's part
   searches, and so possibly its result. Nothing after search_topscores()
   depends on the index, so a query whose part heaps no target enters has
   exactly the result it would have had with the targets indexed. */
auto result_is_stale(struct chimera_query_result_s const & result,
                     std::vector<struct batch_target_s> const & targets,
                     unsigned int const minwordmatches) noexcept -> bool
{
  return std::any_of(result.part_searches.cbegin(), result.part_searches.cend(),
                     [&](struct chimera_part_search_s const & part) -> bool {
                       return std::any_of(targets.cbegin(), targets.cend(),
                                          [&](struct batch_target_s const & target) -> bool {
                                            return enters_part_heap(part, target, minwordmatches);
                                          });
                     });
}


/* What one thread needs to run the detection core, as chimera_thread_core
   keeps it: its chimera_info_s, the hit buffer and the linear memory
   aligner */
struct denovo_worker_s {
  struct chimera_info_s ci;
  std::vector<struct hit> allhits_list = std::vector<struct hit>(maxcandidates);
  std::unique_ptr<LinearMemoryAligner> lma;
};


/* Queries per batch, per thread. Larger batches amortise the two
   synchronisations per batch and even out the threads' loads, but
   more queries are then detected without the ones before them, and have
   to be detected again when one of those would have reached them. The
   best value depends on the input: 1 to 2 on 61k 16S V4 amplicons (35 %
   non-singletons), 2 to 4 on 219k 18S V9 ones (29 %), measured with 8, 16
   and 24 threads (2026-09-24). 2 stays within 20 % of the best on both. */
constexpr auto batch_size_per_thread = 2U;


/* Denovo detection with more than one thread, as batch speculation with
   in-order validation. Each batch of queries is detected in parallel
   against the index as it stood before the batch. The main thread then
   commits the batch in query order: a query that a sequence found not
   chimeric earlier in the same batch would have reached is detected again,
   serially, against the index that now holds that sequence; then the query
   is written, counted and, if not chimeric, indexed, exactly as by the
   single-threaded loop. Every query is thus either shown to be unaffected
   by the queries it was detected without, or detected with them, and the
   output is identical to --threads 1, in the same order. */
auto chimera_denovo_batches(struct chimera_cli_state_s & state) -> void
{
  auto const & database = state.db;
  auto const thread_count = static_cast<std::size_t>(state.detection_parameters.opt_threads);
  assert(thread_count > 1);
  auto const batch_size = static_cast<unsigned int>(batch_size_per_thread * thread_count);
  int const tophits = static_cast<int>(state.detection_parameters.opt_maxaccepts +
                                       state.detection_parameters.opt_maxrejects);
  struct Scoring const scoring = scoring_from_options(state.parameters);

  /* one per worker thread, and one for the main thread's recomputations */
  std::vector<struct denovo_worker_s> workers(thread_count + 1);
  for (auto & worker : workers)
    {
      chimera_thread_init(&worker.ci, tophits, state.detection_parameters,
                          state.dbindex, database, state.mode);
      worker.lma = make_unique<LinearMemoryAligner>(scoring);
    }
  auto & committer = workers.back();

  std::vector<struct chimera_query_result_s> results(batch_size);
  std::vector<struct batch_target_s> batch_targets;
  Uniquer target_kmers;  /* the committer's own: the index's is not to be shared */

  /* the batch being detected, set by the main thread between two runs */
  unsigned int batch_start = 0;
  unsigned int batch_end = 0;
  unsigned int next_query = 0;
  std::mutex mutex_input;

  auto const detect = [&database](struct denovo_worker_s & worker,
                            unsigned int const seqno,
                            struct chimera_query_result_s & result) -> void {
    load_denovo_query(&worker.ci, database, seqno);
    auto const status = chimera_process_query(&worker.ci, worker.allhits_list,
                                              *worker.lma, database);
    collect_query_result(worker.ci, status, 0, result);
  };

  ThreadRunner runner(thread_count, [&](uint64_t const nth_thread) -> void {
    auto & worker = workers[nth_thread];
    unsigned int seqno = 0;
    run_worker_loop(mutex_input,
                    [&]() -> bool {
                      if (next_query >= batch_end)
                        {
                          return false;
                        }
                      seqno = next_query;
                      ++next_query;
                      return true;
                    },
                    [&]() -> void {
                      auto & result = results[seqno - batch_start];
                      detect(worker, seqno, result);
                      collect_part_searches(worker.ci, result);
                    });
  });

  auto const query_count = static_cast<unsigned int>(database.getsequencecount());
  auto const minwordmatches = state.dbindex.minwordmatches;
  for (batch_start = 0; batch_start < query_count; batch_start = batch_end)
    {
      batch_end = batch_start + std::min(batch_size, query_count - batch_start);
      next_query = batch_start;
      runner.run();  // C++20 refactoring: a std::barrier could replace the per-batch run()

      batch_targets.clear();
      for (auto seqno = batch_start; seqno < batch_end; ++seqno)
        {
          auto & result = results[seqno - batch_start];
          if (result_is_stale(result, batch_targets, minwordmatches))
            {
              detect(committer, seqno, result);
            }

          {
            std::lock_guard<std::mutex> const output_lock(state.mutex_output);
            output_query_result(state, result, database);
          }
          ++state.seqno;  /* the --log summary counts the queries with it */

          if (result.status >= Status::suspicious)
            {
              continue;
            }
          /* output_query_result indexed it: the rest of the batch was
             detected without it */
          batch_targets.emplace_back();
          auto & target = batch_targets.back();
          target.seqno = seqno;
          target.length = static_cast<unsigned int>(database.getsequencelen(seqno));
          auto const kmers = target_kmers.count(static_cast<int>(state.dbindex.wordlength),
                                                database.sequence_view(seqno),
                                                state.parameters.opt_qmask);
          target.set_kmers(kmers);
        }
    }

  for (auto & worker : workers)
    {
      chimera_thread_exit(&worker.ci);
    }
}
}  // anonymous namespace


/* Defined below (next to the library detection entry that also uses it). */
static auto chimera_detection_parameters(struct Parameters const & parameters,
                                         ChimeraMode mode) -> struct Parameters;


auto chimera(ChimeraMode const mode, struct Parameters const & parameters) -> void
{
  /* Per-invocation CLI state, owned here and threaded through the worker pool
     (E4). It also holds the report output handles/mutex, used by the report
     writers that process_query calls after detection (E6). */
  struct chimera_cli_state_s state(parameters, mode);

  OutputFileHandle chimeras_handle = open_optional_output_file(parameters.opt_chimeras, OutputOption{"--chimeras"});
  state.fp_chimeras = chimeras_handle.get();
  OutputFileHandle nonchimeras_handle = open_optional_output_file(parameters.opt_nonchimeras, OutputOption{"--nonchimeras"});
  state.fp_nonchimeras = nonchimeras_handle.get();
  OutputFileHandle borderline_handle = open_optional_output_file(parameters.opt_borderline, OutputOption{"--borderline"});
  state.fp_borderline = borderline_handle.get();

  OutputFileHandle uchimealns_handle;
  OutputFileHandle uchimeout_handle;
  if (mode == ChimeraMode::chimeras_denovo)
    {
      uchimealns_handle = open_optional_output_file(parameters.opt_alnout, OutputOption{"--alnout"});
      uchimeout_handle = open_optional_output_file(parameters.opt_tabbedout, OutputOption{"--tabbedout"});
    }
  else
    {
      uchimealns_handle = open_optional_output_file(parameters.opt_uchimealns, OutputOption{"--uchimealns"});
      uchimeout_handle = open_optional_output_file(parameters.opt_uchimeout, OutputOption{"--uchimeout"});
    }
  state.fp_uchimealns = uchimealns_handle.get();
  state.fp_uchimeout = uchimeout_handle.get();


  /* Build the detection configuration (maxaccepts/maxrejects/id/weak_id, and in
     denovo mode self/selfid/maxsizeratio) through the shared builder so the CLI
     and library detection paths stay identical, rather than mutating the shared
     opt_* config globals (E1). The private copy is threaded to the detection
     core (via si->parameters) and to chimera_threads_run (pool size). */
  state.detection_parameters = chimera_detection_parameters(parameters, mode);

  if (parameters.opt_strand)
    {
      fatal("Only --strand plus is allowed with uchime_ref.");
    }

  /* CLI-only: denovo detection is order-dependent (each query is compared
     against previously processed sequences). With more than one thread it
     runs in batches that reproduce the single-threaded result
     (chimera_denovo_batches), which keep their own per-thread state. */
  auto const runs_in_batches =
    chimera_is_denovo(mode) and (state.detection_parameters.opt_threads > 1);

  uint64_t progress_total = 0;
  state.chimera_count = 0;
  state.nonchimera_count = 0;
  state.progress = 0;
  state.seqno = 0;

  /* prepare per-thread chimera detection state */
  if (not runs_in_batches)
    {
      state.cia.resize(static_cast<size_t>(state.detection_parameters.opt_threads));
    }

  /* prepare queries / database */
  if (mode == ChimeraMode::uchime_ref)
    {
      /* check if the reference database may be an UDB file */

      auto const is_udb = udb_detect_isudb(parameters.opt_db);

      if (is_udb)
        {
          udb_read(parameters.opt_db, UdbUse::search, state.dbindex, state.db, parameters);
        }
      else
        {
          state.db.read(parameters.opt_db, 0, parameters);
          apply_masking(state.db, parameters.opt_dbmask, parameters);
          state.dbindex.prepare(parameters.opt_dbmask, state.db, parameters);
          state.dbindex.add_all_sequences(parameters.opt_dbmask, state.db, parameters);
        }

      state.query_fasta_h = fastx_open(parameters.input_filename, parameters);
      progress_total = state.query_fasta_h->get_size();

      /* The query file is parsed inside the worker threads
         (chimera_thread_core). Defer parse errors so a malformed query
         stops the pool cooperatively instead of calling fatal()/std::exit()
         from a worker while siblings are writing output (CC3); reported
         after the pool joins, below, from the main thread. */
      state.query_fasta_h->enable_deferred_errors();
    }
  else
    {

      state.db.read(parameters.input_filename, 0, parameters);

      apply_masking(state.db, parameters.opt_qmask, parameters);

      state.db.sortbyabundance(parameters);
      state.dbindex.prepare(parameters.opt_qmask, state.db, parameters);
      progress_total = state.db.getnucleotidecount();
    }

  if (parameters.fp_log != nullptr)
    {
      if ((mode == ChimeraMode::uchime_ref) or (mode == ChimeraMode::uchime_denovo))
        {
          std::fprintf(parameters.fp_log, "%8.2f", parameters.opt_minh);
          fprint(parameters.fp_log, "  minh\n");
        }
      /* the four --uchime* commands, i.e. everything but --chimeras_denovo */
      auto const is_a_uchime_command = (mode != ChimeraMode::chimeras_denovo);
      if (is_a_uchime_command)
        {
          std::fprintf(parameters.fp_log, "%8.2f", parameters.opt_xn);
          fprint(parameters.fp_log, "  xn\n");
          std::fprintf(parameters.fp_log, "%8.2f", parameters.opt_dn);
          fprint(parameters.fp_log, "  dn\n");
          std::fprintf(parameters.fp_log, "%8.2f", 1.0);
          fprint(parameters.fp_log, "  xa\n");
        }

      if ((mode == ChimeraMode::uchime_ref) or (mode == ChimeraMode::uchime_denovo))
        {
          std::fprintf(parameters.fp_log, "%8.2f", parameters.opt_mindiv);
          fprint(parameters.fp_log, "  mindiv\n");
        }

      std::fprintf(parameters.fp_log, "%8.2f", state.detection_parameters.opt_id);
      fprint(parameters.fp_log, "  id\n");

      if (is_a_uchime_command)
        {
          fprint_integer(parameters.fp_log, 2, 8);
          fprint(parameters.fp_log, "  maxp\n");
        }

      fprint(parameters.fp_log, '\n');
    }


  {
    Progress progress_bar("Detecting chimeras", progress_total, parameters);
    state.progress_bar = &progress_bar;
    if (runs_in_batches)
      {
        chimera_denovo_batches(state);
      }
    else
      {
        chimera_threads_run(state);
      }
  }

  /* all workers joined; report a deferred query parse error (CC3, uchime_ref
     only) from the main thread so it does not race a worker's output */
  if ((mode == ChimeraMode::uchime_ref) and state.query_fasta_h->get_error())
    {
      fatal(state.query_fasta_h->get_errmsg());
    }

  if (not parameters.opt_quiet)
    {
      if (state.total_count > 0)
        {
          if (mode == ChimeraMode::chimeras_denovo)
            {
              fprint(stderr, "Found ");
              fprint_integer(stderr, state.chimera_count);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * state.chimera_count / state.total_count);
              fprint(stderr, "%) chimeras and ");
              fprint_integer(stderr, state.nonchimera_count);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * state.nonchimera_count / state.total_count);
              fprint(stderr, "%) non-chimeras in ");
              fprint_integer(stderr, state.total_count);
              fprint(stderr, " unique sequences.\n");
            }
          else
            {
              fprint(stderr, "Found ");
              fprint_integer(stderr, state.chimera_count);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * state.chimera_count / state.total_count);
              fprint(stderr, "%) chimeras, ");
              fprint_integer(stderr, state.nonchimera_count);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * state.nonchimera_count / state.total_count);
              fprint(stderr, "%) non-chimeras,\nand ");
              fprint_integer(stderr, state.borderline_count);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * state.borderline_count / state.total_count);
              fprint(stderr, "%) borderline sequences in ");
              fprint_integer(stderr, state.total_count);
              fprint(stderr, " unique sequences.\n");
            }
        }
      else
        {
          if (mode == ChimeraMode::chimeras_denovo)
            {
              fprint(stderr, "Found ");
              fprint_integer(stderr, state.chimera_count);
              fprint(stderr, " chimeras and ");
              fprint_integer(stderr, state.nonchimera_count);
              fprint(stderr, " non-chimeras in ");
              fprint_integer(stderr, state.total_count);
              fprint(stderr, " unique sequences.\n");
            }
          else
            {
              fprint(stderr, "Found ");
              fprint_integer(stderr, state.chimera_count);
              fprint(stderr, " chimeras, ");
              fprint_integer(stderr, state.nonchimera_count);
              fprint(stderr, " non-chimeras,\nand ");
              fprint_integer(stderr, state.borderline_count);
              fprint(stderr, " borderline sequences in ");
              fprint_integer(stderr, state.total_count);
              fprint(stderr, " unique sequences.\n");
            }
        }

      if (state.total_abundance > 0)
        {
          if (mode == ChimeraMode::chimeras_denovo)
            {
              fprint(stderr, "Taking abundance information into account, this corresponds to\n");
              fprint_integer(stderr, state.chimera_abundance);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * static_cast<double>(state.chimera_abundance) / static_cast<double>(state.total_abundance));
              fprint(stderr, "%) chimeras and ");
              fprint_integer(stderr, state.nonchimera_abundance);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * static_cast<double>(state.nonchimera_abundance) / static_cast<double>(state.total_abundance));
              fprint(stderr, "%) non-chimeras in ");
              fprint_integer(stderr, state.total_abundance);
              fprint(stderr, " total sequences.\n");
            }
          else
            {
              fprint(stderr, "Taking abundance information into account, this corresponds to\n");
              fprint_integer(stderr, state.chimera_abundance);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * static_cast<double>(state.chimera_abundance) / static_cast<double>(state.total_abundance));
              fprint(stderr, "%) chimeras, ");
              fprint_integer(stderr, state.nonchimera_abundance);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * static_cast<double>(state.nonchimera_abundance) / static_cast<double>(state.total_abundance));
              fprint(stderr, "%) non-chimeras,\nand ");
              fprint_integer(stderr, state.borderline_abundance);
              fprint(stderr, " (");
              std::fprintf(stderr, "%.1f", 100.0 * static_cast<double>(state.borderline_abundance) / static_cast<double>(state.total_abundance));
              fprint(stderr, "%) borderline sequences in ");
              fprint_integer(stderr, state.total_abundance);
              fprint(stderr, " total sequences.\n");
            }
        }
      else
        {
          if (mode == ChimeraMode::chimeras_denovo)
            {
              fprint(stderr, "Taking abundance information into account, this corresponds to\n");
              fprint_integer(stderr, state.chimera_abundance);
              fprint(stderr, " chimeras, ");
              fprint_integer(stderr, state.nonchimera_abundance);
              fprint(stderr, " non-chimeras in ");
              fprint_integer(stderr, state.total_abundance);
              fprint(stderr, " total sequences.\n");
            }
          else
            {
              fprint(stderr, "Taking abundance information into account, this corresponds to\n");
              fprint_integer(stderr, state.chimera_abundance);
              fprint(stderr, " chimeras, ");
              fprint_integer(stderr, state.nonchimera_abundance);
              fprint(stderr, " non-chimeras,\nand ");
              fprint_integer(stderr, state.borderline_abundance);
              fprint(stderr, " borderline sequences in ");
              fprint_integer(stderr, state.total_abundance);
              fprint(stderr, " total sequences.\n");
            }
        }
    }

  if (parameters.fp_log != nullptr)
    {
      std::fputs(parameters.input_filename, parameters.fp_log);

      if (state.seqno > 0)
        {
          fprint(parameters.fp_log, ": ");
          fprint_integer(parameters.fp_log, state.chimera_count);
          fprint(parameters.fp_log, '/');
          fprint_integer(parameters.fp_log, state.seqno);
          fprint(parameters.fp_log, " chimeras (");
          std::fprintf(parameters.fp_log, "%.1f", 100.0 * state.chimera_count / state.seqno);
          fprint(parameters.fp_log, "%)\n");
        }
      else
        {
          fprint(parameters.fp_log, ": ");
          fprint_integer(parameters.fp_log, state.chimera_count);
          fprint(parameters.fp_log, '/');
          fprint_integer(parameters.fp_log, state.seqno);
          fprint(parameters.fp_log, " chimeras\n");
        }
    }


  if (mode == ChimeraMode::uchime_ref)
    {
      state.query_fasta_h->report_stripped_warning(parameters);
    }

  state.dbindex.clear();
  state.db.clear();

}


/* === Library API implementation === */

auto chimera_info_alloc() -> struct chimera_info_s *
{
  return new chimera_info_s {};
}

auto chimera_info_free(struct chimera_info_s * ci) -> void
{
  delete ci;
}

/* Build the chimera-detection configuration shared by the CLI (chimera()) and
   the library (chimera_detect_thread_init): a copy of the caller's parameters
   with the detection knobs applied — maxaccepts/maxrejects/id, and in denovo
   mode self/selfid/maxsizeratio.

   Chimera detection has no weak band: the candidate parents are exactly the hits
   accepted at chimera_id (chimera_process_query keeps only hit.accepted and
   discards the weak hits search_acceptable_aligned would otherwise retain), and
   the weak/rejected distinction does not affect the maxaccepts/maxrejects search
   termination (both increment `rejects`). We therefore set weak_id == id (an
   empty [id, id) band): this is deterministic and identical on the CLI and
   library paths, and it never inherits the opt_weak_id 10.0 sentinel — which,
   left unclamped, would make search_acceptable_aligned's gate reject everything.

   opt_threads is left at the caller's value; the CLI overrides it to 1 for the
   order-dependent denovo pool (see chimera()). */
static auto chimera_detection_parameters(struct Parameters const & parameters,
                                         ChimeraMode const mode) -> struct Parameters
{
  struct Parameters detection = parameters;
  detection.opt_maxaccepts = few;
  detection.opt_maxrejects = rejects;
  detection.opt_id = chimera_id;
  detection.opt_weak_id = detection.opt_id;

  /* For denovo mode, set opt_self/opt_selfid so sequences don't match
     themselves as candidate parents, and set opt_maxsizeratio for
     abundance skew filtering. */
  /* Denovo detection is what the caller asked for, not "opt_uchime_ref is
     null": on the library path no command option was ever set at all, and the
     documented default there is reference-based detection (LIBRARY_API.md), so
     absence must not select the denovo knobs. The mode says which it is, and
     the library entry point states uchime_ref explicitly. */
  if (chimera_is_denovo(mode))
    {
      detection.opt_self = 1;
      detection.opt_selfid = 1;
      detection.opt_maxsizeratio = 1.0 / parameters.opt_abskew;
    }
  return detection;
}


auto chimera_session_init(struct Parameters const & /*parameters*/) -> void
{
  /* Session-level initialization for chimera detection. Formerly overwrote the
     search-shaping opt_* globals to the detection defaults; those are now built
     per-thread in chimera_detect_thread_init and read through ci->parameters,
     so this no longer mutates any global (E1). Kept as a stable API symbol.

     The detection core takes no output lock on the library path (it receives a
     null CLI-output sink), so there is no shared mutex to initialize here. */
}


auto chimera_session_cleanup() -> void
{
  /* nothing to release: the detection core holds no file-static output
     state. Kept as a stable API symbol. */
}


auto chimera_detect_thread_init(struct chimera_info_s * ci, struct Parameters const & parameters,
                                struct Dbindex const & dbindex, struct Database const & db,
                                ChimeraMode const mode) -> void
{
  /* Per-thread initialization: SIMD aligners, k-mer finders, working
     buffers. Safe to call concurrently for different ci instances.
     Build this thread's chimera-detection configuration (formerly the globals
     set by chimera_session_init) and store it in ci so the detection core reads
     it through ci->parameters (E1). tophits sizes the per-part minheaps; it is
     derived from those options here rather than read from a shared file-static
     (E4). */
  ci->detection_parameters = chimera_detection_parameters(parameters, mode);
  int const tophits = static_cast<int>(ci->detection_parameters.opt_maxaccepts +
                                       ci->detection_parameters.opt_maxrejects);

  chimera_thread_init(ci, tophits, ci->detection_parameters, dbindex, db, mode);

  /* Allocate per-thread working state for chimera_process_query.
     These mirror the locals in chimera_thread_core but persist
     across calls to chimera_detect_single. */
  ci->api_allhits_list.resize(maxcandidates);

  struct Scoring const scoring = scoring_from_options(parameters);
  ci->api_lma_ptr = make_unique<LinearMemoryAligner>(scoring);
}


auto chimera_detect_init(struct chimera_info_s * ci, struct Parameters const & parameters,
                         struct Dbindex const & dbindex, struct Database const & db,
                         ChimeraMode const mode) -> void
{
  /* Convenience wrapper: session init + per-thread init in one call.
     Use for single-threaded detection (one chimera_info_s per session).
     For multi-threaded detection, call chimera_session_init() once then
     chimera_detect_thread_init() per thread. */
  chimera_session_init(parameters);
  chimera_detect_thread_init(ci, parameters, dbindex, db, mode);
}

auto chimera_detect_single(struct chimera_info_s * ci,
                           struct query_record_s const & query,
                           struct chimera_result_s * result) -> int
{
  /* Validate the caller-supplied query. The S18 hazard the checks here used to
     guard against is now unrepresentable: the buffers were sized from a
     query_len the caller passed separately, then filled from the C-string it
     pointed at, so a query_len shorter than strlen(query_seq) overflowed
     ci->query_seq on the heap. A View carries the bytes and their count as one
     value, which cannot disagree with itself, so that check has nothing left to
     compare. Likewise the null-header check: a header of length zero is an
     empty header, which is legal input and was already spelled "" before.

     What survives is the empty query, which the detection core does not handle,
     and a null data pointer on a non-empty view -- View's constructor asserts
     against that, but asserts are stripped under NDEBUG and this is a library
     entry point. There is no recoverable error channel (the function always
     returns 0), so an invalid call is fatal, as elsewhere in vsearch. */
  if (query.sequence.empty())
    {
      fatal("chimera_detect_single: query sequence must not be empty");
    }
  if (query.sequence.data() == nullptr)
    {
      fatal("chimera_detect_single: query sequence must not be null");
    }

  /* Populate query in the chimera_info_s.
     ci is per-thread state — must NOT be shared across threads. */
  ci->query_no = 0;
  ci->query_len = static_cast<int>(query.sequence.size());
  ci->query_size = query.abundance;

  realloc_arrays(ci, *ci->db, query.header.size());

  ci->query_head = copy_into_scratch(query.header, ci->query_head_v);
  copy_into_scratch(query.sequence, ci->query_seq);

  /* Clear result. Non-chimeric results will have only query_label and
     flag='N' populated; all other fields remain zero. */
  *result = {};
  ci->result_out = result;

  /* Use the SAME processing code as the CLI path: the detection core
     populates ci->result_out, writes no files and takes no output lock (the
     CLI writes its reports afterwards, from process_query). */
  auto const status = chimera_process_query(ci, ci->api_allhits_list,
                                            *ci->api_lma_ptr, *ci->db);

  if (status == Status::no_parents)
    {
      /* Populate result for no-parents case */
      copy_label(result->query_label, ci->query_head);
      result->flag = 'N';
    }

  ci->result_out = nullptr;
  return 0;
}

auto chimera_detect_thread_cleanup(struct chimera_info_s * ci) -> void
{
  /* Per-thread cleanup: frees all resources allocated by
     chimera_detect_thread_init (SIMD aligners, unique k-mer finders,
     minheaps, CIGAR strings, linear memory aligner). */
  ci->nwcigar.clear();
  ci->nwcigar.shrink_to_fit();
  chimera_thread_exit(ci);

  /* Release API working state */
  ci->api_lma_ptr.reset();
  ci->api_allhits_list.clear();
  ci->api_allhits_list.shrink_to_fit();
}


auto chimera_detect_cleanup(struct chimera_info_s * ci) -> void
{
  /* Convenience wrapper: per-thread cleanup + session cleanup in one call.
     Use for single-threaded detection (one chimera_info_s per session).
     For multi-threaded detection, call chimera_detect_thread_cleanup()
     per thread, then chimera_session_cleanup() once. */
  chimera_detect_thread_cleanup(ci);
  chimera_session_cleanup();
}


/* === Batch chimera detection API === */


/* Deleter for a per-thread chimera_info_s owned via unique_ptr: runs the same
   teardown as the normal path (chimera_detect_thread_cleanup + chimera_info_free)
   and is safe on a handle that was allocated but only partially initialised, so
   a fatal() unwinding mid-init frees every element created so far. noexcept: the
   teardown frees already-allocated buffers and never fatal()s. */
struct chimera_info_thread_deleter {
  auto operator()(struct chimera_info_s * ci) const noexcept -> void
  {
    chimera_detect_thread_cleanup(ci);
    chimera_info_free(ci);
  }
};

struct chimera_batch_context_s {
  View<struct query_record_s> queries;
  Span<struct chimera_result_s> results;

  /* per-thread chimera state arrays (sized to opt_threads). Owned unique_ptrs so
     a fatal() during per-thread init unwinds them, freeing every element built
     so far and the array itself. */
  std::vector<std::unique_ptr<struct chimera_info_s, chimera_info_thread_deleter>> ci_array;

  /* work-stealing counter */
  std::mutex mutex;
  int next_query;
};


static auto chimera_batch_worker_fn(struct chimera_batch_context_s & ctx,
                                    uint64_t const tid) -> void
{
  struct chimera_info_s * ci = ctx.ci_array[tid].get();

  int qi {0};

  auto const has_work_to_claim = [&]() -> bool {
    qi = ctx.next_query++;
    return static_cast<std::size_t>(qi) < ctx.queries.size();
  };

  auto const process_query = [&]() -> void {
    auto const index = static_cast<std::size_t>(qi);
    chimera_detect_single(ci, ctx.queries[index], &ctx.results[index]);
  };

  run_worker_loop(ctx.mutex, has_work_to_claim, process_query);
}


auto chimera_detect_batch(struct Parameters const & parameters,
                          struct Dbindex const & dbindex,
                          struct Database const & db,
                          View<struct query_record_s> const queries,
                          Span<struct chimera_result_s> const results,
                          ChimeraMode const mode) -> void
{
  assert(results.size() == queries.size());
  if (queries.empty())
    {
      return;
    }

  int const nthreads = std::max(1, static_cast<int>(parameters.opt_threads));

  /* Session-level init (no longer mutates globals; the per-thread detection
     configuration is built in chimera_detect_thread_init) */
  chimera_session_init(parameters);

  /* Allocate per-thread chimera state */
  struct chimera_batch_context_s ctx;
  ctx.queries = queries;
  ctx.results = results;
  ctx.next_query = 0;

  ctx.ci_array.reserve(static_cast<size_t>(nthreads));

  for (int t = 0; t < nthreads; t++)
    {
      /* own the handle before initialising it, so a fatal() in
         chimera_detect_thread_init frees this element and all prior ones. */
      ctx.ci_array.emplace_back(chimera_info_alloc());
      chimera_detect_thread_init(ctx.ci_array.back().get(), parameters, dbindex, db,
                                 mode);
    }

  /* run all queries through the worker pool (work-stealing on next_query) */
  {
    ThreadRunner threadrunner(static_cast<std::size_t>(nthreads),
                              [&ctx](uint64_t const tid) -> void {
                                chimera_batch_worker_fn(ctx, tid);
                              });
    threadrunner.run();
  }

  /* Cleanup per-thread state: clearing the vector runs the deleter
     (chimera_detect_thread_cleanup + chimera_info_free) on each element. */
  ctx.ci_array.clear();

  /* Session-level cleanup */
  chimera_session_cleanup();
}
