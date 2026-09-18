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

/*
  fastq_denoise: model-based correction of amplicon sequencing errors.

  This is an implementation of the divisive amplicon denoising algorithm (DADA)
  as it is *actually implemented* in the dada2 R package (Callahan et al. 2016,
  Nat Methods 13:581), which differs in several ways from the description in
  Rosen et al. 2012 (BMC Bioinformatics 13:283). In particular there is no
  expectation-maximisation over a mixture model; instead:

    1. Reads are dereplicated into "uniques", each with an abundance and a
       per-position mean quality score.
    2. All uniques start in one partition whose center is the most abundant
       unique. Every unique is compared with the center, giving lambda: the
       probability that a read of the center is observed as this unique,
         lambda = prod_j  err[center_nt(j) -> unique_nt(j)][quality(j)]
    3. Abundance p-value: given lambda and n reads in the partition, the number
       of reads of the unique expected from errors alone is Poisson(n * lambda).
       The p-value is P(X >= observed | X >= 1).
    4. The unique with the smallest p-value, if Bonferroni-significant
       (p * n_uniques < omega_a), "buds" off as the center of a new partition.
    5. All uniques are compared with the new center, then "shuffled": each
       unique joins the partition that maximises  lambda * partition_reads.
       Repeat from 3 until nothing is significant.
    6. Each unique is corrected to the center of its partition.

  The error matrix err[16 transitions][quality] is not known in advance. It is
  learnt by alternating steps 2-6 with re-estimation of err from the observed
  center -> read transitions (loess fit of log10 error rate against quality,
  weighted by number of observations), starting from the most pessimistic
  matrix possible (all ones, a single partition), until the matrix repeats
  itself ("self-consistency") or --denoise_maxconsist rounds have passed.

  Indels. dada2's error model knows substitutions only: gaps are skipped when
  lambda is computed, so a unique that differs from its center by indels alone
  is at distance zero, is never tested, and is silently merged into its
  center. That is --denoise_indels ignore, the default, for parity with dada2.
  With --denoise_indels model, every interior gap position is an error event
  too, with its own rate (one for insertions, one for deletions, independent
  of quality since a deleted base has none), learnt in the same
  self-consistency loop as the substitution rates:
    lambda = prod_j err[...] * p_ins^n_ins * p_del^n_del * (1-p_ins-p_del)^L
  The distance becomes substitutions + gap positions, so that indel-only
  variants are tested like any other, and are split off when they are too
  abundant to be explained by indel errors. Terminal gaps stay free in both
  modes: a length difference at the end of a read cannot be told from
  trimming. This matters for markers such as 12S and 16S rRNA, in which
  closely related species are often separated by indels only.

  Deliberate deviations from dada2, all documented where they occur:
    - probabilities are kept in log space, so p-values below 1e-308 do not
      collapse to zero and tie;
    - the loess fit is evaluated directly at each quality score, where R's
      default interpolates between k-d tree vertices;
    - the gapless shortcut is not taken when mismatches cluster at an end of
      the reads, where the k-mer test cannot see an indel;
    - pairs that need a true alignment go through vsearch's 16-channel SIMD
      aligner (search16), one center against many uniques at a time, where
      dada2 aligns pair by pair within a band of 16. The alignment is
      therefore unbanded: an indel longer than 16 is found, not missed.
*/

#include "commands/fastq_denoise.hpp"
#include "core/align_simd.hpp"  // search16, CELL
#include "core/attributes.hpp"  // OutputAnnotations
#include "core/db.hpp"  // struct Database
#include "core/eestats.hpp"  // vsearch::QualityScoreTable
#include "core/fasta.hpp"  // fasta_print_general
#include "core/fastq.hpp"  // fastq_open, fastq_print_general
#include "core/fastx.hpp"
#include "core/linmemalign.hpp"  // LinearMemoryAligner, struct Scoring
#include "core/quality_range.hpp"  // vsearch::check_quality_score
#include "core/seq_record.hpp"  // SeqRecord
#include "utils/base_mapping.hpp"  // Mapping
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"  // open_optional_output_file
#include "utils/progress.hpp"
#include "utils/span.hpp"  // Span
#include "utils/threads.hpp"  // ThreadRunner
#include "utils/view.hpp"  // View<char>
#include "vsearch.hpp"
#include <algorithm>  // std::find, std::max, std::min, std::nth_element, std::stable_sort
#include <array>
#include <cmath>  // std::log, std::exp, std::expm1, std::log1p, std::lgamma
#include <cstddef>  // std::size_t
#include <cstdint>  // uint8_t, uint16_t, uint32_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <fstream>
#include <iterator>  // std::next
#include <limits>
#include <memory>  // std::unique_ptr
#include <numeric>  // std::iota
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>  // std::move, std::swap
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* constants taken from dada2 (R/dada.R, src/dada.h) */
  constexpr auto kmer_length = std::size_t{5};
  constexpr auto n_kmers = std::size_t{1} << (2U * kmer_length);  // 4^5
  constexpr auto kmer_distance_cutoff = 0.42;  // KDIST_CUTOFF
  constexpr auto terminal_window = 2 * kmer_length;  // see has_terminal_mismatch_cluster()
  constexpr auto max_shuffles = 10;  // MAX_SHUFFLE
  constexpr auto min_error_rate = 1e-7;  // MIN_ERROR_RATE (loessErrfun)
  constexpr auto max_error_rate = 0.25;  // MAX_ERROR_RATE (loessErrfun)
  constexpr auto loess_span = 0.75;  // default span of R's loess()
  constexpr auto min_n_qualities = std::size_t{41};
  constexpr auto alignment_match = 5;  // MATCH
  constexpr auto alignment_mismatch = -4;  // MISMATCH
  constexpr auto alignment_gap = 8;  // GAP_PENALTY, linear

  constexpr auto n_nucleotides = std::size_t{4};
  constexpr auto n_transitions = n_nucleotides * n_nucleotides;
  constexpr auto unaligned = uint8_t{4};  // read position facing a gap
  constexpr auto no_unique = std::numeric_limits<uint32_t>::max();
  constexpr auto minus_infinity = -std::numeric_limits<double>::infinity();

  constexpr std::array<char, n_nucleotides> nucleotide_symbols = {{'A', 'C', 'G', 'T'}};


  /* one row of the error matrix per (true nucleotide, observed nucleotide),
     in dada2's order: A2A, A2C, A2G, A2T, C2A, ... */
  auto transition_index(uint8_t const from, uint8_t const to) -> std::size_t {
    return (n_nucleotides * from) + to;
  }


  /* -------------------------------------------------------------------- */
  /*  data structures                                                     */
  /* -------------------------------------------------------------------- */

  /* What dada2 calls a 'Comparison': the outcome of aligning one unique with
     the center of one partition. */
  struct Comparison {
    uint32_t partition = 0;
    uint32_t unique = 0;
    double log_lambda = minus_infinity;
    uint32_t distance = 0;  // substitutions (+ interior gap positions)
  };


  /* What dada2 calls a 'Raw': a dereplicated sequence. */
  struct Unique {
    /* The header (of the first read seen with this sequence) and the sequence
       itself are moved into a Database once dereplication is over, under the
       index of the unique: search16 aligns database entries. */
    std::string header;
    std::string sequence;  // A, C, G, T
    std::vector<uint8_t> nucleotides;  // 0, 1, 2, 3
    std::vector<double> quality_sums;  // summed during dereplication ...
    std::vector<uint8_t> qualities;  // ... then rounded mean, per position
    std::vector<uint8_t> kmer_counts;  // 4^k saturating counters
    std::vector<uint16_t> kmer_order;  // the k-mer starting at each position
    uint64_t reads = 0;

    /* state of the partitioning, reset at each round */
    uint32_t partition = 0;
    Comparison comparison;  // with the center of its current partition
    double log_pvalue = 0.0;
    double log_expected_minmax = minus_infinity;  // E_minmax
    bool is_locked = false;
    bool is_corrected = true;
  };


  /* What dada2 calls a 'Bi'. */
  struct Partition {
    uint32_t center = 0;
    std::vector<uint32_t> members;
    uint64_t reads = 0;
    bool pvalues_are_stale = true;  // update_e
    bool locks_are_unchecked = true;  // check_locks
    std::vector<Comparison> comparisons;  // uniques this partition could attract
  };


  /* The error model: P(observed nucleotide | true nucleotide, quality).
     Stored as natural logs, since it is only ever used in products over a
     whole read. */
  class ErrorModel {
  public:
    explicit ErrorModel(std::size_t const n_qualities)
      : n_qualities_(n_qualities),
        rates_(n_transitions * n_qualities, 1.0),
        log_rates_(n_transitions * n_qualities, 0.0) {}

    auto n_qualities() const -> std::size_t { return n_qualities_; }

    auto rate(std::size_t const transition, std::size_t const quality) const -> double {
      return rates_[(transition * n_qualities_) + quality];
    }

    auto log_rate(std::size_t const transition, std::size_t const quality) const -> double {
      return log_rates_[(transition * n_qualities_) + quality];
    }

    auto set_rate(std::size_t const transition, std::size_t const quality, double const value) -> void {
      rates_[(transition * n_qualities_) + quality] = value;
      log_rates_[(transition * n_qualities_) + quality] = std::log(value);
    }

    /* Indel rates, per position, used with --denoise_indels model only. The
       probability of no indel at a position is set apart, as the
       self-transitions of the substitution matrix are: it must be forced to
       one in the starting model, where all rates are one. */
    auto insertion_rate() const -> double { return insertion_rate_; }
    auto deletion_rate() const -> double { return deletion_rate_; }
    auto log_insertion_rate() const -> double { return std::log(insertion_rate_); }
    auto log_deletion_rate() const -> double { return std::log(deletion_rate_); }
    auto log_no_indel_rate() const -> double { return log_no_indel_rate_; }

    auto set_indel_rates(double const insertion, double const deletion) -> void {
      insertion_rate_ = insertion;
      deletion_rate_ = deletion;
      log_no_indel_rate_ = std::log(std::max(1.0 - insertion - deletion, min_error_rate));
    }

    auto set_no_indel_rate_to_one() -> void { log_no_indel_rate_ = 0.0; }

  private:
    std::size_t n_qualities_;
    std::vector<double> rates_;
    std::vector<double> log_rates_;
    double insertion_rate_ = 1.0;
    double deletion_rate_ = 1.0;
    double log_no_indel_rate_ = 0.0;
  };


  /* What the reads, once assigned to their centers, say about errors. */
  struct ErrorCounts {
    std::vector<uint64_t> substitutions;  // n_transitions x n_qualities
    uint64_t insertions = 0;  // read bases facing an interior gap
    uint64_t deletions = 0;  // center bases facing an interior gap
    uint64_t positions = 0;  // aligned pairs of bases

    auto operator==(ErrorCounts const & rhs) const -> bool {
      return (substitutions == rhs.substitutions) and (insertions == rhs.insertions)
        and (deletions == rhs.deletions) and (positions == rhs.positions);
    }
  };


  /* -------------------------------------------------------------------- */
  /*  Poisson abundance p-value                                           */
  /* -------------------------------------------------------------------- */

  /* log P(X >= reads), X ~ Poisson(mean), with the mean given as its log.

     dada2 calls R's ppois(); vsearch has no Rmath. Above the mean the tail is
     summed upwards from its first term,
       P(X >= x) = pmf(x) * (1 + m/(x+1) + m^2/((x+1)(x+2)) + ...),
     which converges geometrically and loses no precision however small the
     result. At or below the mean the p-value is large and uninteresting, so
     the complement of the lower tail is accurate enough. */
  auto log_poisson_upper_tail(uint64_t const reads, double const log_mean) -> double {
    if (reads == 0) {
      return 0.0;
    }
    if (log_mean == minus_infinity) {
      return minus_infinity;
    }
    auto const mean = std::exp(log_mean);
    auto const observed = static_cast<double>(reads);
    if (observed > mean) {
      auto const log_pmf = -mean + (observed * log_mean) - std::lgamma(observed + 1.0);
      auto sum = 1.0;
      auto term = 1.0;
      for (auto i = 1.0; term > sum * 1e-17; i += 1.0) {
        term *= mean / (observed + i);
        sum += term;
      }
      return log_pmf + std::log(sum);
    }
    /* lower tail P(X <= reads - 1), summed downwards from pmf(reads - 1) */
    auto const last = observed - 1.0;
    auto const log_pmf = -mean + (last * log_mean) - std::lgamma(last + 1.0);
    auto sum = 1.0;
    auto term = 1.0;
    for (auto i = last; (i > 0.0) and (term > sum * 1e-17); i -= 1.0) {
      term *= i / mean;
      sum += term;
    }
    auto const lower_tail = std::min(1.0, std::exp(log_pmf) * sum);
    return std::log1p(-lower_tail);
  }


  /* calc_pA() in dada2. When 'is_conditional', the p-value is conditioned on
     the unique having been observed at all (X >= 1): that is the test used to
     create partitions. The unconditional form is used for the final decision
     to correct a unique or not (omega_c). */
  auto log_abundance_pvalue(uint64_t const reads, double const log_expected,
                            bool const is_conditional) -> double {
    auto log_pvalue = log_poisson_upper_tail(reads, log_expected);
    if (is_conditional and (log_expected != minus_infinity)) {
      auto const expected = std::exp(log_expected);
      /* log(1 - exp(-E)); for tiny E that is log(E - E^2/2) */
      auto const log_norm = (expected < 1e-5)
        ? log_expected + std::log1p(-0.5 * expected)
        : std::log(-std::expm1(-expected));
      log_pvalue -= log_norm;
    }
    return std::min(log_pvalue, 0.0);
  }


  /* -------------------------------------------------------------------- */
  /*  k-mer screens                                                       */
  /* -------------------------------------------------------------------- */

  auto index_kmers(Unique & unique) -> void {
    auto const length = unique.nucleotides.size();
    auto const n_positions = length - kmer_length + 1;
    unique.kmer_counts.assign(n_kmers, 0);
    unique.kmer_order.assign(n_positions, 0);
    for (auto i = std::size_t{0}; i < n_positions; ++i) {
      auto kmer = std::size_t{0};
      for (auto j = i; j < i + kmer_length; ++j) {
        kmer = (n_nucleotides * kmer) + unique.nucleotides[j];
      }
      unique.kmer_order[i] = static_cast<uint16_t>(kmer);
      if (unique.kmer_counts[kmer] < std::numeric_limits<uint8_t>::max()) {
        ++unique.kmer_counts[kmer];
      }
    }
  }


  /* number of k-mers shared, regardless of their position (kmer_dist) */
  auto count_shared_kmers(Unique const & lhs, Unique const & rhs) -> std::size_t {
    auto shared = std::size_t{0};
    for (auto i = std::size_t{0}; i < n_kmers; ++i) {
      shared += std::min(lhs.kmer_counts[i], rhs.kmer_counts[i]);
    }
    return shared;
  }


  /* number of k-mers shared at the same position (kord_dist) */
  auto count_collinear_kmers(Unique const & lhs, Unique const & rhs) -> std::size_t {
    auto shared = std::size_t{0};
    for (auto i = std::size_t{0}; i < lhs.kmer_order.size(); ++i) {
      shared += (lhs.kmer_order[i] == rhs.kmer_order[i]) ? 1U : 0U;
    }
    return shared;
  }


  /* -------------------------------------------------------------------- */
  /*  comparison of a unique with a center                                */
  /* -------------------------------------------------------------------- */

  /* An alignment is reduced to what the error model needs: for each position
     of the unique, the center nucleotide facing it (or 'unaligned'), and the
     number of interior gap positions on each side. Terminal gaps are not
     counted: the alignment is ends-free, as a length difference at either end
     of a read says nothing about the molecule it was read from. */
  struct Alignment {
    std::vector<uint8_t> facing;
    uint32_t insertions = 0;  // bases of the unique facing an interior gap
    uint32_t deletions = 0;  // bases of the center facing an interior gap
  };


  enum struct IndelPolicy : unsigned char {
    ignore,  // dada2: gaps are invisible to the error model
    model,  // gap positions are errors with a rate of their own
  };


  auto ends_free_scoring() -> struct Scoring {
    struct Scoring scoring;
    scoring.match = alignment_match;
    scoring.mismatch = alignment_mismatch;
    scoring.gap_open_query_interior = 0;
    scoring.gap_open_target_interior = 0;
    scoring.gap_extension_query_interior = alignment_gap;
    scoring.gap_extension_target_interior = alignment_gap;
    /* terminal gaps: all eight penalties stay at zero */
    return scoring;
  }


  /* The k-mer test for the gapless shortcut is blind near the ends of the
     reads: no k-mer fits between an indel and an end less than k bases away,
     so the shifted tail of a read truncated to a fixed length passes for a
     few substitutions. Their joint probability is minute, and at high depth
     such reads were inferred as new sequences where they are single indel
     errors (or, with --denoise_indels ignore, nothing at all). Two or more
     mismatches within 2k bases of an end send the pair to the aligner, which
     settles the matter. dada2 has no such guard. */
  auto has_terminal_mismatch_cluster(Unique const & center, Unique const & unique) -> bool {
    auto const length = unique.nucleotides.size();
    auto const window = std::min(terminal_window, length);
    auto n_leading = 0U;
    auto n_trailing = 0U;
    for (auto i = std::size_t{0}; i < window; ++i) {
      n_leading += (center.nucleotides[i] != unique.nucleotides[i]) ? 1U : 0U;
      n_trailing += (center.nucleotides[length - 1 - i] != unique.nucleotides[length - 1 - i]) ? 1U : 0U;
    }
    return (n_leading >= 2) or (n_trailing >= 2);
  }


  /* The first half of raw_align() in dada2: what the k-mers say of a pair. */
  enum struct Screen : unsigned char {
    rejected,  // too few k-mers in common: lambda is taken to be zero
    gapless,  // no indel: the alignment is the trivial one
    gapped,  // needs a true alignment
  };


  auto screen(Unique const & center, Unique const & unique, double const kmer_cutoff) -> Screen {
    auto const center_length = center.nucleotides.size();
    auto const unique_length = unique.nucleotides.size();
    auto const n_positions = std::min(center_length, unique_length) - kmer_length + 1;
    auto const shared = count_shared_kmers(center, unique);
    auto const kmer_distance = 1.0 - (static_cast<double>(shared) / static_cast<double>(n_positions));
    if (kmer_distance > kmer_cutoff) {
      return Screen::rejected;
    }
    /* Gapless shortcut: if every shared k-mer is shared at the same position,
       no indel separates the two sequences. This is the case for nearly all
       pairs of Illumina amplicon reads, and it turns the comparison into a
       linear scan. */
    if ((center_length == unique_length) and (count_collinear_kmers(center, unique) == shared)
        and (not has_terminal_mismatch_cluster(center, unique))) {
      return Screen::gapless;
    }
    return Screen::gapped;
  }


  auto set_gapless(Unique const & center, Alignment & alignment) -> void {
    alignment.facing = center.nucleotides;
    alignment.insertions = 0;
    alignment.deletions = 0;
  }


  /* al2subs() in dada2. The cigar describes the alignment of the center (the
     query: 'D' consumes it) with the unique (the target: 'I' consumes it). */
  auto set_from_cigar(char const * cigar, Unique const & center, Unique const & unique,
                      Alignment & alignment) -> void {
    alignment.facing.assign(unique.nucleotides.size(), unaligned);
    alignment.insertions = 0;
    alignment.deletions = 0;
    auto center_pos = std::size_t{0};
    auto unique_pos = std::size_t{0};
    auto run = std::size_t{0};
    /* gap positions seen since the last aligned pair: they become interior
       once another aligned pair follows, and stay terminal otherwise */
    auto pending_insertions = uint32_t{0};
    auto pending_deletions = uint32_t{0};
    auto has_aligned_pair = false;
    for (auto const * p = cigar; *p != '\0'; ++p) {
      if ((*p >= '0') and (*p <= '9')) {
        run = (10 * run) + static_cast<std::size_t>(*p - '0');
        continue;
      }
      run = std::max(run, std::size_t{1});  // an omitted run length means 1
      switch (*p) {
      case 'M':
        for (auto i = std::size_t{0}; i < run; ++i) {
          alignment.facing[unique_pos + i] = center.nucleotides[center_pos + i];
        }
        if (has_aligned_pair) {  // otherwise these are leading gaps
          alignment.insertions += pending_insertions;
          alignment.deletions += pending_deletions;
        }
        has_aligned_pair = true;
        pending_insertions = 0;
        pending_deletions = 0;
        center_pos += run;
        unique_pos += run;
        break;
      case 'I':
        pending_insertions += static_cast<uint32_t>(run);
        unique_pos += run;
        break;
      default:  // 'D'
        pending_deletions += static_cast<uint32_t>(run);
        center_pos += run;
        break;
      }
      run = 0;
    }
  }


  /* One center against many uniques: the shape of problem search16 was written
     for. The query profile of the center is computed once, then the uniques
     stream through the 16 channels of the SIMD kernel, a finished alignment
     making room for the next unique. One CenterAligner per thread.

     search16 computes 16-bit scores. A pair whose score could overflow comes
     back with the maximal score as a sentinel and no cigar; it is then
     re-aligned with the scalar, linear-memory aligner, as vsearch does
     everywhere else. At 5 per match that takes sequences of 6.5 kb. */
  class CenterAligner {
  public:
    CenterAligner(struct Scoring const & scoring, Database const & database)
      : database_(database),
        simd_(search16_init(scoring)),
        scalar_(scoring) {}

    ~CenterAligner() { search16_exit(simd_); }
    CenterAligner(CenterAligner const &) = delete;
    auto operator=(CenterAligner const &) -> CenterAligner & = delete;
    CenterAligner(CenterAligner &&) = delete;
    auto operator=(CenterAligner &&) -> CenterAligner & = delete;

    /* Align the center with each of 'targets', and hand each alignment to
       'consume' (unique index, alignment). */
    template <typename Consumer>
    auto align(uint32_t const center_index, std::vector<Unique> const & uniques,
               std::vector<unsigned int> const & targets, Consumer consume) -> void {
      if (targets.empty()) {
        return;
      }
      auto const & center = uniques[center_index];
      auto const center_view = database_.sequence_view(center_index);
      search16_qprep(*simd_, center_view);
      auto const n_targets = targets.size();
      scores_.resize(n_targets);
      for (auto * counts : {&aligned_, &matches_, &mismatches_, &gaps_}) {
        counts->resize(n_targets);
      }
      cigars_.resize(n_targets);
      search16(*simd_,
               View<unsigned int>{targets.data(), n_targets},
               Span<CELL>{scores_.data(), n_targets},
               Span<unsigned short>{aligned_.data(), n_targets},
               Span<unsigned short>{matches_.data(), n_targets},
               Span<unsigned short>{mismatches_.data(), n_targets},
               Span<unsigned short>{gaps_.data(), n_targets},
               Span<std::string>{cigars_.data(), n_targets},
               database_);
      for (auto i = std::size_t{0}; i < n_targets; ++i) {
        auto const & unique = uniques[targets[i]];
        auto const * cigar = cigars_[i].c_str();
        if (scores_[i] == std::numeric_limits<CELL>::max()) {
          cigar = scalar_.align(center_view, database_.sequence_view(targets[i]));
        }
        set_from_cigar(cigar, center, unique, alignment_);
        consume(targets[i], alignment_);
      }
    }

  private:
    Database const & database_;
    s16info_s * simd_;
    LinearMemoryAligner scalar_;
    Alignment alignment_;
    std::vector<CELL> scores_;
    std::vector<unsigned short> aligned_;
    std::vector<unsigned short> matches_;
    std::vector<unsigned short> mismatches_;
    std::vector<unsigned short> gaps_;
    std::vector<std::string> cigars_;
  };


  /* compute_lambda() in dada2, plus the indel terms */
  auto compare(Unique const & unique, Alignment const & alignment,
               ErrorModel const & error_model, IndelPolicy const indel_policy,
               Comparison & comparison) -> void {
    auto log_lambda = 0.0;
    auto distance = uint32_t{0};
    for (auto i = std::size_t{0}; i < unique.nucleotides.size(); ++i) {
      auto const observed = unique.nucleotides[i];
      auto const truth = (alignment.facing[i] == unaligned) ? observed : alignment.facing[i];
      distance += (truth != observed) ? 1U : 0U;
      log_lambda += error_model.log_rate(transition_index(truth, observed), unique.qualities[i]);
    }
    if (indel_policy == IndelPolicy::model) {
      distance += alignment.insertions + alignment.deletions;
      log_lambda += static_cast<double>(unique.nucleotides.size()) * error_model.log_no_indel_rate();
      /* guarded: a rate may be one (log 0) or, read from a file, tiny */
      if (alignment.insertions != 0) {
        log_lambda += alignment.insertions * error_model.log_insertion_rate();
      }
      if (alignment.deletions != 0) {
        log_lambda += alignment.deletions * error_model.log_deletion_rate();
      }
    }
    comparison.log_lambda = log_lambda;
    comparison.distance = distance;
  }


  /* -------------------------------------------------------------------- */
  /*  the divisive partitioning                                           */
  /* -------------------------------------------------------------------- */

  class Denoiser {
  public:
    Denoiser(std::vector<Unique> & uniques, Database const & database,
             struct Parameters const & parameters)
      : uniques_(uniques),
        n_threads_(static_cast<std::size_t>(std::max<int64_t>(parameters.opt_threads, 1))),
        log_omega_a_(std::log(parameters.opt_denoise_omega_a)),
        log_omega_c_(std::log(parameters.opt_denoise_omega_c)),
        indel_policy_(parameters.opt_denoise_indels_model ? IndelPolicy::model : IndelPolicy::ignore),
        results_(uniques.size()),
        alignments_(n_threads_),
        batches_(n_threads_),
        thread_runner_(n_threads_, [this](uint64_t const nth_thread) -> void { worker(nth_thread); }) {
      auto const scoring = ends_free_scoring();
      for (auto i = std::size_t{0}; i < n_threads_; ++i) {
        aligners_.emplace_back(new CenterAligner(scoring, database));
      }
      for (auto const & unique : uniques_) {
        total_reads_ += unique.reads;
      }
    }

    /* run_dada() in dada2. Returns the observed transition counts. */
    auto run(ErrorModel const & error_model, std::size_t max_partitions) -> ErrorCounts {
      error_model_ = &error_model;
      initialize();
      compare_all(0, 1.0);  // no k-mer screen for the first partition
      update_pvalues();
      if (max_partitions == 0) {
        max_partitions = uniques_.size();
      }
      while (partitions_.size() < max_partitions) {
        auto const new_partition = bud();
        if (new_partition == 0) {
          break;
        }
        compare_all(new_partition, kmer_distance_cutoff);
        for (auto i = 0; (i < max_shuffles) and shuffle(); ++i) {
          /* keep shuffling until nobody moves */
        }
        update_pvalues();
      }
      decide_corrections();
      return count_errors();
    }

    auto partitions() const -> std::vector<Partition> const & { return partitions_; }

  private:
    auto initialize() -> void {
      partitions_.clear();
      partitions_.emplace_back();
      auto & first = partitions_.front();
      first.members.resize(uniques_.size());
      std::iota(first.members.begin(), first.members.end(), uint32_t{0});
      first.reads = total_reads_;
      first.center = 0;  // uniques are sorted by decreasing abundance
      for (auto & unique : uniques_) {
        unique.partition = 0;
        unique.comparison = Comparison{};
        unique.log_pvalue = 0.0;
        unique.log_expected_minmax = minus_infinity;
        unique.is_locked = false;
        unique.is_corrected = true;
      }
    }

    /* One thread's share of the uniques: a contiguous block would give the
       first thread all the abundant uniques, so the work is dealt out like
       cards instead. Each thread writes only to its own slots of results_. */
    auto worker(uint64_t const nth_thread) -> void {
      auto const center_index = partitions_[job_partition_].center;
      auto const & center = uniques_[center_index];
      auto & alignment = alignments_[nth_thread];
      auto & batch = batches_[nth_thread];
      batch.clear();
      for (auto index = static_cast<std::size_t>(nth_thread); index < uniques_.size(); index += n_threads_) {
        auto const & unique = uniques_[index];
        auto & result = results_[index];
        result = Comparison{};
        result.partition = static_cast<uint32_t>(job_partition_);
        result.unique = static_cast<uint32_t>(index);
        /* "greedy" mode of dada2: a unique more abundant than the center
           cannot be an error of it, and a locked unique is already explained
           by its own center alone; neither needs aligning. */
        if ((unique.reads > center.reads) or unique.is_locked) {
          continue;
        }
        switch (screen(center, unique, job_kmer_cutoff_)) {
        case Screen::rejected:
          break;
        case Screen::gapless:
          set_gapless(center, alignment);
          compare(unique, alignment, *error_model_, indel_policy_, result);
          break;
        case Screen::gapped:
          batch.push_back(static_cast<unsigned int>(index));
          break;
        }
      }
      /* the pairs that need a true alignment, all at once */
      aligners_[nth_thread]->align(center_index, uniques_, batch,
        [this](unsigned int const index, Alignment const & gapped_alignment) -> void {
          compare(uniques_[index], gapped_alignment, *error_model_, indel_policy_, results_[index]);
        });
    }

    /* b_compare() in dada2 */
    auto compare_all(std::size_t const partition_index, double const kmer_cutoff) -> void {
      job_partition_ = partition_index;
      job_kmer_cutoff_ = kmer_cutoff;
      thread_runner_.run();

      auto & partition = partitions_[partition_index];
      auto const log_total_reads = std::log(static_cast<double>(total_reads_));
      auto const log_center_reads = std::log(static_cast<double>(uniques_[partition.center].reads));
      for (auto index = std::size_t{0}; index < uniques_.size(); ++index) {
        auto & unique = uniques_[index];
        auto const & result = results_[index];
        /* Keep the comparison only if this partition could ever attract the
           unique, that is if it would do so even if it held every read. The
           first partition keeps them all: shuffle() relies on that. */
        auto const is_useful = (result.log_lambda + log_total_reads) > unique.log_expected_minmax;
        if ((partition_index != 0) and (not is_useful)) {
          continue;
        }
        unique.log_expected_minmax = std::max(unique.log_expected_minmax,
                                              result.log_lambda + log_center_reads);
        partition.comparisons.push_back(result);
        if ((partition_index == 0) or (index == partition.center)) {
          unique.comparison = result;
        }
      }
    }

    /* b_p_update() in dada2 */
    auto update_pvalues() -> void {
      for (auto & partition : partitions_) {
        auto const log_reads = std::log(static_cast<double>(partition.reads));
        if (partition.pvalues_are_stale) {
          for (auto const index : partition.members) {
            auto & unique = uniques_[index];
            /* Singletons are never significant, nor is a unique at distance
               zero from its center: under --denoise_indels ignore that
               includes the variants that differ by indels only. */
            if ((unique.reads == 1) or (unique.comparison.distance == 0)) {
              unique.log_pvalue = 0.0;
            } else {
              unique.log_pvalue = log_abundance_pvalue(unique.reads,
                                                       unique.comparison.log_lambda + log_reads,
                                                       true);
            }
          }
          partition.pvalues_are_stale = false;
        }
        if (partition.locks_are_unchecked) {
          auto const log_center_reads = std::log(static_cast<double>(uniques_[partition.center].reads));
          for (auto const index : partition.members) {
            auto & unique = uniques_[index];
            auto const log_expected = unique.comparison.log_lambda + log_center_reads;
            if ((index == partition.center) or (log_expected > std::log(static_cast<double>(unique.reads)))) {
              unique.is_locked = true;
            }
          }
          partition.locks_are_unchecked = false;
        }
      }
    }

    /* b_bud() in dada2. Returns the index of the new partition, or zero. */
    auto bud() -> std::size_t {
      auto best = no_unique;
      auto best_log_pvalue = 0.0;
      for (auto const & partition : partitions_) {
        for (auto const index : partition.members) {
          auto const & unique = uniques_[index];
          if ((index == partition.center) or (unique.comparison.distance == 0)) {
            continue;
          }
          auto const is_better = (unique.log_pvalue < best_log_pvalue)
            or ((best != no_unique) and (unique.log_pvalue == best_log_pvalue)
                and (unique.reads > uniques_[best].reads));
          if (is_better) {
            best = index;
            best_log_pvalue = unique.log_pvalue;
          }
        }
      }
      /* Bonferroni correction for the number of uniques tested */
      auto const log_n_tests = std::log(static_cast<double>(uniques_.size()));
      if ((best == no_unique) or ((best_log_pvalue + log_n_tests) >= log_omega_a_)) {
        return 0;
      }
      move(best, partitions_.size());
      auto & partition = partitions_.back();
      partition.center = best;
      uniques_[best].is_locked = false;
      return partitions_.size() - 1;
    }

    /* b_shuffle2() in dada2 */
    auto shuffle() -> bool {
      std::vector<Comparison const *> best(uniques_.size(), nullptr);
      std::vector<double> best_log_expected(uniques_.size(), minus_infinity);
      for (auto const & partition : partitions_) {
        auto const log_reads = std::log(static_cast<double>(partition.reads));
        for (auto const & comparison : partition.comparisons) {
          auto const log_expected = comparison.log_lambda + log_reads;
          if ((best[comparison.unique] == nullptr) or (log_expected > best_log_expected[comparison.unique])) {
            best[comparison.unique] = &comparison;
            best_log_expected[comparison.unique] = log_expected;
          }
        }
      }
      auto has_shuffled = false;
      for (auto index = std::size_t{0}; index < uniques_.size(); ++index) {
        auto & unique = uniques_[index];
        auto const destination = best[index]->partition;
        if ((destination == unique.partition) or (index == partitions_[unique.partition].center)) {
          continue;  // the center of a partition cannot leave it
        }
        move(static_cast<uint32_t>(index), destination);
        unique.comparison = *best[index];
        has_shuffled = true;
      }
      return has_shuffled;
    }

    auto move(uint32_t const index, std::size_t const destination) -> void {
      if (destination == partitions_.size()) {
        partitions_.emplace_back();  // before any reference is taken: may reallocate
      }
      auto & unique = uniques_[index];
      auto & source = partitions_[unique.partition];
      source.members.erase(std::find(source.members.begin(), source.members.end(), index));
      source.reads -= unique.reads;
      source.pvalues_are_stale = true;
      auto & target = partitions_[destination];
      target.members.push_back(index);
      target.reads += unique.reads;
      target.pvalues_are_stale = true;
      unique.partition = static_cast<uint32_t>(destination);
    }

    /* A unique that is significantly too abundant to be an error of its
       center, yet was not significant enough to found a partition, is left
       uncorrected (omega_c). This time the p-value is not conditional. */
    auto decide_corrections() -> void {
      for (auto const & partition : partitions_) {
        auto const log_reads = std::log(static_cast<double>(partition.reads));
        for (auto const index : partition.members) {
          auto & unique = uniques_[index];
          if (index == partition.center) {
            continue;
          }
          auto const log_pvalue = log_abundance_pvalue(unique.reads,
                                                       unique.comparison.log_lambda + log_reads,
                                                       false);
          unique.is_corrected = (log_pvalue >= log_omega_c_);
        }
      }
    }

    /* b_make_transition_by_quality_matrix() in dada2: every corrected read
       is, position by position, an observation of center -> read at the
       quality of the read. Interior gap positions are tallied as well; they
       feed the indel rates of --denoise_indels model. */
    auto count_errors() -> ErrorCounts {
      ErrorCounts counts;
      counts.substitutions.assign(n_transitions * error_model_->n_qualities(), 0);
      auto & alignment = alignments_.front();
      auto & batch = batches_.front();
      for (auto const & partition : partitions_) {
        auto const & center = uniques_[partition.center];
        batch.clear();
        for (auto const index : partition.members) {
          if (not uniques_[index].is_corrected) {
            continue;
          }
          if (screen(center, uniques_[index], 1.0) == Screen::gapless) {
            set_gapless(center, alignment);
            tally(uniques_[index], alignment, counts);
          } else {
            batch.push_back(index);
          }
        }
        aligners_.front()->align(partition.center, uniques_, batch,
          [this, &counts](unsigned int const index, Alignment const & gapped_alignment) -> void {
            tally(uniques_[index], gapped_alignment, counts);
          });
      }
      return counts;
    }

    auto tally(Unique const & unique, Alignment const & alignment, ErrorCounts & counts) const -> void {
      auto const n_qualities = error_model_->n_qualities();
      for (auto i = std::size_t{0}; i < unique.nucleotides.size(); ++i) {
        if (alignment.facing[i] == unaligned) {
          continue;
        }
        auto const transition = transition_index(alignment.facing[i], unique.nucleotides[i]);
        counts.substitutions[(transition * n_qualities) + unique.qualities[i]] += unique.reads;
        counts.positions += unique.reads;
      }
      counts.insertions += alignment.insertions * unique.reads;
      counts.deletions += alignment.deletions * unique.reads;
    }

    std::vector<Unique> & uniques_;
    std::size_t n_threads_;
    double log_omega_a_;
    double log_omega_c_;
    IndelPolicy indel_policy_;
    uint64_t total_reads_ = 0;
    ErrorModel const * error_model_ = nullptr;
    std::vector<Partition> partitions_;
    std::vector<Comparison> results_;  // one slot per unique, filled by the workers
    std::vector<Alignment> alignments_;  // one scratch buffer per thread
    std::vector<std::vector<unsigned int>> batches_;  // per thread: uniques to align
    std::vector<std::unique_ptr<CenterAligner>> aligners_;  // one per thread
    std::size_t job_partition_ = 0;
    double job_kmer_cutoff_ = 1.0;
    ThreadRunner thread_runner_;  // last member: threads start once all else is built
  };


  /* -------------------------------------------------------------------- */
  /*  error model estimation (loessErrfun in dada2)                       */
  /* -------------------------------------------------------------------- */

  /* Weighted local polynomial regression, evaluated at x0: the value at zero
     of the polynomial of the given degree fitted to (x - x0, y) with weights
     prior weight x tricube(distance / bandwidth). Solved through the normal
     equations, which are at most 3 x 3. Returns false if they are singular. */
  auto fit_local_polynomial(std::vector<double> const & xs, std::vector<double> const & ys,
                            std::vector<double> const & weights, double const x0,
                            double const bandwidth, std::size_t const degree,
                            double & prediction) -> bool {
    auto const n_terms = degree + 1;
    std::array<std::array<double, 4>, 3> system = {{}};  // augmented matrix
    for (auto i = std::size_t{0}; i < xs.size(); ++i) {
      auto const distance = std::abs(xs[i] - x0) / bandwidth;
      if (distance >= 1.0) {
        continue;
      }
      auto const tricube = std::pow(1.0 - std::pow(distance, 3), 3);
      auto const weight = weights[i] * tricube;
      std::array<double, 3> const powers = {{1.0, xs[i] - x0, (xs[i] - x0) * (xs[i] - x0)}};
      for (auto row = std::size_t{0}; row < n_terms; ++row) {
        for (auto col = std::size_t{0}; col < n_terms; ++col) {
          system[row][col] += weight * powers[row] * powers[col];
        }
        system[row][n_terms] += weight * powers[row] * ys[i];
      }
    }
    /* Gaussian elimination with partial pivoting */
    for (auto pivot = std::size_t{0}; pivot < n_terms; ++pivot) {
      auto best_row = pivot;
      for (auto row = pivot + 1; row < n_terms; ++row) {
        if (std::abs(system[row][pivot]) > std::abs(system[best_row][pivot])) {
          best_row = row;
        }
      }
      std::swap(system[pivot], system[best_row]);
      if (std::abs(system[pivot][pivot]) < 1e-10 * std::abs(system[0][0]) or (system[pivot][pivot] == 0.0)) {
        return false;
      }
      for (auto row = pivot + 1; row < n_terms; ++row) {
        auto const factor = system[row][pivot] / system[pivot][pivot];
        for (auto col = pivot; col <= n_terms; ++col) {
          system[row][col] -= factor * system[pivot][col];
        }
      }
    }
    std::array<double, 3> solution = {{}};
    for (auto row = n_terms; row-- > 0; ) {
      auto value = system[row][n_terms];
      for (auto col = row + 1; col < n_terms; ++col) {
        value -= system[row][col] * solution[col];
      }
      solution[row] = value / system[row][row];
    }
    prediction = solution[0];
    return true;
  }


  /* R: predict(loess(y ~ x, weights = w), x0), with loess' defaults of span
     0.75 and degree 2. With few distinct quality scores (binned qualities of
     recent Illumina instruments) a local quadratic is not identifiable: the
     degree is then lowered until the fit succeeds. */
  auto loess(std::vector<double> const & xs, std::vector<double> const & ys,
             std::vector<double> const & weights, double const x0) -> double {
    auto const n_points = xs.size();
    auto const n_neighbours = std::max(std::size_t{1},
      std::min(n_points, static_cast<std::size_t>(std::floor((static_cast<double>(n_points) * loess_span) + 1e-5))));
    std::vector<double> distances(n_points);
    for (auto i = std::size_t{0}; i < n_points; ++i) {
      distances[i] = std::abs(xs[i] - x0);
    }
    std::nth_element(distances.begin(), std::next(distances.begin(), static_cast<std::ptrdiff_t>(n_neighbours - 1)), distances.end());
    auto const bandwidth = std::max(distances[n_neighbours - 1], 1e-9);
    auto prediction = 0.0;
    for (auto degree = std::size_t{2}; ; --degree) {
      if (fit_local_polynomial(xs, ys, weights, x0, bandwidth * (1.0 + 1e-9), degree, prediction) or (degree == 0)) {
        break;
      }
    }
    return prediction;
  }


  auto estimate_error_model(ErrorCounts const & counts, std::size_t const n_qualities) -> ErrorModel {
    ErrorModel error_model(n_qualities);
    auto const count = [&counts, n_qualities](std::size_t const transition, std::size_t const quality) -> double {
      return static_cast<double>(counts.substitutions[(transition * n_qualities) + quality]);
    };

    /* Indel rates: events per aligned position, with the same pseudocount and
       the same bounds as the substitution rates. No loess here: a deleted
       base has no quality to regress on. */
    auto const indel_rate = [&counts](uint64_t const events) -> double {
      auto const rate = (static_cast<double>(events) + 1.0) / static_cast<double>(std::max(counts.positions, uint64_t{1}));
      return std::min(std::max(rate, min_error_rate), max_error_rate);
    };
    error_model.set_indel_rates(indel_rate(counts.insertions), indel_rate(counts.deletions));

    for (auto from = uint8_t{0}; from < n_nucleotides; ++from) {
      std::vector<double> totals(n_qualities, 0.0);
      for (auto quality = std::size_t{0}; quality < n_qualities; ++quality) {
        for (auto to = uint8_t{0}; to < n_nucleotides; ++to) {
          totals[quality] += count(transition_index(from, to), quality);
        }
      }
      for (auto to = uint8_t{0}; to < n_nucleotides; ++to) {
        if (from == to) {
          continue;
        }
        auto const transition = transition_index(from, to);
        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<double> weights;
        for (auto quality = std::size_t{0}; quality < n_qualities; ++quality) {
          if (totals[quality] > 0.0) {
            xs.push_back(static_cast<double>(quality));
            ys.push_back(std::log10((count(transition, quality) + 1.0) / totals[quality]));  // pseudocount
            weights.push_back(totals[quality]);
          }
        }
        if (xs.empty()) {
          fatal("Error rates could not be estimated (too few reads)");
        }
        for (auto quality = std::size_t{0}; quality < n_qualities; ++quality) {
          /* outside the observed range of qualities: nearest fitted value */
          auto const x0 = std::min(std::max(static_cast<double>(quality), xs.front()), xs.back());
          auto const rate = std::pow(10.0, loess(xs, ys, weights, x0));
          error_model.set_rate(transition, quality, std::min(std::max(rate, min_error_rate), max_error_rate));
        }
      }
      for (auto quality = std::size_t{0}; quality < n_qualities; ++quality) {
        auto self_rate = 1.0;
        for (auto to = uint8_t{0}; to < n_nucleotides; ++to) {
          if (to != from) {
            self_rate -= error_model.rate(transition_index(from, to), quality);
          }
        }
        error_model.set_rate(transition_index(from, from), quality, self_rate);
      }
    }
    return error_model;
  }


  /* learnErrors() / dada(selfConsist = TRUE) in dada2 */
  auto learn_error_model(Denoiser & denoiser, std::size_t const n_qualities,
                         struct Parameters const & parameters) -> ErrorModel {
    /* Round 0: a matrix of ones and a single partition. Every difference with
       the most abundant unique is counted as an error, which gives an upper
       bound on the error rates to start from. */
    ErrorModel error_model(n_qualities);
    auto counts = denoiser.run(error_model, 1);
    error_model = estimate_error_model(counts, n_qualities);
    for (auto nucleotide = uint8_t{0}; nucleotide < n_nucleotides; ++nucleotide) {
      for (auto quality = std::size_t{0}; quality < n_qualities; ++quality) {
        error_model.set_rate(transition_index(nucleotide, nucleotide), quality, 1.0);
      }
    }
    error_model.set_no_indel_rate_to_one();
    /* The model is a deterministic function of the counts, so it repeats
       itself exactly when the counts do. */
    std::vector<ErrorCounts> history;
    for (auto round = int64_t{1}; round <= parameters.opt_denoise_maxconsist; ++round) {
      counts = denoiser.run(error_model, 0);
      if (not parameters.opt_quiet) {
        std::fprintf(stderr, "Self-consistency round %lld: %zu partitions\n",
                     static_cast<long long>(round), denoiser.partitions().size());
      }
      error_model = estimate_error_model(counts, n_qualities);
      if (std::find(history.begin(), history.end(), counts) != history.end()) {
        return error_model;
      }
      history.push_back(counts);
    }
    if (not parameters.opt_quiet) {
      std::fprintf(stderr, "Warning: the error model did not converge in %lld rounds\n",
                   static_cast<long long>(parameters.opt_denoise_maxconsist));
    }
    return error_model;
  }


  /* -------------------------------------------------------------------- */
  /*  input and output                                                    */
  /* -------------------------------------------------------------------- */

  auto transition_name(std::size_t const transition) -> std::string {
    return std::string{nucleotide_symbols[transition / n_nucleotides]} + "2"
      + nucleotide_symbols[transition % n_nucleotides];
  }


  /* Same layout as the matrix returned by dada2::getErrors(): one row per
     transition, one column per quality score. Under --denoise_indels model
     two rows follow, 'ins' and 'del'; the indel rates do not depend on
     quality, so each row repeats one value, which keeps the table
     rectangular (and its first 16 rows usable as a dada2 matrix). */
  auto write_error_model(std::FILE * output_handle, ErrorModel const & error_model,
                         bool const has_indel_rates) -> void {
    std::fprintf(output_handle, "transition");
    for (auto quality = std::size_t{0}; quality < error_model.n_qualities(); ++quality) {
      std::fprintf(output_handle, "\t%zu", quality);
    }
    std::fprintf(output_handle, "\n");
    for (auto transition = std::size_t{0}; transition < n_transitions; ++transition) {
      std::fprintf(output_handle, "%s", transition_name(transition).c_str());
      for (auto quality = std::size_t{0}; quality < error_model.n_qualities(); ++quality) {
        std::fprintf(output_handle, "\t%.6e", error_model.rate(transition, quality));
      }
      std::fprintf(output_handle, "\n");
    }
    if (not has_indel_rates) {
      return;
    }
    for (auto const is_insertion : {true, false}) {
      std::fprintf(output_handle, is_insertion ? "ins" : "del");
      for (auto quality = std::size_t{0}; quality < error_model.n_qualities(); ++quality) {
        std::fprintf(output_handle, "\t%.6e", is_insertion ? error_model.insertion_rate() : error_model.deletion_rate());
      }
      std::fprintf(output_handle, "\n");
    }
  }


  auto read_error_model(char const * filename, std::size_t const n_qualities,
                        bool const needs_indel_rates) -> ErrorModel {
    std::ifstream input(filename);
    if (not input) {
      fatal("Unable to open the error model file given with --denoise_errin");
    }
    ErrorModel error_model(n_qualities);
    std::string line;
    std::getline(input, line);  // column names
    for (auto transition = std::size_t{0}; transition < n_transitions; ++transition) {
      if (not std::getline(input, line)) {
        fatal("The error model file must have one row per transition (16)");
      }
      std::istringstream fields(line);
      std::string name;
      fields >> name;
      if (name != transition_name(transition)) {
        fatal("Unexpected row name in the error model file");
      }
      auto rate = 0.0;
      auto last_rate = 0.0;
      for (auto quality = std::size_t{0}; quality < n_qualities; ++quality) {
        /* as dada2 does: a matrix narrower than the data is extended by
           repeating its last column */
        last_rate = (fields >> rate) ? rate : last_rate;
        if ((last_rate <= 0.0) or (last_rate > 1.0)) {
          fatal("Error rates must be in the interval ]0, 1]");
        }
        error_model.set_rate(transition, quality, last_rate);
      }
    }
    if (not needs_indel_rates) {
      return error_model;
    }
    std::array<double, 2> indel_rates = {{0.0, 0.0}};
    for (auto const * expected_name : {"ins", "del"}) {
      std::string name;
      auto rate = 0.0;
      if (std::getline(input, line)) {
        std::istringstream fields(line);
        fields >> name >> rate;
      }
      if ((name != expected_name) or (rate <= 0.0) or (rate > 1.0)) {
        fatal("--denoise_indels model needs an error model file with 'ins' and 'del' rows, "
              "as written by --denoise_errout under --denoise_indels model");
      }
      indel_rates[(name == "ins") ? 0 : 1] = rate;
    }
    error_model.set_indel_rates(indel_rates[0], indel_rates[1]);
    return error_model;
  }


  auto encode(std::string const & sequence, std::vector<uint8_t> & nucleotides) -> bool {
    nucleotides.resize(sequence.size());
    for (auto i = std::size_t{0}; i < sequence.size(); ++i) {
      switch (sequence[i]) {
      case 'A': nucleotides[i] = 0; break;
      case 'C': nucleotides[i] = 1; break;
      case 'G': nucleotides[i] = 2; break;
      case 'T': nucleotides[i] = 3; break;
      default: return false;  // N and other ambiguous symbols
      }
    }
    return true;
  }


  /* derepFastq() in dada2: full-length dereplication, keeping for each unique
     the mean quality score at each position. 'read_to_unique' remembers, for
     every read of the file, which unique it is, so that the reads can be
     written back in their original order with their original headers. */
  auto dereplicate(struct Parameters const & parameters, std::vector<Unique> & uniques,
                   Database & database, std::vector<uint32_t> & read_to_unique) -> std::size_t {
    auto input_handle = fastq_open(parameters.input_filename, parameters);
    if (input_handle->is_pipe_input()) {
      fatal("--fastq_denoise reads its input twice and cannot read from a pipe");
    }
    vsearch::QualityScoreTable const score_table(parameters);
    std::unordered_map<std::string, uint32_t> index_of;
    std::vector<uint8_t> nucleotides;
    auto max_quality = 0;
    {
      Progress progress("Dereplicating", input_handle->get_size(), parameters);
      while (input_handle->next(HeaderTruncation::keep_whole, Mapping::upcase)) {
        auto const sequence_view = input_handle->sequence_view();
        auto const quality_view = input_handle->quality_view();
        std::string const sequence(sequence_view.data(), sequence_view.size());
        progress.update(input_handle->get_position());
        if ((sequence.size() <= kmer_length) or (not encode(sequence, nucleotides))) {
          read_to_unique.push_back(no_unique);
          continue;
        }
        auto const insertion = index_of.emplace(sequence, static_cast<uint32_t>(uniques.size()));
        if (insertion.second) {
          uniques.emplace_back();
          auto & unique = uniques.back();
          auto const header_view = input_handle->header_view();
          unique.header.assign(header_view.data(), header_view.size());
          unique.sequence = sequence;
          unique.nucleotides = nucleotides;
          unique.quality_sums.assign(sequence.size(), 0.0);
        }
        auto & unique = uniques[insertion.first->second];
        ++unique.reads;
        for (auto i = std::size_t{0}; i < quality_view.size(); ++i) {
          if (not score_table.accepts(quality_view[i])) {
            vsearch::check_quality_score(quality_view[i] - parameters.opt_fastq_ascii,
                                         parameters, input_handle->quality_location());
          }
          auto const score = score_table.score(quality_view[i]);
          unique.quality_sums[i] += score;
          max_quality = std::max(max_quality, score);
        }
        read_to_unique.push_back(insertion.first->second);
      }
    }

    /* sort by decreasing abundance (ties: order of appearance), and remap */
    std::vector<uint32_t> order(uniques.size());
    std::iota(order.begin(), order.end(), uint32_t{0});
    std::stable_sort(order.begin(), order.end(), [&uniques](uint32_t const lhs, uint32_t const rhs) -> bool {
      return uniques[lhs].reads > uniques[rhs].reads;
    });
    std::vector<uint32_t> new_index(uniques.size());
    std::vector<Unique> sorted;
    sorted.reserve(uniques.size());
    for (auto rank = std::size_t{0}; rank < order.size(); ++rank) {
      new_index[order[rank]] = static_cast<uint32_t>(rank);
      sorted.push_back(std::move(uniques[order[rank]]));
    }
    uniques = std::move(sorted);
    for (auto & index : read_to_unique) {
      index = (index == no_unique) ? no_unique : new_index[index];
    }

    for (auto & unique : uniques) {
      unique.qualities.resize(unique.quality_sums.size());
      for (auto i = std::size_t{0}; i < unique.quality_sums.size(); ++i) {
        unique.qualities[i] = static_cast<uint8_t>(std::lround(unique.quality_sums[i] / static_cast<double>(unique.reads)));
      }
      std::vector<double>().swap(unique.quality_sums);  // no longer needed
      index_kmers(unique);
      /* entry number = index of the unique, the order being final by now */
      database.add(false,
                   SeqRecord{View<char>{unique.header.data(), unique.header.size()},
                             View<char>{unique.sequence.data(), unique.sequence.size()},
                             View<char>{}},
                   static_cast<int64_t>(unique.reads));
      std::string().swap(unique.header);  // the database has them now
      std::string().swap(unique.sequence);
    }
    return std::max(min_n_qualities, static_cast<std::size_t>(max_quality) + 1);
  }

}  // end of anonymous namespace


auto fastq_denoise(struct Parameters const & parameters) -> void
{
  if ((parameters.opt_fastqout == nullptr) and (parameters.opt_fastaout == nullptr)
      and (parameters.opt_denoise_errout == nullptr)) {
    fatal("No output files specified");
  }
  if ((parameters.opt_denoise_omega_a < 0.0) or (parameters.opt_denoise_omega_a >= 1.0)) {
    fatal("The argument to --denoise_omega_a must be in the interval [0, 1[");
  }
  if ((parameters.opt_denoise_omega_c < 0.0) or (parameters.opt_denoise_omega_c >= 1.0)) {
    fatal("The argument to --denoise_omega_c must be in the interval [0, 1[");
  }
  if (parameters.opt_denoise_maxconsist < 1) {
    fatal("The argument to --denoise_maxconsist must be at least 1");
  }

  std::vector<Unique> uniques;
  std::vector<uint32_t> read_to_unique;
  Database database;
  database.init();
  auto const n_qualities = dereplicate(parameters, uniques, database, read_to_unique);
  if (uniques.empty()) {
    fatal("No valid sequences found (sequences with N are ignored)");
  }

  Denoiser denoiser(uniques, database, parameters);

  /* step 1: the error model, given or learnt from the reads themselves */
  auto const error_model = (parameters.opt_denoise_errin != nullptr)
    ? read_error_model(parameters.opt_denoise_errin, n_qualities, parameters.opt_denoise_indels_model)
    : learn_error_model(denoiser, n_qualities, parameters);

  /* step 2: the final partition under that model */
  denoiser.run(error_model, 0);
  auto const & partitions = denoiser.partitions();

  /* step 3: output */
  if (parameters.opt_denoise_errout != nullptr) {
    auto const handle = open_optional_output_file(parameters.opt_denoise_errout, OutputOption{"--denoise_errout"});
    write_error_model(handle.get(), error_model, parameters.opt_denoise_indels_model);
  }

  /* each partition's size, counting only the reads that were corrected */
  std::vector<uint64_t> corrected_reads(partitions.size(), 0);
  for (auto i = std::size_t{0}; i < partitions.size(); ++i) {
    for (auto const index : partitions[i].members) {
      corrected_reads[i] += uniques[index].is_corrected ? uniques[index].reads : 0;
    }
  }

  if (parameters.opt_fastaout != nullptr) {
    auto const handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
    for (auto i = std::size_t{0}; i < partitions.size(); ++i) {
      fasta_print_general(handle.get(),
                          database.sequence_view(partitions[i].center),
                          database.header_view(partitions[i].center),
                          OutputAnnotations{corrected_reads[i], static_cast<int64_t>(i + 1)},
                          parameters);
    }
  }

  auto n_corrected = uint64_t{0};
  auto n_unchanged = uint64_t{0};
  auto n_discarded = uint64_t{0};
  if ((parameters.opt_fastqout != nullptr) or (parameters.opt_fastqout_discarded != nullptr)) {
    /* The corrected reads carry the mean quality of the reads of their
       center: their own quality string describes base calls that have just
       been replaced, and may not even have the right length. */
    std::vector<std::string> center_qualities(partitions.size());
    for (auto i = std::size_t{0}; i < partitions.size(); ++i) {
      for (auto const quality : uniques[partitions[i].center].qualities) {
        center_qualities[i].push_back(static_cast<char>(quality + parameters.opt_fastq_ascii));
      }
    }
    auto const kept_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});
    auto const discarded_handle = open_optional_output_file(parameters.opt_fastqout_discarded, OutputOption{"--fastqout_discarded"});
    auto input_handle = fastq_open(parameters.input_filename, parameters);
    Progress progress("Writing corrected reads", input_handle->get_size(), parameters);
    auto ordinal = int64_t{0};
    for (auto nth_read = std::size_t{0}; input_handle->next(HeaderTruncation::keep_whole, Mapping::upcase); ++nth_read) {
      progress.update(input_handle->get_position());
      auto const index = read_to_unique[nth_read];
      if ((index == no_unique) or (not uniques[index].is_corrected)) {
        ++n_discarded;
        if (discarded_handle.get() != nullptr) {
          fastq_print_general(discarded_handle.get(), input_handle->record(),
                              OutputAnnotations{0, static_cast<int64_t>(n_discarded)}, parameters);
        }
        continue;
      }
      auto const nth_partition = uniques[index].partition;
      auto const center_index = partitions[nth_partition].center;
      (index == center_index) ? ++n_unchanged : ++n_corrected;
      if (kept_handle.get() != nullptr) {
        ++ordinal;
        fastq_print_general(kept_handle.get(),
                            database.sequence_view(center_index),
                            input_handle->header_view(),
                            View<char>{center_qualities[nth_partition].data(), center_qualities[nth_partition].size()},
                            OutputAnnotations{0, ordinal},
                            parameters);
      }
    }
  }

  for (auto * stream : {parameters.opt_quiet ? nullptr : stderr, parameters.fp_log}) {
    if (stream == nullptr) {
      continue;
    }
    std::fprintf(stream, "%zu reads, %zu unique sequences, %zu denoised sequences\n",
                 read_to_unique.size(), uniques.size(), partitions.size());
    if ((parameters.opt_fastqout != nullptr) or (parameters.opt_fastqout_discarded != nullptr)) {
      std::fprintf(stream, "%llu reads unchanged, %llu corrected, %llu discarded\n",
                   static_cast<unsigned long long>(n_unchanged),
                   static_cast<unsigned long long>(n_corrected),
                   static_cast<unsigned long long>(n_discarded));
    }
  }
}
