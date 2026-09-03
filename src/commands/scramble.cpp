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


#include "commands/scramble.hpp"
#include "utils/span.hpp"
#include "utils/view.hpp"
#include "vsearch.hpp"
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"
#include "core/fastq.hpp"
#include "core/fastx.hpp"
#include "utils/progress.hpp"
#include "utils/fatal.hpp"
#include "utils/grow_to_fit.hpp"  // vsearch::grow_to_fit
#include "utils/maps.hpp"
#include "utils/open_file.hpp"
#include "utils/hash_table_size.hpp"  // vsearch::table_size_half
#include "utils/random.hpp"
#include <algorithm>  // std::copy, std::find
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint32_t, uint64_t
#include <cstdio>  // std::FILE, std::size_t
#include <iterator>  // std::next
#include <limits>  // std::numeric_limits
#include <utility>  // std::swap
#include <vector>


constexpr auto initial_memory_allocation = 512;


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  /* the number of bits in a byte, without the CHAR_BIT macro: for
     unsigned char, digits is exactly that by definition (widened to
     an unsigned type, as it only ever feeds shifts and products) */
  constexpr auto bits_per_byte =
    static_cast<std::size_t>(std::numeric_limits<unsigned char>::digits);

  /* Flat open-addressing vertex table for every (k-1)-mer width (1 to
     8 packed bytes), replacing std::unordered_map: one find-or-insert
     probe sequence per position instead of emplace (which allocates a
     node even for a duplicate key -- one malloc/free per input base,
     60-70% of every k >= 2 run in the 2026-08-31 profile), contiguous
     slots instead of node chains. Slots are epoch-stamped:
     next_record() invalidates the whole table in O(1), so nothing is
     cleared between records, and the stamp doubles as the empty-slot
     marker -- no key value needs to be reserved, so the full 64-bit
     key space of width 8 stays valid. The epoch is 64 bits on
     purpose: a 32-bit counter would let a stale slot alias a live one
     after exactly 2^32 records. next_record() sizes a power-of-two probe
     region from the record's own position count: the load factor
     stays at or below one half, so linear probing terminates and
     probe chains stay short, and no rehash can happen mid-record --
     ids and slots stay put until the next record. The region is
     per-record so a short record after a long one probes a small,
     cache-resident region of the high-water allocation. */
  class FlatVertexTable {
  public:
    /* prepare for one record inserting at most n_keys distinct keys:
       grow the slot vector if needed (fresh slots carry epoch 0,
       stale by construction) and invalidate every slot in O(1) */
    auto next_record(std::size_t const n_keys) -> void {
      std::size_t const region = vsearch::table_size_half(n_keys);
      if (region > slots_.size()) {
        slots_.resize(region);
      }
      mask_ = region - 1;
      ++epoch_;
    }

    /* return the dense vertex id of the (k-1)-mer packed in 'key',
       assigning 'next_id' on its first occurrence this record; a
       result equal to next_id signals a first occurrence */
    /* the bucket index comes from vsearch::splitmix64_mix(): keys are
       raw sequence bytes packed into a word, so without mixing they
       cluster in the low bits and linear probing would degenerate into
       long chains */
    auto find_or_insert(uint64_t const key, uint32_t const next_id) noexcept -> uint32_t {
      std::size_t index = vsearch::splitmix64_mix(key) & mask_;
      while (true) {
        Slot & slot = slots_[index];
        if (slot.epoch != epoch_) {  // empty this record: claim it
          slot.epoch = epoch_;
          slot.key = key;
          slot.vertex_id = next_id;
          return next_id;
        }
        if (slot.key == key) {
          return slot.vertex_id;
        }
        index = (index + 1) & mask_;
      }
    }

  private:
    struct Slot {
      uint64_t epoch = 0;  // 0 = never used; epoch_ is bumped above it before any lookup
      uint64_t key = 0;    // packed (k-1)-mer, meaningful only when epoch matches
      uint32_t vertex_id = 0;
    };
    std::vector<Slot> slots_;
    std::size_t mask_ = 0;  // current probe region size minus one (region is a power of two)
    uint64_t epoch_ = 0;  // the current record's stamp, see next_record()
  };


  /* De Bruijn multigraph of one record, in CSR (compressed sparse row)
     form: one vertex per distinct (k-1)-mer, one directed edge per
     position (the k-mer starting there), edges grouped by source
     vertex. Reading the record left to right traverses every edge
     exactly once -- an Eulerian path from start_vertex() to
     end_vertex() -- and conversely any Eulerian path between those two
     vertices spells a sequence of the same length with exactly the
     same j-mer counts for every j <= k. The class owns the CSR data
     and its consistency only: no randomness, no walk state, so it is
     buildable and checkable deterministically. All buffers are reused
     across rebuild() calls (clear() keeps capacity). */
  class DeBruijnGraph {
  public:
    auto rebuild(View<char> const sequence, int64_t const kmer) -> void {
      assert(kmer >= 2);
      auto const width = static_cast<std::size_t>(kmer - 1);  // (k-1)-mer length
      assert(width <= sizeof(uint64_t));  // ids are (k-1)-mers packed into 8 bytes
      assert(sequence.size() > static_cast<std::size_t>(kmer));  // at least two edges
      auto const n_positions = sequence.size() - width + 1;  // (k-1)-mer occurrences
      if (n_positions > std::numeric_limits<uint32_t>::max()) {
        fatal("--scramble_kmer 2 or more does not support sequences longer than 2^32 nucleotides");
      }

      /* pass 1: map each position's (k-1)-mer to a dense vertex id,
         with a rolling key (shift in one byte, mask); ids are assigned
         in order of first occurrence, so the mapping is deterministic */
      position_vertices_.clear();
      last_chars_.clear();
      map_vertices_hashed(sequence, width);
      assert(position_vertices_.size() == n_positions);
      start_vertex_ = position_vertices_.front();
      end_vertex_ = position_vertices_.back();

      /* pass 2: count out-degrees, prefix-sum into offsets_ */
      auto const n_edges = static_cast<uint32_t>(n_positions - 1);
      offsets_.assign(static_cast<std::size_t>(vertex_count()) + 1, 0);
      for (std::size_t pos = 0; pos + 1 < n_positions; ++pos) {
        ++offsets_[static_cast<std::size_t>(position_vertices_[pos]) + 1];
      }
      for (std::size_t vertex = 1; vertex < offsets_.size(); ++vertex) {
        offsets_[vertex] += offsets_[vertex - 1];
      }
      assert(offsets_.back() == n_edges);

      /* pass 3: scatter successor ids, grouped by source vertex */
      adjacency_.assign(n_edges, 0);
      scatter_cursors_.assign(offsets_.begin(), std::prev(offsets_.end()));
      for (std::size_t pos = 0; pos + 1 < n_positions; ++pos) {
        auto const source = position_vertices_[pos];
        adjacency_[scatter_cursors_[source]] = position_vertices_[pos + 1];
        ++scatter_cursors_[source];
      }
    }

    auto vertex_count() const noexcept -> uint32_t {
      return static_cast<uint32_t>(last_chars_.size());
    }
    auto start_vertex() const noexcept -> uint32_t { return start_vertex_; }
    auto end_vertex() const noexcept -> uint32_t { return end_vertex_; }
    auto last_char(uint32_t const vertex) const noexcept -> char {
      assert(vertex < vertex_count());
      return last_chars_[vertex];
    }
    /* mutable on purpose, not a leak: the class invariant is the edge
       *multiset* per source vertex; the order within a slice is
       explicitly free, so a caller reordering a slice cannot break
       anything the class promises. */
    auto out_edges(uint32_t const vertex) noexcept -> Span<uint32_t> {
      assert(vertex < vertex_count());
      auto const begin_index = offsets_[vertex];
      auto const end_index = offsets_[static_cast<std::size_t>(vertex) + 1];
      return Span<uint32_t>{
        std::next(adjacency_.data(), static_cast<std::ptrdiff_t>(begin_index)),
        static_cast<std::size_t>(end_index - begin_index)};
    }

  private:
    /* the general vertex mapper: any width up to 8 packed bytes,
       through the flat open-addressing table */
    auto map_vertices_hashed(View<char> const sequence, std::size_t const width) -> void {
      vertex_ids_.next_record(sequence.size() - width + 1);
      uint64_t key = 0;
      auto const mask = (width == sizeof(uint64_t))
        ? std::numeric_limits<uint64_t>::max()
        : ((uint64_t{1} << (bits_per_byte * width)) - 1);
      for (std::size_t pos = 0; pos < sequence.size(); ++pos) {
        key = ((key << bits_per_byte)
               | static_cast<uint64_t>(static_cast<unsigned char>(sequence[pos]))) & mask;
        if (pos + 1 < width) { continue; }  // key does not hold a full (k-1)-mer yet
        /* key holds the (k-1)-mer ending at pos */
        auto const next_id = static_cast<uint32_t>(last_chars_.size());
        auto const vertex = vertex_ids_.find_or_insert(key, next_id);
        if (vertex == next_id) {  // first occurrence this record
          last_chars_.push_back(sequence[pos]);
        }
        position_vertices_.push_back(vertex);
      }
    }

    std::vector<uint32_t> offsets_;      // vertex_count() + 1 prefix sums into adjacency_
    std::vector<uint32_t> adjacency_;    // one successor id per edge, grouped by source
    std::vector<char> last_chars_;       // per-vertex last byte of its (k-1)-mer
    std::vector<uint32_t> position_vertices_;  // build scratch: vertex id at each position
    std::vector<uint32_t> scatter_cursors_;    // build scratch for pass 3
    FlatVertexTable vertex_ids_;  // packed (k-1)-mer -> dense id
    uint32_t start_vertex_ = 0;
    uint32_t end_vertex_ = 0;
  };


  /* One Scrambler per run owns the scramble operation. k = 1 is a
     plain portable Fisher-Yates over the record's bytes. For k >= 2 it
     samples a uniformly random Eulerian path of the record's de Bruijn
     multigraph -- uShuffle's method (Jiang et al. 2008, BMC
     Bioinformatics 9:192), generalizing the Altschul-Erickson (1985)
     dinucleotide shuffle: draw a random arborescence oriented toward
     the end vertex, make each vertex's arborescence edge its *last*
     exit, shuffle every vertex's remaining out-edges, then walk from
     the start vertex taking unused edges in slice order. The
     arborescence guarantees the walk consumes every edge, and by the
     BEST theorem the resulting path is uniform among all sequences
     with the input's j-mer counts, j <= k. Every draw goes through
     random_bounded(), so a given --randseed yields the same output on
     any platform. All buffers, the graph included, are reused across
     records. */
  class Scrambler {
  public:
    template <typename URBG>
    auto scramble(Span<char> const sequence, int64_t const kmer, URBG & generator) -> void {
      assert(kmer >= 1);
      if (kmer == 1) {
        random_shuffle(sequence, generator);
        return;
      }
      /* a record holding at most one k-mer admits only one
         arrangement: pass it through unchanged, consuming zero draws
         (uShuffle's convention) */
      if (sequence.size() <= static_cast<std::size_t>(kmer)) {
        return;
      }
      graph_.rebuild(static_cast<View<char>>(sequence), kmer);
      select_last_exits(generator);
      shuffle_out_edges(generator);
      rewrite(sequence, kmer);
    }

  private:
    /* Wilson's algorithm: from each vertex, a loop-erased random walk
       to the growing tree rooted at the end vertex (overwriting
       last_exits_[u] on every revisit is the loop erasure). Picking a
       uniformly random slot of the adjacency slice weights each
       successor by its edge multiplicity, so the arborescence is
       uniform over the *multigraph*'s arborescences -- exactly the
       BEST-theorem weighting that makes the final path uniform. The
       walk terminates with probability 1 because every vertex of the
       record's graph can reach the end vertex (the record's own
       suffix is such a route). */
    template <typename URBG>
    auto select_last_exits(URBG & generator) -> void {
      auto const n_vertices = graph_.vertex_count();
      last_exits_.assign(n_vertices, 0);
      in_tree_.assign(n_vertices, 0);
      in_tree_[graph_.end_vertex()] = 1;
      for (uint32_t vertex = 0; vertex < n_vertices; ++vertex) {
        auto walker = vertex;
        while (in_tree_[walker] == 0) {
          auto const edges = graph_.out_edges(walker);
          assert(not edges.empty());  // only the end vertex may lack out-edges
          auto const pick = static_cast<std::size_t>(random_bounded(generator, edges.size()));
          last_exits_[walker] = edges[pick];
          walker = edges[pick];
        }
        walker = vertex;
        while (in_tree_[walker] == 0) {
          in_tree_[walker] = 1;
          walker = last_exits_[walker];
        }
      }
    }

    /* reserve one edge toward each vertex's arborescence target as its
       last exit (swapped to the back of the slice; parallel edges are
       interchangeable, they spell the same k-mer), then Fisher-Yates
       the rest of the slice; the end vertex reserves nothing */
    template <typename URBG>
    auto shuffle_out_edges(URBG & generator) -> void {
      auto const n_vertices = graph_.vertex_count();
      auto const root = graph_.end_vertex();
      for (uint32_t vertex = 0; vertex < n_vertices; ++vertex) {
        auto edges = graph_.out_edges(vertex);
        if (vertex != root) {
          auto * const reserved = std::find(edges.begin(), edges.end(),
                                            last_exits_[vertex]);
          assert(reserved != edges.end());
          std::swap(*reserved, edges.back());
          edges = edges.first(edges.size() - 1);
        }
        random_shuffle(edges, generator);
      }
    }

    /* walk from the start vertex, consuming each slice front to back:
       the first k-1 bytes spell the start vertex and are already in
       place, and each traversed edge appends its target's last byte */
    auto rewrite(Span<char> const sequence, int64_t const kmer) -> void {
      cursors_.assign(graph_.vertex_count(), 0);
      auto current = graph_.start_vertex();
      for (auto position = static_cast<std::size_t>(kmer) - 1;
           position < sequence.size(); ++position) {
        auto const edges = graph_.out_edges(current);
        assert(cursors_[current] < edges.size());  // the arborescence forbids getting stuck
        auto const target = edges[cursors_[current]];
        ++cursors_[current];
        sequence[position] = graph_.last_char(target);
        current = target;
      }
      assert(current == graph_.end_vertex());
    }

    DeBruijnGraph graph_;
    std::vector<uint32_t> last_exits_;  // per-vertex arborescence target vertex
    std::vector<char> in_tree_;         // Wilson state; not vector<bool>, index-heavy
    std::vector<uint32_t> cursors_;     // walk state: next unused edge per vertex
  };

}  // end of anonymous namespace


auto scramble(struct Parameters const & parameters) -> void
{
  /* deliberately declared outside the record loop: the buffer grows to
     the longest record seen and is reused, so after the first few
     records no allocation happens at all (cppcheck's variableScope
     hint would reallocate per record) */
  std::vector<char> seq_buffer(initial_memory_allocation);

  if ((parameters.opt_fastaout == nullptr) && (parameters.opt_fastqout == nullptr)) {
    fatal("No output files specified");
  }

  auto input_handle = fastx_open(parameters.opt_scramble, parameters);

  if ((parameters.opt_fastqout != nullptr) && not input_handle->is_fastq_input())
    {
      fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");
    }

  auto const filesize = input_handle->get_size();

  /* the RandomSeed carries the full 64-bit --randseed (or an OS value
     when 0); every draw goes through random_bounded(), so a given seed
     yields the same output on any platform (see utils/random.hpp).
     RandomSeed's constructor may throw (std::random_device on the
     --randseed 0 path), so scramble() cannot be noexcept. */
  RandomSeed const seed(parameters);
  Scrambler scrambler;

  {
    /* declared here, ahead of the progress bar, so that leaving this scope
       closes the outputs -- and reports any deferred write error -- after the
       final progress line and before the stripped-character warning below,
       which is the order the explicit close used to give. Which of the two is
       closed first is no longer observable: outputs naming the same target
       share one std::FILE (see utils/open_file.hpp). */
    auto fastaout_handle = open_optional_output_file(parameters.opt_fastaout, OutputOption{"--fastaout"});
    auto fastqout_handle = open_optional_output_file(parameters.opt_fastqout, OutputOption{"--fastqout"});

    int64_t count = 0;  // the ordinal fed to --relabel; int would wrap at 2^31 records
    Progress progress("Scrambling", filesize, parameters);
    while (input_handle->next(false, chrmap_no_change()))
      {
        ++count;

        /* header */

        auto const header = input_handle->header_view();
        auto const abundance = static_cast<uint64_t>(input_handle->get_abundance());


        /* sequence */

        auto const sequence = input_handle->sequence_view();
        auto const length = sequence.size();

        vsearch::grow_to_fit(seq_buffer, length);

        std::copy(sequence.cbegin(), sequence.cend(), seq_buffer.begin());

        /* per-record substream (sintax-style): each record's scramble
           depends only on the base seed, its ordinal, and its length,
           never on how many draws earlier records consumed -- so fasta
           and fastq runs over the same records scramble identically */
        SplitMix64 generator(seed.substream(static_cast<uint64_t>(count)));
        scrambler.scramble(make_span(seq_buffer).first(length),
                           parameters.opt_scramble_kmer, generator);


        /* quality values */

        /* the quality string is passed through untouched: quality is
           never scrambled, so the positional quality profile of each
           record is preserved exactly while the base<->quality pairing
           is deliberately broken (the view is empty for fasta input) */

        if (parameters.opt_fastaout != nullptr)
          {
            fasta_print_general(fastaout_handle.get(),
                                nullptr,
                                make_view(seq_buffer).first(length),
                                header,
                                OutputAnnotations{abundance, count},
                                parameters);
          }

        if (parameters.opt_fastqout != nullptr)
          {
            fastq_print_general(fastqout_handle.get(),
                                make_view(seq_buffer).first(length),
                                header,
                                input_handle->quality_view(),
                                OutputAnnotations{abundance, count},
                                parameters);
          }

        progress.update(input_handle->get_position());
      }
  }

  input_handle->report_stripped_warning(parameters);
}
