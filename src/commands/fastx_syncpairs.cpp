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
#include <memory>  // std::unique_ptr
#include "core/attributes.hpp"  // struct OutputAnnotations
#include "core/fasta.hpp"  // fasta_print_general
#include "core/fastq.hpp"  // fastq_print_general
#include "core/fastx.hpp"  // fastx_handle
#include "core/seq_record.hpp"  // struct SeqRecord
#include "utils/cityhash.hpp"  // hash_cityhash64
#include "utils/print_view.hpp"  // fprint
#include "utils/progress.hpp"
#include "utils/base_mapping.hpp"  // Mapping::none
#include "utils/fatal.hpp"
#include "utils/open_file.hpp"
#include "utils/view.hpp"  // View<char>
#include <algorithm>  // std::equal, std::find_if, std::max
#include <cassert>  // assert
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE
#include <iterator>  // std::distance, std::next
#include <limits>  // std::numeric_limits
#include <string>
#include <utility>  // std::swap
#include <vector>


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // default set of characters that may introduce the mate number 1 or 2
  // (used when the user does not provide --read_separators). A whitespace
  // is always treated as a separator, regardless of this set.
  char const * const default_read_separators = "/";

  struct output_file {
    char * name = nullptr;
    OutputFileHandle handle;
  };

  // a destination for a category of reads (synced or orphaned), which may
  // be written in fasta and/or fastq format depending on user options
  struct output_pair {
    output_file fasta;
    output_file fastq;
  };

  struct output_files {
    output_pair synced_fwd;
    output_pair synced_rev;
    output_pair orphans_fwd;
    output_pair orphans_rev;
  };

  // a single read kept in memory while indexing the reverse file: where its
  // bytes are in the store's arena, and how they are divided
  struct read_record {
    uint32_t chunk = 0;  // which of the store's chunks holds the bytes
    uint32_t offset = 0;  // where the record starts inside that chunk
    uint32_t header_len = 0;
    uint32_t sequence_len = 0;
    uint32_t quality_len = 0;  // 0 when the input is in fasta format
    uint32_t key_length = 0;  // the matching key is a prefix of the header
    int64_t abundance = 1;
  };


  // Size of one arena chunk (see RecordStore below). Big enough that the
  // per-chunk waste is immaterial, small enough that overshooting the last
  // one costs nothing worth counting.
  //
  // At namespace scope rather than a static constexpr member of RecordStore
  // because std::max() binds it by reference, which odr-uses it: a C++11
  // static constexpr member would then need an out-of-class definition.
  // C++17 refactoring: static constexpr members are implicitly inline, so
  // this can move inside the class.
  constexpr std::size_t arena_chunk_size = 8UL * 1024 * 1024;
  static_assert(arena_chunk_size <= std::numeric_limits<uint32_t>::max(),
                "a record's offset within its chunk is a 32-bit field");

  // The reverse reads kept in memory: one fixed-size descriptor each, and
  // their bytes end to end in a chunked arena, header first.
  //
  // The three std::strings a record used to be were three separate
  // allocations, all past the small-string threshold at typical read lengths
  // (a ~50-byte header, a 200 nt sequence, a 200-byte quality line), so a
  // single record's own three fields sat in three unrelated places and half a
  // million reverse reads cost one and a half million mallocs. Locality
  // matters more here than it usually would, because the forward pass visits
  // the stored records in FORWARD file order, which is unrelated to the order
  // they were stored in: an out-of-order reverse file -- the input this
  // command exists for -- makes every one of those lookups a random access,
  // and a fully shuffled reverse file cost 30% more than an already
  // synchronized one for that reason alone.
  //
  // Chunks, rather than one buffer for the whole file. A single buffer has to
  // be sized in advance and the only figure available is the input file's own
  // size, which is the *compressed* size for a gzip input -- about nine times
  // short on the measured input. It then doubles repeatedly, copying
  // everything stored so far each time, and ends up holding nearly twice the
  // bytes it needs: 448 MB peak against 256 MB for the chunked form, and
  // slower than the three separate strings on exactly the input that matters
  // most. Chunks take the estimate out of the problem: nothing is ever
  // copied, the overshoot is bounded by one chunk instead of by the total,
  // and the peak no longer depends on how well the input happened to
  // compress.
  //
  // A record never straddles a chunk boundary, so a chunk's trailing bytes go
  // unused -- at most one record's worth against 8 MB. Because a record is
  // located by a chunk index rather than by arithmetic over one buffer,
  // chunks may differ in size, which is what lets a record too big for an
  // empty chunk have a chunk of its own.
  class RecordStore {
  public:
    // Store the record the handle is currently on. Its matching key is the
    // first key.size() bytes of its header, and its abundance is whatever the
    // caller decided an output could print (see printable_abundance).
    auto add(fastx_handle handle, bool const is_fastq,
             View<char> const key, int64_t const abundance) -> void {
      auto const stored = handle->record();
      assert(key.size() <= stored.header.size());
      auto const quality_size =
        is_fastq ? stored.quality.size() : std::size_t{0};
      // a std::size_t deliberately: each of the three lengths fits a 32-bit
      // field, but their sum need not, so the span is never narrowed
      auto const span =
        stored.header.size() + stored.sequence.size() + quality_size;

      open_chunk_for(span);
      read_record record;
      record.chunk = to_field(chunks_.size() - 1);
      record.offset = to_field(chunks_.back().size());
      append(stored.header);
      append(stored.sequence);
      if (is_fastq) {
        append(stored.quality);
      }
      room_ -= span;

      record.header_len = to_field(stored.header.size());
      record.sequence_len = to_field(stored.sequence.size());
      record.quality_len = to_field(quality_size);
      record.key_length = to_field(key.size());
      record.abundance = abundance;
      records_.push_back(record);
    }

    auto size() const -> std::size_t { return records_.size(); }

    // the stored record's three fields, as views over its chunk
    auto seq_record(std::size_t const index) const -> SeqRecord {
      auto const & record = records_[at(index)];
      auto const bytes = record_bytes(record);
      return SeqRecord{bytes.first(record.header_len),
                       bytes.subspan(record.header_len, record.sequence_len),
                       bytes.last(record.quality_len),};
    }

    // the bytes a hash match is verified against: the header sits at the
    // start of the record's span, and the matching key is a prefix of it
    auto key(std::size_t const index) const -> View<char> {
      auto const & record = records_[at(index)];
      return record_bytes(record).first(record.key_length);
    }

    auto abundance(std::size_t const index) const -> int64_t {
      return records_[at(index)].abundance;
    }

  private:
    // Narrow a byte count to a descriptor field. Every count narrowed here is
    // already bounded: fastx_filter_header() and
    // fastx_filter_sequence_length() reject a header or a sequence longer
    // than INT_MAX minus the buffer headroom, fatally, at the single point
    // every FASTA/FASTQ read passes through, and a FASTQ quality line is
    // exactly as long as its sequence. An offset is bounded by arena_chunk_size:
    // only the first record of a chunk may exceed that size, and it exhausts
    // its chunk, so every non-zero offset is an offset into an ordinary one.
    static auto to_field(std::size_t const count) -> uint32_t {
      assert(count <= std::numeric_limits<uint32_t>::max());
      return static_cast<uint32_t>(count);
    }

    // an index into records_, checked
    auto at(std::size_t const index) const -> std::size_t {
      assert(index < records_.size());
      return index;
    }

    // Make room for a record of 'span' bytes, opening a chunk when the
    // current one cannot hold the whole of it. room_ starts at zero, so the
    // first record always opens one and append() always has a chunk to write
    // to; the emptiness test is what keeps that true for a zero-byte record.
    auto open_chunk_for(std::size_t const span) -> void {
      if ((not chunks_.empty()) and (span <= room_)) {
        return;
      }
      chunks_.emplace_back();
      room_ = std::max(span, arena_chunk_size);
      chunks_.back().reserve(room_);
    }

    // never reallocates: open_chunk_for() reserved the whole record
    auto append(View<char> const bytes) -> void {
      chunks_.back().insert(chunks_.back().end(), bytes.begin(), bytes.end());
    }

    auto record_bytes(read_record const & record) const -> View<char> {
      assert(record.chunk < chunks_.size());
      // widened before the sum, for the reason given at 'span' above
      auto const span = static_cast<std::size_t>(record.header_len)
        + record.sequence_len + record.quality_len;
      // subspan() asserts that the span is inside the chunk
      return make_view(chunks_[record.chunk]).subspan(record.offset, span);
    }

    std::vector<read_record> records_;
    std::vector<std::vector<char>> chunks_;
    std::size_t room_ = 0;  // bytes still free in the current chunk
  };

  // Positions of the reverse reads, keyed by their matching key: a flat
  // open-addressing table (linear probing, power-of-two capacity, doubled at
  // 50% occupancy) of {key hash, position} slots. The keys themselves are not
  // stored again -- a probe that matches on the 64-bit hash is verified
  // against the key-length prefix of the record's stored header -- so a hash
  // collision costs one extra comparison but can never mispair reads, and
  // each lookup chases at most one record instead of a chain of map nodes,
  // each holding a separately allocated key string.
  class KeyIndex {
  public:
    static auto npos() -> std::size_t {
      return std::numeric_limits<std::size_t>::max();
    }

    // index the key of the record about to be stored at 'position' in
    // 'store'; false when an equal key is already indexed (a duplicate
    // read label)
    auto insert(View<char> const key, std::size_t const position,
                RecordStore const & store) -> bool {
      grow_if_needed();
      auto const hash = hash_cityhash64(key);
      auto slot_number = static_cast<std::size_t>(hash) & mask();
      while (slots_[slot_number].position_plus_one != 0) {
        if (matches(slots_[slot_number], hash, key, store)) {
          return false;
        }
        slot_number = (slot_number + 1) & mask();
      }
      slots_[slot_number].hash = hash;
      slots_[slot_number].position_plus_one = position + 1;
      ++n_entries_;
      return true;
    }

    // position in 'store' of the reverse record with an equal key, or npos()
    auto find(View<char> const key,
              RecordStore const & store) const -> std::size_t {
      auto const hash = hash_cityhash64(key);
      auto slot_number = static_cast<std::size_t>(hash) & mask();
      while (slots_[slot_number].position_plus_one != 0) {
        if (matches(slots_[slot_number], hash, key, store)) {
          return slots_[slot_number].position_plus_one - 1;
        }
        slot_number = (slot_number + 1) & mask();
      }
      return npos();
    }

  private:
    struct Slot {
      uint64_t hash;
      std::size_t position_plus_one;  // 0 marks an empty slot
    };

    static constexpr std::size_t initial_n_slots = 1024;  // any power of two

    auto mask() const -> std::size_t { return slots_.size() - 1; }

    static auto matches(Slot const & slot, uint64_t const hash,
                        View<char> const key,
                        RecordStore const & store) -> bool {
      if (slot.hash != hash) {
        return false;
      }
      auto const stored_key = store.key(slot.position_plus_one - 1);
      return (stored_key.size() == key.size()) and
        std::equal(key.begin(), key.end(), stored_key.begin());
    }

    auto grow_if_needed() -> void {
      if (2 * (n_entries_ + 1) <= slots_.size()) {
        return;
      }
      // stored hashes make rehashing a plain redistribution: the entries are
      // pairwise distinct already, so no equality checks are needed
      std::vector<Slot> old_slots(2 * slots_.size());
      std::swap(old_slots, slots_);
      for (auto const & slot : old_slots) {
        if (slot.position_plus_one == 0) {
          continue;
        }
        auto slot_number = static_cast<std::size_t>(slot.hash) & mask();
        while (slots_[slot_number].position_plus_one != 0) {
          slot_number = (slot_number + 1) & mask();
        }
        slots_[slot_number] = slot;
      }
    }

    std::vector<Slot> slots_ = std::vector<Slot>(initial_n_slots);
    std::size_t n_entries_ = 0;
  };


  auto check_parameters(struct Parameters const & parameters) -> void {
    if (parameters.opt_reverse == nullptr) {
      fatal("No reverse reads file specified with --reverse");
    }

    if ((parameters.opt_fastaout == nullptr) and
        (parameters.opt_fastqout == nullptr) and
        (parameters.opt_fastaout_rev == nullptr) and
        (parameters.opt_fastqout_rev == nullptr) and
        (parameters.opt_fastaout_orphans == nullptr) and
        (parameters.opt_fastqout_orphans == nullptr) and
        (parameters.opt_fastaout_orphans_rev == nullptr) and
        (parameters.opt_fastqout_orphans_rev == nullptr)) {
      fatal("No output files specified");
    }
  }


  auto requests_fastq_output(struct Parameters const & parameters) -> bool {
    return (parameters.opt_fastqout != nullptr) or
      (parameters.opt_fastqout_rev != nullptr) or
      (parameters.opt_fastqout_orphans != nullptr) or
      (parameters.opt_fastqout_orphans_rev != nullptr);
  }


  auto open_output(char * name, char const * option) -> output_file {
    output_file outfile;
    outfile.name = name;
    outfile.handle = open_optional_output_file(name, OutputOption{option});
    return outfile;
  }


  auto open_output_files(struct Parameters const & parameters) -> output_files {
    output_files outfiles;
    outfiles.synced_fwd.fasta = open_output(parameters.opt_fastaout, "--fastaout");
    outfiles.synced_fwd.fastq = open_output(parameters.opt_fastqout, "--fastqout");
    outfiles.synced_rev.fasta = open_output(parameters.opt_fastaout_rev, "--fastaout_rev");
    outfiles.synced_rev.fastq = open_output(parameters.opt_fastqout_rev, "--fastqout_rev");
    outfiles.orphans_fwd.fasta = open_output(parameters.opt_fastaout_orphans, "--fastaout_orphans");
    outfiles.orphans_fwd.fastq = open_output(parameters.opt_fastqout_orphans, "--fastqout_orphans");
    outfiles.orphans_rev.fasta = open_output(parameters.opt_fastaout_orphans_rev, "--fastaout_orphans_rev");
    outfiles.orphans_rev.fastq = open_output(parameters.opt_fastqout_orphans_rev, "--fastqout_orphans_rev");
    return outfiles;
  }


  auto close_output_files(output_files & outfiles) -> void {
    /* called before the stripped-character warnings, so that a deferred write
       error is fatal ahead of them. The order among the eight is not
       significant: outputs naming the same target share one std::FILE (see
       utils/open_file.hpp). */
    for (auto * pair : {& outfiles.synced_fwd, & outfiles.synced_rev,
                        & outfiles.orphans_fwd, & outfiles.orphans_rev,}) {
      pair->fasta.handle.reset();
      pair->fastq.handle.reset();
    }
  }


  // Derive the key shared by the two mates of a pair. Casava 1.8+ headers
  // ("instrument... 1:N:0:" and "... 2:N:0:") already share the substring
  // before the first whitespace, so truncating there is enough. Older
  // headers carry the mate number as a "/1" or "/2" suffix, removed here
  // when its separator belongs to the configured set.
  auto matching_key(View<char> const header,
                    std::string const & separators) -> View<char> {
    auto const * const blank =
      std::find_if(header.begin(), header.end(),
                   [](char const symbol) -> bool {
                     return (symbol == ' ') or (symbol == '\t');
                   });
    auto key = header.first(
      static_cast<std::size_t>(std::distance(header.begin(), blank)));

    if (key.size() >= 2) {
      auto const last = *std::next(key.begin(), static_cast<std::ptrdiff_t>(key.size()) - 1);
      auto const separator = *std::next(key.begin(), static_cast<std::ptrdiff_t>(key.size()) - 2);
      if (((last == '1') or (last == '2')) and
          (separators.find(separator) != std::string::npos)) {
        key = key.first(key.size() - 2);
      }
    }

    return key;
  }


  // The abundance an output of this command could print. Asking the reader
  // for it means scanning the header for a ";size=" annotation
  // (header_find_attribute(), ~375 instructions a record and 5.6% of the run
  // when done for every record of both files), and the value reaches an
  // output through exactly two places, both in fprint_header_annotations()
  // and both gated on --sizeout. Since the CLI rejects --sizeout, --xsize and
  // --relabel* for fastx_syncpairs, the parse currently cannot reach a single
  // byte of any output; gating on the option rather than hard-coding the
  // fallback keeps that true, and honest, if the option table ever changes.
  // 1 is what get_abundance() itself returns when the annotation is absent.
  auto printable_abundance(fastx_handle handle,
                           struct Parameters const & parameters) -> int64_t {
    return parameters.opt_sizeout ? handle->get_abundance() : 1;
  }


  // write a record straight from its views (a reader's record() or a stored
  // record's fields), so a forward read needs no intermediate copy
  auto write_record(output_pair const & destination,
                    SeqRecord const & record,
                    OutputAnnotations const & annotations,
                    struct Parameters const & parameters) -> void {
    if (destination.fastq.handle != nullptr) {
      fastq_print_general(destination.fastq.handle.get(),
                          record, annotations, parameters);
    }
    if (destination.fasta.handle != nullptr) {
      fasta_print_general(destination.fasta.handle.get(), record, annotations, parameters);
    }
  }


  // write a stored reverse read, by position: its bytes live in the store, so
  // a descriptor on its own cannot name them
  auto write_stored(output_pair const & destination,
                    RecordStore const & store,
                    std::size_t const position,
                    int64_t const ordinal,
                    struct Parameters const & parameters) -> void {
    write_record(destination, store.seq_record(position),
                 OutputAnnotations{static_cast<uint64_t>(store.abundance(position)),
                                   ordinal},
                 parameters);
  }


  // Read the reverse file once, keeping every record in memory (in file
  // order, for orphan output) and mapping its matching key to its position.
  auto index_reverse(fastx_handle reverse_handle,
                     bool const is_fastq,
                     std::string const & separators,
                     RecordStore & store,
                     KeyIndex & index,
                     struct Parameters const & parameters) -> void {
    Progress progress("Indexing reverse reads", reverse_handle->get_size(), parameters);
    while (reverse_handle->next(false, Mapping::none)) {
      auto const key = matching_key(reverse_handle->header_view(), separators);
      auto const position = store.size();
      if (not index.insert(key, position, store)) {
        fatal("Duplicate read label in reverse file");
      }
      store.add(reverse_handle, is_fastq, key,
                printable_abundance(reverse_handle, parameters));
      progress.update(reverse_handle->get_position());
    }
  }


  auto stats_message(std::FILE * output_stream,
                     uint64_t const pairs,
                     uint64_t const orphans_fwd,
                     uint64_t const orphans_rev) -> void {
    fprint_integer(output_stream, pairs);
    fprint(output_stream, " pairs synchronized, ");
    fprint_integer(output_stream, orphans_fwd);
    fprint(output_stream, " forward and ");
    fprint_integer(output_stream, orphans_rev);
    fprint(output_stream, " reverse orphan reads\n");
  }

}  // end of anonymous namespace


auto fastx_syncpairs(struct Parameters const & parameters) -> void
{
  /* check parameters */

  check_parameters(parameters);

  /* open and check input files */

  auto forward_handle = fastx_open(parameters.input_filename, parameters);
  auto reverse_handle = fastx_open(parameters.opt_reverse, parameters);

  auto const forward_empty = forward_handle->is_empty_input();
  auto const reverse_empty = reverse_handle->is_empty_input();
  auto const forward_is_fastq = forward_handle->is_fastq_input();
  auto const reverse_is_fastq = reverse_handle->is_fastq_input();

  if ((not forward_empty) and (not reverse_empty) and
      (forward_is_fastq != reverse_is_fastq)) {
    fatal("Forward and reverse files must both be FASTA or both FASTQ");
  }

  /* the effective format follows the non-empty input file */
  auto const is_fastq = forward_empty ? reverse_is_fastq : forward_is_fastq;

  if ((not is_fastq) and requests_fastq_output(parameters)) {
    fatal("Cannot write FASTQ output from FASTA input (no quality scores)");
  }

  /* open output files */

  auto outfiles = open_output_files(parameters);

  std::string const separators =
    (parameters.opt_read_separators != nullptr) ?
    parameters.opt_read_separators : default_read_separators;

  /* index the reverse file (read once, kept in memory) */

  RecordStore reverse_records;
  KeyIndex reverse_index;
  index_reverse(reverse_handle.get(), is_fastq, separators, reverse_records, reverse_index, parameters);

  /* stream the forward file, emitting synced pairs in forward order */

  std::vector<bool> reverse_used(reverse_records.size(), false);
  uint64_t pairs = 0;
  uint64_t orphans_fwd = 0;

  {
    Progress progress("Synchronizing reads", forward_handle->get_size(), parameters);
    while (forward_handle->next(false, Mapping::none)) {
      auto const key = matching_key(forward_handle->header_view(), separators);
      auto const position = reverse_index.find(key, reverse_records);
      if (position == KeyIndex::npos()) {
        write_record(outfiles.orphans_fwd, forward_handle->record(),
                     OutputAnnotations{static_cast<uint64_t>(printable_abundance(forward_handle.get(), parameters)),
                                       static_cast<int64_t>(orphans_fwd + 1)},
                     parameters);
        ++orphans_fwd;
      }
      else {
        // a reverse read already claimed by an earlier forward read means
        // two forward reads share the same matching key: the pairing is
        // ambiguous. Forward orphans that share a key are harmless and are
        // not detected here (they are simply written out twice).
        if (reverse_used[position]) {
          fatal("Duplicate read label in forward file");
        }
        reverse_used[position] = true;
        ++pairs;
        write_record(outfiles.synced_fwd, forward_handle->record(),
                     OutputAnnotations{static_cast<uint64_t>(printable_abundance(forward_handle.get(), parameters)),
                                       static_cast<int64_t>(pairs)},
                     parameters);
        write_stored(outfiles.synced_rev, reverse_records, position,
                     static_cast<int64_t>(pairs), parameters);
      }
      progress.update(forward_handle->get_position());
    }
  }

  /* write the reverse reads that had no forward mate, in reverse order */

  uint64_t orphans_rev = 0;
  for (std::size_t position = 0; position < reverse_records.size(); ++position) {
    if (not reverse_used[position]) {
      write_stored(outfiles.orphans_rev, reverse_records, position,
                   static_cast<int64_t>(orphans_rev + 1), parameters);
      ++orphans_rev;
    }
  }

  /* report */

  if (not parameters.opt_quiet) {
    stats_message(stderr, pairs, orphans_fwd, orphans_rev);
  }
  if (parameters.fp_log != nullptr) {
    stats_message(parameters.fp_log, pairs, orphans_fwd, orphans_rev);
  }

  /* clean up */

  close_output_files(outfiles);
  forward_handle->report_stripped_warning(parameters);
  reverse_handle->report_stripped_warning(parameters);
}
