/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2024, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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
#include "msa.h"
#include <array>
#include <algorithm>  // std::max()
#include <cassert>
#include <cctype>  // std::toupper
#include <cinttypes>  // macro PRId64
#include <climits>  // INT_MAX
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::sscanf
#include <cstdlib>  // str::strtoll
#include <cstring>  // std::memset, std::strlen
#include <iterator> // std::next
#include <numeric> // std::accumulate
#include <vector>


/* Compute multiple sequence alignment (msa), profile, and consensus
   sequence of clustered sequences */

using prof_type = uint64_t;
constexpr auto PROFSIZE = 6;


auto update_profile(char const nucleotide,
                    int const position_in_alignment,
                    prof_type const abundance,
                    std::vector<prof_type>& profile) -> void {
  static constexpr auto A_counter = 0;
  static constexpr auto C_counter = 1;
  static constexpr auto G_counter = 2;
  static constexpr auto U_counter = 3;  // note: T converted to U?
  static constexpr auto N_counter = 4;
  static constexpr auto gap_counter = 5;
  auto const offset = PROFSIZE * position_in_alignment;

  // refactoring: eliminate unused cases
  switch(std::toupper(nucleotide))
    {
    case 'A':
      profile[offset + A_counter] += abundance;
      break;
    case 'C':
      profile[offset + C_counter] += abundance;
      break;
    case 'G':
      profile[offset + G_counter] += abundance;
      break;
    case 'T':
    case 'U':
      profile[offset + U_counter] += abundance;
      break;
    case 'R':
    case 'Y':
    case 'S':
    case 'W':
    case 'K':
    case 'M':
    case 'B':
    case 'D':
    case 'H':
    case 'V':
    case 'N':
      profile[offset + N_counter] += abundance;
      break;
    case '-':
      profile[offset + gap_counter] += abundance;
      break;
    default:
      break;
    }
}


auto update_msa(char const nucleotide, int &position_in_alignment,
                std::vector<char>& alignment) -> void {
  alignment[position_in_alignment] = nucleotide;
  ++position_in_alignment;
}


auto find_runlength_of_leftmost_operation(char * first_character, char ** first_non_digit) -> long long {
  // std::strtoll:
  // - start from the 'first_character' pointed to,
  // - consume as many characters as possible to form a valid integer,
  // - advance pointer to the first non-digit character,
  // - return the valid integer
  // - if there is no valid integer: pointer is not advanced and function returns zero,
  static constexpr auto decimal_base = 10;
  auto runlength = std::strtoll(first_character, first_non_digit, decimal_base);
  assert(runlength <= INT_MAX);

  // in the context of cigar strings, runlength is at least 1
  if (runlength == 0) {
    runlength = 1;
  }
  return runlength;  // is in [1, LLONG_MAX]
}


auto find_max_insertions_per_position(int const target_count,
                                      std::vector<struct msa_target_s> const & target_list_v,
                                      int const centroid_len) -> std::vector<int> {
  std::vector<int> max_insertions(centroid_len + 1);
  for(auto i = 1; i < target_count; ++i)
    {
      char * cigar_start = target_list_v[i].cigar;
      auto const cigar_length = static_cast<long>(std::strlen(cigar_start));
      char * cigar_end = std::next(cigar_start, cigar_length);
      auto * position_in_cigar = cigar_start;
      auto position_in_centroid = 0LL;
      while (position_in_cigar < cigar_end)
        {
          auto** next_operation = &position_in_cigar;  // operations: match (M), insertion (I), or deletion (D)
          auto const runlength = find_runlength_of_leftmost_operation(position_in_cigar, next_operation);
          auto const operation = **next_operation;
          position_in_cigar = std::next(position_in_cigar);
          switch (operation)
            {
            case 'M':
            case 'I':
              position_in_centroid += runlength;
              break;
            case 'D':
              max_insertions[position_in_centroid] = std::max(static_cast<int>(runlength), max_insertions[position_in_centroid]);
              break;
            default:
              break;
            }
        }
  }
  return max_insertions;
}


auto find_total_alignment_length(std::vector<int> const & max_insertions) -> int {
  auto const centroid_len = static_cast<int>(max_insertions.size() - 1);
  return std::accumulate(max_insertions.begin(), max_insertions.end(), centroid_len);
}


auto find_longest_target_on_reverse_strand(int const target_count,
                                           std::vector<struct msa_target_s> const & target_list_v) -> int64_t {
  int64_t longest_reversed = 0;
  for(auto i = 0; i < target_count; ++i)
    {
      auto const & target = target_list_v[i];
      if (target.strand == 0) { continue; }
      auto const len = static_cast<int64_t>(db_getsequencelen(target.seqno));
      longest_reversed = std::max(len, longest_reversed);
    }
  return longest_reversed;
}


auto allocate_buffer_for_reverse_strand_target(int const target_count,
                                               std::vector<struct msa_target_s> const & target_list_v,
                                               std::vector<char> & rc_buffer_v) -> char * {
  /* Find longest target sequence on reverse strand and allocate buffer */
  auto const longest_reversed = find_longest_target_on_reverse_strand(target_count, target_list_v);
  if (longest_reversed > 0)
    {
      rc_buffer_v.resize(longest_reversed + 1);
      return rc_buffer_v.data();
    }
  return nullptr;
}


auto blank_line_before_each_msa(std::FILE * fp_msaout) -> void {
  if (fp_msaout == nullptr) { return ; }
  static_cast<void>(std::fprintf(fp_msaout, "\n"));
}


auto print_header_and_sequence(std::FILE * fp_msaout, char const * header_prefix,
                               int const target_seqno,
                               std::vector<char> & aln_v) -> void {
  // header_prefix == "*" or "", resulting in ">*header" or ">header"
  if (fp_msaout == nullptr) { return ; }

  fasta_print_general(fp_msaout,
                      header_prefix,
                      aln_v.data(),
                      static_cast<int>(aln_v.size() - 1),
                      db_getheader(target_seqno),
                      static_cast<int>(db_getheaderlen(target_seqno)),
                      db_getabundance(target_seqno),
                      0, -1.0, -1, -1, nullptr, 0.0);
}


auto reverse_complement_target_if_need_be(int const strand, int const target_seqno,
                                          char * rc_buffer, char * target_seq) -> char * {
  if (strand == 0) { return target_seq; }
  reverse_complement(rc_buffer, target_seq,
                     static_cast<int64_t>(db_getsequencelen(target_seqno)));
  return rc_buffer;
}


auto process_and_print_centroid(char *rc_buffer,
                                std::vector<struct msa_target_s> const &target_list_v,
                                std::vector<int> const &max_insertions,
                                std::vector<prof_type> &profile,
                                std::vector<char> &aln_v,
                                std::FILE * fp_msaout) -> void {
  auto const centroid_len = static_cast<int>(max_insertions.size() - 1);
  auto const & target = target_list_v.front();
  int const target_seqno = target.seqno;
  char * const target_seq = reverse_complement_target_if_need_be(target.strand, target_seqno, rc_buffer,
                                                                 db_getsequence(target_seqno));
  prof_type const target_abundance = opt_sizein ? db_getabundance(target_seqno) : 1;
  int position_in_alignment = 0;

  for(auto i = 0; i < centroid_len; ++i)
    {
      for(auto j = 0; j < max_insertions[i]; ++j)
        {
          update_profile('-', position_in_alignment, target_abundance, profile);
          update_msa('-', position_in_alignment, aln_v);
        }
      update_profile(*std::next(target_seq, i), position_in_alignment, target_abundance, profile);
      update_msa(*std::next(target_seq, i), position_in_alignment, aln_v);
    }

  // insert
  for(auto j = 0; j < max_insertions[centroid_len]; ++j)
    {
      update_profile('-', position_in_alignment, target_abundance, profile);
      update_msa('-', position_in_alignment, aln_v);
    }

  /* end of sequence string */
  aln_v[position_in_alignment] = '\0';

  /* print header & sequence */
  print_header_and_sequence(fp_msaout, "*", target_seqno, aln_v);
}


auto insert_gaps_in_alignment_and_profile(bool const inserted,
                                          int const max_insertions_at_position,
                                          int & position_in_alignment,
                                          prof_type const target_abundance,
                                          std::vector<prof_type> & profile,
                                          std::vector<char> & aln_v) -> void {
  if (inserted) { return ; }
  for (auto i = 0; i < max_insertions_at_position; ++i) {
    update_profile('-', position_in_alignment, target_abundance, profile);
    update_msa('-', position_in_alignment, aln_v);
  }
}


auto compute_and_print_msa(int const target_count,
                           std::vector<struct msa_target_s> const & target_list_v,
                           std::vector<int> const &max_insertions,
                           std::vector<prof_type> &profile,
                           std::vector<char> &aln_v,
                           std::FILE * fp_msaout) -> void {

  blank_line_before_each_msa(fp_msaout);

  /* Find longest target sequence on reverse strand and allocate buffer */
  std::vector<char> rc_buffer_v;
  char * rc_buffer = allocate_buffer_for_reverse_strand_target(target_count, target_list_v, rc_buffer_v);

  // ------------------------------------------------------- deal with centroid
  process_and_print_centroid(rc_buffer, target_list_v, max_insertions,
                             profile, aln_v, fp_msaout);

  // --------------------------------- deal with other sequences in the cluster
  for(auto i = 1; i < target_count; ++i)
    {
      auto const & target = target_list_v[i];
      auto const target_seqno = target.seqno;
      auto * const target_seq = reverse_complement_target_if_need_be(target.strand, target_seqno,
                                                                     rc_buffer, db_getsequence(target_seqno));
      prof_type const target_abundance = opt_sizein ? db_getabundance(target_seqno) : 1;
      int position_in_alignment = 0;

      auto inserted = false;
      auto qpos = 0;
      auto tpos = 0;

      char * cigar_start = target.cigar;
      auto const cigar_length = static_cast<long>(std::strlen(cigar_start));
      char * cigar_end = std::next(cigar_start, cigar_length);
      auto * position_in_cigar = cigar_start;
      while (position_in_cigar < cigar_end)
        {
          // Consume digits (if any), return the position of the
          // first char (M, D, or I), store it, move cursor to the next byte.
          // Operations: match (M), insertion (I), or deletion (D)
          auto** next_operation = &position_in_cigar;
          auto const runlength = find_runlength_of_leftmost_operation(position_in_cigar, next_operation);
          auto const operation = **next_operation;
          position_in_cigar = std::next(position_in_cigar);

          switch (operation) {
          case 'D':
            for(auto j = 0; j < runlength; ++j)
              {
                update_profile(*std::next(target_seq, tpos), position_in_alignment, target_abundance, profile);
                update_msa(*std::next(target_seq, tpos), position_in_alignment, aln_v);
                ++tpos;
              }
            for(auto j = runlength; j < max_insertions[qpos]; ++j)
              {
                update_profile('-', position_in_alignment, target_abundance, profile);
                update_msa('-', position_in_alignment, aln_v);
              }
            inserted = true;
            break;
          case 'M':
            for(auto j = 0; j < runlength; ++j)
              {
                insert_gaps_in_alignment_and_profile(inserted, max_insertions[qpos],
                                                     position_in_alignment, target_abundance,
                                                     profile, aln_v);
                update_profile(*std::next(target_seq, tpos), position_in_alignment, target_abundance, profile);
                update_msa(*std::next(target_seq, tpos), position_in_alignment, aln_v);
                ++tpos;
                ++qpos;
                inserted = false;
              }
            break;
          case 'I':
            for(auto j = 0; j < runlength; ++j)
              {
                insert_gaps_in_alignment_and_profile(inserted, max_insertions[qpos],
                                                     position_in_alignment, target_abundance,
                                                     profile, aln_v);
                update_profile('-', position_in_alignment, target_abundance, profile);
                update_msa('-', position_in_alignment, aln_v);
                ++qpos;
                inserted = false;
              }
              break;
          default:
            break;
          }
        }

      insert_gaps_in_alignment_and_profile(inserted, max_insertions[qpos],
                                           position_in_alignment, target_abundance,
                                           profile, aln_v);

      /* end of sequence string */
      aln_v[position_in_alignment] = '\0';

      /* print header & sequence */
      print_header_and_sequence(fp_msaout, "", target_seqno, aln_v);
    }
}


auto compute_and_print_consensus(std::vector<int> const &max_insertions,
                                 std::vector<char> &aln_v,
                                 std::vector<char> &cons_v,
                                 std::vector<prof_type> &profile,
                                 std::FILE * fp_msaout) -> void {
  static constexpr char index_of_N = 15;  // 15th char in sym_nt_4bit[] (=> 'N')

  auto const alignment_length = static_cast<int>(aln_v.size() - 1);
  int conslen = 0;

  /* Censor part of the consensus sequence outside the centroid sequence */

  auto const left_censored = max_insertions.front();
  auto const right_censored = max_insertions.back();

  for(auto i = 0; i < alignment_length; ++i)
    {
      if ((i < left_censored) or (i >= alignment_length - right_censored))
        {
          aln_v[i] = '+';
        }
      else
        {
          /* find most common symbol of A, C, G and T */
          char best_sym = 0;
          prof_type best_count = 0;
          for(auto nucleotide = 0U; nucleotide < 4; ++nucleotide)
            {
              auto const count = profile[PROFSIZE * i + nucleotide];
              if (count > best_count)
                {
                  best_count = count;
                  best_sym = static_cast<char>(1U << nucleotide);  // 1, 2, 4, or 8
                }
            }

          /* if no A, C, G, or T, check if there are any N's */
          auto const N_count = profile[PROFSIZE * i + 4];
          if ((best_count == 0) and (N_count > 0))
            {
              best_count = N_count;
              best_sym = index_of_N; // N
            }

          /* compare to the number of gap symbols */
          auto const gap_count = profile[PROFSIZE * i + 5];
          if (best_count >= gap_count)
            {
              auto const index = static_cast<unsigned char>(best_sym);
              auto const sym = sym_nt_4bit[index];  // A, C, G, T, or N
              aln_v[i] = sym;
              cons_v[conslen] = sym;
              ++conslen;
            }
          else
            {
              aln_v[i] = '-';
            }
        }
    }

  aln_v.back() = '\0';
  cons_v[conslen] = '\0';
  cons_v.resize(conslen + 1);

  if (fp_msaout != nullptr)
    {
      fasta_print(fp_msaout, "consensus", aln_v.data(), alignment_length);
    }
}


auto print_consensus_sequence(std::FILE *fp_consout, std::vector<char> & cons_v,
                              int64_t const totalabundance, int const target_count,
                              int const cluster,
                              int const centroid_seqno) -> void {
  if (fp_consout == nullptr) { return ; }
  fasta_print_general(fp_consout,
                      "centroid=",
                      cons_v.data(),
                      static_cast<int>(cons_v.size()),
                      db_getheader(centroid_seqno),
                      static_cast<int>(db_getheaderlen(centroid_seqno)),
                      totalabundance,
                      cluster + 1,
                      -1.0,
                      target_count,
                      opt_clusterout_id ? cluster : -1,
                      nullptr, 0.0);
}


auto print_alignment_profile(std::FILE *fp_profile, std::vector<char> &aln_v,
                             std::vector<prof_type> const &profile,
                             int64_t const totalabundance, int const target_count,
                             int const cluster,
                             int const centroid_seqno) -> void {
  if (fp_profile == nullptr) { return ; }

  // Note: gaps before Ns in profile output
  // 0 = A, 1 = C, 2 = G, 3 = T, 4 = N, 5 = '-' (gap)
  static const std::array<int, 6> symbol_indexes = {0, 1, 2, 3, 5, 4};
  fasta_print_general(fp_profile,
                      "centroid=",
                      nullptr,
                      0,
                      db_getheader(centroid_seqno),
                      static_cast<int>(db_getheaderlen(centroid_seqno)),
                      totalabundance,
                      cluster + 1,
                      -1.0,
                      target_count,
                      opt_clusterout_id ? cluster : -1,
                      nullptr, 0.0);

  aln_v.pop_back(); // remove last element ('\0')
  auto counter = 0;
  for (auto const nucleotide: aln_v) {
    static_cast<void>(std::fprintf(fp_profile, "%d\t%c", counter, nucleotide));
      // A, C, G and T, then gap '-', then N
      for (auto const symbol_index : symbol_indexes) {
        static_cast<void>(std::fprintf(fp_profile, "\t%" PRId64, profile[PROFSIZE * counter + symbol_index]));
      }
      static_cast<void>(std::fprintf(fp_profile, "\n"));
      ++counter;
    }
  static_cast<void>(std::fprintf(fp_profile, "\n"));
}


auto msa(std::FILE * fp_msaout, std::FILE * fp_consout, std::FILE * fp_profile,
         int cluster,
         int const target_count, std::vector<struct msa_target_s> const & target_list_v,
         int64_t totalabundance) -> void
{
  int const centroid_seqno = target_list_v[0].seqno;
  auto const centroid_length = static_cast<int>(db_getsequencelen(centroid_seqno));

  /* find max insertions in front of each position in the centroid sequence */
  auto const max_insertions = find_max_insertions_per_position(target_count, target_list_v, centroid_length);
  auto const alignment_length = find_total_alignment_length(max_insertions);

  /* allocate memory for profile (for consensus) and aligned seq */
  std::vector<prof_type> profile(static_cast<unsigned long>(PROFSIZE) * alignment_length);  // refactoring: std::vector<std::array<prof_type, PROFSIZE>>(alnlen);??
  std::vector<char> aln_v(alignment_length + 1);
  std::vector<char> cons_v(alignment_length + 1);

  /* msaout: multiple sequence alignment ... */
  compute_and_print_msa(target_count, target_list_v, max_insertions,
                        profile, aln_v,
                        fp_msaout);

  /* msaout: ... and consensus sequence at the end */
  compute_and_print_consensus(max_insertions,
                              aln_v,
                              cons_v,
                              profile,
                              fp_msaout);

  /* consout: consensus sequence (dedicated input) */
  print_consensus_sequence(fp_consout, cons_v,
                           totalabundance, target_count,
                           cluster,
                           centroid_seqno);

  /* profile: multiple sequence alignment profile (dedicated input) */
  print_alignment_profile(fp_profile, aln_v,
                          profile,
                          totalabundance, target_count,
                          cluster,
                          centroid_seqno);
}
