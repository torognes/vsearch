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
#include <algorithm>  // std::max()
#include <cctype>  // std::toupper
#include <cinttypes>  // macro PRId64
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::sscanf
#include <cstring>  // std::memset, std::strlen
#include <iterator> // std::next
#include <numeric> // std::accumulate
#include <vector>


/* Compute consensus sequence and msa of clustered sequences */

using prof_type = uint64_t;
constexpr auto PROFSIZE = 6;

static char * aln;
static int alnpos;


auto update_profile(char const nucleotide, prof_type const abundance,
             std::vector<prof_type>& profile) -> void {
  static constexpr auto A_counter = 0;
  static constexpr auto C_counter = 1;
  static constexpr auto G_counter = 2;
  static constexpr auto U_counter = 3;  // note: T converted to U?
  static constexpr auto N_counter = 4;
  static constexpr auto gap_counter = 5;
  auto const offset = PROFSIZE * alnpos;

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


auto update_msa(char const nucleotide, std::vector<char>& alignment) -> void {
  alignment[alnpos] = nucleotide;
  ++alnpos;  // refactoring: pass as copy
}


auto find_max_insertions_per_position(int const target_count,
                                      struct msa_target_s * target_list,
                                      int const centroid_len) -> std::vector<int> {
  std::vector<int> max_insertions(centroid_len + 1);
  for(auto i = 1; i < target_count; ++i)
    {
      char * position = target_list[i].cigar;
      char * end = position + std::strlen(position);
      int pos = 0;  // refactoring: rename?
      while (position < end)
        {
          int64_t run = 1;
          int scanlength = 0;
          std::sscanf(position, "%" PRId64 "%n", &run, &scanlength);
          position += scanlength;
          char const operation = *position;
          ++position;
          switch (operation)
            {
            case 'M':
            case 'I':
              pos += run;
              break;
            case 'D':
              max_insertions[pos] = std::max(static_cast<int>(run), max_insertions[pos]);
              break;
            default:
              break;
            }
        }
  }
  return max_insertions;
}


auto find_total_alignment_length(std::vector<int> const & max_insertions) -> int {
  int const centroid_len = max_insertions.size() - 1;
  return std::accumulate(max_insertions.begin(), max_insertions.end(), centroid_len);
}


auto find_longest_target_on_reverse_strand(int const target_count,
                                           struct msa_target_s * target_list) -> int64_t {
  int64_t longest_reversed = 0;
  for(auto i = 0; i < target_count; ++i)
    {
      if (target_list[i].strand != 0)
        {
          auto const len = static_cast<int64_t>(db_getsequencelen(target_list[i].seqno));
          longest_reversed = std::max(len, longest_reversed);
        }
    }
  return longest_reversed;
}


auto compute_and_print_msa(int const target_count,
                           struct msa_target_s *target_list,
                           std::vector<int> const &max_insertions,
                           std::vector<prof_type> &profile,
                           std::vector<char> &aln_v,
                           std::FILE * fp_msaout) -> void {
  /* Find longest target sequence on reverse strand and allocate buffer */
  auto const longest_reversed = find_longest_target_on_reverse_strand(target_count, target_list);
  char * rc_buffer = nullptr;
  if (longest_reversed > 0)
    {
      std::vector<char> rc_buffer_v(longest_reversed + 1);
      rc_buffer = rc_buffer_v.data();
    }

  int const centroid_len = max_insertions.size() - 1;
  for(auto i = 0; i < target_count; ++i)
    {
      int const target_seqno = target_list[i].seqno;
      char * target_seq = db_getsequence(target_seqno);
      prof_type const target_abundance = opt_sizein ?
        db_getabundance(target_seqno) : 1;

      if (target_list[i].strand != 0)
        {
          reverse_complement(rc_buffer, target_seq,
                             db_getsequencelen(target_seqno));
          target_seq = rc_buffer;
        }

      auto inserted = false;
      int qpos = 0;
      int tpos = 0;
      alnpos = 0;

      if (i == 0)
        {
          for(auto j = 0; j < centroid_len; ++j)
            {
              for(auto k = 0; k < max_insertions[qpos]; ++k)
                {
                  update_profile('-', target_abundance, profile);
                  update_msa('-', aln_v);
                }
              update_profile(target_seq[tpos], target_abundance, profile);
              update_msa(target_seq[tpos], aln_v);
              ++tpos;
              ++qpos;
            }
        }
      else
        {
          char * position = target_list[i].cigar;
          char * end = position + std::strlen(position);
          while (position < end)
            {
              int64_t run = 1;
              int scanlength = 0;
              std::sscanf(position, "%" PRId64 "%n", &run, &scanlength);
              position += scanlength;
              char const operation = *position;
              ++position;

              if (operation == 'D')
                {
                  for(auto j = 0; j < max_insertions[qpos]; ++j)
                    {
                      if (j < run)
                        {
                          update_profile(target_seq[tpos], target_abundance, profile);
                          update_msa(target_seq[tpos], aln_v);
                          ++tpos;
                        }
                      else
                        {
                          update_profile('-', target_abundance, profile);
                          update_msa('-', aln_v);
                        }
                    }
                  inserted = true;
                }
              else
                {
                  for(auto j = 0; j < run; ++j)
                    {
                      if (not inserted)
                        {
                          for(auto k = 0; k < max_insertions[qpos]; ++k)
                            {
                              update_profile('-', target_abundance, profile);
                              update_msa('-', aln_v);
                            }
                        }

                      if (operation == 'M')
                        {
                          update_profile(target_seq[tpos], target_abundance, profile);
                          update_msa(target_seq[tpos], aln_v);
                          ++tpos;
                        }
                      else
                        {
                          update_profile('-', target_abundance, profile);
                          update_msa('-', aln_v);
                        }

                      ++qpos;
                      inserted = false;
                    }
                }
            }
        }

      if (not inserted)
        {
          for(auto j = 0; j < max_insertions[qpos]; ++j)
            {
              update_profile('-', target_abundance, profile);
              update_msa('-', aln_v);
            }
        }

      /* end of sequence string */
      aln_v[alnpos] = 0;

      /* print header & sequence */
      if (fp_msaout != nullptr)
        {
          fasta_print_general(fp_msaout,
                              i != 0 ? "" : "*",
                              aln_v.data(),
                              aln_v.size() - 1,
                              db_getheader(target_seqno),
                              static_cast<int>(db_getheaderlen(target_seqno)),
                              db_getabundance(target_seqno),
                              0, -1.0, -1, -1, nullptr, 0.0);
        }
    }
}


auto compute_and_print_consensus(std::vector<int> const &max_insertions,
                                 std::vector<char> &aln_v,
                                 std::vector<char> &cons_v,
                                 std::vector<prof_type> &profile,
                                 std::FILE * fp_msaout) -> void {
  int const alignment_length = aln_v.size() - 1;
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
          for(auto nucleotide = 0; nucleotide < 4; ++nucleotide)
            {
              auto const count = profile[PROFSIZE * i + nucleotide];
              if (count > best_count)
                {
                  best_count = count;
                  best_sym = 1 << nucleotide;
                }
            }

          /* if no A, C, G, or T, check if there are any N's */
          auto const n_count = profile[PROFSIZE * i + 4];
          if ((best_count == 0) and (n_count > 0))
            {
              best_count = n_count;
              best_sym = 15; // N
            }

          /* compare to the number of gap symbols */
          auto const gap_count = profile[PROFSIZE * i + 5];
          if (best_count >= gap_count)
            {
              auto const sym = sym_nt_4bit[static_cast<unsigned char>(best_sym)];
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


auto msa(std::FILE * fp_msaout, std::FILE * fp_consout, std::FILE * fp_profile,
         int cluster,
         int const target_count, struct msa_target_s * target_list,
         int64_t totalabundance) -> void
{
  int const centroid_seqno = target_list[0].seqno;
  auto const centroid_len = static_cast<int>(db_getsequencelen(centroid_seqno));

  /* find max insertions in front of each position in the centroid sequence */
  auto const max_insertions = find_max_insertions_per_position(target_count, target_list, centroid_len);
  auto const alnlen = find_total_alignment_length(max_insertions);

  /* allocate memory for profile (for consensus) and aligned seq */
  std::vector<prof_type> profile(static_cast<unsigned long>(PROFSIZE) * alnlen);  // refactoring: std::vector<std::array<prof_type, PROFSIZE>>(alnlen);??
  std::vector<char> aln_v(alnlen + 1);
  aln = aln_v.data();
  std::vector<char> cons_v(alnlen + 1);

  /* blank line before each msa */
  if (fp_msaout != nullptr)
    {
      fprintf(fp_msaout, "\n");
    }

  /* msaout: multiple sequence alignment ... */
  compute_and_print_msa(target_count, target_list, max_insertions,
                        profile, aln_v,
                        fp_msaout);

  /* msaout: ... and consensus sequence at the end */
  compute_and_print_consensus(max_insertions,
                              aln_v,
                              cons_v,
                              profile,
                              fp_msaout);

  /* consout: consensus sequence (dedicated input) */
  if (fp_consout != nullptr)
    {
      fasta_print_general(fp_consout,
                          "centroid=",
                          cons_v.data(),
                          cons_v.size(),
                          db_getheader(centroid_seqno),
                          db_getheaderlen(centroid_seqno),
                          totalabundance,
                          cluster + 1,
                          -1.0,
                          target_count,
                          opt_clusterout_id ? cluster : -1,
                          nullptr, 0.0);
    }

  /* profile: multiple sequence alignment profile (dedicated input) */
  if (fp_profile != nullptr)
    {
      // Note: gaps before Ns in profile output
      // 0 = A, 1 = C, 2 = G, 3 = T, 4 = N, 5 = '-' (gap)
      static const std::array<int, 6> symbol_indexes = {0, 1, 2, 3, 5, 4};
      fasta_print_general(fp_profile,
                          "centroid=",
                          nullptr,
                          0,
                          db_getheader(centroid_seqno),
                          db_getheaderlen(centroid_seqno),
                          totalabundance,
                          cluster + 1,
                          -1.0,
                          target_count,
                          opt_clusterout_id ? cluster : -1,
                          nullptr, 0.0);

      for (auto i = 0; i < alnlen; ++i)
        {
          fprintf(fp_profile, "%d\t%c", i, aln_v[i]);
          // A, C, G and T, then gap '-', then N
          for (auto const symbol_index : symbol_indexes) {
            fprintf(fp_profile, "\t%" PRId64, profile[PROFSIZE * i + symbol_index]);
          }
          fprintf(fp_profile, "\n");
        }
      fprintf(fp_profile, "\n");
    }
}
