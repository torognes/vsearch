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


auto msa_add(char const nucleotide, prof_type const abundance,
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

  *std::next(aln, alnpos) = nucleotide;  // refactoring: extract from function
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

  /* Find longest target sequence on reverse strand and allocate buffer */
  auto const longest_reversed = find_longest_target_on_reverse_strand(target_count, target_list);
  char * rc_buffer = nullptr;
  if (longest_reversed > 0)
    {
      std::vector<char> rc_buffer_v(longest_reversed + 1);
      rc_buffer = rc_buffer_v.data();
    }

  /* blank line before each msa */
  if (fp_msaout != nullptr)
    {
      fprintf(fp_msaout, "\n");
    }

  for(auto j = 0; j < target_count; ++j)
    {
      int const target_seqno = target_list[j].seqno;
      char * target_seq = db_getsequence(target_seqno);
      prof_type const target_abundance = opt_sizein ?
        db_getabundance(target_seqno) : 1;

      if (target_list[j].strand != 0)
        {
          reverse_complement(rc_buffer, target_seq,
                             db_getsequencelen(target_seqno));
          target_seq = rc_buffer;
        }

      auto inserted = false;
      int qpos = 0;
      int tpos = 0;
      alnpos = 0;

      if (j == 0)
        {
          for(auto x = 0; x < centroid_len; ++x)
            {
              for(auto y = 0; y < max_insertions[qpos]; ++y)
                {
                  msa_add('-', target_abundance, profile);
                }
              msa_add(target_seq[tpos], target_abundance, profile);
              ++tpos;
              ++qpos;
            }
        }
      else
        {
          char * p = target_list[j].cigar;
          char * end = p + std::strlen(p);
          while (p < end)
            {
              int64_t run = 1;
              int scanlength = 0;
              std::sscanf(p, "%" PRId64 "%n", &run, &scanlength);
              p += scanlength;
              char const operation = *p;
              ++p;

              if (operation == 'D')
                {
                  for(auto x = 0; x < max_insertions[qpos]; ++x)
                    {
                      if (x < run)
                        {
                          msa_add(target_seq[tpos], target_abundance, profile);
                          ++tpos;
                        }
                      else
                        {
                          msa_add('-', target_abundance, profile);
                        }
                    }
                  inserted = true;
                }
              else
                {
                  for(auto x = 0; x < run; ++x)
                    {
                      if (not inserted)
                        {
                          for(auto y = 0; y < max_insertions[qpos]; ++y)
                            {
                              msa_add('-', target_abundance, profile);
                            }
                        }

                      if (operation == 'M')
                        {
                          msa_add(target_seq[tpos], target_abundance, profile);
                          ++tpos;
                        }
                      else
                        {
                          msa_add('-', target_abundance, profile);
                        }

                      ++qpos;
                      inserted = false;
                    }
                }
            }
        }

      if (not inserted)
        {
          for(auto x = 0; x < max_insertions[qpos]; ++x)
            {
              msa_add('-', target_abundance, profile);
            }
        }

      /* end of sequence string */
      aln_v[alnpos] = 0;

      /* print header & sequence */
      if (fp_msaout != nullptr)
        {
          fasta_print_general(fp_msaout,
                              j != 0 ? "" : "*",
                              aln_v.data(),
                              alnlen,
                              db_getheader(target_seqno),
                              static_cast<int>(db_getheaderlen(target_seqno)),
                              db_getabundance(target_seqno),
                              0, -1.0, -1, -1, nullptr, 0.0);
        }
    }


  /* consensus */

  int conslen = 0;

  /* Censor part of the consensus sequence outside the centroid sequence */

  auto const left_censored = max_insertions.front();
  auto const right_censored = max_insertions.back();

  for(auto i = 0; i < alnlen; ++i)
    {
      if ((i < left_censored) or (i >= alnlen - right_censored))
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

  if (fp_msaout != nullptr)
    {
      fasta_print(fp_msaout, "consensus", aln_v.data(), alnlen);
    }

  if (fp_consout != nullptr)
    {
      fasta_print_general(fp_consout,
                          "centroid=",
                          cons_v.data(),
                          conslen,
                          db_getheader(centroid_seqno),
                          db_getheaderlen(centroid_seqno),
                          totalabundance,
                          cluster+1,
                          -1.0,
                          target_count,
                          opt_clusterout_id ? cluster : -1,
                          nullptr, 0.0);
    }

  if (fp_profile != nullptr)
    {
      fasta_print_general(fp_profile,
                          "centroid=",
                          nullptr,
                          0,
                          db_getheader(centroid_seqno),
                          db_getheaderlen(centroid_seqno),
                          totalabundance,
                          cluster+1,
                          -1.0,
                          target_count,
                          opt_clusterout_id ? cluster : -1,
                          nullptr, 0.0);

      for (auto i = 0; i < alnlen; ++i)
        {
          fprintf(fp_profile, "%d\t%c", i, aln_v[i]);
          // A, C, G and T
          for (auto c = 0; c < 4; ++c)
            {
              fprintf(fp_profile, "\t%" PRId64, profile[PROFSIZE * i + c]);
            }
          // Gap symbol
          fprintf(fp_profile, "\t%" PRId64, profile[PROFSIZE * i + 5]);
          // Ambiguous nucleotide (Ns and others)
          fprintf(fp_profile, "\t%" PRId64, profile[PROFSIZE * i + 4]);
          fprintf(fp_profile, "\n");
        }
      fprintf(fp_profile, "\n");
    }
}
