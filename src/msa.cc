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
#include <cstddef>  // std::ptrdiff_t
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
static prof_type * profile;


auto msa_add(char const nucleotide, prof_type const abundance) -> void
{
  static constexpr auto A_counter = 0;
  static constexpr auto C_counter = 1;
  static constexpr auto G_counter = 2;
  static constexpr auto U_counter = 3;  // note: T converted to U?
  static constexpr auto N_counter = 4;
  static constexpr auto gap_counter = 5;
  auto * const position_profile = std::next(profile, static_cast<std::ptrdiff_t>(PROFSIZE) * alnpos);

  switch(std::toupper(nucleotide))
    {
    case 'A':
      *std::next(position_profile, A_counter) += abundance;
      break;
    case 'C':
      *std::next(position_profile, C_counter) += abundance;
      break;
    case 'G':
      *std::next(position_profile, G_counter) += abundance;
      break;
    case 'T':
    case 'U':
      *std::next(position_profile, U_counter) += abundance;
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
      *std::next(position_profile, N_counter) += abundance;
      break;
    case '-':
      *std::next(position_profile, gap_counter) += abundance;
      break;
    default:
      break;
    }

  *std::next(aln, alnpos) = nucleotide;
  ++alnpos;
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
              if (run > max_insertions[pos])
                {
                  max_insertions[pos] = run;
                }
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


auto msa(std::FILE * fp_msaout, std::FILE * fp_consout, std::FILE * fp_profile,
         int cluster,
         int const target_count, struct msa_target_s * target_list,
         int64_t totalabundance) -> void
{
  int const centroid_seqno = target_list[0].seqno;
  int const centroid_len = db_getsequencelen(centroid_seqno);

  /* find max insertions in front of each position in the centroid sequence */
  auto const max_insertions = find_max_insertions_per_position(target_count, target_list, centroid_len);
  auto const alnlen = find_total_alignment_length(max_insertions);

  /* allocate memory for profile (for consensus) and aligned seq */
  profile = static_cast<prof_type *>(xmalloc(PROFSIZE * sizeof(prof_type) * alnlen));
  for (auto i = 0; i < PROFSIZE * alnlen; ++i)
    {
      profile[i] = 0;
    }
  aln = static_cast<char *>(xmalloc(alnlen + 1));
  char * cons = static_cast<char *>(xmalloc(alnlen + 1));

  /* Find longest target sequence on reverse strand and allocate buffer */
  int64_t longest_reversed = 0;
  for(auto i = 0; i < target_count; ++i)
    {
      if (target_list[i].strand != 0)
        {
          int64_t const len = db_getsequencelen(target_list[i].seqno);
          longest_reversed = std::max(len, longest_reversed);
        }
    }
  char * rc_buffer = nullptr;
  if (longest_reversed > 0)
    {
      rc_buffer = static_cast<char *>(xmalloc(longest_reversed + 1));
    }

  /* blank line before each msa */
  if (fp_msaout != nullptr)
    {
      fprintf(fp_msaout, "\n");
    }

  for(auto j = 0; j < target_count; ++j)
    {
      int target_seqno = target_list[j].seqno;
      char * target_seq = db_getsequence(target_seqno);
      prof_type target_abundance = opt_sizein ?
        db_getabundance(target_seqno) : 1;

      if (target_list[j].strand)
        {
          reverse_complement(rc_buffer, target_seq,
                             db_getsequencelen(target_seqno));
          target_seq = rc_buffer;
        }

      int inserted = 0;
      int qpos = 0;
      int tpos = 0;
      alnpos = 0;

      if (not j)
        {
          for(int x = 0; x < centroid_len; ++x)
            {
              for(int y = 0; y < max_insertions[qpos]; ++y)
                {
                  msa_add('-', target_abundance);
                }
              msa_add(target_seq[tpos++], target_abundance);
              ++qpos;
            }
        }
      else
        {
          char * p = target_list[j].cigar;
          char * e = p + std::strlen(p);
          while (p < e)
            {
              int64_t run = 1;
              int scanlength = 0;
              std::sscanf(p, "%" PRId64 "%n", &run, &scanlength);
              p += scanlength;
              char const op = *p;
              ++p;

              if (op == 'D')
                {
                  for(int x = 0; x < max_insertions[qpos]; ++x)
                    {
                      if (x < run)
                        {
                          msa_add(target_seq[tpos++], target_abundance);
                        }
                      else
                        {
                          msa_add('-', target_abundance);
                        }
                    }
                  inserted = 1;
                }
              else
                {
                  for(int x = 0; x < run; ++x)
                    {
                      if (not inserted)
                        {
                          for(int y = 0; y < max_insertions[qpos]; ++y)
                            {
                              msa_add('-', target_abundance);
                            }
                        }

                      if (op == 'M')
                        {
                          msa_add(target_seq[tpos++], target_abundance);
                        }
                      else
                        {
                          msa_add('-', target_abundance);
                        }

                      ++qpos;
                      inserted = 0;
                    }
                }
            }
        }

      if (not inserted)
        {
          for(int x = 0; x < max_insertions[qpos]; ++x)
            {
              msa_add('-', target_abundance);
            }
        }

      /* end of sequence string */
      aln[alnpos] = 0;

      /* print header & sequence */
      if (fp_msaout)
        {
          fasta_print_general(fp_msaout,
                              j ? "" : "*",
                              aln,
                              alnlen,
                              db_getheader(target_seqno),
                              db_getheaderlen(target_seqno),
                              db_getabundance(target_seqno),
                              0, -1.0, -1, -1, nullptr, 0.0);
        }
    }

  if (rc_buffer)
    {
      xfree(rc_buffer);
    }

  /* consensus */

  int conslen = 0;

  /* Censor part of the consensus sequence outside the centroid sequence */

  auto const left_censored = max_insertions.front();
  auto const right_censored = max_insertions.back();

  for(int i = 0; i < alnlen; i++)
    {
      if ((i < left_censored) || (i >= alnlen - right_censored))
        {
          aln[i] = '+';
        }
      else
        {
          /* find most common symbol of A, C, G and T */
          char best_sym = 0;
          prof_type best_count = 0;
          for(int c = 0; c < 4; c++)
            {
              prof_type count = profile[PROFSIZE * i + c];
              if (count > best_count)
                {
                  best_count = count;
                  best_sym = 1 << c;
                }
            }

          /* if no A, C, G, or T, check if there are any N's */
          prof_type n_count = profile[PROFSIZE * i + 4];
          if ((best_count == 0) && (n_count > 0))
            {
              best_count = n_count;
              best_sym = 15; // N
            }

          /* compare to the number of gap symbols */
          prof_type gap_count = profile[PROFSIZE * i + 5];
          if (best_count >= gap_count)
            {
              char sym = sym_nt_4bit[(int)best_sym];
              aln[i] = sym;
              cons[conslen++] = sym;
            }
          else
            {
              aln[i] = '-';
            }
        }
    }

  aln[alnlen] = 0;
  cons[conslen] = 0;

  if (fp_msaout)
    {
      fasta_print(fp_msaout, "consensus", aln, alnlen);
    }

  if (fp_consout)
    {
      fasta_print_general(fp_consout,
                          "centroid=",
                          cons,
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

  if (fp_profile)
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

      for (int i = 0; i < alnlen; i++)
        {
          fprintf(fp_profile, "%d\t%c", i, aln[i]);
          // A, C, G and T
          for (int c = 0; c < 4; c++)
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

  xfree(aln);
  xfree(cons);
  xfree(profile);
}
