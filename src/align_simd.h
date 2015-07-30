/*
    Copyright (C) 2014 Torbjorn Rognes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

typedef signed short CELL;
typedef unsigned short WORD;
typedef unsigned char BYTE;

struct s16info_s;

extern unsigned int * chrmap; // TODO should be exchanged with a global field that is set in the beginning

void backtrack16(s16info_s * s,
                 char * dseq,
                 unsigned long dlen,
                 unsigned long offset,
                 unsigned long channel,
                 unsigned short * paligned,
                 unsigned short * pmatches,
                 unsigned short * pmismatches,
                 unsigned short * pgaps);

/*
 * Initialises the alignment process.
 *
 * Requires the score matrix (class ScoreMatrix) to already be intialised.
 */
struct s16info_s *
search16_init_2(CELL penalty_gap_open_query_left,
              CELL penalty_gap_open_target_left,
              CELL penalty_gap_open_query_interior,
              CELL penalty_gap_open_target_interior,
              CELL penalty_gap_open_query_right,
              CELL penalty_gap_open_target_right,
              CELL penalty_gap_extension_query_left,
              CELL penalty_gap_extension_target_left,
              CELL penalty_gap_extension_query_interior,
              CELL penalty_gap_extension_target_interior,
              CELL penalty_gap_extension_query_right,
              CELL penalty_gap_extension_target_right);

/*
 * Deprecated function kept for compatibility reasons. Should be
 * exchanged for search16_init_2.
 *
 * This function additionally initialises the score matrix, while
 * this should be done at a different point in the program.
 */
struct s16info_s *
search16_init(int score_match,
              int score_mismatch,
              CELL penalty_gap_open_query_left,
              CELL penalty_gap_open_target_left,
              CELL penalty_gap_open_query_interior,
              CELL penalty_gap_open_target_interior,
              CELL penalty_gap_open_query_right,
              CELL penalty_gap_open_target_right,
              CELL penalty_gap_extension_query_left,
              CELL penalty_gap_extension_target_left,
              CELL penalty_gap_extension_query_interior,
              CELL penalty_gap_extension_target_interior,
              CELL penalty_gap_extension_query_right,
              CELL penalty_gap_extension_target_right);

void
search16_exit(s16info_s * s);

void
search16_qprep(s16info_s * s, char * qseq, int qlen);

void
search16(s16info_s * s,
         unsigned int sequences,
         unsigned int * seqnos,
         CELL * pscores,
         unsigned short * paligned,
         unsigned short * pmatches,
         unsigned short * pmismatches,
         unsigned short * pgaps,
         char * * pcigar);
