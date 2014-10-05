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

struct nwinfo_s;

struct nwinfo_s * nw_init();

void nw_exit(struct nwinfo_s * nw);

void nw_align(char * dseq,
              char * dend,
              char * qseq,
              char * qend,
              long * score_matrix,
              long gapopen_q_left,
              long gapopen_q_internal,
              long gapopen_q_right,
              long gapopen_t_left,
              long gapopen_t_internal,
              long gapopen_t_right,
              long gapextend_q_left,
              long gapextend_q_internal,
              long gapextend_q_right,
              long gapextend_t_left,
              long gapextend_t_internal,
              long gapextend_t_right,
              long * nwscore,
              long * nwdiff,
              long * nwgaps,
              long * nwindels,
              long * nwalignmentlength,
              char ** nwalignment,
              long queryno,
              long dbseqno,
              struct nwinfo_s * nw);
