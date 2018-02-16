
/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2017, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

void rereplicate()
{
  opt_xsize = 1;

  FILE * fp_output = 0;

  if (opt_output)
    {
      fp_output = fopen_output(opt_output);
      if (!fp_output)
        fatal("Unable to open FASTA output file for writing");
    }

  fastx_handle fh = fasta_open(opt_rereplicate);
  int64_t filesize = fasta_get_size(fh);

  progress_init("Rereplicating", filesize);

  int64_t i = 0;
  int64_t n = 0;
  while (fasta_next(fh, ! opt_notrunclabels, chrmap_no_change))
    {
      n++;
      int64_t abundance = fasta_get_abundance(fh);
      
      for(int64_t j=0; j<abundance; j++)
        {
          i++;
          if (opt_output)
            fasta_print_general(fp_output,
                                0,
                                fasta_get_sequence(fh),
                                fasta_get_sequence_length(fh),
                                fasta_get_header(fh),
                                fasta_get_header_length(fh),
                                1,
                                i,
                                -1, -1, 0, 0.0);
        }

      progress_update(fasta_get_position(fh));
    }
  progress_done();
  
  fprintf(stderr, "Rereplicated %" PRId64 " reads from %" PRId64 " amplicons\n", i, n);

  fasta_close(fh);

  fclose(fp_output);
}
