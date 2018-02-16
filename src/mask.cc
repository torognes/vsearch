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

static int dust_word = 3;
static int dust_level = 20;
static int dust_window = 64;
static int dust_window2 = dust_window >> 1;
static int word_count = 1 << (dust_word << 1);
static int bitmask = word_count - 1;

int wo(int len, const char *s, int *beg, int *end)
{
  int l1 = len - dust_word + 1 - 5; /* smallest possible region is 8 */
  if (l1 < 0)
    return 0;
  
  int bestv = 0;
  int besti = 0;
  int bestj = 0;
  int counts[word_count];
  int words[dust_window];
  int word = 0;

  for (int j = 0; j < len; j++)
    {
      word <<= 2;
      word |= chrmap_2bit[(int)(s[j])];
      words[j] = word & bitmask;
    }

  for (int i=0; i < l1; i++)
    {
      memset(counts, 0, sizeof(counts));
      
      int sum = 0;
      
      for (int j = dust_word-1; j<len-i; j++)
        {
          word = words[i+j];
          int c = counts[word];
          if (c)
            {
              sum += c;
              int v = 10 * sum / j;
              
              if (v > bestv)
                {
                  bestv = v;
                  besti = i;
                  bestj = j;
                }
            }
          counts[word]++;
        }
    }
  
  *beg = besti;
  *end = besti + bestj;
  
  return bestv;
}

void dust(char * m, int len)
{
  int a, b;

  /* make a local copy of the original sequence */
  char * s = (char*) xmalloc(len+1);
  strcpy(s, m);

  if (! opt_hardmask)
    {
      /* convert sequence to upper case unless hardmask in effect */
      for(int i=0; i < len; i++)
        m[i] = toupper(m[i]);
      m[len] = 0;
    }

  for (int i=0; i < len; i += dust_window2)
    {
      int l = (len > i + dust_window) ? dust_window : len-i;
      int v = wo(l, s+i, &a, &b);
      
      if (v > dust_level)
        {
          if (opt_hardmask)
            for(int j=a+i; j<=b+i; j++)
              m[j] = 'N';
          else
            for(int j=a+i; j<=b+i; j++)
              m[j] = s[j] | 0x20;
          
          if (b < dust_window2)
            i += dust_window2 - b;
        }
    }

  xfree(s);
}

static pthread_t * pthread;
static pthread_attr_t attr;
static pthread_mutex_t mutex;
static int nextseq = 0;
static int seqcount = 0;

void * dust_all_worker(void * vp)
{
  while(1)
    {
      pthread_mutex_lock(&mutex);
      int seqno = nextseq;
      if (seqno < seqcount)
        {
          nextseq++;
          progress_update(seqno);
          pthread_mutex_unlock(&mutex);
          dust(db_getsequence(seqno), db_getsequencelen(seqno));
        }
      else
        {
          pthread_mutex_unlock(&mutex);
          break;
        }
    }
  return 0;
}

void dust_all()
{
  nextseq = 0;
  seqcount = db_getsequencecount();
  progress_init("Masking", seqcount);

  pthread_mutex_init(&mutex, NULL);

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  pthread = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));

  for(int t=0; t<opt_threads; t++)
    if (pthread_create(pthread+t, &attr, dust_all_worker, (void*)(int64_t)t))
      fatal("Cannot create thread");

  for(int t=0; t<opt_threads; t++)
    if (pthread_join(pthread[t], NULL))
      fatal("Cannot join thread");

  xfree(pthread);

  pthread_attr_destroy(&attr);

  pthread_mutex_destroy(&mutex);

  progress_done();
}

void hardmask(char * seq, int len)
{
  /* convert all lower case letters in seq to N */
  
  for(int j=0; j<len; j++)
    if (seq[j] & 0x20)
      seq[j] = 'N';
}

void hardmask_all()
{
  for(uint64_t i=0; i<db_getsequencecount(); i++)
    hardmask(db_getsequence(i), db_getsequencelen(i));
}

void maskfasta()
{
  FILE * fp_output = fopen_output(opt_output);
  if (!fp_output)
    fatal("Unable to open mask output file for writing");

  db_read(opt_maskfasta, 0);
  show_rusage();

  seqcount = db_getsequencecount();

  if (opt_qmask == MASK_DUST)
    dust_all();
  else if ((opt_qmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();
  show_rusage();

  progress_init("Writing output", seqcount);
  for(int i=0; i<seqcount; i++)
    {
      fasta_print_db_relabel(fp_output, i, i+1);
      progress_update(i);
    }
  progress_done();
  show_rusage();

  db_free();
  fclose(fp_output);
}

void fastx_mask()
{
  FILE * fp_fastaout = 0;
  FILE * fp_fastqout = 0;

  if (opt_fastaout)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (!fp_fastaout)
        fatal("Unable to open mask output FASTA file for writing");
    }

  if (opt_fastqout)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (!fp_fastqout)
        fatal("Unable to open mask output FASTQ file for writing");
    }

  db_read(opt_fastx_mask, 0);
  show_rusage();

  if (fp_fastqout && ! db_is_fastq())
    fatal("Cannot write FASTQ output with a FASTA input file, lacking quality scores");

  seqcount = db_getsequencecount();

  if (opt_qmask == MASK_DUST)
    dust_all();
  else if ((opt_qmask == MASK_SOFT) && (opt_hardmask))
    hardmask_all();
  show_rusage();

  int kept = 0;
  int discarded_less = 0;
  int discarded_more = 0;
  progress_init("Writing output", seqcount);
  for(int i=0; i<seqcount; i++)
    {
      int unmasked = 0;
      char * seq = db_getsequence(i);
      int len = db_getsequencelen(i);
      if (opt_qmask == MASK_NONE)
        {
          unmasked = len;
        }
      else if (opt_hardmask)
        {
          for(int j=0; j<len; j++)
            if (seq[j] != 'N')
              unmasked++;
        }
      else
        {
          for(int j=0; j<len; j++)
            if (isupper(seq[j]))
              unmasked++;
        }
      double unmasked_pct = 100.0 * unmasked / len;

      if (unmasked_pct < opt_min_unmasked_pct)
        discarded_less++;
      else if (unmasked_pct >  opt_max_unmasked_pct)
        discarded_more++;
      else
        {
          kept++;

          if (opt_fastaout)
            fasta_print_general(fp_fastaout,
                                0,
                                seq,
                                len,
                                db_getheader(i),
                                db_getheaderlen(i),
                                db_getabundance(i),
                                kept,
                                -1, -1, 0, 0.0);

          if (opt_fastqout)
            fastq_print_general(fp_fastqout,
                                seq,
                                len,
                                db_getheader(i),
                                db_getheaderlen(i),
                                db_getquality(i),
                                db_getabundance(i),
                                kept,
                                0, 0.0);
        }

      progress_update(i);
    }
  progress_done();

  if (!opt_quiet)
    {
      if (opt_min_unmasked_pct > 0.0)
        fprintf(stderr, "%d sequences with less than %.1lf%% unmasked residues discarded\n", discarded_less, opt_min_unmasked_pct);
      if (opt_max_unmasked_pct < 100.0)
        fprintf(stderr, "%d sequences with more than %.1lf%% unmasked residues discarded\n", discarded_more, opt_max_unmasked_pct);
      fprintf(stderr, "%d sequences kept\n", kept);
    }

  if (opt_log)
    {
      if (opt_min_unmasked_pct > 0.0)
        fprintf(fp_log, "%d sequences with less than %.1lf%% unmasked residues discarded\n", discarded_less, opt_min_unmasked_pct);
      if (opt_max_unmasked_pct < 100.0)
        fprintf(fp_log, "%d sequences with more than %.1lf%% unmasked residues discarded\n", discarded_more, opt_max_unmasked_pct);
      fprintf(fp_log, "%d sequences kept\n", kept);
    }   

  show_rusage();
  db_free();

  if (fp_fastaout)
    fclose(fp_fastaout);
  if (fp_fastqout)
    fclose(fp_fastqout);
}
