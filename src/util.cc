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

#include "vsearch.h"

static const char * progress_prompt;
static unsigned long progress_next;
static unsigned long progress_size;
static unsigned long progress_chunk;
static const unsigned long progress_granularity = 200;

#ifdef HAVE_BZLIB
static unsigned char magic_bzip[] = "\x42\x5a";
#endif
#ifdef HAVE_ZLIB
static unsigned char magic_gzip[] = "\x1f\x8b";
#endif

void progress_init(const char * prompt, unsigned long size)
{
  progress_prompt = prompt;
  progress_size = size;
  progress_chunk = size < progress_granularity ? 
    1 : size / progress_granularity;
  progress_next = 0;
  fprintf(stderr, "%s %.0f%%", prompt, 0.0);
}

void progress_update(unsigned long progress)
{
  if (progress >= progress_next)
    {
      fprintf(stderr, "  \r%s %.0f%%", progress_prompt,
	      100.0 * progress / progress_size);
      progress_next = progress + progress_chunk;
    }
}

void progress_done()
{
  fprintf(stderr, "  \r%s %.0f%%\n", progress_prompt, 100.0);
}

long gcd(long a, long b)
{
  if (b == 0)
  {
    return a;
  }
  else
  {
    return gcd(b, a % b);
  }
}

void  __attribute__((noreturn)) fatal(const char * msg)
{
  fprintf(stderr, "Error: %s\n", msg);
  exit(1);
}

void  __attribute__((noreturn)) fatal(const char * format, 
				      const char * message)
{
  fprintf(stderr, format, message);
  fprintf(stderr, "\n");
  exit(1);
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t;
  posix_memalign(& t, alignment, size);

  if (t==NULL)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

char * xstrchrnul(char *s, int c)
{
  char * r = strchr(s, c);

  if (r)
    return r;
  else
    return (char *)s + strlen(s);
}

unsigned long hash_fnv_1a_64(char * s, unsigned long n)
{
  /*
    This is the Fowler Noll Vo (FNV) hash function,
    version 1a (FNV-1a), 64 bit
    https://en.wikipedia.org/wiki/Fowler-Noll-Vo_hash_function
    It is in the public domain.
  */

  const unsigned long fnv_offset = 14695981039346656037UL;
  const unsigned long fnv_prime = 1099511628211;

  unsigned long hash = fnv_offset;

  for(unsigned long i = 0; i < n; i++)
    {
      unsigned char c = (unsigned char) *s++;
      hash = (hash ^ c) * fnv_prime;
    }

  return hash;
}

unsigned long hash_fnv_1a_64_uc(char * s, unsigned long n)
{
  /*
    compute hash of upper case string
  */

  const unsigned long fnv_offset = 14695981039346656037UL;
  const unsigned long fnv_prime = 1099511628211;

  unsigned long hash = fnv_offset;

  for(unsigned long i = 0; i < n; i++)
    {
      unsigned char c = toupper(*s++);
      hash = (hash ^ c) * fnv_prime;
    }

  return hash;
}

unsigned long hash_cityhash64(char * s, unsigned long n)
{
  return CityHash64((const char*)s, n);
}

long getusec(void)
{
  struct timeval tv;
  if(gettimeofday(&tv,0) != 0) return 0;
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

void show_rusage()
{
#if 0
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);
  
  fprintf(stderr, "Time: %.3fs (user)", r_usage.ru_utime.tv_sec * 1.0 + (double) r_usage.ru_utime.tv_usec * 1.0e-6);
  fprintf(stderr, " %.3fs (sys)", r_usage.ru_stime.tv_sec * 1.0 + r_usage.ru_stime.tv_usec * 1.0e-6);

  fprintf(stderr, " Memory: %luMB\n",  arch_get_memused() / 1024 / 1024);
#endif
}

void fprint_fasta_hdr_only(FILE * fp, char * hdr)
{
  fprintf(fp, ">%s\n", hdr);
}

void fprint_fasta_seq_only(FILE * fp, char * seq, unsigned long len, int width)
{
  /* 
     The actual length of the sequence may be longer than "len", but only
     "len" characters are printed.
     
     Specify width of lines - zero (or <1)  means linearize (all on one line).
  */

  if (width < 1)
    fprintf(fp, "%.*s\n", (int)(len), seq);
  else
    {
      long rest = len;
      for(unsigned long i=0; i<len; i += width)
	{
	  fprintf(fp, "%.*s\n", (int)(MIN(rest,width)), seq+i);
	  rest -= width;
	}
    }
}

void db_fprint_fasta(FILE * fp, unsigned long seqno)
{
  char * hdr = db_getheader(seqno);
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  
  fprint_fasta_hdr_only(fp, hdr);
  fprint_fasta_seq_only(fp, seq, seqlen, opt_fasta_width);
}

void db_fprint_fasta_with_size(FILE * fp, unsigned long seqno, unsigned long size)
{
  char * hdr = db_getheader(seqno);
  int hdrlen = db_getheaderlen(seqno);
  char * seq = db_getsequence(seqno);
  long seqlen = db_getsequencelen(seqno);
  
  /* remove any previous size annotation */
  /* regexp search for "(^|;)(\d+)(;|$)" */
  /* replace by ';' if not at either end */

  regmatch_t pmatch[1];
  if (!regexec(&db_regexp, hdr, 1, pmatch, 0))
    {
      int pat_start = pmatch[0].rm_so;
      int pat_end = pmatch[0].rm_eo;

      /* print initial part */
      fprintf(fp, ">%.*s", pat_start, hdr);

      /* replace old size with ";" if not at any end */
      if ((pat_start > 0) && (pat_end < hdrlen))
	fprintf(fp, ";");

      /* print remaining part */
      fprintf(fp, "%.*s", hdrlen - pat_end, hdr + pat_end);

      /* print size annotation */
      fprintf(fp, ";size=%lu;\n", size);
    }
  else
    {
      /* print entire header + size annotation */
      fprintf(fp, ">%s;size=%lu;\n", hdr, size);
    }

  fprint_fasta_seq_only(fp, seq, seqlen, opt_fasta_width);
}

void reverse_complement(char * rc, char * seq, long len)
{
  /* Write the reverse complementary sequence to rc.
     The memory for rc must be long enough for the rc of the sequence
     (identical to the length of seq + 1. */

#if 0
  printf("Sequence:\n");
  fprint_fasta_seq_only(stdout, seq, len, 60);
#endif
  for(long i=0; i<len; i++)
    {
      rc[i] = chrmap_complement[(int)(seq[len-1-i])];
    }
  rc[len] = 0;
#if 0
  printf("RC:\n");
  fprint_fasta_seq_only(stdout, rc, len, 60);
#endif
}

#ifdef HAVE_BZLIB
char * bz_fgets (char * s, int size, BZFILE * stream, long linealloc,
                 int * bz_error_ptr, char * buf_internal, long * buf_internal_len_ptr)
{
  long linelen    = 0;
  long buf_internal_len = *buf_internal_len_ptr;

  /* fill the buffer */
  long bytes_read = BZ2_bzRead(bz_error_ptr, stream, buf_internal + buf_internal_len,
                               linealloc - buf_internal_len - 1);

  buf_internal_len              += bytes_read;
  buf_internal[buf_internal_len] = 0;

  if (buf_internal_len)
   {
     linelen = xstrchrnul(buf_internal, '\n') - buf_internal;
     memcpy (s, buf_internal, linelen);
     s[linelen] = 0;
   }

  /* if newline found */
  if (buf_internal[linelen])
   memmove (buf_internal, buf_internal + linelen + 1, buf_internal_len - linelen); 

  if (buf_internal_len - linelen > 0)
    buf_internal_len = buf_internal_len - linelen - 1;

  *buf_internal_len_ptr = buf_internal_len;
  return s;
}
#endif

int detect_compress_format (const char * filename)
{
  /* check for magic numbers to detect file type */
  unsigned char magic[3];
  int cnt;
  FILE * fp;

  if (!(fp = fopen(filename, "r")))
    fatal("Error: Unable to open database file (%s)", filename);

  cnt = fread(magic, sizeof(unsigned char), 2, fp);
  fclose(fp);
  if (cnt < 2) return (0); 

#ifdef HAVE_BZLIB
  if (!memcmp(magic, magic_bzip, 2)) return FORMAT_BZIP;
#endif

#ifdef HAVE_ZLIB
  if (!memcmp(magic, magic_gzip, 2)) return FORMAT_GZIP;
#endif

  return FORMAT_PLAIN;
}

