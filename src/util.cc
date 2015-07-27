/*
    Copyright (C) 2014-2015 Torbjorn Rognes & Tomas Flouri

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
  if (! opt_quiet)
    {
      progress_prompt = prompt;
      progress_size = size;
      progress_chunk = size < progress_granularity ? 
        1 : size / progress_granularity;
      progress_next = 0;
      fprintf(stderr, "%s %.0f%%", prompt, 0.0);
    }
}

void progress_update(unsigned long progress)
{
  if (! opt_quiet)
    if (progress >= progress_next)
      {
        fprintf(stderr, "  \r%s %.0f%%", progress_prompt,
                100.0 * progress / progress_size);
        progress_next = progress + progress_chunk;
      }
}

void progress_done()
{
  if (! opt_quiet)
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
  fprintf(stderr, "\n\n");
  fprintf(stderr, "Fatal error: %s\n", msg);

  if (opt_log)
    {
      fprintf(fp_log, "\n\n");
      fprintf(fp_log, "Fatal error: %s\n", msg);
    }

  exit(EXIT_FAILURE);
}

void  __attribute__((noreturn)) fatal(const char * format, 
                                      const char * message)
{
  fprintf(stderr, "\n\n");
  fprintf(stderr, format, message);
  fprintf(stderr, "\n");

  if (opt_log)
    {
      fprintf(fp_log, "\n\n");
      fprintf(fp_log, format, message);
      fprintf(fp_log, "\n");
    }

  exit(EXIT_FAILURE);
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t;
  (void) posix_memalign(& t, alignment, size);
  if (t==0)
    fatal("Unable to allocate enough memory.");
  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  if (size == 0)
    size = 1;
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

char * xstrdup(const char *s)
{
  size_t len = strlen(s);
  char * p = (char*) xmalloc(len+1);
  return strcpy(p, s);
}


char * xstrchrnul(char *s, int c)
{
  char * r = strchr(s, c);

  if (r)
    return r;
  else
    return (char *)s + strlen(s);
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

  if (opt_log)
    {
      fprintf(fp_log, "Time: %.3fs (user)", r_usage.ru_utime.tv_sec * 1.0 + (double) r_usage.ru_utime.tv_usec * 1.0e-6);
      fprintf(fp_log, " %.3fs (sys)", r_usage.ru_stime.tv_sec * 1.0 + r_usage.ru_stime.tv_usec * 1.0e-6);
      
      fprintf(fp_log, " Memory: %luMB\n",  arch_get_memused() / 1024 / 1024);
    }
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

void reverse_complement(char * rc, char * seq, long len)
{
  /* Write the reverse complementary sequence to rc.
     The memory for rc must be long enough for the rc of the sequence
     (identical to the length of seq + 1. */

  for(long i=0; i<len; i++)
    rc[i] = chrmap_complement[(int)(seq[len-1-i])];
  rc[len] = 0;
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

void random_init()
{
  /* initialize pseudo-random number generator */
  unsigned int seed = opt_randseed;
  if (seed == 0)
    {
      int fd = open("/dev/urandom", O_RDONLY);
      if (fd < 0)
        fatal("Unable to open /dev/urandom");
      if (read(fd, & seed, sizeof(seed)) < 0)
        fatal("Unable to read from /dev/urandom");
      close(fd);
    }
  srandom(seed);
}

long random_int(long n)
{
  /*
    Generate a random integer in the range 0 to n-1, inclusive.
    n must be > 0
    The random() function returns a random number in the range
    0 to 2147483647 (=2^31-1=RAND_MAX), inclusive.
    We should avoid some of the upper generated numbers to
    avoid modulo bias.
  */

  long random_max = RAND_MAX;
  long limit = random_max - (random_max + 1) % n;
  long r = random();
  while (r > limit)
    r = random();
  return r % n;
}

unsigned long random_ulong(unsigned long n)
{
  /*
    Generate a random integer in the range 0 to n-1, inclusive,
    n must be > 0
  */

  unsigned long random_max = ULONG_MAX;
  unsigned long limit = random_max - (random_max - n + 1) % n;
  unsigned long r = ((random() << 48) ^ (random() << 32) ^
                     (random() << 16) ^ (random()));
  while (r > limit)
    r = ((random() << 48) ^ (random() << 32) ^
         (random() << 16) ^ (random()));
  return r % n;
}
