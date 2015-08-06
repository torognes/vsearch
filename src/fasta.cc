/*

Simple FASTA reader capable of indefinite long lines and compressed files

1. Start where the initial ">" is.
2. If there is no ">": abort
3. Find header until end of line.
4. Read sequence lines
5. Stop when line starts with ">" or NULL.

RETURN

1. pointer to buffer with header
2. length of header
3. pointer to buffer with sequence
4. length of sequence


Use these functions:

int BZ2_bzRead(int * bzerror, BZFILE * b, void * buf, int len);
+ BZ2_bzReadOpen and BZ2_vzReadClose

ZEXTERN int ZEXPORT gzread OF((gzFile file, voidp buf, unsigned len));
+ gzopen and gzclose

size_t fread(void *restrict ptr, size_t size, size_t nitems, FILE *restrict stream);
+ fopen and fclose


One buffer for reading, perhaps 4k, fixed
One buffer for the header, perhaps 4k, growing
One buffer for the sequence, perhaps 4k, growing

These functions do not filter or analyse the sequences in any way.

Buffers are reused and caller must copy data elsewhere between calls.
*/

#include "vsearch.h"

#define FASTA_BUFFER_ALLOC 8192
static unsigned char MAGIC_GZIP[] = "\x1f\x8b";
static unsigned char MAGIC_BZIP[] = "BZ";

fasta_handle fasta_file_open(const char * filename)
{
  fasta_handle h = (fasta_handle) xmalloc(sizeof(struct fasta_s));
  
  h->fp = NULL;
  h->fp = fopen(filename, "rb");
  if(!h->fp)
    fatal("Error: Unable to open fasta file for reading (%s)", filename);
  
  /* detect compression (plain, gzipped or bzipped) */
  
  unsigned char magic[2];
  if (fread(&magic, 1, 2, h->fp) < 2)
    fatal("Error: File too small (less than 2 characters)");

  if (!memcmp(magic, MAGIC_GZIP, 2))
    h->format = FORMAT_GZIP;
  else if (!memcmp(magic, MAGIC_BZIP, 2))
    h->format = FORMAT_BZIP;
  else
    h->format = FORMAT_PLAIN;

  if (fseek(h->fp, 0, SEEK_END))
    fatal("Error: Unable to seek in fasta file (%s)", filename);
  h->file_size = ftell(h->fp);
  rewind(h->fp);

#ifdef HAVE_ZLIB
  /* GZIP: Close ordinary file and open again as gzipped file */
  if (h->format == FORMAT_GZIP)
    {
      fclose(h->fp);
      if (! (h->fp_gz = gzopen(filename, "rb")))
        fatal("Unable to open gzip compressed fasta file (%s)", filename);
    }
#endif

#ifdef HAVE_ZLIB
  /* BZIP: Keep original file open, then open as bzipped file as well */
  if (h->format == FORMAT_BZIP)
    {
      int bzError;
      if (! (h->fp_bz = BZ2_bzReadOpen(& bzError, h->fp,
                                       BZ_VERBOSE_0, BZ_MORE_MEM, NULL, 0)))
        fatal("Unable to open bzip2 compressed fasta file (%s)", filename);
    }
#endif

  h->file_position = 0;
  h->file_buffer_position = 0;

  h->file_buffer_alloc = FASTA_BUFFER_ALLOC;
  h->header_buffer_alloc = FASTA_BUFFER_ALLOC;
  h->sequence_buffer_alloc = FASTA_BUFFER_ALLOC;
  
  h->file_buffer = (char*) xmalloc(h->file_buffer_alloc);
  h->header_buffer = (char*) xmalloc(h->header_buffer_alloc);
  h->sequence_buffer = (char*) xmalloc(h->sequence_buffer_alloc);

  h->header_buffer[0] = 0;
  h->sequence_buffer[0] = 0;

  h->file_buffer_length = 0;
  h->header_buffer_length = 0;
  h->sequence_buffer_length = 0;

  return h;
}

void fasta_file_close(fasta_handle h)
{
#ifdef HAVE_ZLIB
  if (h->format == FORMAT_GZIP)
    {
      gzclose(h->fp_gz);
      h->fp_gz = 0;
    }
#endif
  
#ifdef HAVE_BZLIB
  if (h->format == FORMAT_BZIP)
    {
      int bz_error;
      BZ2_bzReadClose(&bz_error, h->fp_bz);
      h->fp_bz = 0;
    }
#endif

  if (h->format != FORMAT_GZIP)
    {
      fclose(h->fp);
      h->fp = 0;
    }

  if (h->file_buffer)
    free(h->file_buffer);
  if (h->header_buffer)
    free(h->header_buffer);
  if (h->sequence_buffer)
    free(h->sequence_buffer);
  
  h->file_buffer = 0;
  h->header_buffer = 0;
  h->sequence_buffer = 0;

  h->file_buffer_length = 0;
  h->header_buffer_length = 0;
  h->sequence_buffer_length = 0;

  h->file_buffer_alloc = 0;
  h->header_buffer_alloc = 0;
  h->sequence_buffer_alloc = 0;

  h->file_size = 0;
  h->file_position = 0;
  h->file_buffer_position = 0;
}


size_t fasta_file_fill_buffer(fasta_handle h)
{
  /* read more data if necessary */
  size_t rest = h->file_buffer_length - h->file_buffer_position;
  
  if (rest > 0)
    return rest;
  else
    {
      size_t space = h->file_buffer_alloc - h->file_buffer_length;

      if (space == 0)
        {
          /* back to beginning of buffer */
          h->file_buffer_position = 0;
          h->file_buffer_length = 0;
          space = h->file_buffer_alloc;
        }
      
      size_t bytes_read_plain = 0;
      int bytes_read_gz = 0;
      int bytes_read_bz = 0;
      int bzError = 0;

      switch(h->format)
        {
        case FORMAT_GZIP:
#ifdef HAVE_ZLIB
          bytes_read_gz = gzread(h->fp_gz,
                                 h->file_buffer + h->file_buffer_position,
                                 space);
          if (bytes_read_gz < 0)
            fatal("Error reading gzip compressed fasta file");
          h->file_buffer_length += bytes_read_gz;
          return bytes_read_gz;
#else
          fatal("gzip compressed fasta files not supported.");
          return 0;
#endif
          
        case FORMAT_BZIP:
#ifdef HAVE_BZLIB
          bytes_read_bz = BZ2_bzRead(& bzError,
                                     h->fp_bz,
                                     h->file_buffer + h->file_buffer_position,
                                     space);
          if ((bytes_read_bz < 0) || 
              ! ((bzError == BZ_OK) || (bzError == BZ_STREAM_END) || (bzError == BZ_SEQUENCE_ERROR)))
            {
              //              printf("bzError = %d bytes_read = %d\n", bzError, bytes_read_bz);
              fatal("Error reading bzip2 compressed fasta file");
            }
          h->file_buffer_length += bytes_read_bz;
          return bytes_read_bz;
#else
          fatal("bzip2 compressed fasta files not supported.");
          return 0;
#endif

        case FORMAT_PLAIN:
          bytes_read_plain = fread(h->file_buffer + h->file_buffer_position,
                                   1,
                                   space,
                                   h->fp);
          h->file_buffer_length += bytes_read_plain;
          return bytes_read_plain;
          
        default:
          fatal("Internal error");
          return 0;
        }
    }
}

void fasta_header_extend(fasta_handle h, char * buf, size_t len)
{
  if (h->header_buffer_length + len + 1 > h->header_buffer_alloc)
    {
      /* alloc space for len more characters + terminating zero,
         but round up to nearest block size */
      h->header_buffer_alloc = 
        (FASTA_BUFFER_ALLOC * 
         ((h->header_buffer_length + len) / FASTA_BUFFER_ALLOC) + 1);
      h->header_buffer = (char*) xrealloc(h->header_buffer, 
                                          h->header_buffer_alloc);
    }

  /* copy string */
  memcpy(h->header_buffer + h->header_buffer_length,
         buf, len);
  h->header_buffer_length += len;

  /* add terminator */
  h->header_buffer[h->header_buffer_length] = 0;
}

void fasta_sequence_extend(fasta_handle h, char * buf, size_t len)
{
  if (h->sequence_buffer_length + len + 1 > h->sequence_buffer_alloc)
    {
      /* alloc space for len more characters + terminating zero,
         but round up to nearest block size */
      h->sequence_buffer_alloc = 
        (FASTA_BUFFER_ALLOC * 
         ((h->sequence_buffer_length + len) / FASTA_BUFFER_ALLOC) + 1);
      h->sequence_buffer = (char*) xrealloc(h->sequence_buffer, 
                                          h->sequence_buffer_alloc);
    }

  /* copy string */
  memcpy(h->sequence_buffer + h->sequence_buffer_length,
         buf, len);
  h->sequence_buffer_length += len;

  /* add terminator */
  h->sequence_buffer[h->sequence_buffer_length] = 0;
}

bool fasta_file_read(fasta_handle h, 
                     char * * header, size_t * header_length,
                     char * * sequence, size_t * sequence_length)
{
  h->header_buffer_length = 0;
  h->sequence_buffer_length = 0;

  size_t rest = fasta_file_fill_buffer(h);

  if (rest == 0)
    {
      header = 0;
      header_length = 0;
      sequence = 0;
      sequence_length = 0;
      return 0;
    }

  /* read header */

  /* check initial > character */
  
  if (h->file_buffer[h->file_buffer_position] != '>')
    fatal("Invalid FASTA - header must start with > character");
  h->file_buffer_position++;
  rest--;

  char * lf = 0;
  while (lf == 0)
    {
      /* get more data */
      rest = fasta_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTA - header must be terminated with newline");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer + h->file_buffer_position,
                           '\n', rest);

      /* copy to header buffer */
      size_t len;
      if (lf == 0)
        {
          /* no LF in buffer, copy all */
          len = rest;
        }
      else
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer + h->file_buffer_position) + 1;
        }
      fasta_header_extend(h, h->file_buffer + h->file_buffer_position, len);
      h->file_buffer_position += len;
      rest -= len;
    }

  /* read one or more sequence lines */

  while (1)
    {
      /* get more data, if necessary */
      rest = fasta_file_fill_buffer(h);
      if (rest == 0)
        break;

      if (lf && (h->file_buffer[h->file_buffer_position] == '>'))
        break;

      /* find LF */
      lf = (char *) memchr(h->file_buffer + h->file_buffer_position,
                           '\n', rest);

      size_t len;
      if (lf == 0)
        {
          /* no LF in buffer, copy all */
          len = rest;
        }
      else
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer + h->file_buffer_position) + 1;
        }
      fasta_sequence_extend(h, h->file_buffer + h->file_buffer_position, len);
      h->file_buffer_position += len;
      rest -= len;
    }

  * header = h->header_buffer;
  * header_length = h->header_buffer_length;
  * sequence = h->sequence_buffer;
  * sequence_length = h->sequence_buffer_length;

  return 1;
}


size_t fasta_file_getpos(fasta_handle h)
{
  return h->file_position;
}

size_t fasta_file_getsize(fasta_handle h)
{
  return h->file_size;
}
