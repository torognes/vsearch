/*

Simple FASTA reader capable of indefinite long lines and compressed files

*/

#include "vsearch.h"

#define FASTA_BUFFER_ALLOC 8192
static unsigned char MAGIC_GZIP[] = "\x1f\x8b";
static unsigned char MAGIC_BZIP[] = "BZ";

void buffer_init(struct buffer_s * buffer)
{
  buffer->alloc = FASTA_BUFFER_ALLOC;
  buffer->data = (char*) xmalloc(buffer->alloc);
  buffer->data[0] = 0;
  buffer->length = 0;
  buffer->position = 0;
}

void buffer_free(struct buffer_s * buffer)
{
  if (buffer->data)
    free(buffer->data);
  buffer->data = 0;
  buffer->alloc = 0;
  buffer->length = 0;
  buffer->position = 0;
}

fasta_handle fasta_open(const char * filename)
{
  fasta_handle h = (fasta_handle) xmalloc(sizeof(struct fasta_s));
  
  h->fp = NULL;
  h->fp = fopen(filename, "rb");
  if(!h->fp)
    fatal("Error: Unable to open fasta file for reading (%s)", filename);
  
  /* detect compression (plain, gzipped or bzipped) */
  
  unsigned char magic[2];
  h->format = FORMAT_PLAIN;
  if (fread(&magic, 1, 2, h->fp) >= 2)
    {
      if (!memcmp(magic, MAGIC_GZIP, 2))
        h->format = FORMAT_GZIP;
      else if (!memcmp(magic, MAGIC_BZIP, 2))
        h->format = FORMAT_BZIP;
    }

  if (fseek(h->fp, 0, SEEK_END))
    fatal("Error: Unable to seek in fasta file (%s)", filename);
  h->file_size = ftell(h->fp);
  rewind(h->fp);

  /* GZIP: Close ordinary file and open again as gzipped file */
  if (h->format == FORMAT_GZIP)
    {
#ifdef HAVE_ZLIB
      fclose(h->fp);
      if (! (h->fp_gz = gzopen(filename, "rb")))
        fatal("Unable to open gzip compressed fasta file (%s)", filename);
#else
      fatal("Files compressed with gzip are not supported");
#endif
    }

  /* BZIP: Keep original file open, then open as bzipped file as well */
  if (h->format == FORMAT_BZIP)
    {
#ifdef HAVE_ZLIB
      int bzError;
      if (! (h->fp_bz = BZ2_bzReadOpen(& bzError, h->fp,
                                       BZ_VERBOSE_0, BZ_MORE_MEM, NULL, 0)))
        fatal("Unable to open bzip2 compressed fasta file (%s)", filename);
#else
      fatal("Files compressed with bzip2 are not supported");
#endif
    }

  h->stripped_all = 0;

  for(int i=0; i<256; i++)
    h->stripped[i] = 0;

  if (regcomp(& h->size_regexp, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
    fatal("Compilation of regular expression for abundance annotation failed");
  
  h->file_position = 0;

  buffer_init(& h->file_buffer);
  buffer_init(& h->header_buffer);
  buffer_init(& h->sequence_buffer);

  h->lineno = 1;
  h->seqno = -1;

  return h;
}

void fasta_close(fasta_handle h)
{
  /* Warn about stripped chars */

  if (h->stripped_all)
    {
      fprintf(stderr, "WARNING: invalid characters stripped from fasta file:");
      for (int i=0; i<256;i++)
        if (h->stripped[i])
          fprintf(stderr, " %c(%lu)", i, h->stripped[i]);
      fprintf(stderr, "\n");

      if (opt_log)
        {
          fprintf(fp_log, "WARNING: invalid characters stripped from fasta file:");
          for (int i=0; i<256;i++)
            if (h->stripped[i])
              fprintf(fp_log, " %c(%lu)", i, h->stripped[i]);
          fprintf(fp_log, "\n");
        }
    }

#ifdef HAVE_BZLIB
  int bz_error;
#endif
  
  switch(h->format)
    {
    case FORMAT_PLAIN:
      fclose(h->fp);
      h->fp = 0;
      break;

    case FORMAT_GZIP:
#ifdef HAVE_ZLIB
      gzclose(h->fp_gz);
      h->fp_gz = 0;
      break;
#endif
      
    case FORMAT_BZIP:
#ifdef HAVE_BZLIB
      BZ2_bzReadClose(&bz_error, h->fp_bz);
      h->fp_bz = 0;
      break;
#endif

    default:
      fatal("Internal error");
    }
  
  regfree(& h->size_regexp);

  buffer_free(& h->file_buffer);
  buffer_free(& h->header_buffer);
  buffer_free(& h->sequence_buffer);

  h->file_size = 0;
  h->file_position = 0;

  h->lineno = 0;
  h->seqno = -1;
}


size_t fasta_file_fill_buffer(fasta_handle h)
{
  /* read more data if necessary */
  size_t rest = h->file_buffer.length - h->file_buffer.position;
  
  if (rest > 0)
    return rest;
  else
    {
      size_t space = h->file_buffer.alloc - h->file_buffer.length;

      if (space == 0)
        {
          /* back to beginning of buffer */
          h->file_buffer.position = 0;
          h->file_buffer.length = 0;
          space = h->file_buffer.alloc;
        }
      
      size_t bytes_read_plain = 0;

#ifdef HAVE_ZLIB
      int bytes_read_gz = 0;
#endif

#ifdef HAVE_BZLIB
      int bytes_read_bz = 0;
      int bzError = 0;
#endif

      switch(h->format)
        {
        case FORMAT_PLAIN:
          bytes_read_plain = fread(h->file_buffer.data
                                   + h->file_buffer.position,
                                   1,
                                   space,
                                   h->fp);
          h->file_buffer.length += bytes_read_plain;
          return bytes_read_plain;

        case FORMAT_GZIP:
#ifdef HAVE_ZLIB
          bytes_read_gz = gzread(h->fp_gz,
                                 h->file_buffer.data + h->file_buffer.position,
                                 space);
          if (bytes_read_gz < 0)
            fatal("Error reading gzip compressed fasta file");
          h->file_buffer.length += bytes_read_gz;
          return bytes_read_gz;
#endif
          
        case FORMAT_BZIP:
#ifdef HAVE_BZLIB
          bytes_read_bz = BZ2_bzRead(& bzError,
                                     h->fp_bz,
                                     h->file_buffer.data 
                                     + h->file_buffer.position,
                                     space);
          if ((bytes_read_bz < 0) ||
              ! ((bzError == BZ_OK) ||
                 (bzError == BZ_STREAM_END) ||
                 (bzError == BZ_SEQUENCE_ERROR)))
            fatal("Error reading bzip2 compressed fasta file");
          h->file_buffer.length += bytes_read_bz;
          return bytes_read_bz;
#endif
          
        default:
          fatal("Internal error");
          return 0;
        }
    }
}

void buffer_extend(struct buffer_s * buffer, char * buf, size_t len)
{
  if (buffer->length + len + 1 > buffer->alloc)
    {
      /* alloc space for len more characters + terminating zero,
         but round up to nearest block size */
      buffer->alloc = 
        (FASTA_BUFFER_ALLOC * 
         ((buffer->length + len) / FASTA_BUFFER_ALLOC) + 1);
      buffer->data = (char*) xrealloc(buffer->data, buffer->alloc);
    }

  /* copy string */
  memcpy(buffer->data + buffer->length, buf, len);
  buffer->length += len;

  /* add terminator */
  buffer->data[buffer->length] = 0;
}

void fasta_truncate_header(fasta_handle h, bool truncateatspace)
{
  /* Truncate the zero-terminated header string by inserting a new
     terminator (zero byte) at the first space (if truncateatspace)
     or first linefeed. */
  
  if (truncateatspace)
    h->header_buffer.length = strcspn(h->header_buffer.data, " \n");
  else
    h->header_buffer.length = strcspn(h->header_buffer.data, "\n");
  
  h->header_buffer.data[h->header_buffer.length] = 0;
}

void fasta_filter_sequence(fasta_handle h,
                           unsigned int * char_action,
                           char * char_mapping)
{
  /* Strip unwanted characters from the sequence and raise warnings or
     errors on certain characters. */

  char * p = h->sequence_buffer.data;
  char * q = p;
  char c;
  char msg[200];

  while ((c = *p++))
    {
      char m = char_action[(int)c];

      switch(m)
        {
        case 0:
          /* stripped */
          h->stripped_all++;
          h->stripped[(int)c]++;
          break;
              
        case 1:
          /* legal character */
          *q++ = char_mapping[(int)(c)];
          break;
          
        case 2:
          /* fatal character */
          if (c>=32)
            snprintf(msg,
                     200,
                     "illegal character '%c' on line %lu in fasta file",
                     c,
                     h->lineno);
          else
            snprintf(msg,
                     200,
                     "illegal character %#.2x (hexadecimal) on line %lu in fasta file",
                     c,
                     h->lineno);
          fatal(msg);
          break;
              
        case 3:
          /* silently stripped chars (whitespace) */
          break;

        case 4:
          /* newline (silently stripped) */
          h->lineno++;
          break;
        }
    }

  /* add zero after sequence */
  *q = 0;
  h->sequence_buffer.length = q - h->sequence_buffer.data;
}

bool fasta_next(fasta_handle h,
                bool truncateatspace,
                unsigned int * char_action,
                char * char_mapping)
{
  h->header_buffer.length = 0;
  h->sequence_buffer.length = 0;

  size_t rest = fasta_file_fill_buffer(h);

  if (rest == 0)
    return 0;

  /* read header */

  /* check initial > character */
  
  if (h->file_buffer.data[h->file_buffer.position] != '>')
    fatal("Invalid FASTA - header must start with > character");
  h->file_buffer.position++;
  rest--;

  char * lf = 0;
  while (lf == 0)
    {
      /* get more data if buffer empty*/
      rest = fasta_file_fill_buffer(h);
      if (rest == 0)
        fatal("Invalid FASTA - header must be terminated with newline");
      
      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n',
                           rest);

      /* copy to header buffer */
      size_t len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
          h->lineno++;
        }
      buffer_extend(& h->header_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  /* read one or more sequence lines */

  while (1)
    {
      /* get more data, if necessary */
      rest = fasta_file_fill_buffer(h);

      /* end if no more data */
      if (rest == 0)
        break;

      /* end if new sequence starts */
      if (lf && (h->file_buffer.data[h->file_buffer.position] == '>'))
        break;

      /* find LF */
      lf = (char *) memchr(h->file_buffer.data + h->file_buffer.position,
                           '\n', rest);

      size_t len = rest;
      if (lf)
        {
          /* LF found, copy up to and including LF */
          len = lf - (h->file_buffer.data + h->file_buffer.position) + 1;
        }
      buffer_extend(& h->sequence_buffer,
                    h->file_buffer.data + h->file_buffer.position,
                    len);
      h->file_buffer.position += len;
      rest -= len;
    }

  h->seqno++;

  fasta_truncate_header(h, truncateatspace);
  fasta_filter_sequence(h, char_action, char_mapping);


  /* Need to change if compressed !!! */
  h->file_position = ftell(h->fp);

  return 1;
}

long fasta_get_abundance(fasta_handle h)
{
  /* read size/abundance annotation */
  
  unsigned long abundance = 1;
  regmatch_t pmatch[4];
  
  if (!regexec(& h->size_regexp, h->header_buffer.data, 4, pmatch, 0))
    {
      long number = atol(h->header_buffer.data + pmatch[2].rm_so);
      if (number > 0)
        abundance = number;
      else
        fatal("Invalid (zero) abundance annotation in fasta header");
    }
  return abundance;
}

size_t fasta_get_position(fasta_handle h)
{
  return h->file_position;
}

size_t fasta_get_size(fasta_handle h)
{
  return h->file_size;
}

size_t fasta_get_lineno(fasta_handle h)
{
  return h->lineno;
}

size_t fasta_get_seqno(fasta_handle h)
{
  return h->seqno;
}

size_t fasta_get_header_length(fasta_handle h)
{
  return h->header_buffer.length;
}

size_t fasta_get_sequence_length(fasta_handle h)
{
  return h->sequence_buffer.length;
}

char * fasta_get_header(fasta_handle h)
{
  return h->header_buffer.data;
}

char * fasta_get_sequence(fasta_handle h)
{
  return h->sequence_buffer.data;
}

