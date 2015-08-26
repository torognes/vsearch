struct buffer_s
{
  char * data;
  size_t length;
  size_t alloc;
  size_t position;
};

struct fasta_s
{
  FILE * fp;

#ifdef HAVE_ZLIB
  gzFile fp_gz;
#endif

#ifdef HAVE_BZLIB
  BZFILE * fp_bz;
#endif

  int format;

  regex_t size_regexp;

  struct buffer_s file_buffer;
  struct buffer_s header_buffer;
  struct buffer_s sequence_buffer;

  size_t file_size;
  size_t file_position;

  size_t lineno;
  long seqno;

  size_t stripped_all;
  size_t stripped[256];
};

typedef struct fasta_s * fasta_handle;

fasta_handle fasta_open(const char * filename);
void fasta_close(fasta_handle h);
bool fasta_next(fasta_handle h,
                bool truncateatspace,
                unsigned int * char_action,
                char * char_mapping);
size_t fasta_get_position(fasta_handle h);
size_t fasta_get_size(fasta_handle h);
size_t fasta_get_lineno(fasta_handle h);
size_t fasta_get_seqno(fasta_handle h);
char * fasta_get_header(fasta_handle h);
char * fasta_get_sequence(fasta_handle h);
size_t fasta_get_header_length(fasta_handle h);
size_t fasta_get_sequence_length(fasta_handle h);
long fasta_get_abundance(fasta_handle h);
