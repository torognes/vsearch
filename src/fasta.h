struct fasta_s
{
  FILE * fp;

  int format;

#ifdef HAVE_ZLIB
  gzFile fp_gz;
#endif

#ifdef HAVE_BZLIB
  BZFILE * fp_bz;
#endif

  char * file_buffer;
  char * header_buffer;
  char * sequence_buffer;

  size_t file_buffer_length;
  size_t header_buffer_length;
  size_t sequence_buffer_length;

  size_t file_buffer_alloc;
  size_t header_buffer_alloc;
  size_t sequence_buffer_alloc;

  size_t file_buffer_position;
  size_t file_size;
  size_t file_position;

};

typedef struct fasta_s * fasta_handle;

fasta_handle fasta_file_open(const char * filename);
void fasta_file_close(fasta_handle h);
bool fasta_file_read(fasta_handle h, 
                     char * * header, size_t * header_length,
                     char * * sequence, size_t * sequence_length);
size_t fasta_file_getpos(fasta_handle h);
size_t fasta_file_getsize(fasta_handle h);
