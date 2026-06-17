`--db` *filename*
: Read reference sequences from *filename*, in fasta, fastq, or UDB
  format. UDB files must have been created with `--wordlength` 12 (see
  [`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md)). This
  option is mandatory. Unlike most input options, `--db` does not accept
  `-` to mean standard input: *filename* must be a path that can be
  opened directly. To read the database from a stream, give an explicit
  path such as `/dev/stdin`, a named pipe, or a process substitution.
