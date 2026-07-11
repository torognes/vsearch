`--db` *filename*
: Compare query sequences against the target sequences in *filename*, in
  fasta or fastq format. Alternatively, the name of a preformatted UDB
  database created with `--makeudb_usearch` may be specified (see
  [`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md) and
  [`vsearch-udb(5)`](../formats/vsearch-udb.5.md)). `--db` accepts `-` to
  read the database from standard input, as well as an explicit stream
  path such as `/dev/stdin`, a named pipe, or a process substitution. The
  query and the database cannot both be `-`, however, as they would
  compete for the same standard input; give at least one of them an
  explicit path.
