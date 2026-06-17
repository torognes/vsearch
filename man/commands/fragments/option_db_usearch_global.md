`--db` *filename*
: Compare query sequences against the target sequences in *filename*, in
  fasta or fastq format. Alternatively, the name of a preformatted UDB
  database created with `--makeudb_usearch` may be specified (see
  [`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md) and
  [`vsearch-udb(5)`](../formats/vsearch-udb.5.md)). Unlike most input
  options, `--db` does not accept `-` to mean standard input: *filename*
  must be a path that can be opened directly. To read the database from
  a stream, give an explicit path such as `/dev/stdin`, a named pipe, or
  a process substitution.
