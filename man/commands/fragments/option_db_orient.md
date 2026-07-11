`--db` *filename*
: Read reference sequences from *filename*, in fasta, fastq, or UDB
  format. UDB files must have been created with `--wordlength` 12 (see
  [`vsearch-makeudb_usearch(1)`](./vsearch-makeudb_usearch.1.md)). This
  option is mandatory. `--db` accepts `-` to read the database from
  standard input, as well as an explicit stream path such as
  `/dev/stdin`, a named pipe, or a process substitution. The query and
  the database cannot both be `-`, however, as they would compete for the
  same standard input; give at least one of them an explicit path.
