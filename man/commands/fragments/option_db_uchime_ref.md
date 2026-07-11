`--db` *filename*
: Read chimera-free reference sequences from *filename*. Chimeras
  cannot be detected if their parents, or sufficiently close
  relatives, are not present in the database. *filename* must refer
  to a fasta file or to a UDB file. If a UDB file is used, it should
  be created using the `--makeudb_usearch` command with the `--dbmask
  dust` option. This option is mandatory for `--uchime_ref`. `--db`
  accepts `-` to read the database from standard input, as well as an
  explicit stream path such as `/dev/stdin`, a named pipe, or a process
  substitution. The query and the database cannot both be `-`, however,
  as they would compete for the same standard input; give at least one of
  them an explicit path.
