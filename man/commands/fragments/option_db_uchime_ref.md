`--db` *filename*
: Read chimera-free reference sequences from *filename*. Chimeras
  cannot be detected if their parents, or sufficiently close
  relatives, are not present in the database. *filename* must refer
  to a fasta file or to a UDB file. If a UDB file is used, it should
  be created using the `--makeudb_usearch` command with the `--dbmask
  dust` option. This option is mandatory for `--uchime_ref`.
