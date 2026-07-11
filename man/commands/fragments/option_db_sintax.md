`--db` *filename*
: Read reference sequences from *filename*, in fasta, fastq, or UDB format
  (see [`vsearch-udb(5)`](../formats/vsearch-udb.5.md)). Reference sequence
  headers must carry taxonomic annotations (see DESCRIPTION). This option
  is mandatory. `--db` accepts `-` to read the database from standard
  input, as well as an explicit stream path such as `/dev/stdin`, a named
  pipe, or a process substitution. The query and the database cannot both
  be `-`, however, as they would compete for the same standard input; give
  at least one of them an explicit path.
