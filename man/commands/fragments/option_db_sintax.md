`--db` *filename*
: Read reference sequences from *filename*, in fasta, fastq, or UDB format
  (see [`vsearch-udb(5)`](../formats/vsearch-udb.5.md)). Reference sequence
  headers must carry taxonomic annotations (see DESCRIPTION). This option
  is mandatory. Unlike most input options, `--db` does not accept `-` to
  mean standard input: *filename* must be a path that can be opened
  directly. To read the database from a stream, give an explicit path
  such as `/dev/stdin`, a named pipe, or a process substitution.
