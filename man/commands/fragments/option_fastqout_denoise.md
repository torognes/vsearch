`--fastqout` *filename*
: Write the corrected reads to *filename*, in fastq format (see
  [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). Reads are
  written in their input order, with their original header. The
  sequence of each read is replaced with the center of its partition,
  and its quality string with the mean quality scores of the reads of
  that center. Reads that are not corrected are omitted (see
  `--fastqout_discarded`).
