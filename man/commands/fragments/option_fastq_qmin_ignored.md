`--fastq_qmin` *positive integer*
: Option is ignored: quality scores are not checked when reading fastq
  sequences, so no score is rejected for being too low. The argument
  itself is still validated, as the offset (see `--fastq_ascii`) plus
  the minimal score must be at least 33, the first printable ASCII
  character.
