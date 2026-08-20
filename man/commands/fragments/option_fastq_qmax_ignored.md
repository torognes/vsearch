`--fastq_qmax` *positive integer*
: Option is ignored: quality scores are not checked when reading fastq
  sequences, so no score is rejected for being too high. The argument
  itself is still validated, as the offset (see `--fastq_ascii`) plus
  the maximal score may not exceed 126, the last printable ASCII
  character.
