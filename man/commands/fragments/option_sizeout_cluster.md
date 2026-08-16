`--sizeout`
: Add abundance annotations to sequence headers when writing fasta or
  fastq files. Add the pattern `;size=integer`. Centroid records
  receive the total abundance of their cluster: the sum of the
  members' abundance annotations when `--sizein` is used, or the
  number of member sequences otherwise.
