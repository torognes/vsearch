`--sizeout`
: Add abundance annotations to sequence headers when writing fasta or
  fastq files. Add the pattern `;size=integer`, where the integer is
  the total abundance of the cluster: the sum of the members'
  abundance annotations when `--sizein` is used, or the number of
  merged sequences otherwise.
