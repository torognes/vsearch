`--sizeout`
: Add abundance annotations to sequence headers when writing fasta or
  fastq files. Add the pattern `;size=integer`. Existing `;size=`
  annotations are reported unchanged; entries without one receive
  `;size=1`. For this command `--sizein` is not needed: abundance
  annotations are always parsed from the input headers.
