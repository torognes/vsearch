`--sizeout`
: Add abundance annotations to sequence headers when writing fasta or
  fastq files. Add the pattern `;size=integer`, where the integer is
  the abundance after subsampling. Without `--sizeout`, existing
  abundance annotations are written unchanged, even when only a part
  of a sequence's abundance was sampled.
