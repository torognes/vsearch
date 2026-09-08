`--sizeout`
: Add abundance annotations to sequence headers when writing fasta
  files. Add the pattern `;size=integer`, where the integer is the
  abundance read from the header stored in the UDB file when `--sizein`
  is used, or `1` otherwise. Any existing `;size=` annotation in the
  stored header is replaced, so `--sizeout` without `--sizein` reduces
  every entry to `size=1`.
