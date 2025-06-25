`--sizein`
: Use the abundance annotations present in sequence headers when
  reading fasta or fastq file. Search for the pattern
  `[>@;]size=integer[;]`. Entries without abundance annotations are
  silently assumed to be of `size=1`. With *de novo* clustering
  commands, option `--sizein` and `--sizeout` are always implied, and
  do not need to be specified.
