`--fastq_allowmergestagger`
: Allow merging of staggered read pairs. Staggered pairs arise when a
  fragment shorter than a read is sequenced: the 3' end of one read
  then extends beyond the 5' end of the other. The overhanging portion
  is discarded. By default, staggered pairs are not merged (see
  `--fastq_nostagger`).
