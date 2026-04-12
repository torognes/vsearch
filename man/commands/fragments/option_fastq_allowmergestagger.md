`--fastq_allowmergestagger`
: Allow merging of staggered read pairs. Staggered pairs arise when a
  very short fragment is sequenced: the 3' end of the reverse read
  extends beyond the 5' end of the forward read. The overhanging
  portion of the reverse read is discarded. By default, staggered
  pairs are not merged (see `--fastq_nostagger`).
