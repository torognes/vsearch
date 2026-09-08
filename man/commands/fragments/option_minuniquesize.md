`--minuniquesize` *positive integer*
: Discard sequences with a post-dereplication abundance smaller than
  *positive integer*. When that value is larger than the
  `--maxuniquesize` value, no sequence can satisfy the abundance
  filter: vsearch issues a warning and writes no sequence.
