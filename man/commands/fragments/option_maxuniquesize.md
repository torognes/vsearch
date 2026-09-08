`--maxuniquesize` *positive integer*
: Discard sequences with a post-dereplication abundance greater than
  *positive integer*. When that value is smaller than the
  `--minuniquesize` value, no sequence can satisfy the abundance
  filter: vsearch issues a warning and writes no sequence.
