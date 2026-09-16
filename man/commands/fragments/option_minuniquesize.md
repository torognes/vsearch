`--minuniquesize` *positive integer*
: Discard sequences with a post-dereplication abundance smaller than
  *positive integer*. When that value is larger than the
  `--maxuniquesize` value, no sequence can satisfy the abundance
  filter: vsearch issues a warning and writes no sequence.

    The filter applies to the sequence output only. A cluster summary
    written with `--uc`, where the command offers it, still lists every
    cluster, discarded ones included. This matches usearch.
