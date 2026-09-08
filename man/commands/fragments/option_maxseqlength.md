`--maxseqlength` *positive integer*
: Discard sequences longer than *positive integer* (50,000 nucleotides
  by default). The value must not exceed 2,147,481,646 (`INT_MAX` minus
  2,001). When that value is smaller than the effective
  `--minseqlength` value (whose default is command-specific), no
  sequence can pass the length filter: vsearch issues a warning.
