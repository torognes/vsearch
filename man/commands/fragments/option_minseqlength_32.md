`--minseqlength` *non-negative integer*
: Discard sequences shorter than *non-negative integer* (32
  nucleotides by default). A value of 0 retains empty sequences. When
  that value is larger than the `--maxseqlength` value, no sequence
  can pass the length filter: vsearch issues a warning.
