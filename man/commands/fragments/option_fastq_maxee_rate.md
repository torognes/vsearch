`--fastq_maxee_rate` *non-negative real*
: Discard sequences with an average expected error per base greater
  than *real*. The average expected error per base is the total
  expected error divided by the sequence length, and is always in
  the [0.0, 1.0] range; negative arguments are rejected, and values
  greater than 1.0 effectively disable the filter. Applied after
  trimming. See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).
