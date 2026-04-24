`--fastq_maxee` *positive real*
: Discard sequences with an expected error greater than *real*. The
  expected error is the sum of error probabilities for all positions
  in the sequence, and is strictly positive (zero or negative
  arguments are rejected, as they would discard every sequence).
  Applied after trimming. See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).
