`--fastq_truncee` *non-negative real*
: Truncate reads so that their cumulative expected error does not
  exceed *real*. Truncation occurs at the first position where the
  threshold would be exceeded. Negative arguments are rejected. See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).
