`--fastq_truncee_rate` *real*
: Truncate reads so that their average expected error per base does
  not exceed *real*. Truncation occurs at the first position where the
  threshold would be exceeded. The average expected error per base is
  the total expected error divided by the length of the truncated
  sequence. See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).
