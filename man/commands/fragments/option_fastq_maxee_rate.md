`--fastq_maxee_rate` *real*
: Discard sequences with an average expected error per base greater
  than *real* (values from 0.0 to 1.0). The average expected error
  per base is the total expected error divided by the sequence length.
  Applied after trimming. See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).
