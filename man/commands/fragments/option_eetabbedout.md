`--eetabbedout` *filename*
: Write expected error statistics for each successfully *merged* pair
  to *filename* (pairs that could not be merged are not reported), in a
  tab-separated format with four columns: the expected errors in the
  forward read, the expected errors in the reverse read, the observed
  differences in the forward read within the overlap region, and the
  observed differences in the reverse read within the overlap region.
  See
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md).
