`--query_cov` *real*
: Reject the sequence match if the fraction of the query aligned to
  the target is lower than *real* (value ranging from 0.0 to 1.0
  included). Query coverage is computed as (matches + mismatches) /
  query sequence length, not counting internal or terminal gaps.
