`--target_cov` *real*
: Reject the sequence match if the fraction of the target sequence
  aligned to the query is lower than *real*. Target coverage is
  computed as (matches + mismatches) / target sequence length, not
  counting internal or terminal gaps.
