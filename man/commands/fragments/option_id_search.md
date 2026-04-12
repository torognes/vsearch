`--id` *real*
: Reject the match if the pairwise identity with the target sequence is
  lower than *real* (value ranging from 0.0 to 1.0 included). The
  pairwise identity is defined by default as (matching columns) /
  (alignment length - terminal gaps). That definition can be modified
  with `--iddef`.
