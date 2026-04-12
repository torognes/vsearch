`--id` *real*
: Reject a sequence if its pairwise identity with the cluster centroid
  is lower than *real* (value ranging from 0.0 to 1.0 included). The
  pairwise identity is defined by default as (matching columns) /
  (alignment length - terminal gaps). That definition can be modified
  with `--iddef`.
