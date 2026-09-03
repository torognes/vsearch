`--id` *real*
: Reject a sequence if its pairwise identity with the cluster centroid
  is lower than *real* (value ranging from 0.0 to 1.0 included). The
  pairwise identity is defined by default as (matching columns) /
  (alignment length - terminal gaps). That definition can be modified
  with `--iddef`.

    Which pairs reach the alignment stage where `--id` is applied is
    decided beforehand by a *k*-mer pre-filter (see `--minwordmatches`
    and `--wordlength`): a pair sharing too few words is never aligned
    and never reported, whatever its identity. Short or heavily masked
    sequences share few words, so lowering `--id` alone does not make
    them match; lower `--minwordmatches` too. Below an `--id` of about
    0.5, it is the pre-filter rather than `--id` that decides the
    outcome.
