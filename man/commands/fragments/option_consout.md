`--consout` *filename*
: Write cluster consensus sequences to *filename*. For each cluster, a
  center-star multiple sequence alignment is computed with the centroid
  as the center, using a fast algorithm. The consensus record's header
  is `>centroid=<label>;seqs=<n>`, where *n* is the number of
  sequences in the cluster. Each alignment column within the
  centroid's span contributes its most frequent nucleotide, or nothing
  when gaps outnumber every nucleotide; columns outside the centroid's
  span (terminal extensions contributed by longer members) are always
  omitted. If `--sizein` is specified, sequence abundances are taken
  into account.
