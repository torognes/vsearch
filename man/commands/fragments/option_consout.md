`--consout` *filename*
: Write cluster consensus sequences to *filename*. For each cluster, a
  center-star multiple sequence alignment is computed with the centroid
  as the center, using a fast algorithm. A consensus sequence is
  constructed by taking the majority symbol (nucleotide or gap) from
  each column of the alignment. Columns containing a majority of gaps
  are skipped, except for terminal gaps. If `--sizein` is specified,
  sequence abundances are taken into account.
