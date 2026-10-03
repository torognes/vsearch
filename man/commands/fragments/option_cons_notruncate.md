`--cons_notruncate`
: Compute the consensus sequences (`--consout`, `--msaout`) over the
  positions covered by each sequence. By default, the terminal gaps of
  a sequence (gaps before its first or after its last nucleotide in the
  multiple alignment) vote like internal gaps, and columns outside the
  centroid's span are left out of the consensus, so flanking regions
  present in only part of a cluster are lost. With this option,
  terminal gaps do not vote, columns outside the centroid's span are
  part of the consensus, and terminal gaps are not counted in the gap
  column of `--profile`. Internal gaps still vote. In the multiple
  alignment, the part of a sequence extending beyond the centroid's
  start is placed immediately before the centroid's first nucleotide
  (so that extensions of different lengths line up by their end). For
  example, a centroid `ACTGACCCAGT` and five members `ACCCAGT` give
  the consensus `ACCCAGT` by default, and `ACTGACCCAGT` with
  `--cons_notruncate`.
