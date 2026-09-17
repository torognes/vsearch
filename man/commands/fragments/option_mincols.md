`--mincols` *positive integer*
: Reject the sequence match if the alignment length is shorter than
  *integer* columns. Terminal gaps are not counted here, so the length
  compared is the internal alignment length, the one `--iddef 2`
  divides by. An alignment of 79 columns whose 78 terminal gap columns
  leave a single aligned column therefore counts as one column, and is
  rejected by `--mincols 2`.
