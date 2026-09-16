`--maxdiffs` *positive integer*
: Reject the sequence match if the alignment contains more than
  *integer* substitutions, insertions, or deletions. Terminal gaps are
  not counted, so a query is never rejected for being shorter than its
  target; only internal differences count. The `diffs` userfield
  reports the quantity compared against.
