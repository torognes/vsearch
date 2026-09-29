`--maxdiffs` *positive integer*
: Report an occurrence only if its alignment has at most *integer*
  differences: substitutions plus gap columns (inserted or deleted
  nucleotides). The default is 2. Any value up to the length of the
  oligo is accepted, but the number of random occurrences grows very
  fast with it (see DESCRIPTION). Oligo nucleotides hanging off a read
  end are not differences (see `--target_cov`). The `diffs` userfield
  reports the quantity compared against.
