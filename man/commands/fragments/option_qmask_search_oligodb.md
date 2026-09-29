`--qmask` *none*|*dust*|*soft*
: Mask regions in the query sequences using the *dust* method or the
  *soft* method, or *none* to suppress masking. Values are
  case-insensitive. The default is *none*, unlike the other search
  commands: the search ignores case, so soft masking has no effect on
  it, and only matters with `--hardmask`, which turns the masked
  nucleotides into Ns. Beware that an N matches any oligo nucleotide
  (see `--n_mismatch`), so hard masking creates occurrences in the
  masked regions rather than hiding them.
