`--qmask` *none*|*dust*|*soft*
: Mask regions in query sequences using the *dust* method or the
  *soft* method, or *none* to suppress masking. Values are
  case-insensitive, so `DUST`, `Dust`, and `dust` are all accepted.
  See [`vsearch-fastx_mask(1)`](./vsearch-fastx_mask.1.md) for more
  details. Warning, when using *soft* masking, search commands become
  case sensitive: masking excludes masked regions from the *k*-mer
  pre-filter that selects candidate targets (the pairwise alignment
  itself always ignores case). A query with no unmasked stretch of at
  least the word length samples no *k*-mers and is therefore compared
  against every database sequence. The default is to mask using
  *dust*.
