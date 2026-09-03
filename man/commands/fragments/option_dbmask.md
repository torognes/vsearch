`--dbmask` *none*|*dust*|*soft*
: Mask regions in the database sequences using the *dust* method or
  the *soft* method, or *none* to suppress masking. See
  [`vsearch-fastx_mask(1)`](./vsearch-fastx_mask.1.md) for more details. Warning,
  when using *soft* masking, search commands become case sensitive:
  masking excludes masked regions from the *k*-mer index used to
  select candidate targets (the pairwise alignment itself always
  ignores case). A database sequence with no unmasked stretch of at
  least the word length contributes no *k*-mer to the index and is
  never selected as a candidate target; it can still be reported when
  the query itself samples no *k*-mer, as such a query is compared
  against every database sequence (see `--qmask`). vsearch warns when
  there are such sequences, giving how many. The default is to mask
  using *dust*.
