`--dbmask` *none*|*dust*|*soft*
: Mask regions in the database sequences using the *dust* method or
  the *soft* method, or *none* to suppress masking. See
  [`vsearch-fastx_mask(1)`](./vsearch-fastx_mask.1.md) for more details. Warning,
  when using *soft* masking, search commands become case sensitive:
  masking excludes masked regions from the *k*-mer index used to
  select candidate targets (the pairwise alignment itself always
  ignores case). The default is to mask using *dust*.
