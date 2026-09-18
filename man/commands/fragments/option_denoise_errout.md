`--denoise_errout` *filename*
: Write the error model to *filename*, as a tab-separated table with
  a header line, one row per transition (`A2A`, `A2C`, ..., `T2T`: true
  nucleotide, then observed nucleotide) and one column per quality
  score, starting at zero. The first sixteen rows have the layout of
  the error matrix of DADA2. With `--denoise_indels model`, two rows
  named `ins` and `del` follow, holding the per-position insertion and
  deletion rates (repeated in each column, as these rates do not
  depend on quality).
