`--denoise_errin` *filename*
: Read the error model from *filename* instead of learning it from
  the input reads. The expected format is the one written by
  `--denoise_errout`. A model with fewer quality columns than the
  input requires is extended by repeating its last column. With
  `--denoise_indels model`, the file must contain the `ins` and `del`
  rows. Typical usage is to learn the model once, on a subsample of a
  sequencing run (see
  [`vsearch-fastx_subsample(1)`](./vsearch-fastx_subsample.1.md)),
  and to apply it to each sample of that run.
