`--denoise_maxconsist` *positive integer*
: Set the maximum number of self-consistency rounds used to learn the
  error model. Learning stops earlier if the model repeats itself. The
  default is 10, as in DADA2 (`MAX_CONSIST`). Ignored when
  `--denoise_errin` is used.
