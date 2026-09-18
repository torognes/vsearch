`--denoise_indels` ignore|model
: Choose how insertions and deletions are treated by the error model.
  With `ignore` (the default), the model knows substitutions only, as
  in DADA2: gaps in the alignment of a sequence with the center of its
  partition are skipped, so a sequence that differs from its center by
  indels only is at distance zero, is never tested, and is merged into
  its center. With `model`, each interior gap position is an error
  event with a rate of its own (one rate for insertions, one for
  deletions, independent of quality), learnt from the reads along with
  the substitution rates. Sequences differing by indels only are then
  tested like any other, and are kept as distinct sequences when they
  are too abundant to be explained by indel errors. Use `model` for
  markers in which indels are diagnostic, such as 12S and 16S rRNA.
  Terminal gaps are free and are never counted, in either mode.
