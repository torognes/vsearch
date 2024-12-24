`--xn` *strictly positive real number*
: Set the weight of 'no' votes, corresponding to the parameter *beta*
  in the chimera scoring function. Default value is 8.0. Increasing
  `--xn` reduces the likelihood of tagging a sequence as a chimera
  (less false positives, but also more false negatives). Decreasing
  `--xn` reduces false negative, but increases false positives.
