`--min_centroid_size` *positive integer*
: Set the minimum abundance a sequence must have to become the
  centroid of a new cluster. A sequence with an abundance smaller than
  *positive integer* can still join an existing cluster, but if it
  matches no centroid it is discarded: it is assigned to no cluster and
  is absent from all output files. The number of sequences discarded in
  this way is reported. Abundances are read from the headers whether or
  not `--sizein` is given. The default value is 1, which lets any
  sequence become a centroid.
