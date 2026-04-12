`--sizeorder`
: When a sequence is close to two or more centroids within the
  distance specified by `--id`, resolve the ambiguity by assigning it
  to the centroid with the highest abundance, rather than the closest
  one. Only takes effect when `--maxaccepts` is greater than one. This
  option enables abundance-based greedy clustering (AGC), as opposed
  to the default distance-based greedy clustering (DGC).
