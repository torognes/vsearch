`--scramble_kmer` *positive non-null integer*
: Length *k* of the words whose counts are preserved by the
  scrambling (default is 1). With the default value, each sequence
  becomes a uniformly random permutation of its own nucleotides
  (mononucleotide scrambling). Values greater than 1 --- which would
  preserve the counts of all *j*-mers for *j* <= *k* --- are not
  supported yet and are rejected.
