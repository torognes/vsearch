`--scramble_kmer` *positive non-null integer*
: Length *k* of the words whose counts are preserved by the
  scrambling (range: 1 to 9; default is 1). With the default value,
  each sequence becomes a uniformly random permutation of its own
  nucleotides (mononucleotide scrambling). With *k* = 2 or more, each
  sequence becomes a uniformly random sequence with exactly the same
  counts of all *j*-mers for every *j* <= *k* (dinucleotide
  scrambling, trinucleotide scrambling, ...), obtained by sampling a
  random Eulerian path of the sequence's de Bruijn graph, as in the
  classical Altschul-Erickson dinucleotide shuffle and its *k*-mer
  generalization uShuffle. Three consequences of *k* >= 2 to expect:
  the first and the last *k* - 1 nucleotides always keep their
  positions; a sequence containing at most one *k*-mer (length <=
  *k*) passes through unchanged; and a short or low-complexity
  sequence may admit exactly one valid arrangement, in which case the
  output equals the input --- correct behaviour, not a failure to
  scramble. Increasing *k* preserves more of the local structure and
  therefore scrambles less. For *k* >= 2, memory usage is
  proportional to the sequence length (a few bytes per position,
  transiently).
