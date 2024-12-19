`--wordlength` *positive integer*
: Set the length of words (i.e. *k*-mers) used for sequence
  comparisons. The range of possible values goes from 3 to 15 (default
  value is 12). Longer words may reduce the sensitivity/recall for
  weak similarities, but can increase precision. On the other hand,
  shorter words may increase sensitivity or recall, but may reduce
  precision.

    Computation time generally increases with shorter words and
    decreases with longer words, but it increases again for very long
    words. Memory requirements for a part of the index increase by a
    factor of 4 each time word length increases by one nucleotide, and
    this generally becomes significant for words longer than 12.
