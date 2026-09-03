`--minwordmatches` *non-negative integer*
: Set the minimum number of shared words (i.e. *k*-mers) required for a
  target sequence to be considered further. The default value is 12 for
  the default word length of 8 (see `--wordlength`); for word lengths 3
  to 15 the default values are 18, 17, 16, 15, 14, 12, 11, 10, 9, 8, 7,
  5, and 3, respectively. Neither sequence can share more words than it
  contains, so the requirement is capped by the number of distinct words
  of the query and by that of the target: when either contains fewer
  distinct words than the value above, all the words of that sequence
  must match. If the argument is 0, no word match is required and every
  target sequence is compared to the query.

    Short or heavily masked sequences yield few distinct words, so they
    may share fewer words than required even when their pairwise
    identity is high, and such a match is then never reported. When
    searching or clustering short sequences, lower `--minwordmatches` (1
    is usually enough, and is as sensitive as 0 while much faster), or
    lower `--wordlength`.
