`--join_padgapq` *string*
: Use *string* as the quality padding inserted between the forward and
  reverse reads. The default is a string of `I`'s of the same length
  as the sequence padding string (see `--join_padgap`). The letter `I`
  corresponds to a base quality score of 40, indicating a very high
  quality base with error probability of 0.0001 (see
  [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)).
  Accepts an empty
  string, or any other combination of ASCII characters.
