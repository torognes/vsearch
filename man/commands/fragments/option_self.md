`--self`
: Reject the sequence match if the query and target sequence labels
  are identical. The label is the header up to the first blank, so two
  records sharing an identifier but carrying different descriptions
  still reject each other; with `--notrunclabels` the whole header is
  the label and they no longer do. Use `--selfid` to reject on
  identical sequences rather than identical labels.
