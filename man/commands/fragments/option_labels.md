`--labels` *filename*
: Read labels to match from *filename*, one per line. Unless
  `--label_substr_match` is given, each label must match the entire
  header. The comparison is not case-sensitive. Lines longer than
  1023 characters are not supported (a warning is emitted); empty
  lines are skipped.
