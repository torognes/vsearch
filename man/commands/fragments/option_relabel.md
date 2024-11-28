`--relabel` *string*
: Replace sequence headers with the prefix *string* and a ticker (1,
  2, 3, etc.). For example, with `--relabel "cluster:"`, the first
  sequence header becomes '>cluster:1', the second sequence header
  becomes '>cluster:2', and so on. To retain annotations, use their
  corresponding options (`--lengthout`, `--eeout`, and
  `--sizeout`). Use `--relabel_keep` to retain old sequence
  identifiers.
