`--relabel_md5`
: Replace each sequence header with the MD5 digest derived from the
  sequence itself. The sequence is converted to upper case, and each
  'U' is replaced with a 'T' before computation of the digest. The MD5
  digest is a 128-bit value, represented using a string of 32 ASCII
  characters. Each pair of characters encodes an hexadecimal value,
  ranging from `x00` to `xff`. See `md5(3)` for more details on MD5,
  and `--relabel_sha1` for an alternative algorithm. To retain
  annotations, use their corresponding options (`--lengthout`,
  `--eeout`, and `--sizeout`). Use `--relabel_keep` to retain old
  sequence identifiers.
