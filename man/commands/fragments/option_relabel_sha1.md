`--relabel_sha1`
: Replace each sequence header with the SHA1 digest derived from the
  sequence itself. The sequence is converted to upper case, and each
  'U' is replaced with a 'T' before computation of the digest. The
  SHA1 digest is a 160-bit value, represented using a string of 40
  ascii characters. Each pair of characters encodes an hexadecimal
  value, ranging from `x00` to `xff`. See `sha1(3)` for more details
  on SHA1, and `--relabel_md5` for an alternative algorithm. To retain
  annotations, use their corresponding options (`--lengthout`,
  `--eeout`, and `--sizeout`). Use `--relabel_keep` to retain old
  sequence identifiers.
