`--relabel_md5`
: Replace each sequence header with the MD5 digest derived from the
  sequence itself. The sequence is converted to upper case, and each
  'U' is replaced with a 'T' before computation of the digest. The MD5
  digest is a 128-bit value, represented by a string of 32 ascii
  characters. Each pair of characters encodes an hexadecimal value,
  ranging from `x00` to `xff`. See `md5(3)` for more details on MD5,
  and `--relabel_sha1` for an alternative algorithm. Former headers
  are discarded. Annotations can be conserved with the corresponding
  options (`--lengthout`, `--eeout`, and `--sizeout`).
