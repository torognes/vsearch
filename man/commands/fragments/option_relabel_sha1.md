`--relabel_sha1`
: Relabel sequences using the SHA1 message digest algorithm applied to
  each sequence. It is similar to the `--relabel_md5` option but uses
  the SHA1 algorithm instead of the MD5 algorithm. SHA1 generates a
  160-bit (20-byte) digest that is represented by 20 hexadecimal
  numbers (40 symbols). The probability of a collision (two
  non-identical sequences resulting in the same digest) is smaller for
  the SHA1 algorithm than it is for the MD5 algorithm. Use `--sizeout`
  to conserve the abundance annotations.
