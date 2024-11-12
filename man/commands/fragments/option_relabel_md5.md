`--relabel_md5`
: Relabel sequences using the MD5 message digest algorithm applied to
  each sequence. Former sequence headers are discarded. The sequence
  is converted to upper case and each 'U' is replaced by a 'T' before
  computation of the digest.

: The MD5 digest is a cryptographic hash function designed to minimize
  the probability that two different inputs give the same output, even
  for very similar, but non-identical inputs. Still, there is a very
  small, but non-zero, probability that two different inputs give the
  same digest (i.e. a collision). MD5 generates a 128-bit (16-byte)
  digest that is represented by 16 hexadecimal numbers (using 32
  symbols among 0123456789abcdef). Use `--sizeout` to conserve the
  abundance annotations.
