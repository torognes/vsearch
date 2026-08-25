`--gzip_decompress`
: Specify that the *input pipe* is streaming data compressed using
  Lempel-Ziv coding. See `gzip(1)` for more details. This option is
  required when compressed data arrives on standard input through a
  pipe ('-'), where the format cannot be detected without consuming
  the stream. It is not needed when reading from a regular file
  compressed with gzip, nor when such a file is redirected to standard
  input: compression is then detected automatically, and a
  contradicting option is ignored (with a warning when the input is
  standard input). Pipes other than standard input, such as shell
  process substitutions and named FIFOs, are always read as
  uncompressed data; compressed data must arrive on standard input or
  as a named file.
