`--relabel` *string*
: Replace sequence headers with the prefix *string* and a ticker (1,
  2, 3, etc.). For example, with `--relabel "cluster:"`, the first
  sequence header becomes '>cluster:1', the second sequence header
  becomes '>cluster:2', and so on. To retain annotations, use their
  corresponding options (`--lengthout`, `--eeout`, and
  `--sizeout`). Use `--relabel_keep` to also retain old sequence
  identifiers.

  As a special case, when *string* is exactly `@`, it is not used as a
  literal prefix: the prefix becomes a *sample identifier* derived from
  the name of the input file, and a period is inserted before the
  ticker. With an input file named `sampleA_R1.fastq`, headers become
  '>sampleA.1', '>sampleA.2', and so on. Only the whole argument is
  special: `--relabel "@x"`, `--relabel "x@"` and `--relabel "@@"`
  remain literal prefixes.

  The identifier is derived from the *base name* of the input file, any
  directory part being discarded first: if the base name contains an
  underscore, everything from the first underscore is dropped;
  otherwise, everything from the first period is dropped. The
  underscore takes precedence over the period, so `a.b_c.fastq` yields
  'a.b' and not 'a'. The result is then truncated at the first ';' or
  blank character (space, tab, newline, carriage return, vertical tab
  or form feed), for the reason given for `--sample`, so
  `my sample.fastq` yields 'my' and `a;b.fastq` yields 'a'.

  A base name beginning with an underscore or a period leaves nothing:
  headers become '>.1', '>.2', and so on, and vsearch issues a warning
  saying so. Standard input (`-`) carries no file name to derive from
  and is rejected with a fatal error. Other ways of naming a stream do
  carry a name and are used as given: a process substitution is usually
  seen as `/dev/fd/63`, yielding '>63.1', whereas a named pipe called
  `runA_1.fifo` yields '>runA.1'. A pipeline that needs a stable label
  should pass an explicit prefix rather than `@`.

  This option is provided for compatibility with usearch, where it is
  spelled `-relabel @`. vsearch accepts it wherever `--relabel` is
  accepted, whereas usearch honours it for two commands only and treats
  it as a literal prefix elsewhere (see
  [`vsearch-usearch(7)`](../misc/vsearch-usearch.7.md)). To record a
  sample identifier, `--sample` remains the recommended option: it adds
  a `;sample=` annotation and leaves the sequence identifier intact.
