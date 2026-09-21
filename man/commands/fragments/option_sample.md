`--sample` *string*
: Add the given sample identifier *string* to sequence headers when
  writing fasta or fastq files. For instance, if *string* is 'ABC',
  the text `;sample=ABC` will be added to the headers. *string* is
  silently truncated at the first ';' or whitespace character
  (space, tab, newline, carriage return, vertical tab or form feed),
  so such characters should not be used in *string*. Other
  characters (alphabetical, numerical and punctuations) are
  accepted. When nothing is left after truncation --- an empty
  *string*, or one starting with ';' or a blank character --- vsearch
  issues a warning and writes a bare `;sample=` annotation.

  This is the recommended way to record a sample identifier, as it
  leaves the sequence identifier intact. The alternative, `--relabel @`,
  replaces the identifier with one derived from the input file name; it
  exists for compatibility with usearch (see `--relabel`).
