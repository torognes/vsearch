`--sample` *string*
: Add the given sample identifier *string* to sequence headers when
  writing fasta or fastq files. For instance, if *string* is 'ABC',
  the text `;sample=ABC` will be added to the headers. *string* is
  silently truncated at the first ';' or whitespace character
  (space, tab, newline, carriage return, vertical tab or form feed),
  so such characters should not be used in *string*. Other
  characters (alphabetical, numerical and punctuations) are
  accepted.
