`--sample` *string*
: Add the the given sample identifier *string* to sequence headers
  when writing fasta or fastq files. For instance, if *string* is
  'ABC', the text `;sample=ABC` will be added to the headers. Note
  that *string* will be truncated at the first ';' or blank
  character. Other characters (alphabetical, numerical and
  punctuations) are accepted.
