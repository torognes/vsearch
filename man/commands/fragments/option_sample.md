`--sample` *string*
: When writing fasta or fastq files, add the the given sample
  identifier *string* to sequence headers. For instance, if *string*
  is 'ABC', the text `;sample=ABC` will be added to the header. Note
  that *string* will be truncated at the first ';' or blank
  character. Other characters (alphabetical, numerical and
  punctuations) are accepted.
