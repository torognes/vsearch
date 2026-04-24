`--sample` *string*
: Add the the given sample identifier *string* to sequence headers
  when writing fasta or fastq files. For instance, if *string* is
  'ABC', the text `;sample=ABC` will be added to the headers. The
  *string* is appended verbatim: it is not truncated at the first
  ';' or blank character, so avoid embedding separator characters in
  *string* if downstream tools parse `;sample=...` up to the first
  separator.
