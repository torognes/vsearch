`--db` *filename*
: Search for the oligos (primers, tags, barcodes...) in *filename*, in
  fasta format. Oligos can contain IUPAC ambiguity codes (see
  [`vsearch-nucleotides(7)`](../misc/vsearch-nucleotides.7.md)), and
  must be 1 to 64 nucleotides long: an empty oligo, or an oligo longer
  than 64 nucleotides, is a fatal error. A UDB database is not accepted.
