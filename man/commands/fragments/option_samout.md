`--samout` *filename*
: Write alignment results to *filename* in the SAM format, see
  [`vsearch-sam(5)`](../formats/vsearch-sam.5.md). Use `--samheader`
  to include header lines. Each non-header line is a SAM record
  representing either a query-target alignment or the absence of a
  match. The alignment column of each record uses the CIGAR format,
  see [`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md).
