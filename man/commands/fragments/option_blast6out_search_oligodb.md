`--blast6out` *filename*
: Write the occurrences to *filename* using a BLAST-like tab-separated
  format, one line per occurrence, with twelve fields: query label,
  oligo label, percentage identity, alignment length, mismatches, gap
  openings, first and last query positions of the occurrence (swapped
  on the minus strand), first and last oligo positions aligned,
  expectation value (always -1), and bit score (always 0). The
  positions are those of the `qlo`, `qhi`, `tlo` and `thi` userfields.
