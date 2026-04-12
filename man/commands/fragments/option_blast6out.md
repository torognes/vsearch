`--blast6out` *filename*
: Write results to *filename* using a BLAST-like tab-separated format
  with twelve fields per query-target match: query label, target
  label, percentage identity, alignment length, mismatches, gap
  openings, query start, query end, target start, target end,
  expectation value (always -1), and bit score (always 0). If
  `--output_no_hits` is used, non-matching queries are also written.
  Note that vsearch uses global pairwise alignments, not BLAST's
  seed-and-extend algorithm.
