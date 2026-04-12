`--iddef` *0|1|2|3|4*
: Change the pairwise identity definition used with `--id`. Accepted
  values are:

  0. CD-HIT definition: (matching columns) / (shortest sequence length).
  1. edit distance: (matching columns) / (alignment length).
  2. edit distance excluding terminal gaps (default definition for
     `--id`).
  3. Marine Biological Lab definition, counting each gap opening
     (internal or terminal) as a single mismatch, whether or not the
     gap was extended: 1.0 - [(mismatches + gap openings)/(longest
     sequence length)].
  4. BLAST definition, equivalent to `--iddef 1` for global pairwise
     alignments.
