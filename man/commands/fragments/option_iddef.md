`--iddef` *0|1|2|3|4|5*
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
  5. score-based definition: 1.0 - [(`--match` * *L* - score) /
     ((`--match` - `--mismatch`) * *L*)], where score is the alignment
     score (userfield `raw`) and *L* the shortest sequence length,
     clamped to the range 0.0 to 1.0. Requires `--match` to be greater
     than `--mismatch`.

  In definitions 0 to 4, a column holding an ambiguous symbol is a
  matching column whenever the two symbols share at least one of the
  nucleotides they represent; an N is thus a match against anything,
  unless `--n_mismatch` is given.

  Definitions 0 to 4 count columns of the chosen alignment; none of
  them reads the alignment score. The scoring options (`--match`,
  `--mismatch`, `--gapopen`, `--gapext`) therefore act on identity only
  indirectly, by changing which alignment is optimal, and a pair whose
  optimal alignment does not change keeps the identity it had. See
  [`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md).

  Definition 5 reads the score instead: *L* * `--match` is the highest
  score the pair can reach, and the amount by which the alignment falls
  short of it is counted in mismatch equivalents, one mismatch costing
  `--match` - `--mismatch`. A substitution thus counts as one mismatch,
  a gap as its penalty divided by `--match` - `--mismatch` (with the
  default scores, an internal gap of length 1 counts as 20 / 6 = 3.3
  mismatches), and a column holding an ambiguous symbol, which scores
  zero, as `--match` / (`--match` - `--mismatch`), a third of a
  mismatch by default. Terminal gaps count at their own penalties: with
  the default scores a sequence that is a prefix of the other is below
  100%, and with `--gapopen 0E --gapext 0E` it is at 100%. For a
  gapless alignment without ambiguous symbols, definition 5 equals
  definitions 0, 1 and 2.
