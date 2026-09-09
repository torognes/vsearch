`--n_mismatch`
: Count alignments of nucleotides against Ns as mismatches. By
  default, an alignment column holding an N is neutral: it scores
  zero, and it is counted as a *matching* column when the identity
  percentage is computed. With `--n_mismatch`, any column where at
  least one of the two symbols is an N (regardless of case) is scored
  and counted as a mismatch instead, N against N included. Both the
  alignment score and the post-alignment count of matches and
  mismatches are affected, so the identity percentage compared against
  `--id` (and reported in the output files) changes too.

    The option targets N only: the other ambiguous symbols
    (BDHKMRSVWY) keep their default behaviour, and still count as
    matching any symbol they share a nucleotide with. See
    [`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md)
    for the default treatment of ambiguous symbols, and for why a query
    can align to a long run of Ns with 100% identity.
