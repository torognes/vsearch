`--id` *real*
: Reject the match if the pairwise identity with the target sequence is
  lower than *real* (value ranging from 0.0 to 1.0 included). The
  pairwise identity is defined by default as (matching columns) /
  (alignment length - terminal gaps). That definition can be modified
  with `--iddef`.

    A column holding an ambiguous symbol counts as a matching column
    whenever the two symbols share at least one of the nucleotides they
    represent, so an N matches anything: a query aligned over a run of
    Ns is reported at 100% identity. Use `--n_mismatch` to count these
    columns as mismatches instead.

    Unlike `--usearch_global`, no *k*-mer pre-filter decides which pairs
    reach the alignment stage: every target is aligned, so `--id` alone
    decides the outcome, at any value. That is the purpose of this
    command, and the reason low values behave as asked rather than being
    overruled by the pre-filter. Bear in mind what a low value admits:
    two unrelated sequences of similar length routinely align at 30% to
    45% identity (see the note in the description above).
