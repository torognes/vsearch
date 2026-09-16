`--hardmask`
: Replace masked nucleotides with Ns, rather than lowercasing them.

    This is also the only masking option that reaches the alignment.
    Soft and dust masking (see `--qmask` and `--dbmask`) only keep
    masked words out of the *k*-mer pre-filter that selects candidates;
    the masked region is still aligned and scored like any other. An N,
    on the other hand, scores zero and counts as a matching column, so
    hard masking makes a masked region match whatever it is aligned
    against, and the reported identity can only go up. Add
    `--n_mismatch` to count those columns as mismatches instead.
