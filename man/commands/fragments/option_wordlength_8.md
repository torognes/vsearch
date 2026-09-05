`--wordlength` *positive integer*
: Set the length of words (i.e. *k*-mers) used for sequence indexing and
  comparisons. Valid values range from 3 to 15. The default is 8. Note
  that the default `--minwordmatches` is derived from this value, so
  changing one changes both (see `--minwordmatches`).

    Longer words make the *k*-mer index more selective, so fewer targets
    are offered as candidates and the search itself gets faster: from
    word length 5 to 11 the search phase shrank by a factor of 2.4 on
    130-nucleotide amplicons and 5.8 on full-length reference sequences.
    Working against that, the index has 4^*wordlength* slots, so the
    memory it needs and the time spent building it both quadruple with
    each added nucleotide. On a 400 000-sequence database the whole run
    needed 0.3 GB at word length 10, 1.3 GB at 12 and 16 GB at 15.

    The best setting balances the two, and depends on how many queries
    are searched against a given database: with few queries the index
    build dominates and a shorter word is cheaper overall, while with
    many queries the search dominates and a longer word repays its index.
    Changing the word length is not output-neutral, so it should be
    chosen for a workload rather than tuned per run.
