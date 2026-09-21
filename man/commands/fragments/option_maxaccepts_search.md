`--maxaccepts` *positive integer*
: Set the maximum number of matching target sequences to accept before
  stopping the search for a given query. The default value is 1. Use
  together with `--maxrejects`. If both `--maxaccepts` and
  `--maxrejects` are set to 0, the complete database is searched, save
  for the targets the word pre-filter removes beforehand (see
  `--minwordmatches`).

    Target sequences are considered in order of decreasing number of
    words shared with the query, a proxy for similarity, and each is
    aligned and then accepted or rejected according to `--id` and the
    other criteria.

    Raising `--maxaccepts` widens the set of targets the query can be
    matched against. All accepted targets are reported, sorted by
    decreasing identity, unless `--maxhits` or `--top_hits_only`
    restricts them. With the default value of 1, the search stops at
    the first acceptable target, which is not necessarily the closest
    one.
