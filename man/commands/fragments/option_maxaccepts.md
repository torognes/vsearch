`--maxaccepts` *positive integer*
: Set the maximum number of matching target sequences to accept before
  stopping the search for a given query. The default value is 1. Use
  together with `--maxrejects`. If both `--maxaccepts` and
  `--maxrejects` are set to 0, the complete database is searched, save
  for the targets the word pre-filter removes beforehand (see
  `--minwordmatches`).

    Target sequences (when clustering, the centroids of the clusters
    created so far) are considered in order of decreasing number of
    words shared with the query, a proxy for similarity, and each is
    aligned and then accepted or rejected according to `--id` and the
    other criteria. Raising `--maxaccepts` widens the set of targets
    the outcome is chosen from. When clustering, it never places a
    query in several clusters: the query joins the accepted centroid
    with the highest identity, or the most abundant one with
    `--sizeorder`. When searching, all accepted targets are reported,
    sorted by decreasing identity, unless `--maxhits` or
    `--top_hits_only` restricts them.
