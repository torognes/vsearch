`--maxaccepts` *positive integer*
: Set the maximum number of matching target sequences to accept before
  stopping the search for a given query. The default value is 1. Use
  together with `--maxrejects`. If both `--maxaccepts` and
  `--maxrejects` are set to 0, the complete database is searched, save
  for the targets the word pre-filter removes beforehand (see
  `--minwordmatches`).

    The target sequences are the centroids of the clusters created so
    far. They are considered in order of decreasing number of words
    shared with the query, a proxy for similarity, and each is aligned
    and then accepted or rejected according to `--id` and the other
    criteria.

    Raising `--maxaccepts` never places a sequence in several clusters:
    it widens the set of centroids the outcome is chosen from. The
    sequence joins the accepted centroid with the highest identity, or
    the most abundant one with `--sizeorder`. With the default value of
    1, the search stops at the first acceptable centroid, which is not
    necessarily the closest one.
