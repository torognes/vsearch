`--maxclustersize` *positive integer*
: Discard clusters larger than *positive integer* from the outputs
  written once clustering is complete: `--clusters`, `--centroids`,
  `--consout`, `--msaout` and `--profile`. The size of a cluster is
  measured as for `--minclustersize`, and the same rules apply: `--uc`
  and the OTU tables are not filtered, cluster numbers are kept, and a
  line on stderr and in the `--log` file reports how many clusters
  were discarded. By default, no cluster is discarded.
