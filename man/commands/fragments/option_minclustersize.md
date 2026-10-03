`--minclustersize` *positive integer*
: Discard clusters smaller than *positive integer* from the outputs
  written once clustering is complete: `--clusters` (no file is
  created for a discarded cluster), `--centroids`, `--consout`,
  `--msaout` and `--profile`. The size of a cluster is the one
  reported in `--uc` cluster records: its number of sequences, or the
  sum of their abundances if `--sizein` is specified. Default is 1 (no
  cluster is discarded). The remaining clusters keep their cluster
  numbers, so `--clusters` file names and `clusterid` values may have
  gaps. A line on stderr and in the `--log` file reports how many
  clusters were discarded.

    `--uc` and the OTU tables (`--otutabout`, `--biomout`,
    `--mothur_shared_out`) are not filtered: they still list every
    cluster. When the value is larger than the `--maxclustersize`
    value, no cluster can satisfy the filter: vsearch issues a warning.
    vsearch also issues a warning when none of the filtered outputs is
    specified. Not to be confused with `--minsize`, which filters input
    sequences by abundance (among the clustering commands,
    `--cluster_unoise` only).
