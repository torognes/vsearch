`--lca_cutoff` *real*
: Set the fraction of hits that must agree at a given taxonomic rank for
  that rank to be included in the last common ancestor (LCA) output (see
  `--lcaout`). The default value is 1.0, requiring all hits to agree. Lower
  values (e.g., 0.95) tolerate a small fraction of disagreeing hits. The
  value must be greater than 0.5 and at most 1.0.
