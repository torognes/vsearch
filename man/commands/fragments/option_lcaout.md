`--lcaout` *filename*
: Write last common ancestor (LCA) information to *filename* in a
  tab-separated format. The first column contains the query identifier and
  the second column contains the deepest taxonomic lineage shared by a
  sufficient fraction of hits (see `--lca_cutoff`). Database sequence headers
  must carry taxonomic annotations in the format used by `--sintax`
  (e.g., `tax=k:Archaea,p:Euryarchaeota,c:Halobacteria`). Set `--maxaccepts`
  to a value greater than 1 for the LCA to be meaningful.
