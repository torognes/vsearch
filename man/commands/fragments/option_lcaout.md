`--lcaout` *filename*
: Write last common ancestor (LCA) information to *filename* in a
  tab-separated format. The first column contains the query identifier and
  the second column contains the deepest taxonomic lineage shared by a
  sufficient fraction of hits (see `--lca_cutoff`). Database sequence headers
  must carry taxonomic annotations in the format used by `--sintax`
  (e.g., `tax=k:Archaea,p:Euryarchaeota,c:Halobacteria`). Set `--maxaccepts`
  to a value greater than 1 for the LCA to be meaningful. The vote runs
  over the reported hits (after `--maxhits`/`--top_hits_only`,
  including weak hits from `--weak_id`). A line is written for every
  query: queries without hits produce a line with an empty second
  column, whether or not `--output_no_hits` is given.
